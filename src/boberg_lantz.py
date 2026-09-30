"""
Boberg-Lantz cyclic steam stimulation model.

This module implements the workflow used in Example 8.6 of:

    Green, D. W. and Willhite, G. P.
    Enhanced Oil Recovery, 2nd ed.

The implementation calculates the heated-zone geometry, temperature-
dependent oil properties, productivity enhancement, heat removal,
water-oil ratio evolution, and production response after steam injection.

Notes
-----
The published worked example obtains radial dimensionless temperatures
from a graphical Boberg-Lantz solution. This implementation evaluates
the stated analytical approximation directly so that the workflow is
fully computational. Small differences from the published tabulated
trajectory are therefore expected.
"""

from dataclasses import dataclass
import math
from typing import Iterable, Tuple

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class ModelInputs:
    """Input parameters for the Boberg-Lantz calculation."""

    h: float = 41.0
    T_s: float = 381.82
    T_r: float = 93.0
    Ms: float = 41.0
    Mr: float = 38.5
    api: float = 13.5
    mu_oc: float = 2026.0
    B_oc: float = 1.01
    r_e: float = 233.5
    r_w: float = 0.292
    W_is: float = 1000.0
    t_inj: float = 14.0
    fsd: float = 0.7
    alpha: float = 0.731
    H_wT: float = 355.4
    L_vdh: float = 843.3
    H_wr: float = 61.0
    qoi: float = 5.79
    Pres: float = 125.0
    BHP: float = 25.0

    # ASTM oil-property correlation constants
    a: float = 16.1368
    b: float = 5.6934
    beta: float = 3.5e-4

    # Produced-fluid properties
    Fwoc: float = 0.25
    C_w: float = 1.016
    rho_w: float = 54.4
    C_o: float = 0.5


DEFAULT_TIMESTEPS = (
    4, 60, 90, 120, 150, 180, 210, 240, 260,
    290, 320, 350, 380, 410, 440, 470, 500, 530,
)


def calculate_initial_oil_density(api: float) -> float:
    """Return initial oil specific gravity from API gravity."""
    return 141.5 / (api + 131.5)


def calculate_oil_density(rho_oi: float, T: float, beta: float = 3.5e-4) -> float:
    """Return oil specific gravity at temperature T in degF."""
    return rho_oi * math.exp(-beta * (T - 60.0))


def calculate_oil_viscosity(
    rho_o: float,
    T: float,
    a: float = 16.1368,
    b: float = 5.6934,
) -> float:
    """Return oil viscosity in cp from the ASTM correlation."""
    lhs = a - b * math.log10(T + 460.0)
    return rho_o * (10.0 ** (10.0 ** lhs) - 0.6)


def G(tD: float) -> float:
    """Marx-Langenheim dimensionless heat function."""
    if tD < 0.0:
        raise ValueError("tD must be non-negative.")
    return (
        math.exp(tD) * math.erfc(math.sqrt(tD))
        + 2.0 * math.sqrt(tD / math.pi)
        - 1.0
    )


def calculate_TDr(tDr: float) -> float:
    """
    Calculate radial dimensionless temperature.

    This is the truncated analytical approximation used by the
    computational implementation. The original worked example reads
    radial temperature values from a graphical solution.
    """
    if tDr <= 0.0:
        raise ValueError("tDr must be positive.")

    return 1.0 - math.sqrt(tDr / math.pi) * (
        2.0
        - tDr / 2.0
        - 3.0 * tDr**2 / 16.0
        - 15.0 * tDr**3 / 64.0
        - 525.0 * tDr**4 / 1024.0
    )


def calculate_TDz(tDz: float) -> float:
    """Calculate vertical dimensionless temperature."""
    if tDz <= 0.0:
        raise ValueError("tDz must be positive.")

    return (
        math.erf(1.0 / math.sqrt(tDz))
        - math.sqrt(tDz / math.pi) * (1.0 - math.exp(-1.0 / tDz))
    )


def calculate_Fwo(Wp: float, W_is_total: float, Fwoc: float = 0.25) -> float:
    """Calculate water-oil ratio after stimulation."""
    if W_is_total <= 0.0:
        raise ValueError("W_is_total must be positive.")
    return Fwoc + 1.83 * math.exp(-7.38 * (Wp / W_is_total))


def calculate_geometry(inputs: ModelInputs) -> dict:
    """Calculate steam heat content and heated-zone geometry."""
    M = (inputs.Ms + inputs.Mr) / 2.0
    H_s = inputs.H_wT + inputs.fsd * inputs.L_vdh - inputs.H_wr
    tD = 4.0 * inputs.alpha * inputs.t_inj / inputs.h**2
    G_tD = G(tD)

    # Steam/water mass injection rate, lbm/day
    m_s_rate = inputs.W_is * 5.615 * 62.4

    A_h = (
        (m_s_rate * H_s * inputs.h)
        / (4.0 * (inputs.T_s - inputs.T_r) * inputs.alpha * M)
        * G_tD
    )

    r_h = math.sqrt(A_h / math.pi)

    z = (
        (m_s_rate * H_s * inputs.t_inj)
        / (
            math.pi
            * r_h**2
            * (inputs.T_s - inputs.T_r)
            * M
        )
        - inputs.h
    )

    return {
        "M": M,
        "H_s": H_s,
        "tD_injection": tD,
        "G_tD": G_tD,
        "m_s_rate": m_s_rate,
        "A_h": A_h,
        "r_h": r_h,
        "z": z,
    }


def calculate_dimensionless_times(
    t: float,
    ti: float,
    alpha: float,
    r_h: float,
    h: float,
    z: float,
) -> Tuple[float, float]:
    """Return radial and vertical dimensionless times."""
    elapsed = t - ti
    if elapsed <= 0.0:
        raise ValueError("t - ti must be positive.")

    tD_r = alpha * elapsed / r_h**2
    tD_z = alpha * elapsed / (((h + z) / 2.0) ** 2)
    return tD_r, tD_z


def run_model(
    inputs: ModelInputs = ModelInputs(),
    timesteps: Iterable[float] = DEFAULT_TIMESTEPS,
) -> tuple[pd.DataFrame, dict]:
    """
    Run the cyclic-steam-stimulation production-response calculation.

    Returns
    -------
    results : pandas.DataFrame
        One row per production time.
    geometry : dict
        Heated-zone geometry and steam-heat quantities.
    """
    timesteps = list(timesteps)
    if len(timesteps) < 2:
        raise ValueError("At least two timesteps are required.")
    if any(t2 <= t1 for t1, t2 in zip(timesteps[:-1], timesteps[1:])):
        raise ValueError("Timesteps must be strictly increasing.")

    geometry = calculate_geometry(inputs)

    H_s = geometry["H_s"]
    m_s_rate = geometry["m_s_rate"]
    r_h = geometry["r_h"]
    z = geometry["z"]

    delta_N_p_o = 0.0
    delta_N_p_w = 0.0
    Wp = 0.0
    Np = 0.0
    CumQ_pi = 0.0

    deltaP = inputs.Pres - inputs.BHP
    Jc = inputs.qoi / deltaP

    delta = 0.0
    T_p = inputs.T_s
    T_Dp = 1.0

    rows = []

    for i, timestep in enumerate(timesteps):

        if i == 0:
            # The worked example assumes Tp = Ts at the initial 4-day
            # production point because the dimensionless times are below
            # the graphical range used in the source.
            T_D = 1.0
            T_Dp = T_D
            T_p = inputs.T_s

            tD_r = np.nan
            TD_r = np.nan
            tD_z = np.nan
            TD_z = np.nan

        else:
            # Add production from the preceding interval.
            Np += delta_N_p_o
            Wp += delta_N_p_w

            # Temperature entering the present production interval comes
            # from the corrected dimensionless temperature calculated at
            # the preceding row.
            T_p = inputs.T_r + (inputs.T_s - inputs.T_r) * T_Dp

        F_wo = calculate_Fwo(
            Wp,
            inputs.W_is * inputs.t_inj,
            inputs.Fwoc,
        )

        rho_o_i = calculate_initial_oil_density(inputs.api)
        rho_o = calculate_oil_density(rho_o_i, T_p, inputs.beta)
        mu_oh = calculate_oil_viscosity(rho_o, T_p, inputs.a, inputs.b)
        B_oh = (rho_o_i / rho_o) * inputs.B_oc

        J = (
            (inputs.B_oc / B_oh)
            * np.log(inputs.r_e / inputs.r_w)
            / (
                (mu_oh / inputs.mu_oc) * np.log(r_h / inputs.r_w)
                + np.log(inputs.r_e / r_h)
            )
        )

        # Reservoir-volume and surface-volume oil rates.
        q_oh = J * Jc * deltaP * B_oh
        q_o = q_oh / B_oh

        # Heat-removal rate.
        Q_pi = (
            5.615
            * q_oh
            * (
                rho_o * 62.4 * inputs.C_o
                + F_wo * inputs.rho_w * inputs.C_w
            )
            * (T_p - inputs.T_r)
        )

        # Duration represented by this row.
        if i < len(timesteps) - 1:
            delta_t = timesteps[i + 1] - timestep
        else:
            delta_t = timesteps[i] - timesteps[i - 1]

        CumQ_pi += Q_pi * delta_t

        delta = CumQ_pi / (
            2.0 * m_s_rate * H_s * inputs.t_inj
        )

        delta_N_p_o = q_o * delta_t
        delta_N_p_w = delta_N_p_o * F_wo

        if i == 0:
            T_Dp = T_D * (1.0 - delta) - delta
        else:
            tD_r, tD_z = calculate_dimensionless_times(
                timestep,
                0.0,
                inputs.alpha,
                r_h,
                inputs.h,
                z,
            )

            TD_r = calculate_TDr(tD_r)
            TD_z = calculate_TDz(tD_z)
            T_D = TD_r * TD_z
            T_Dp = T_D * (1.0 - delta) - delta

        rows.append(
            {
                "t - ti (days)": timestep,
                "Tp (F)": T_p,
                "mu_oh (cp)": mu_oh,
                "qo (STB/D)": q_o,
                "qoh (RB/D)": q_oh,
                "Fwo (bbl/bbl)": F_wo,
                "Qpi (Btu/D)": Q_pi,
                "CumQpi (Btu)": CumQ_pi,
                "delta": delta,
                "tDr": tD_r,
                "TDr": TD_r,
                "tDz": tD_z,
                "TDz": TD_z,
                "TD": T_D,
                "TDp": T_Dp,
                "dNp (bbl)": delta_N_p_o,
                "dWp (bbl)": delta_N_p_w,
                "Wp begin (bbl)": Wp,
                "Np begin (bbl)": Np,
            }
        )

    return pd.DataFrame(rows), geometry
