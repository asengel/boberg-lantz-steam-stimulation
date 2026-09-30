# Model formulation

This repository implements a computational version of the Boberg-Lantz
cyclic steam stimulation workflow used in Example 8.6 of Green and
Willhite, *Enhanced Oil Recovery*, Second Edition.

## 1. Steam heat content

The effective heat carried into the reservoir is

\[
H_s = H_{wT} + f_{sd} L_{vdh} - H_{wr}.
\]

## 2. Heated-zone geometry

The injection-period dimensionless time is

\[
t_D = \frac{4\alpha t_{\rm inj}}{h^2}.
\]

The Marx-Langenheim heat function is evaluated as

\[
G(t_D)
=
e^{t_D}\operatorname{erfc}(\sqrt{t_D})
+
2\sqrt{\frac{t_D}{\pi}}
-1.
\]

The heated area is

\[
A_h =
\frac{\dot m_s H_s h}
{4(T_s-T_r)\alpha M}
G(t_D),
\]

and the equivalent heated radius is

\[
r_h=\sqrt{\frac{A_h}{\pi}}.
\]

An incremental thickness \(z\) accounts for heat transferred outside the
reservoir interval during steam injection.

## 3. Dimensionless cooling

Radial and vertical dimensionless times are

\[
t_{Dr} = \frac{\alpha(t-t_i)}{r_h^2},
\]

\[
t_{Dz} =
\frac{\alpha(t-t_i)}
{\left[(h+z)/2\right]^2}.
\]

The computational model evaluates \(T_{Dr}\) using the implemented
truncated analytical approximation and evaluates \(T_{Dz}\) directly.
The combined dimensionless temperature is

\[
T_D=T_{Dr}T_{Dz}.
\]

## 4. Oil density and viscosity

Initial oil specific gravity is calculated from API gravity:

\[
\rho_{oi}=\frac{141.5}{131.5+API}.
\]

Temperature-dependent oil density is

\[
\rho_o
=
\rho_{oi}
\exp[-\beta(T-60)].
\]

The ASTM double-log viscosity correlation is then evaluated at the
current produced-fluid temperature.

## 5. Productivity response

The ratio of stimulated to cold productivity is evaluated from the
heated and unheated radial-flow regions. The model calculates both:

- surface oil rate, \(q_o\), in STB/D;
- hot reservoir-volume oil rate, \(q_{oh}\), in RB/D.

## 6. Heat removal by produced fluids

At each production timestep, the energy removed with produced oil and
water is accumulated. The heat-removal fraction is

\[
\delta
=
\frac{\sum \dot Q_{pi}\Delta t_i}
{2\dot m_s H_s t_{\rm inj}}.
\]

The corrected dimensionless produced-fluid temperature is

\[
T_{Dp}=T_D(1-\delta)-\delta,
\]

and the next produced-fluid temperature is

\[
T_p=T_r+(T_s-T_r)T_{Dp}.
\]

## 7. Water-oil ratio

The post-stimulation WOR correlation is updated from cumulative produced
water. This affects the sensible heat removed by the produced water.

## Numerical note

The published Example 8.6 uses graphically obtained radial
dimensionless-temperature values in its tabulated calculation. This
repository intentionally calculates the radial response rather than
hard-coding those tabulated values. As a result, the first worked
interval is closely reproduced, while the later temperature trajectory
can differ modestly from the printed table.

This distinction is documented rather than tuned away so that the code
remains computational and reusable.
