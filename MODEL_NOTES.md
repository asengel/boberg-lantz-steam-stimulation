# Model formulation

This repository implements a computational version of the **Boberg–Lantz cyclic steam stimulation model** used in Example 8.6 of Green and Willhite, *Enhanced Oil Recovery*, Second Edition.

The equations below use GitHub-supported Markdown math syntax so they render directly on GitHub.

## 1. Steam heat content

The effective heat carried into the reservoir is

$$
H_s = H_{wT} + f_{sd}L_{vdh} - H_{wr}
$$

where:

- $H_s$ = effective steam heat input, Btu/lbm
- $H_{wT}$ = water enthalpy at steam-injection temperature, Btu/lbm
- $f_{sd}$ = fraction of steam condensed in the reservoir
- $L_{vdh}$ = latent heat of vaporization, Btu/lbm
- $H_{wr}$ = water enthalpy at reservoir temperature, Btu/lbm

For Example 8.6,

$$
H_s \approx 884.7 \ \mathrm{Btu/lbm}
$$

## 2. Heated-zone geometry

The average volumetric heat capacity used in the Boberg–Lantz model is

$$
M = \frac{M_s + M_r}{2}
$$

The dimensionless injection time is

$$
t_D = \frac{4\alpha t_{\mathrm{inj}}}{h^2}
$$

The Marx–Langenheim heat function is

$$
G(t_D)
=
e^{t_D}\operatorname{erfc}\left(\sqrt{t_D}\right)
+
2\sqrt{\frac{t_D}{\pi}}
-1
$$

The steam/water mass-injection rate is

$$
\dot{m}_s
=
W_{is}(5.615)(62.4)
$$

where $W_{is}$ is in bbl/day and $\dot{m}_s$ is in lbm/day.

The heated area is

$$
A_h
=
\frac{\dot{m}_s H_s h}
{4(T_s-T_r)\alpha M}
G(t_D)
$$

and the equivalent heated radius is

$$
r_h
=
\sqrt{\frac{A_h}{\pi}}
$$

The incremental thickness used to account for heat transferred outside the reservoir interval during injection is

$$
z
=
\frac{\dot{m}_s H_s t_{\mathrm{inj}}}
{\pi r_h^2 (T_s-T_r)M}
-h
$$

For Example 8.6, the implementation gives approximately

$$
r_h \approx 51.21 \ \mathrm{ft}
$$

## 3. Dimensionless cooling

The radial dimensionless time is

$$
t_{Dr}
=
\frac{\alpha(t-t_i)}{r_h^2}
$$

and the vertical dimensionless time is

$$
t_{Dz}
=
\frac{\alpha(t-t_i)}
{\left[(h+z)/2\right]^2}
$$

The radial dimensionless temperature is evaluated with the implemented truncated analytical approximation:

$$
T_{Dr}
=
1
-
\sqrt{\frac{t_{Dr}}{\pi}}
\left(
2
-\frac{t_{Dr}}{2}
-\frac{3t_{Dr}^2}{16}
-\frac{15t_{Dr}^3}{64}
-\frac{525t_{Dr}^4}{1024}
\right)
$$

The vertical dimensionless temperature is

$$
T_{Dz}
=
\operatorname{erf}
\left(
\frac{1}{\sqrt{t_{Dz}}}
\right)
-
\sqrt{\frac{t_{Dz}}{\pi}}
\left(
1-e^{-1/t_{Dz}}
\right)
$$

The combined dimensionless temperature is

$$
T_D = T_{Dr}T_{Dz}
$$

At the first production point in Example 8.6, the published solution assumes

$$
T_D = 1
$$

because the dimensionless times are below the graphical range used in the worked example.

## 4. Oil density and viscosity

Initial oil specific gravity is calculated from API gravity:

$$
\rho_{oi}
=
\frac{141.5}{131.5 + API}
$$

The temperature-dependent oil specific gravity is

$$
\rho_o
=
\rho_{oi}
\exp\left[-\beta(T-60)\right]
$$

The ASTM viscosity correlation is written as

$$
\log_{10}
\left[
\log_{10}
\left(
\frac{\mu_o}{\rho_o}+0.6
\right)
\right]
=
a-b\log_{10}(T+460)
$$

Solving explicitly for viscosity gives

$$
\mu_o
=
\rho_o
\left[
10^{\,10^{\,a-b\log_{10}(T+460)}}-0.6
\right]
$$

where $\mu_o$ is in cp.

## 5. Oil formation volume factor

The hot-oil formation volume factor is calculated from the density change:

$$
B_{oh}
=
\frac{\rho_{oi}}{\rho_o}B_{oc}
$$

## 6. Productivity response

The cold productivity index is

$$
J_c
=
\frac{q_{oi}}{P_{\mathrm{res}}-P_{\mathrm{bh}}}
$$

The ratio of stimulated to cold productivity is

$$
\frac{J_h}{J_c}
=
\frac{
\left(B_{oc}/B_{oh}\right)
\ln\left(r_e/r_w\right)
}{
\left(\mu_{oh}/\mu_{oc}\right)
\ln\left(r_h/r_w\right)
+
\ln\left(r_e/r_h\right)
}
$$

In the code this ratio is stored as `J`.

The hot reservoir-volume oil rate is

$$
q_{oh}
=
\left(\frac{J_h}{J_c}\right)
J_c
\Delta P
B_{oh}
$$

and the surface oil rate is

$$
q_o
=
\frac{q_{oh}}{B_{oh}}
$$

The model therefore reports both:

- $q_o$ in STB/D
- $q_{oh}$ in RB/D

## 7. Water-oil ratio

The post-stimulation water-oil ratio is calculated from cumulative produced water:

$$
F_{wo}
=
F_{woc}
+
1.83
\exp
\left(
-7.38\frac{W_p}{W_{is,\mathrm{total}}}
\right)
$$

where

$$
W_{is,\mathrm{total}}
=
W_{is}t_{\mathrm{inj}}
$$

and $W_p$ is cumulative produced water.

## 8. Heat removal by produced fluids

The produced-fluid heat-removal rate is

$$
\dot{Q}_{pi}
=
5.615 \,
q_{oh}
\left(
\rho_o C_o
+
F_{wo}\rho_w C_w
\right)
(T_p-T_r)
$$

In the implementation, oil specific gravity is converted to lbm/ft³ using the factor 62.4 before entering this equation.

For a timestep $\Delta t_i$, cumulative heat removal is approximated by

$$
\sum_i
\dot{Q}_{pi}\Delta t_i
$$

The fraction of the initially injected heat removed by produced fluids is

$$
\delta
=
\frac{
\sum_i \dot{Q}_{pi}\Delta t_i
}{
2\dot{m}_s H_s t_{\mathrm{inj}}
}
$$

## 9. Produced-fluid temperature correction

The corrected dimensionless produced-fluid temperature is

$$
T_{Dp}
=
T_D(1-\delta)-\delta
$$

The produced-fluid temperature used for the next timestep is

$$
T_p
=
T_r
+
(T_s-T_r)T_{Dp}
$$

This corrected temperature controls the oil density, viscosity, formation volume factor, productivity response, and subsequent heat-removal calculation.

## 10. Time-stepping sequence

For each production timestep, the implementation:

1. adds oil and water produced during the preceding interval to the cumulative totals;
2. updates $F_{wo}$ from cumulative produced water;
3. calculates the current produced-fluid temperature $T_p$;
4. calculates oil density, viscosity, and $B_{oh}$;
5. evaluates the stimulated productivity ratio;
6. calculates $q_o$ and $q_{oh}$;
7. calculates heat removal $\dot{Q}_{pi}$;
8. updates cumulative heat removal and $\delta$;
9. calculates $t_{Dr}$ and $t_{Dz}$;
10. calculates $T_{Dr}$, $T_{Dz}$, and $T_D$;
11. calculates the corrected $T_{Dp}$ for the next timestep.

## Numerical note

The published Example 8.6 obtains radial dimensionless-temperature values from the Boberg–Lantz graphical solution.

This repository intentionally **calculates** the radial dimensionless-temperature response rather than hard-coding the values from the published table. Consequently, the first worked interval is reproduced closely, while later temperature, viscosity, and production values can differ modestly from the printed table.

The implemented $T_{Dr}$ expression is a truncated analytical approximation. Its range and accuracy should be checked before applying the code substantially outside the conditions of the validation example.

## Reference

Green, Don W., and G. Paul Willhite. *Enhanced Oil Recovery*. Second Edition. Example 8.6, “Estimation of Steam Stimulation Response With Boberg and Lantz Model,” pp. 538–542.
