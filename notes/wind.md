# Solar wind: how global models drive it (survey) and options for CZTR_3D

Survey of 2026-11-09 (a Claude subagent with web search; my summary). Tags:
[V] read in the source, [S] read in a citing paper or abstract, [R] recalled,
not verified. Full texts not obtained: Groth 2000, Manchester 2004, Roussev
2003, Cohen 2007, Usmanov 2000/2011, Leer & Holzer 1980, Hansteen & Leer
1995, Withbroe 1988, Hammer 1982, Mikic 1999, Lionello 2001/2009.

Why: atm_relax (2026-11-09) with a closed wall at 2.5 R: the outer corona
cools and falls (-86 km/s at 0.43 R up). The user: the missing wind is the
cause, and at 2.5 R the wind is still subsonic and sub-Alfvenic (sonic point
~4.6 R at 1.5 MK, ~7 R at 1 MK), so the top must move to ~20 R, where every
characteristic leaves, with a coarse radial grid out there.

## 1. MAS / Predictive Science (thermodynamic MHD)

Torok et al. 2018, App. A (arXiv:1801.05903) [V]:
- WKB Alfven waves (undamped, acceleration only; Jacques 1977, Mikic 1999):
  d eps+-/dt = -div[eps+- (v +- v_A b)] - (eps+-/2) div v; momentum gets
  -(1/rho) grad[(eps+ + eps-)/2].
- Empirical heating (sum of Lionello 2009 terms); Spitzer conduction switched
  off by f_c = 0.5 (1 - tanh[(r - 10 R)/0.5 R]) and replaced by a
  collisionless flux ~ f_nc(r) T v beyond; the TR broadened: kappa and Q(T)
  ~ (T/T_cut)^(5/2) below T_cut = 5e5 K.
- 1-20 R, dr 4.5e-4 R in the TR to 0.4 R at 20 R; base T 2e4 K,
  n 2e12 cm^-3; characteristic BCs (Linker & Mikic 1997) at both ends;
  relaxed 160 tau_A ~ 64 h. The wind speed "is predominantly determined by
  the coronal heating and the Alfven-wave pressure" (came out a bit fast).
- Empirical heating (Lionello 2009, as in Downs 2010) [V]: quiet Sun
  H0 exp(-(r - R)/lambda), 7.28e-5 erg cm^-3 s^-1, lambda 40 Mm; coronal
  hole/extended 5e-7 erg cm^-3 s^-1 with lambda 0.7 R (everywhere); AR
  4e-5 |B|; open cells (PFSS tracing) get only the long term. Downs 2010:
  outer boundary 24 R, zero gradient.
- Polytropic MAS (gamma = 1.05) [S]: right morphology, wrong temperatures.
- Wave-turbulence-driven MAS, Lionello 2014 (arXiv:1402.4188, 1D) [V]:
  H = rho (|z-| z+^2 + |z+| z-^2)/(4 lambda_perp), lambda_perp = lambda_sun
  sqrt(A/A_sun); dz+-/dt = -(U +- V_a) dz+-/ds + R1 z+- + R2 z-+ - |z-+| z+-/
  (2 lambda_perp); base n0 2e12, T0 2e4 K (lower n: evaporation risk);
  z_sun 23-31 km/s with lambda_sun 0.03-0.07 R fit observations (25 km/s,
  0.05 R: 753 km/s, n 2.06 cm^-3 at 1 AU); relaxation 80 + 40 + 200 h.
  Downs 2016 [V]: n_e0 6e12, T0 1.75e4 K, waves injected at n_e = 2e11 with
  z0 ~10-15 km/s; closed loops show thermal non-equilibrium (e.g. 4.6 h).

## 2. SWMF / BATS-R-US before AWSoM

- Groth et al. 2000 [V via Gombosi 2000]: gamma = 5/3 from a hot base
  (2.85 MK, 1.5e8 cm^-3) stays hydrostatic until the base T > ~10-12 MK; fix:
  heating ~ (T0 - T), T0 = 1.75 T in the coronal hole (> 72.5 deg), T0 = T
  outside, scale height 4.5 R (equator) to 9 R (poles). Explicit form
  Q = rho q0 (T0 - gamma p/rho) exp[-(r - R)^2/sigma0^2] [R] (Manchester 2004).
  The slow wind from numerical reconnection at the cusp.
- Roussev 2003 [S], Cohen 2007 [S]: variable polytropic index (WSA through
  the Bernoulli integral, gamma -> ~1.5 by ~12.5 R [R]).
- Cohen 2017 (arXiv:1701.01146) [V]: polytropic gamma 1.05/1.1 never above
  600 km/s; 1.01 fast but too much mass loss (2-3e-14 Msun/yr); AWSoM > 600
  km/s with 1.3-1.6e-14 Msun/yr.

## 3. Classical theory

- Leer & Holzer 1980 [S]: energy added below the sonic point raises the
  scale height and the mass flux (little speed); above it, the asymptotic
  speed. Momentum addition likewise [R].
- Hansteen & Leer 1995 [S]: chromosphere-TR-corona-wind with conduction:
  most energy deposited on open field goes into the wind; the mass flux is
  fairly insensitive to where heat goes; the exponential T dependence of
  thermally driven mass fluxes is an artifact of a fixed TR pressure.
- Hammer 1982 [R]: the TR balance (conduction down = radiation + enthalpy)
  sets the base pressure and the mass flux. Withbroe 1988 [S]: ~10% of the
  energy into enthalpy; coronal holes need ~8e5 erg cm^-2 s^-1 [R].
- Parker/Chapman [R]: conduction-only winds at 1-2 MK are slow; a fast wind
  needs a non-thermal flux (Cranmer 2010, LRSP 2010-4) [S].

## 4. Wave pressure (WKB) and w+- models

- WKB [V]: dw/dt + div[(u + V_A) w] + (w/2) div u = 0; force -grad P_A,
  P_A = w/2, work -u.grad P_A in the energy equation; damping adds Gamma w to
  the heating. Steady 1D: w (u + V_A)^2 A/V_A = const (Jacques 1977) [R].
- AWSoM, van der Holst 2014 [V]: P_A = (w+ + w-)/2; dw+-/dt + div[(u +- V_A)
  w+-] + (w+-/2) div u = -+ R sqrt(w- w+) - Gamma+- w+-, Gamma+- = (2/L_perp)
  sqrt(w-+/rho); S_A/B = 1.1e6 W m^-2 T^-1, L_perp sqrt(B) = 1.5e5 m T^1/2,
  N 2e17 m^-3 (deliberately ~10x high against evaporation), T 5e4 K,
  dr_min 1e-3 R; top 24 R; Spitzer blended to the Hollweg flux
  q = (3/2) alpha p_e u (alpha 1.05) with f_S = 1/(1 + (r/5 R)^2). App. A: a
  local-dissipation simplification (Gamma = the reflection rate).
- Suzuki 2006 [V]: energetics estimate of v at 1 AU from the wave flux.
- Usmanov 2000 [S]: only the solutions with waves match the Ulysses fast wind.

## 5. Practicalities and pitfalls

- Outer boundary: 20 R (MAS, characteristic) or 24 R (SWMF, zero gradient);
  dr ~0.4 R out there is used in production.
- Inner boundary: upper chromosphere (2e4-5e4 K, 2e12-6e12 cm^-3) with TR
  broadening, or a hot base (polytropic).
- Pitfalls: base density too low -> evaporation; Spitzer conduction to 20 R
  is unphysical and stiff (switch it off by 5-10 R, collisionless flux
  beyond); heating too strong or too low -> too much mass flux, too little
  speed; polytropic models get the mass flux wrong; stratified heating on
  closed loops -> TNE; relaxation tens of hours (3D) to hundreds (1D).

## Recommendation (ranked)

1. WKB w+ transport (no damping at first, the outward wave only) plus the
   extended open-field heating (H0 ~5e-7 erg cm^-3 s^-1, lambda 0.7 R) on top
   of the existing heating; -grad P_A in momentum, -u.grad P_A in energy;
   Spitzer blended to a collisionless (Hollweg) flux by 5-10 R; outflow at
   20 R. Without a field yet: V_A from a reference open field B0 (R/r)^2 in
   the wave equation only; the injected flux from S_A/B = 1.1e6 W m^-2 T^-1.
   AWSoM later = switching on Gamma+- and reflection.
2. AWSoM App. A local dissipation (L_perp sqrt(B)) as the minimal w+-
   heating, replacing the empirical open-field term.
3. Not: variable-gamma/WSA or target-temperature heating (they bypass our
   thermodynamics and mis-set the mass flux).

## Minimal 1D tests

a. Isothermal Parker wind (1.5 MK) on the coarse grid to 20 R: transonic,
   matches the analytic solution, no reflection at the top.
b. Thermal wind, chromosphere to 20 R (Kopp-Holzer tube f_max 12.5,
   r_f 1.31 R, sigma 0.51 R), conduction + losses + heating, no waves: slow
   wind, mass flux ~2-5e8 cm^-2 s^-1 scaled to 1 AU.
c. Undamped WKB w+: wave action conserved; speed at 20 R up, mass flux
   nearly unchanged.
d. Damped: reproduce Lionello 2014 (z 25 km/s, lambda 0.05 R).
e. Later: Downs 2016's closed loop (TNE).

Main links: arXiv 1801.05903 (Torok 2018), 1311.4093 (AWSoM), 1402.4188
(Lionello 2014), 0912.2647 (Downs 2010), 1701.01146 (Cohen 2017),
astro-ph/0602062 (Suzuki 2006); doi 10.3847/0004-637X/832/2/180 (Downs 2016),
10.1029/2000JA900093 (Groth 2000), 10.1088/0004-637X/690/1/902 (Lionello
2009), 10.1029/JA085iA09p04681 (Leer & Holzer 1980), 10.1029/95JA02300
(Hansteen & Leer 1995), 10.1029/1999JA000233 (Usmanov 2000); Cranmer 2010,
LRSP 2010-4.
