# Coronal heating: options for CZTR_3D

Literature survey of 2026-11-07 (done by a Claude subagent with web search;
my summary). Tags: [V] read in the source (full text or abstract), [S] read
in a later paper quoting it, [R] recalled, not verified. Units cgs unless noted.

Constraint: at 0.47 deg (5.7 Mm cells at the surface) granulation is not
resolved, so the resolved Poynting flux from braiding is far too small; the
heating must be a parametrization fed by quantities the run has (|B|, rho,
T, v, the MLT velocity, the surface field).

## 1. Empirical laws of global thermodynamic models

**Downs et al. 2010** (ApJ 712, 1219; arXiv:0912.2647) [V], the SWMF lower
corona (after MAS/Lionello):
- Exponential: Q = H0 exp(-(r - R)/lambda), H0 = 7.28e-5 erg cm^-3 s^-1,
  lambda = 40 Mm (a surface flux of 3.3e5 erg cm^-2 s^-1).
- |B|-weighted (Abbett 2007, Bercik 2005 form): effectively Q = H_B |B| with
  H_B ~ 4e-5 erg cm^-3 s^-1 G^-1, times the 40 Mm exponential.
- Wind/coronal-hole term: the same form, H0 = 5e-7, lambda = 0.7 R.
- The short-scale-height term only on closed field (from PFSS).
- Best match to EIT/SXT ("Run D"): exponential, blended to 4e-5 |B| above
  30 G, plus the open-field cutoff. Uniform exponential heating already gives
  the mean low corona; active regions need the |B| term.
- Also used: kappa0 = 1.23e-6, transition-region broadening below 3e5 K.

**Lionello, Linker & Mikic 2009** (ApJ 690, 902) [V abstract; formula S via
COCONUT, Baratashvili et al. 2024, arXiv:2407.17903]:
- Q = Q_exp + Q_QS + Q_AR
- Q_exp = H0 exp(-(r - R)/0.7 R)
- Q_QS = 1.18e-5 f(r) B_t^2/[B (|B_r| + 0.55 G)]
- Q_AR = 1.87e-5 g(B) (B/1 G)^1.2
- f(r) = 1/2 [1 + tanh((1.7 - r/R)/0.1)] exp(-(r/R - 1)/0.2),
  g(B) = 1/2 [1 + tanh((B - 18.1 G)/3.97 G)]
- The thresholds are case-tuned (COCONUT moved g's 18.1 G to 4 G); expect to
  retune at our resolution.

**Schrijver et al. 2004** (ApJ 615, 512) [S]: loop ensemble, flux F_H ~ B L^-1
(volumetric ~ B L^-2). Needs each cell's field-line length: tracing on
Yin-Yang every step is awkward. Not for now.

**Mok et al. 2005** (ApJ 621, 1098) [V abstract]: active region, Q ~ |B|,
synthetic EUV/SXR. Footpoint-concentrated steady heating gives thermal
non-equilibrium cycles (Mok 2016, Mikic 2013) [S].

## 2. Braiding / Poynting-flux heating (boxes)

Gudiksen & Nordlund 2005 [V citation]; Bingert & Peter 2011 (A&A 530, A112,
Q = eta j^2, eta = 1e10 m^2/s) [V abstract]; van Wettum, Bingert & Peter 2013
(A&A 554, A39) [V]: compared with van Ballegooijen 2011 waves
Q ~ B^0.55 L^-0.92 v_rms^1.65 and Rappazzo 2008 turbulence
Q ~ B^1.75 L^-1.75 rho^0.125 v_rms^1.25 l^0.75 (v_rms = footpoint velocity:
our MLT velocity could supply it). Warnecke & Peter 2019 (A&A 624, L12) [V
abstract]. Rempel 2017 (ApJ 834, 10) [V]: no explicit heating, numerical
dissipation of the self-consistent Poynting flux; coronal losses 5e5 (quiet),
2e6 (arcade), 4e6 (AR) erg cm^-2 s^-1. Chen, Rempel & Fan 2022 (ApJ,
arXiv:2106.14055) [V]: quiet-Sun S_z = 9e5 erg cm^-2 s^-1, mostly viscous.
All need granular-scale driving (<~400 km): eta j^2 on our grid would badly
underheat. Fine as an extra term (large-scale stresses, emergence), not as
the main heating.

## 3. Alfven-wave turbulence (AWSoM, MAS WTD)

**AWSoM, van der Holst et al. 2014** (ApJ 782, 81; arXiv:1311.4093) [V];
Sokolov et al. 2013 (ApJ 778, 176):
- d w+-/dt + div[(u +- V_A) w+-] + (w+-/2) div u = -+ R sqrt(w- w+) - Gamma+- w+-
- Gamma+- = (2/L_perp) sqrt(w-+/rho); Q = Gamma+ w+ + Gamma- w-
- L_perp sqrt(B) = 1.5e5 m T^1/2 (~1.5e9 cm G^1/2: L_perp ~ 1.5 Mm at 100 G)
- reflection R = min[R_imb, max(Gamma+-)] x imbalance factor,
  R_imb = sqrt[((V_A.grad) log V_A)^2 + (b.curl u)^2]
- inner boundary: outgoing S_A/B = 1.1e6 W m^-2 T^-1 (~1.1e5 erg cm^-2 s^-1
  G^-1); best fits over the years 0.9e6 (2010) ... 1.5e6 (2017) [S].
- ARs brighten in synthetic EUV (strong B: short L_perp, more dissipation).

**MAS wave-turbulence-driven (WTD), Downs et al. 2016** (ApJ 832, 180;
arXiv:1610.02113) [V], after Lionello et al. 2014 [R]:
- Q_w = rho (|z-| z+^2 + |z+| z-^2)/(4 lambda_perp),
  lambda_perp = lambda0 sqrt(6.09 G/B), lambda0 = 0.01-0.07 R
- waves injected at 5-15 km/s, WKB (no reflection/dissipation) until
  n_e = 2e11 cm^-3; an extra exponential (1.1 Mm) heating balances the
  chromospheric losses. Used for the 2017 eclipse prediction (Mikic et al.
  2018, Nat. Astron. 2, 913).

## 4. CZ-to-corona models

Abbett 2007, Abbett & Fisher 2012 (RADMHD, boxes): the |B|-weighted heating
(Downs model 1). Warnecke et al. 2010-2012: spherical wedges to ~1.5 R with a
polytropic "corona", no heating physics. No global spherical radiative
CZ-to-corona model with a physical heating law was found: ours appears new.

## Recommendation (for the first demo run)

1. **Q = (H_B |B| + H0) exp(-h/lambda) + a weak wind term** (Downs 2010
   Run D / Mok 2005): H_B ~ 4e-5 erg cm^-3 s^-1 G^-1, H0 ~ 7e-5, lambda ~
   40 Mm; wind 5e-7 with 0.7 R. Local, cheap, ARs brighten by themselves.
   For the coupled model: h from the local photosphere, and gate it on
   density or temperature (e.g. only where n_e < 2e11 cm^-3, as Downs 2016),
   or the kG field in the CZ/photosphere gets heated. Footpoint-heavy heating
   can drive thermal non-equilibrium: TR broadening tames it.
2. If ARs come out too faint: the Lionello 2009 composite (B^1.2 AR term,
   retuned threshold).
3. Stepping stone to AWSoM: normalize H_B so the flux per tube matches
   S_A/B = 1.1e6 W m^-2 T^-1 (Downs' H|B| over 40 Mm gives ~1.6e5 B erg
   cm^-2 s^-1 vs AWSoM's 1.1e5 B); then evolve w+- with
   Gamma+- = (2/L_perp) sqrt(w-+/rho), L_perp sqrt(B) = 1.5e5 m T^1/2,
   injected at a density threshold; finally replace the fixed S_A/B by
   S_A = rho <dv^2> V_A with dv from the MLT velocity (the user's idea,
   CZTR_2D/notes/alfven_wave_heating.md).

Full-text extracts of the papers (this session's scratchpad, temporary):
/private/tmp/claude-502/-Users-xianyu-Documents-FALCON-EOS/d51ba730-a454-464d-b78e-af444a7d5cac/scratchpad/txt/
