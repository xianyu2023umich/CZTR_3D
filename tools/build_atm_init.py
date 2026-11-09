#!/usr/bin/env python3
"""CZTR_3D init file: the CZ background plus an atmosphere to the top of the
extended grid (tools/grid/extend_grid.py), in the format cztr3d.x reads
("n k", then 'i r z rho T v rho0 T0 g' for cells -1..n+2).

Cells <= k (the CZ, k = the join cell of the CZ init file): state and
background exactly those of the CZ init file (static_init_015a1.dat), so the
CZ perturbations of CZ_3D restarts keep their meaning.

Initial atmosphere (cells > k):
  - z <= z_cut: CZTR's 1D step-3 state (chromosphere, TR, low corona),
    ln rho and ln T interpolated in z, rho scaled by rho_CZ/rho_CZTR at the
    CZ top (a uniform pressure scale keeps it hydrostatic);
  - z > z_cut: T rises smoothly (C1) from T(z_cut) to T_corona,
    T = T_c - (T_c - T_cut) exp(-(z - z_cut)/w_T), and rho follows from
    hydrostatic equilibrium dp/dr = -rho g with p = rho k T/(mu m_H)
    (fully ionized, mu = mu_cor), integrated from z_cut.
Background above the join (CZTR's cztr_init_mod rules: smooth (C1), not
following the TR, not above the state):
  rho0 = (rho_k - rho_fl) exp(-int dr/H) + rho_fl,
  H(r) = H_k + (H_inf - H_k)(1 - exp(-((r - r_k)/w)^2)),
  T0   = T_k + T'_k w (1 - exp(-(r - r_k)/w)),
  H_k and T'_k from least-squares slopes of the top n_fit CZ cells of the
  background (C1 at the join), w = blend_km, H_inf the largest value with
  rho0 <= rho in every cell more than w above the join (bisection), rho_fl =
  min(rho_floor, 0.99 min rho).
g = G M_sun / r^2 above the join (the CZ file's g below).

Optional (--parker): an isothermal Parker wind at T_corona (sound speed
c = sqrt(k T/(mu m_H)), sonic point r_s = G M/(2 c^2), the transonic
solution u^2 - ln u^2 = 4 ln(r/r_s) + 4 r_s/r - 3) joined to the static
atmosphere: ln rho blends from the hydrostatic profile to the Parker density
(normalized to the hydrostatic one at r1) between r1 and r2 (a cos^2 step);
the mass flux F = rho u c r^2 of the Parker solution then gives
v = F/(rho r^2) everywhere above the TR (switched on by a cos^2 step between
v_z0 and v_z1), so the start conserves mass. The background gets a floor
rho_fl = floor_frac x min(rho) (>= --rho0-min) so that it stays within the
EOS table far out and below the state.

Optional (--state-profile file): the atmosphere state (cells above the
join) taken from a profile file's <rho>, <T>, <v_r> (a relaxed run on the
same grid, e.g. the steady 1D wind), the ghosts above the top extrapolated
(rho ~ r^-2, T and v copied); the background is then built from that state
(with --parker's coronal background part).

Usage: build_atm_init.py grid.dat cz_init.dat cztr_step3_init.dat out.dat
       [--z-cut-km 3000] [--T-corona 1.5e6] [--wT-km 10000] [--mu 0.593]
       [--blend-km 300] [--n-fit 4] [--Tb-frac 0.8] [--fB 0.1]
       [--parker] [--parker-r1 1.2] [--parker-r2 2.0] [--v-z0-km 2000] [--v-z1-km 5000]
       [--floor-frac 0.3] [--rho0-min 2e-21]
"""
import argparse
import numpy as np

G, MSUN, RSUN = 6.67430e-8, 1.98892e33, 6.957e10
KB, MH = 1.380649e-16, 1.6735575e-24


def read_grid(path):
    rows = [l.split() for l in open(path) if not l.startswith('#')]
    n = int(rows[0][0])
    rc = np.array([float(r[1]) for r in rows[1:n + 5]])
    rf = np.array([float(r[1]) for r in rows[n + 5:2 * n + 6]])
    return n, rc, rf


def read_init(path):
    rows = [l.split() for l in open(path) if not l.startswith('#')]
    n, k = int(rows[0][0]), int(rows[0][1])
    d = np.array(rows[1:n + 5], dtype=float)   # cells -1..n+2
    return n, k, d


def ls_slope(x, y):
    x = np.asarray(x); y = np.asarray(y)
    xm, ym = x.mean(), y.mean()
    return np.sum((x - xm) * (y - ym)) / np.sum((x - xm) ** 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('grid'); ap.add_argument('cz_init'); ap.add_argument('cztr_init'); ap.add_argument('out')
    ap.add_argument('--z-cut-km', type=float, default=3000.0)
    ap.add_argument('--T-corona', type=float, default=1.5e6)
    ap.add_argument('--wT-km', type=float, default=10000.0)
    ap.add_argument('--mu', type=float, default=0.593)
    ap.add_argument('--blend-km', type=float, default=300.0)
    ap.add_argument('--n-fit', type=int, default=4)
    ap.add_argument('--Tb-frac', type=float, default=0.8, help='coronal background part: T_b = Tb_frac T_corona')
    ap.add_argument('--fB', type=float, default=0.1, help='coronal background part: fraction of the state at z_cut')
    ap.add_argument('--parker', action='store_true', help='join an isothermal Parker wind at T_corona')
    ap.add_argument('--parker-r1', type=float, default=1.2, help='blend start [R_sun]')
    ap.add_argument('--parker-r2', type=float, default=2.0, help='blend end [R_sun]')
    ap.add_argument('--v-z0-km', type=float, default=2000.0, help='the wind velocity switched on from here ...')
    ap.add_argument('--v-z1-km', type=float, default=5000.0, help='... to here (cos^2 step)')
    ap.add_argument('--floor-frac', type=float, default=0.0, help='background floor: this x min(rho) above the join')
    ap.add_argument('--rho0-min', type=float, default=0.0, help='background floor at least this (EOS table)')
    ap.add_argument('--state-profile', default='', help='the atmosphere state from this profile file')
    a = ap.parse_args()

    n, rc, rf = read_grid(a.grid)
    ncz, k, dcz = read_init(a.cz_init)
    _, k3, d3 = read_init(a.cztr_init)
    # the CZ init file's cells -1..k+2 must be the grid's
    if not np.allclose(dcz[:k + 4, 1], rc[:k + 4], rtol=1e-13, atol=0):
        raise SystemExit('the CZ init file and the grid differ in cells -1..k+2')
    z = (rc - RSUN) / 1e5                                 # km, cells -1..n+2 (index i+1)
    ii = np.arange(-1, n + 3)

    rho = np.empty(n + 4); T = np.empty(n + 4); g = np.empty(n + 4)
    rho0 = np.empty(n + 4); T0 = np.empty(n + 4)
    cz = ii <= k
    dc = dcz[:k + 2]                                      # CZ file rows of cells -1..k
    rho[cz], T[cz], rho0[cz], T0[cz], g[cz] = dc[:, 3], dc[:, 4], dc[:, 6], dc[:, 7], dc[:, 8]
    g[~cz] = G * MSUN / rc[~cz] ** 2

    # --- initial atmosphere: CZTR step 3 up to z_cut
    z3, rho3, T3 = d3[:, 2], d3[:, 3], d3[:, 4]
    j3 = k3 + 1                                           # row of CZTR's join cell
    zj = z[k + 1]
    scale = rho[k + 1] / np.exp(np.interp(zj, z3, np.log(rho3)))
    low = (~cz) & (z <= a.z_cut_km)
    rho[low] = scale * np.exp(np.interp(z[low], z3, np.log(rho3)))
    T[low] = np.exp(np.interp(z[low], z3, np.log(T3)))
    # --- above z_cut: smooth T rise and hydrostatic rho (fine sub-steps)
    icut = np.where(low)[0][-1]
    rcut, Tcut = rc[icut], T[icut]
    pcut = rho[icut] * KB * Tcut / (a.mu * MH)
    wT = a.wT_km * 1e5
    Tof = lambda rr: a.T_corona - (a.T_corona - Tcut) * np.exp(-(rr - rcut) / wT)
    # integrate d ln p / dr = -mu m_H g / (k T) on a fine grid up to the last ghost
    rr = np.linspace(rcut, rc[-1], 200001)
    integrand = a.mu * MH * G * MSUN / (rr ** 2 * KB * Tof(rr))
    lnp = np.log(pcut) - np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(rr))])
    up = np.arange(n + 4) > icut
    T[up] = Tof(rc[up])
    rho[up] = np.exp(np.interp(rc[up], rr, lnp)) * a.mu * MH / (KB * T[up])
    v = np.zeros(n + 4)
    parker_info = ''
    if a.parker:
        c = np.sqrt(KB * a.T_corona / (a.mu * MH))
        rs = G * MSUN / (2.0 * c * c)
        def u_of(x):                                     # transonic branch, x = r/r_s
            rhs = 4.0 * np.log(x) + 4.0 / x - 3.0
            lo, hi = (-30.0, 0.0) if x < 1.0 else (0.0, 10.0)   # ln u
            if abs(x - 1.0) < 1e-9:
                return 1.0
            f = lambda lu: np.exp(2 * lu) - 2 * lu - rhs
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                # f is decreasing in lu for u < 1 and increasing for u > 1
                if (f(mid) > 0) == (x < 1.0):
                    lo = mid
                else:
                    hi = mid
            return np.exp(0.5 * (lo + hi))
        r1, r2 = a.parker_r1 * RSUN, a.parker_r2 * RSUN
        uu = np.array([u_of(x) for x in rc / rs])
        i1 = np.argmin(abs(rc - r1))
        C = rho[i1] * rc[i1] ** 2 * uu[i1]               # rho_P r^2 u = C (normalized at r1)
        rhoP = C / (rc ** 2 * uu)
        wb = np.clip((rc - r1) / (r2 - r1), 0.0, 1.0)
        wb = np.sin(0.5 * np.pi * wb) ** 2
        blend = up & (rc > r1)
        rho[blend] = np.exp((1.0 - wb[blend]) * np.log(rho[blend]) + wb[blend] * np.log(rhoP[blend]))
        F = C * c                                         # rho v r^2
        wz = np.clip((z - a.v_z0_km) / (a.v_z1_km - a.v_z0_km), 0.0, 1.0)
        wz = np.sin(0.5 * np.pi * wz) ** 2
        atm_cells = np.arange(n + 4) > k + 1
        v[atm_cells] = wz[atm_cells] * F / (rho[atm_cells] * rc[atm_cells] ** 2)
        parker_info = ('Parker: c %.1f km/s, sonic point %.2f R_sun, mass flux at 1 AU %.3e g/cm2/s '
                       '(n v %.2e cm-2 s-1), blend %.2f-%.2f R_sun'
                       % (c / 1e5, rs / RSUN, F / 1.496e13 ** 2, F / 1.496e13 ** 2 / (a.mu * MH), a.parker_r1, a.parker_r2))

    # --- the atmosphere state from a relaxed profile (same grid)
    if a.state_profile:
        Pf = np.loadtxt(a.state_profile)
        if Pf.shape[0] != n or not np.allclose(Pf[:, 1], rc[2:n + 2], rtol=1e-10, atol=0):
            raise SystemExit('the state profile is not on this grid')
        jj = np.arange(n + 4)
        sel = (jj > k + 1) & (jj <= n + 1)               # array rows of cells k+1..n
        rho[sel] = Pf[jj[sel] - 2, 3];  T[sel] = Pf[jj[sel] - 2, 4];  v[sel] = Pf[jj[sel] - 2, 7]
        for jg in (n + 2, n + 3):                         # ghosts n+1, n+2
            rho[jg] = rho[n + 1]*(rc[n + 1]/rc[jg])**2;  T[jg] = T[n + 1];  v[jg] = v[n + 1]
        parker_info = (parker_info + '; ' if parker_info else '') + 'atmosphere state from ' + a.state_profile

    # --- background above the join (CZTR's construction)
    blend = a.blend_km * 1e5
    ik = k + 1                                            # array index of cell k
    sl = slice(ik - a.n_fit + 1, ik + 1)
    dlnrho = ls_slope(rc[sl], np.log(rho0[sl]))
    dT = ls_slope(rc[sl], T0[sl])
    rho_k, T_k, r_k = rho0[ik], T0[ik], rc[ik]
    atm = np.arange(n + 4) > ik
    # coronal part rho_B: hydrostatic at T_b (falls faster than the state when
    # T_b < T_corona), fB x the state at z_cut
    Tb = a.Tb_frac * a.T_corona
    kb = a.mu * MH * G * MSUN / (KB * Tb)                # rho_B = c exp(kb (1/r - 1/r_cut))
    cB = a.fB * rho[icut]
    if a.parker:
        # with the wind: fB x the (smooth) coronal state above r_cut, below it
        # hydrostatic at the cut temperature (the same slope at r_cut)
        kc = a.mu * MH * G * MSUN / (KB * Tcut)
        lnrc = np.log(rho[icut:])
        def rhoB(x):
            x = np.asarray(x, dtype=float)
            out = np.where(x >= rcut, a.fB * np.exp(np.interp(x, rc[icut:], lnrc)),
                           cB * np.exp(kc * (1.0 / x - 1.0 / rcut)))
            return out if out.ndim else float(out)
        dlnB = -kc / r_k ** 2
    else:
        rhoB = lambda x: cB * np.exp(kb * (1.0 / x - 1.0 / rcut))
        dlnB = -kb / r_k ** 2                            # d ln rho_B / dr at r_k
    # chromospheric part rho_A: C1 with the CZ background at r_k together with rho_B
    rho_fl = max(a.floor_frac * rho[atm].min(), a.rho0_min)
    if rho_fl > 0.0 and rho_fl >= rho[atm].min():
        raise SystemExit('the background floor %.3e is not below the state minimum %.3e' % (rho_fl, rho[atm].min()))
    rhoA_k = rho_k - rhoB(r_k) - rho_fl
    H_k = -rhoA_k / (rho_k * dlnrho - rhoB(r_k) * dlnB)

    def fill(H_inf):
        out = np.empty(atm.sum())
        for m, j in enumerate(np.where(atm)[0]):
            x = np.linspace(r_k, rc[j], 401)
            H = H_k + (H_inf - H_k) * (1.0 - np.exp(-((x - r_k) / blend) ** 2))
            f = 1.0 / H
            out[m] = rhoA_k * np.exp(-np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(x))) + rhoB(rc[j]) + rho_fl
        return out

    beyond = (rc[atm] - r_k) > blend
    lo, hi = 1.0e5, 1.0e12                               # cm: 1 km .. 1e7 km (H_inf may lie below H_k)
    for _ in range(80):
        mid = np.sqrt(lo * hi)
        if np.all(fill(mid)[beyond] <= rho[atm][beyond]):
            lo = mid
        else:
            hi = mid
    H_inf = lo
    if not np.all(fill(H_inf)[beyond] <= rho[atm][beyond]):
        raise SystemExit('no asymptotic scale height keeps rho0 <= rho beyond the blend zone')
    rho0[atm] = fill(H_inf)
    T0[atm] = T_k + dT * blend * (1.0 - np.exp(-(rc[atm] - r_k) / blend))

    with open(a.out, 'w') as f:
        f.write('# CZTR_3D init: tools/build_atm_init.py %s %s %s\n' % (a.grid, a.cz_init, a.cztr_init))
        f.write('# CZ (cells <= %d) from the CZ init file; atmosphere: CZTR step 3 to z = %.0f km (rho x %.4f), '
                'then T -> %.3e K (w %.0f km), hydrostatic, mu %.3f\n' % (k, a.z_cut_km, scale, a.T_corona, a.wT_km, a.mu))
        f.write('# background above the join: rho_A (H_k %.1f km -> H_inf %.1f km, blend %.0f km) + rho_B '
                '(hydrostatic at %.3e K, %.2f x the state at z_cut) + floor %.3e; dT0/dr %.4e K/cm\n'
                % (H_k / 1e5, H_inf / 1e5, a.blend_km, Tb, a.fB, rho_fl, dT))
        if parker_info:
            f.write('# ' + parker_info + '\n')
        f.write('# n, k\n# i  r[cm]  z[km]  rho  T  v[cm/s]  rho0  T0  g[cm/s2]\n')
        f.write('%d %d\n' % (n, k))
        for j in range(n + 4):
            f.write('%6d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n'
                    % (ii[j], rc[j], z[j], rho[j], T[j], v[j], rho0[j], T0[j], g[j]))

    # --- report
    print('n %d, k %d, CZTR rho scale %.4f, H_k %.1f km, H_inf %.1f km, rho_B: T_b %.3e, fB %.2f, floor %.3e'
          % (n, k, scale, H_k / 1e5, H_inf / 1e5, Tb, a.fB, rho_fl))
    if parker_info:
        print(parker_info)
    r_in = slice(1, n + 3)
    ratio = rho[r_in] / rho0[r_in]
    print('min rho/rho0: everywhere %.4f; beyond the blend zone %.4f' % (ratio.min(), ratio[(rc[r_in] - r_k) > blend].min() if np.any((rc[r_in] - r_k) > blend) else np.nan))
    # smoothness on the grid: second difference of ln rho0 and of T0 per cell (in cell index)
    for name, y in [('ln rho0', np.log(rho0[r_in])), ('T0/T0_k', T0[r_in] / T_k)]:
        d1 = np.diff(y); d2 = np.diff(y, 2)
        ja = np.argmax(np.abs(d2))
        print('%s: largest |first difference| per cell %.3f, largest |second difference| %.4f (z %.0f km); '
              'at the join %.2e' % (name, np.abs(d1).max(), abs(d2[ja]), z[r_in][1:-1][ja], abs(d2[k - 1])))
    for zz in [-500, 0, 500, 1000, 2000, 3000, 1e4, 5e4, 2e5, 7e5, 1.04e6, 2.1e6, 4.2e6, 7e6, 1.3e7]:
        j = np.argmin(abs(z - zz))
        if j >= n + 2:
            continue
        print('  z %9.0f km (%5.2f R)  rho %.3e  T %.3e  v %8.2f km/s  rho0 %.3e  T0 %.0f  rho/rho0 %.3g'
              % (z[j], rc[j] / RSUN, rho[j], T[j], v[j] / 1e5, rho0[j], T0[j], rho[j] / rho0[j]))


if __name__ == '__main__':
    main()
