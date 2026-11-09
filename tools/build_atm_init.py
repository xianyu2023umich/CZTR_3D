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

Usage: build_atm_init.py grid.dat cz_init.dat cztr_step3_init.dat out.dat
       [--z-cut-km 3000] [--T-corona 1.5e6] [--wT-km 10000] [--mu 0.593]
       [--blend-km 300] [--n-fit 4] [--rho-floor 1e-17]
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
    rhoB = lambda x: cB * np.exp(kb * (1.0 / x - 1.0 / rcut))
    dlnB = -kb / r_k ** 2                                # d ln rho_B / dr at r_k
    # chromospheric part rho_A: C1 with the CZ background at r_k together with rho_B
    rhoA_k = rho_k - rhoB(r_k)
    H_k = -rhoA_k / (rho_k * dlnrho - rhoB(r_k) * dlnB)

    def fill(H_inf):
        out = np.empty(atm.sum())
        for m, j in enumerate(np.where(atm)[0]):
            x = np.linspace(r_k, rc[j], 401)
            H = H_k + (H_inf - H_k) * (1.0 - np.exp(-((x - r_k) / blend) ** 2))
            f = 1.0 / H
            out[m] = rhoA_k * np.exp(-np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(x))) + rhoB(rc[j])
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
    rho0[atm] = fill(H_inf)
    T0[atm] = T_k + dT * blend * (1.0 - np.exp(-(rc[atm] - r_k) / blend))

    v = np.zeros(n + 4)
    with open(a.out, 'w') as f:
        f.write('# CZTR_3D init: tools/build_atm_init.py %s %s %s\n' % (a.grid, a.cz_init, a.cztr_init))
        f.write('# CZ (cells <= %d) from the CZ init file; atmosphere: CZTR step 3 to z = %.0f km (rho x %.4f), '
                'then T -> %.3e K (w %.0f km), hydrostatic, mu %.3f\n' % (k, a.z_cut_km, scale, a.T_corona, a.wT_km, a.mu))
        f.write('# background above the join: rho_A (H_k %.1f km -> H_inf %.1f km, blend %.0f km) + rho_B '
                '(hydrostatic at %.3e K, %.2f x the state at z_cut); dT0/dr %.4e K/cm\n'
                % (H_k / 1e5, H_inf / 1e5, a.blend_km, Tb, a.fB, dT))
        f.write('# n, k\n# i  r[cm]  z[km]  rho  T  v[cm/s]  rho0  T0  g[cm/s2]\n')
        f.write('%d %d\n' % (n, k))
        for j in range(n + 4):
            f.write('%6d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n'
                    % (ii[j], rc[j], z[j], rho[j], T[j], v[j], rho0[j], T0[j], g[j]))

    # --- report
    print('n %d, k %d, CZTR rho scale %.4f, H_k %.1f km, H_inf %.1f km, rho_B: T_b %.3e, fB %.2f'
          % (n, k, scale, H_k / 1e5, H_inf / 1e5, Tb, a.fB))
    r_in = slice(1, n + 3)
    ratio = rho[r_in] / rho0[r_in]
    print('min rho/rho0: everywhere %.4f; beyond the blend zone %.4f' % (ratio.min(), ratio[(rc[r_in] - r_k) > blend].min() if np.any((rc[r_in] - r_k) > blend) else np.nan))
    # smoothness on the grid: second difference of ln rho0 and of T0 per cell (in cell index)
    for name, y in [('ln rho0', np.log(rho0[r_in])), ('T0/T0_k', T0[r_in] / T_k)]:
        d1 = np.diff(y); d2 = np.diff(y, 2)
        ja = np.argmax(np.abs(d2))
        print('%s: largest |first difference| per cell %.3f, largest |second difference| %.4f (z %.0f km); '
              'at the join %.2e' % (name, np.abs(d1).max(), abs(d2[ja]), z[r_in][1:-1][ja], abs(d2[k - 1])))
    for zz in [-500, 0, 500, 1000, 2000, 3000, 1e4, 5e4, 2e5, 7e5, 1.04e6]:
        j = np.argmin(abs(z - zz))
        print('  z %9.0f km  rho %.3e  T %.3e  rho0 %.3e  T0 %.0f  rho/rho0 %.3g' % (z[j], rho[j], T[j], rho0[j], T0[j], rho[j] / rho0[j]))


if __name__ == '__main__':
    main()
