#!/usr/bin/env python3
"""Extend a CZTR grid file into the corona (CZTR_3D).

Keeps the first k0 cells of a base grid (the CZ and the chromosphere, so the
CZ cells and restarts are unchanged), then continues the map r(xi) with a
spacing that grows geometrically (rate ramping smoothly from the base grid's
own growth at k0 to s_max per cell) until it reaches aspect * r * dtheta,
after which dr = aspect * r * dtheta (smooth minimum of the two, power p).
The number of cells is chosen so the total is a multiple of nb_r, and aspect
is adjusted (bisection) so the last face lands exactly on r_top. Faces at
integer xi, centres at half-integer xi (as CZTR's grids), two ghost centres
beyond each end.

Optionally (--outer-aspect A > 1) the aspect grows smoothly beyond
--outer-start [R_sun] to A over --outer-width [R_sun] (a cos^2 ramp), for a
coarse radial grid far out (the solar wind region, e.g. to 20 R_sun).

Usage:
  extend_grid.py base.dat out.dat [--z0-km 2000] [--r-top 2.5] [--deg 0.46875]
                 [--smax 0.05] [--ramp 10] [--nb-r 28] [--n-total N]
                 [--outer-aspect 1] [--outer-start 2.5] [--outer-width 2.0]
"""
import argparse
import numpy as np

RSUN = 6.957e10


def read_grid(path):
    head, rows = [], []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                head.append(line.rstrip('\n'))
            else:
                rows.append(line.split())
    n = int(rows[0][0])
    rc = np.array([float(r[1]) for r in rows[1:n + 5]])          # cells -1..n+2
    rf = np.array([float(r[1]) for r in rows[n + 5:2 * n + 6]])  # faces 0..n
    return head, n, rc, rf


def aspect_factor(r, A, r1, w):
    """1 below r1, A beyond r1 + w, a cos^2 ramp between (r in cm, r1, w in cm)."""
    if A <= 1.0:
        return 1.0
    x = min(max((r - r1) / w, 0.0), 1.0)
    return 1.0 + (A - 1.0) * np.sin(0.5 * np.pi * x) ** 2


def march(r0, d0, s0, smax, ramp, a, p, n_cells, sub=40, outer=(1.0, 0.0, 1.0)):
    """r(xi) for xi = 0..n_cells (in steps of 1/sub) from dr/dxi =
    smin(d0 exp(int s), a r), s ramping from s0 to smax over `ramp` cells."""
    m = n_cells * sub + 1
    xi = np.linspace(0.0, n_cells, m)
    h = xi[1] - xi[0]
    r = np.empty(m)
    r[0] = r0
    lng = 0.0                                   # ln of the geometric spacing / d0

    def s_of(x):
        w = 0.5 * (1.0 - np.cos(np.pi * min(x / ramp, 1.0)))
        return s0 + (smax - s0) * w

    def rate(x, rr, lg):
        g = d0 * np.exp(lg)
        c = a * rr * aspect_factor(rr, *outer)
        return (g ** -p + c ** -p) ** (-1.0 / p)

    for k in range(m - 1):
        x = xi[k]
        # RK4 on (r, lng): d lng/dxi = s(xi)
        k1r = rate(x, r[k], lng);                         k1g = s_of(x)
        k2r = rate(x + h / 2, r[k] + h / 2 * k1r, lng + h / 2 * k1g); k2g = s_of(x + h / 2)
        k3r = rate(x + h / 2, r[k] + h / 2 * k2r, lng + h / 2 * k2g); k3g = k2g
        k4r = rate(x + h, r[k] + h * k3r, lng + h * k3g); k4g = s_of(x + h)
        r[k + 1] = r[k] + h / 6 * (k1r + 2 * k2r + 2 * k3r + k4r)
        lng += h / 6 * (k1g + 2 * k2g + 2 * k3g + k4g)
    return xi, r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('base'); ap.add_argument('out')
    ap.add_argument('--z0-km', type=float, default=2000.0, help='keep the base grid up to this height')
    ap.add_argument('--r-top', type=float, default=2.5, help='top radius [R_sun]')
    ap.add_argument('--deg', type=float, default=0.46875, help='horizontal spacing for the cap [deg]')
    ap.add_argument('--smax', type=float, default=0.05, help='largest geometric growth per cell')
    ap.add_argument('--ramp', type=float, default=10.0, help='cells over which the growth ramps up')
    ap.add_argument('--p', type=float, default=4.0, help='smooth-minimum power')
    ap.add_argument('--nb-r', type=int, default=28, help='total cells a multiple of this')
    ap.add_argument('--n-cz', type=int, default=196, help='CZ cells (reported only)')
    ap.add_argument('--n-total', type=int, default=0, help='total cells (0: the smallest multiple of nb_r '
                    'for which the cap regime is reached)')
    ap.add_argument('--outer-aspect', type=float, default=1.0, help='aspect factor far out (1: none)')
    ap.add_argument('--outer-start', type=float, default=2.5, help='where the outer ramp starts [R_sun]')
    ap.add_argument('--outer-width', type=float, default=2.0, help='width of the outer ramp [R_sun]')
    args = ap.parse_args()

    head, n0, rc0, rf0 = read_grid(args.base)
    z_f = (rf0 - RSUN) / 1e5
    k0 = int(np.searchsorted(z_f, args.z0_km))            # faces 0..k0 kept
    d0 = rf0[k0] - rf0[k0 - 1]
    s0 = np.log((rf0[k0] - rf0[k0 - 1]) / (rf0[k0 - 1] - rf0[k0 - 2]))
    r_top = args.r_top * RSUN
    dth = np.radians(args.deg)

    outer = (args.outer_aspect, args.outer_start * RSUN, args.outer_width * RSUN)

    def top_for(a, n_new):
        xi, r = march(rf0[k0], d0, s0, args.smax, args.ramp, a * dth, args.p, n_new, outer=outer)
        return xi, r

    # cells above k0 with aspect 1: how many to reach r_top
    xi, r = top_for(1.0, 4000)
    n_new1 = int(np.ceil(np.interp(r_top, r, xi)))
    n_tot = args.n_total or int(np.ceil((k0 + n_new1) / args.nb_r) * args.nb_r)
    n_new = n_tot - k0
    lo, hi = 0.2, 5.0                                       # aspect bisection
    for _ in range(60):
        a = 0.5 * (lo + hi)
        xi, r = top_for(a, n_new)
        if r[-1] > r_top:
            hi = a
        else:
            lo = a
    xi, r = top_for(a, n_new)
    sub = int(round(1.0 / (xi[1] - xi[0])))
    faces_new = r[::sub]                                     # xi = 0..n_new
    faces_new = faces_new * 0 + faces_new                    # copy
    faces_new[-1] = r_top
    centres_new = r[sub // 2::sub]                           # xi = 0.5, 1.5, ...
    # ghost centres beyond the top: continue the last spacing geometrically
    drl = faces_new[-1] - faces_new[-2]
    q = drl / (faces_new[-2] - faces_new[-3])
    g1 = faces_new[-1] + 0.5 * drl * q
    g2 = g1 + drl * q * q

    rf = np.concatenate([rf0[:k0], faces_new])               # faces 0..n_tot
    rc = np.concatenate([rc0[:k0 + 2], centres_new, [g1, g2]])  # cells -1..n_tot+2
    assert len(rf) == n_tot + 1 and len(rc) == n_tot + 4, (len(rf), len(rc), n_tot)
    assert np.all(np.diff(rf) > 0) and np.all(np.diff(rc) > 0)

    dr = np.diff(rf)
    with open(args.out, 'w') as f:
        f.write('# CZTR_3D grid: tools/grid/extend_grid.py %s (base cells 1..%d kept, z <= %.0f km)\n'
                % (args.base, k0, z_f[k0]))
        f.write('# then geometric growth (s %.4f -> %.3f per cell over %.0f cells) to dr = %.3f r dtheta '
                '(%.5f deg), top %.3f R_sun\n' % (s0, args.smax, args.ramp, a, args.deg, args.r_top))
        if args.outer_aspect > 1.0:
            f.write('# outer: the aspect x %.2f beyond %.2f R_sun (ramp width %.2f R_sun)\n'
                    % (args.outer_aspect, args.outer_start, args.outer_width))
        f.write('# base grid header:\n')
        for h in head:
            f.write('#   ' + h.lstrip('#').strip() + '\n')
        f.write('# cells -1..n+2 (i r_c); faces 0..n (i r_f) [cm]\n')
        f.write('%d\n' % n_tot)
        for i, v in enumerate(rc):
            f.write('%6d %25.17E\n' % (i - 1, v))
        for i, v in enumerate(rf):
            f.write('%6d %25.17E\n' % (i, v))

    zc = (0.5 * (rf[1:] + rf[:-1]) - RSUN) / 1e5
    ratio = dr[1:] / dr[:-1]
    print('cells %d (CZ %d, kept base %d, new %d), aspect %.3f, top %.4f R_sun'
          % (n_tot, args.n_cz, k0, n_new, a, rf[-1] / RSUN))
    print('largest spacing ratio between neighbours %.4f (at z %.0f km)' % (ratio.max(), zc[1:][ratio.argmax()]))
    for zz in [0, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1e6, 3e6, 6e6, 1.2e7]:
        i = np.argmin(abs(zc - zz))
        print('  z %10.0f km  dr %9.1f km  r dth %9.1f km' % (zc[i], dr[i] / 1e5,
              (RSUN + zc[i] * 1e5) * np.radians(args.deg) / 1e5))
    k = np.searchsorted(rf, RSUN * np.array([1.0, 1.01, 1.1, 1.3, 2.0, 2.5, 5.0, 10.0]))
    print('  cells above the photosphere: %d; above 1.01/1.1/1.3/2/2.5/5/10 R: %s'
          % (n_tot - k[0], [int(n_tot - kk) for kk in k[1:]]))


if __name__ == '__main__':
    main()
