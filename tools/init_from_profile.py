#!/usr/bin/env python3
"""An init file whose atmosphere is a relaxed 1D state: the cells above the
join (k+1..n) of an init file (state columns rho, T, v) replaced by the
horizontal means <rho>, <T>, <v_r> of a profile file (e.g. the steady solar
wind of a CZTR_2D 1D-mode run, runs/wind1_lts20.nml), on the same radial grid.
The grid, the background (rho0, T0, g), the CZ cells and the ghost rows are
kept. Checks that rho0 <= rho still holds above the join.

Usage: init_from_profile.py base_init.dat profile.dat out.dat
"""
import sys
import numpy as np


def main():
    base, prof, out = sys.argv[1:4]
    lines = open(base).read().split('\n')
    head = [l for l in lines if l.startswith('#')]
    body = [l for l in lines if l.strip() and not l.startswith('#')]
    n, k = map(int, body[0].split())
    rows = np.array([[float(x) for x in l.split()] for l in body[1:n + 5]])   # cells -1..n+2
    P = np.loadtxt(prof)                                                     # cells 1..n
    if P.shape[0] != n:
        raise SystemExit('profile has %d cells, the init file %d' % (P.shape[0], n))
    if not np.allclose(P[:, 1], rows[2:n + 2, 1], rtol=1e-10, atol=0):
        raise SystemExit('profile and init file radii differ')
    atm = np.arange(1, n + 1) > k
    j = np.where(atm)[0] + 2                                                 # rows of cells k+1..n
    rows[j, 3] = P[atm, 3]                                                   # rho
    rows[j, 4] = P[atm, 4]                                                   # T
    rows[j, 5] = P[atm, 7]                                                   # v_r
    ratio = rows[j, 3] / rows[j, 6]
    if ratio.min() < 1.0:
        print('WARNING: rho0 > rho in %d cells (min rho/rho0 %.4f at z %.0f km)'
              % ((ratio < 1).sum(), ratio.min(), rows[j, 2][ratio.argmin()]))
    t = open(prof).readline().strip()
    with open(out, 'w') as f:
        f.write('# CZTR_3D init: tools/init_from_profile.py %s %s\n' % (base, prof))
        f.write('# the atmosphere (cells %d..%d) from the profile (%s); the rest from the base file:\n'
                % (k + 1, n, t.lstrip('# ')))
        for h in head:
            f.write('#   ' + h.lstrip('#').strip() + '\n')
        f.write('%d %d\n' % (n, k))
        for r in rows:
            f.write('%6d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n'
                    % (int(r[0]), *r[1:]))
    print('atmosphere cells %d..%d from %s; min rho/rho0 above the join %.3f; v_r at the top %.1f km/s'
          % (k + 1, n, prof, ratio.min(), rows[n + 1, 5] / 1e5))


if __name__ == '__main__':
    main()
