#!/usr/bin/env python3
"""cz3d.x init file from the static CZ envelope (build/static_bg.x): the state
and the background are both the static model (rho1 = s1 = v = 0), on the
radial grid, in the init format cz3d.x reads (bg_load: i r z rho T v rho0 T0 g
for cells -1..n+2, after a line "n k").

Usage: python3 tools/build_static_init.py <static_profile.dat> <grid_full.dat> <out.dat>

k (the join cell, the top of the CZ domain at 1.0 R) is the last cell with
its centre at or below R_sun (tau = 0.01 of the static model). Cells up to
k + 2 come from the static model: ln rho and ln T interpolated linearly in r
(the model's ~4000 points are far finer than the grid), extrapolated linearly
above R_sun for the two ghost cells of a 1.0 R top; g = G M_r / r^2. Cells
beyond k + 2 are never used by cz3d.x (its top is at most cell k); they get the
values of cell k + 2. z = r - R_sun.
"""
import sys

import numpy as np

G, MSUN, RSUN = 6.67430e-8, 1.98892e33, 6.957e10      # as tools/envelope/const_mod.f90
prof, gridf, out = sys.argv[1:4]

S = np.loadtxt(prof, usecols=range(7))
r_s, T_s, rho_s, Mr_s = S[::-1, 0], S[::-1, 1], S[::-1, 3], S[::-1, 4]   # increasing r

# the grid: comment lines, n, cells -1..n+2 (i r_c), faces 0..n (i r_f)
lines = [l for l in open(gridf).read().split('\n') if l.strip() and not l.startswith('#')]
n = int(lines[0])
cells = np.array([[float(x) for x in l.split()] for l in lines[1:n + 5]])
idx, r = cells[:, 0].astype(int), cells[:, 1]
if idx[0] != -1 or idx[-1] != n + 2:
    raise SystemExit('unexpected grid file layout')
k = int(np.sum((idx >= 1) & (idx <= n) & (r <= RSUN)))
if r[1] < r_s[0]:
    raise SystemExit(f'static model ends at {r_s[0]:.4e} cm, above the grid bottom {r[1]:.4e}')


def interp(rr, y, log=True):
    f = np.log(y) if log else y
    v = np.interp(rr, r_s, f)
    above = rr > r_s[-1]
    if above.any():                                   # linear extrapolation above tau = 0.01
        sl = (f[-1] - f[-6])/(r_s[-1] - r_s[-6])
        v[above] = f[-1] + sl*(rr[above] - r_s[-1])
    return np.exp(v) if log else v


use = idx <= k + 2
ru = np.where(use, r, r[idx == k + 2][0])            # beyond k + 2: the values of cell k + 2
rho, T = interp(ru, rho_s), interp(ru, T_s)
g = G*interp(ru, Mr_s, log=False)*MSUN/ru**2
with open(out, 'w') as f:
    f.write('# static CZ envelope (tools/static_bg.f90) on the radial grid, state = background\n')
    f.write('# n, k\n')
    f.write('# i  r[cm]  z[km]  rho  T  v[cm/s]  rho0  T0  g[cm/s2]\n')
    f.write(f'{n} {k}\n')
    for j in range(len(idx)):
        row = [r[j], (r[j] - RSUN)/1e5, rho[j], T[j], 0.0, rho[j], T[j], g[j]]
        f.write(f'{idx[j]:5d}' + ''.join(f' {x:24.16E}' for x in row) + '\n')
print(f'wrote {out}: {n} cells, join cell k = {k}')
