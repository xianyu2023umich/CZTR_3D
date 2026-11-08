#!/usr/bin/env python3
"""Horizontal kinetic energy budget per shell (budget_kh.dat of a cz3d.x run
with dt_budget > 0; the last record: mean rates in Lsun since t0).

Usage: /usr/bin/python3 scripts/analyze_kh.py <run dir> [depth_max_km]

K_h = rho v_h^2/2. For each shell, over the run:
  dK_h(hydro) = W_h - div F_Kh + R
W_h = -v_h.grad_h p1 (the horizontal pressure work), F_Kh = rho v_r v_h^2/2
(the radial flux of K_h, centre values interpolated to the faces, 0 at the
walls), R the rest: the exchange with the vertical motion through the
advection and the curvature terms, and the scheme's errors (Yin-Yang
interpolation included). The diffusion's change of K_h is separate, and
dK_h(hydro) + dK_h(diffusion) is the shell's net rate. Rows from the top
down to depth_max_km (default 3000), with sums from the top.
"""
import os
import sys

import numpy as np

run = sys.argv[1]
dmax = float(sys.argv[2]) if len(sys.argv) > 2 else 3000.0
recs, cur, head = [], [], None
for line in open(os.path.join(run, 'budget_kh.dat')):
    if line.startswith('# step='):
        head, cur = line.split(), []
    elif line.startswith('#'):
        continue
    elif line.strip():
        cur.append([float(x) for x in line.split()])
    elif cur:
        recs.append((head, np.array(cur)))
        cur = []
head, A = recs[-1]
t, t0 = float(head[4]), float(head[6])
z = -A[:, 1]
dK, Wh, F, dKd, dKdtot, dKall, W, Fall = (A[:, c] for c in range(2, 10))
n = len(z)
Ff = np.zeros(n + 1)
Ff[1:n] = 0.5*(F[:-1] + F[1:])
divF = Ff[1:] - Ff[:-1]
R = dK - Wh + divF
net = dK + dKd
print(f'{run}: K_h budget, steps over t = {t0:.0f}..{t:.0f} s ({(t - t0)/3600:.2f} h), rates in Lsun x 1e6')
print(' depth[km]   dK_h(hydro)   W_h    -divF_Kh     R      dK_h(diff)    net   | cum from top: W_h  -divF   R   diff  net')
c = np.zeros(5)
for k in range(n - 1, -1, -1):
    if z[k] > dmax:
        break
    c += [Wh[k], -divF[k], R[k], dKd[k], net[k]]
    s = 1e6
    print(f'{z[k]:9.0f}  {dK[k]*s:10.2f} {Wh[k]*s:9.2f} {-divF[k]*s:9.2f} {R[k]*s:9.2f}  {dKd[k]*s:10.2f} {net[k]*s:9.2f} |'
          f' {c[0]*s:8.2f} {c[1]*s:7.2f} {c[2]*s:7.2f} {c[3]*s:7.2f} {c[4]*s:7.2f}')
print(f'flux F_Kh at the layer bottom (centre at {z[k]:.0f} km): {F[k]*1e6:.2f} (x 1e-6 Lsun, > 0 upward)')
