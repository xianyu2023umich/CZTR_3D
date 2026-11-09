#!/usr/bin/env python3
"""Where the hydro operator creates or destroys energy, split into the kinetic
and the thermal energy equation (budget_split.dat of a cztr3d.x run with
dt_budget > 0; the last record, mean rates in Lsun since t0).

Usage: /usr/bin/python3 scripts/analyze_split.py <run dir> [out.png]

For the continuous equations, per shell (horizontal fluxes integrate to 0
over the sphere):
  kinetic:  dE_kin = -div L_kin + W + rsst_kin
  thermal:  dE_th  = -div L_th  - W + rsst_th + bg_work      (L_th = L_tot - L_kin)
with W = -v.grad p1 - g rho1 v_r (the pressure and buoyancy work that moves
energy between the two) and bg_work = v_r (dp0/dr + rho0 g) (the dropped
background force). The scheme's error per shell is what is left:
  S_kin = dE_kin + div L_kin - W - rsst_kin
  S_th  = dE_th  + div L_th  + W - rsst_th - bg_work
The radial fluxes are at the cell centres; at the faces they are
interpolated linearly, 0 at the two walls. The cumulative sums from the
bottom are the luminosity each error adds below a radius.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

run = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(run, 'split.png')

recs, cur, head = [], [], None
for line in open(os.path.join(run, 'budget_split.dat')):
    if line.startswith('# step='):
        head = line.split()
        cur = []
    elif line.startswith('#'):
        continue
    elif line.strip():
        cur.append([float(x) for x in line.split()])
    elif cur:
        recs.append((head, np.array(cur)))
        cur = []
if cur:
    recs.append((head, np.array(cur)))
head, A = recs[-1]
t, t0 = float(head[4]), float(head[6])
i, z = A[:, 0], A[:, 1]
dEk, dEt, W, rk, rt_all, bgw, Lt, Lk = (A[:, c] for c in range(2, 10))
rth = rt_all - rk
Lth = Lt - Lk


def div(Lc):
    """shell divergence from centre fluxes: faces by linear interpolation, walls 0"""
    n = len(Lc)
    Lf = np.zeros(n + 1)
    Lf[1:n] = 0.5*(Lc[:-1] + Lc[1:])
    return Lf[1:] - Lf[:-1]


Sk = dEk + div(Lk) - W - rk
St = dEt + div(Lth) + W - rth - bgw
depth = -z/1e3
print(f'{run}: t0 = {t0/86400:.2f} d to t = {t/86400:.2f} d ({len(recs)} records; the last)')
print(f'  whole domain [Lsun]: S_kin {Sk.sum():+.3f}   S_th {St.sum():+.3f}   total {Sk.sum() + St.sum():+.3f}')
print(f'  W (kinetic gain) {W.sum():+.3f}; rsst kin {rk.sum():+.3f} th {rth.sum():+.3f}; bg work {bgw.sum():+.3f}')
print('  shell  depth[Mm]   S_kin    S_th    cum S_kin  cum S_th   W      L_kin(c)  L_th(c)')
ck, ct = np.cumsum(Sk), np.cumsum(St)
for k in list(range(0, len(z), 4)) + [len(z) - 1]:
    print(f'  {int(i[k]):4d}  {depth[k]:8.1f}  {Sk[k]:+7.3f} {St[k]:+7.3f}   {ck[k]:+7.3f}  {ct[k]:+7.3f}  {W[k]:+7.3f}  {Lk[k]:+7.3f}  {Lth[k]:+7.3f}')

fig, ax = plt.subplots(1, 3, figsize=(17, 5.5))
ax[0].plot(depth, Sk, 'o-', ms=3, color='#2ca02c', label='S_kin (kinetic equation)')
ax[0].plot(depth, St, 'o-', ms=3, color='#1f77b4', label='S_th (thermal equation)')
ax[0].plot(depth, Sk + St, 'k-', lw=2, label='total')
ax[0].set_title('scheme error per shell [Lsun]')
ax[1].plot(depth, ck, color='#2ca02c', label='cumulative S_kin')
ax[1].plot(depth, ct, color='#1f77b4', label='cumulative S_th')
ax[1].plot(depth, ck + ct, 'k', lw=2, label='cumulative total')
ax[1].set_title('created below each radius (cumulative from the bottom) [Lsun]')
ax[2].plot(depth, W, color='#d62728', label='W (pressure + buoyancy work)')
ax[2].plot(depth, div(Lk), color='#2ca02c', label='div L_kin')
ax[2].plot(depth, div(Lth), color='#1f77b4', label='div L_th')
ax[2].plot(depth, dEk, ':', color='#2ca02c', label='dE_kin')
ax[2].plot(depth, dEt, ':', color='#1f77b4', label='dE_th')
ax[2].set_title('the terms per shell [Lsun]')
for a in ax:
    a.set_xscale('log')
    a.invert_xaxis()
    a.set_xlabel('depth [Mm]')
    a.axhline(0, color='gray', lw=0.5)
    a.legend(fontsize=8)
fig.suptitle(f'{run}: hydro energy error, days {t0/86400:.1f}-{t/86400:.1f}')
fig.tight_layout()
fig.savefig(out, dpi=100)
print('figure:', out)
