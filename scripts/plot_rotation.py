#!/usr/bin/env python3
"""Differential rotation and meridional flow from rotation_NNNN.dat files
(cz3d.x with omega_rot /= 0): the mean over a range of profiles of
Omega - Omega_0 [nHz] and <v_theta> on the meridional plane, Omega against
radius at several latitudes, and the north-south symmetric part.

Usage: /usr/bin/python3 scripts/plot_rotation.py <run dir> [first] [last] [out.png]

first, last: profile numbers (rotation_NNNN.dat), default the last 30.
"""
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

run = sys.argv[1]
fs = sorted(glob.glob(os.path.join(run, 'rotation_*.dat')))
first = int(sys.argv[2]) if len(sys.argv) > 2 else max(0, len(fs) - 30)
last = int(sys.argv[3]) if len(sys.argv) > 3 else len(fs) - 1
out = sys.argv[4] if len(sys.argv) > 4 else os.path.join(run, f'rotation_{first:04d}-{last:04d}.png')
sel = [f for f in fs if first <= int(f[-8:-4]) <= last]
D = np.mean([np.loadtxt(f) for f in sel], axis=0)
t0 = float(open(sel[0]).readline().split()[3])/86400
t1 = float(open(sel[-1]).readline().split()[3])/86400
ir = D[:, 0].astype(int)
r = D[:, 1]/6.957e10
lat = D[:, 3]
nr, nl = len(set(ir)), len(set(lat))
R = r.reshape(nr, nl)
L = np.radians(lat.reshape(nr, nl))
Om = D[:, 5].reshape(nr, nl)
Vt = D[:, 6].reshape(nr, nl)/100.0
Oms = 0.5*(Om + Om[:, ::-1])                      # north-south symmetric part

fig = plt.figure(figsize=(17, 6))
for n, (Q, name, unit) in enumerate([(Oms, 'Omega - Omega_0 (N-S symmetric)', 'nHz'), (Vt, '<v_theta> (> 0 southward)', 'm/s')]):
    ax = fig.add_subplot(1, 3, n + 1)
    X, Y = R*np.cos(L), R*np.sin(L)
    v = np.nanmax(abs(Q))
    pc = ax.pcolormesh(X, Y, Q, cmap='RdBu_r', norm=TwoSlopeNorm(0, -v, v), shading='gouraud')
    ax.contour(X, Y, Q, 10, colors='k', linewidths=0.4)
    ax.set_aspect('equal')
    ax.set_xlabel('r cos(lat) / R_sun')
    ax.set_ylabel('r sin(lat) / R_sun')
    ax.set_title(name)
    fig.colorbar(pc, ax=ax, label=unit, shrink=0.8)
ax = fig.add_subplot(1, 3, 3)
for la, c in zip([0, 15, 30, 45, 60, 75], plt.cm.viridis(np.linspace(0, 0.9, 6))):
    m = np.argmin(abs(lat.reshape(nr, nl)[0] - la))
    ax.plot(R[:, m], Oms[:, m], color=c, label=f'{lat.reshape(nr, nl)[0, m]:.1f} deg')
ax.axhline(0, color='gray', lw=0.5)
ax.set_xlabel('r / R_sun')
ax.set_ylabel('Omega - Omega_0 [nHz]')
ax.set_title('by latitude (N-S symmetric)')
ax.legend(fontsize=8)
fig.suptitle(f'{run}: mean of {len(sel)} rotation profiles, day {t0:.1f}-{t1:.1f}')
fig.tight_layout()
fig.savefig(out, dpi=100)
print('figure:', out)
print(f'Omega - Omega_0 at the top shell: equator {Oms[-1, nl//2]:+.1f} nHz, lat ~60 {Oms[-1, np.argmin(abs(lat.reshape(nr, nl)[0] - 60))]:+.1f} nHz')
