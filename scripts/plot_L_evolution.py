#!/usr/bin/env python3
"""Evolution of one run's luminosity profiles: L_total and its components
(radiation, MLT, resolved, numerical diffusion) against depth, one curve per
profile (coloured by step).

Usage: /usr/bin/python3 scripts/plot_L_evolution.py <run dir> [every n-th profile] [out.png]
"""
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

run = sys.argv[1]
every = int(sys.argv[2]) if len(sys.argv) > 2 else 1
out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(run, 'L_evolution.png')
fs = sorted(glob.glob(os.path.join(run, 'profile_0*.dat')))
lg = np.loadtxt(os.path.join(run, 'log.txt'), comments='#')
P, st = [], []
for f in fs:
    h = open(f).readline().split()
    # older profiles carry only t: the step from the log
    s = int(h[h.index('step') + 1].rstrip(';')) if 'step' in h else int(round(np.interp(float(h[3]), lg[:, 1], lg[:, 0])))
    if st and s == st[-1]:
        P[-1] = np.loadtxt(f)
        continue
    st.append(s)
    P.append(np.loadtxt(f))
keep = list(range(0, len(P), every))
if keep[-1] != len(P) - 1:
    keep.append(len(P) - 1)
norm = Normalize(min(st), max(st))
cmap = plt.cm.turbo
panels = [('L_total', lambda p: p[:, 10] + p[:, 11] + p[:, 12] + p[:, 14] + p[:, 15]),
          ('radiation', lambda p: p[:, 10]), ('MLT', lambda p: p[:, 11]),
          ('resolved (enthalpy + kinetic)', lambda p: p[:, 12]), ('numerical diffusion', lambda p: p[:, 14])]
fig, ax = plt.subplots(2, 3, figsize=(17, 9))
for a, (name, fn) in zip(ax.flat, panels):
    for i in keep:
        a.plot(-P[i][:, 2]/1e3, fn(P[i]), color=cmap(norm(st[i])), lw=1.3)
    a.set_title(name)
    a.set_xscale('log')
    a.invert_xaxis()
    a.set_xlabel('depth [Mm]')
    a.axhline(0, color='gray', lw=0.5)
    a.axhline(1, color='gray', lw=0.5, ls='--')
ax.flat[-1].axis('off')
fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, label='step', pad=0.01)
fig.suptitle(f'{run}: L / Lsun (horizontal mean)')
fig.savefig(out, dpi=100, bbox_inches='tight')
print('figure:', out)
