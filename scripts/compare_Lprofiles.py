#!/usr/bin/env python3
"""L_total(depth) profiles of two runs at a series of steps, side by side;
the colour is the step on one scale for both (same colour = same step).

Usage: /usr/bin/python3 scripts/compare_Lprofiles.py <run A> <run B> [step interval] [max step] [out.png]
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

a, b = sys.argv[1], sys.argv[2]
every = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
smax = int(sys.argv[4]) if len(sys.argv) > 4 else 40000
out = sys.argv[5] if len(sys.argv) > 5 else os.path.join(a, 'compare_Lprofiles.png')


def load(run):
    lg = np.loadtxt(os.path.join(run, 'log.txt'), comments='#')
    res = {}
    for f in sorted(glob.glob(os.path.join(run, 'profile_0*.dat'))):
        h = open(f).readline().split()
        s = int(h[h.index('step') + 1].rstrip(';')) if 'step' in h else int(round(np.interp(float(h[3]), lg[:, 1], lg[:, 0])))
        res[s] = np.loadtxt(f)
    return res


def L_total(p):
    return p[:, 10] + p[:, 11] + p[:, 12] + p[:, 14] + p[:, 15]


cmap = plt.cm.turbo
norm = Normalize(0, smax)
fig, ax = plt.subplots(1, 2, figsize=(15, 6), sharey=True)
for x, run in [(ax[0], a), (ax[1], b)]:
    P = load(run)
    amax = max(p[:, 17].max() if p.shape[1] > 17 else 1.0 for p in P.values())
    lab = os.path.basename(run.rstrip('/')) + (' (time accurate)' if amax <= 1.0 else f' (alpha up to {amax:.1f})')
    steps = np.array(sorted(P))
    last = steps[steps <= smax].max()
    for target in list(range(0, smax + 1, every)):
        if target > last:
            break
        s = steps[np.abs(steps - target).argmin()]
        if abs(s - target) > every/4:
            continue
        p = P[s]
        x.plot(-p[:, 2]/1e3, L_total(p), color=cmap(norm(s)), lw=1.4)
    x.axhline(1.0, color='gray', lw=0.6)
    x.set_xscale('log')
    x.invert_xaxis()
    x.set_xlabel('depth [Mm]')
    x.set_title(f'{lab}: every {every} steps to step {last}')
ax[0].set_ylabel('L_total / Lsun (horizontal mean)')
ax[0].set_ylim(-0.1, 3.5)
fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, label='step', pad=0.01)
fig.savefig(out, dpi=100, bbox_inches='tight')
print('figure:', out)
