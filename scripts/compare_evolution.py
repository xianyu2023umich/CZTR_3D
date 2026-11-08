#!/usr/bin/env python3
"""Evolution of several runs against the step: log time series and space-time
maps of the horizontally averaged L_total(depth, step) on common axes.

Usage: /usr/bin/python3 scripts/compare_evolution.py <run> <run> [<run> ...] [max=40000] [out=file.png]

Labels: the run directory's name and its dt1 (from the profiles' tsc column:
'time accurate' when tsc = 1 throughout).
"""
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

runs = [a for a in sys.argv[1:] if '=' not in a]
opt = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
smax = int(opt.get('max', 40000))
out = opt.get('out', os.path.join(runs[0], 'compare_evolution.png'))


def load(run):
    lg = np.loadtxt(os.path.join(run, 'log.txt'), comments='#')
    st, prof = [], []
    for f in sorted(glob.glob(os.path.join(run, 'profile_0*.dat'))):
        h = open(f).readline().split()
        s = int(h[h.index('step') + 1].rstrip(';')) if 'step' in h else int(round(np.interp(float(h[3]), lg[:, 1], lg[:, 0])))
        if s <= smax:
            st.append(s)
            prof.append(np.loadtxt(f))
    o = np.argsort(st)
    st = np.array(st)[o]
    prof = [prof[i] for i in o]
    # restarts repeat a step: keep the last profile of each step
    keep = [i for i in range(len(st)) if i == len(st) - 1 or st[i+1] != st[i]]
    return lg[lg[:, 0] <= smax], st[keep], [prof[i] for i in keep]


def L_total(p):
    return p[:, 10] + p[:, 11] + p[:, 12] + p[:, 14] + p[:, 15]


def label(run, P):
    a = max(p[:, 17].max() if p.shape[1] > 17 else 1.0 for p in P)
    return f'{os.path.basename(run.rstrip("/"))} ' + ('(time accurate)' if a <= 1.0 else f'(alpha up to {a:.1f})')


data = [(r, load(r)) for r in runs]
labs = [label(r, d[2]) for r, d in data]
styles = ['--', '-', '-.', ':']
colors = ['k', '#1f77b4', '#d62728', '#2ca02c']
n = len(data)
fig = plt.figure(figsize=(5.5*max(n, 3), 10))
gs = fig.add_gridspec(2, max(n, 3), height_ratios=[1, 1.25])
ax1, ax2, ax3 = (fig.add_subplot(gs[0, i]) for i in range(3))
for (run, (lg, st, P)), lab, ls, c in zip(data, labs, styles, colors):
    ax1.semilogy(lg[:, 0], lg[:, 5], ls, color=c, label=lab)
    ax2.plot(lg[:, 0], lg[:, 15], ls, color=c, label=lab)
    ax3.plot(st, [L_total(p)[0] for p in P], ls, color=c, label=lab)
ax1.set_title('max |v_h| [cm/s]')
ax2.set_title('max over depth of L_total [Lsun]')
ax3.set_title('L_total at the bottom cell [Lsun]')
for x in (ax1, ax2, ax3):
    x.set_xlabel('step')
    x.set_xlim(0, smax)
    x.legend(fontsize=8)
ax2.axhline(1.0, color='gray', lw=0.5)
ax3.axhline(1.0, color='gray', lw=0.5)
norm = TwoSlopeNorm(vmin=0.0, vcenter=1.0, vmax=2.5)
maps = []
for i, ((run, (lg, st, P)), lab) in enumerate(zip(data, labs)):
    x = fig.add_subplot(gs[1, i])
    depth = -P[0][:, 2]/1e3
    M = np.array([L_total(p) for p in P]).T
    im = x.pcolormesh(st, depth, M, cmap='RdBu_r', norm=norm, shading='nearest')
    x.set_yscale('log')
    x.set_ylim(depth.max(), depth.min())
    x.set_xlim(0, smax)
    x.set_xlabel('step')
    x.set_ylabel('depth [Mm]')
    x.set_title(f'L_total: {lab}', fontsize=10)
    maps.append(x)
fig.colorbar(im, ax=maps, label='L_total / Lsun (1 = white)', pad=0.01)
fig.savefig(out, dpi=100, bbox_inches='tight')
print('figure:', out)
