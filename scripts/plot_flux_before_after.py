#!/usr/bin/env python3
"""Flux components before and after a change (e.g. rotation switched on):
the mean profile of run A over days [a0, a1] (solid) and of run B over days
[b0, b1] (dashed), days counted from each run's first profile.

Usage: /usr/bin/python3 scripts/plot_flux_before_after.py out.png <runA> a0 a1 <runB> b0 b1 [title]
"""
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def mean(run, d0, d1):
    fs = sorted(glob.glob(os.path.join(run, 'profile_0*.dat')))
    t0 = float(open(fs[0]).readline().split()[3])
    P = [np.loadtxt(f) for f in fs if d0 <= (float(open(f).readline().split()[3]) - t0)/86400 <= d1]
    return np.mean(P, axis=0), len(P)


out = sys.argv[1]
A, na = mean(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]))
B, nb = mean(sys.argv[5], float(sys.argv[6]), float(sys.argv[7]))
title = sys.argv[8] if len(sys.argv) > 8 else ''
comp = [('total', lambda p: p[:, 10] + p[:, 11] + p[:, 12] + p[:, 14] + p[:, 15], 'k'),
        ('radiation', lambda p: p[:, 10], '#d62728'), ('MLT', lambda p: p[:, 11], '#9467bd'),
        ('resolved', lambda p: p[:, 12], '#1f77b4'), ('numerical diffusion', lambda p: p[:, 14], '#ff7f0e')]
fig, ax = plt.subplots(1, 2, figsize=(15, 5.5))
for name, fn, c in comp:
    if max(abs(fn(A)).max(), abs(fn(B)).max()) < 1e-3:
        continue
    ax[0].plot(-A[:, 2]/1e3, fn(A), '-', color=c, lw=2 if name == 'total' else 1.5, label=f'{name}')
    ax[0].plot(-B[:, 2]/1e3, fn(B), '--', color=c, lw=2 if name == 'total' else 1.5)
ax[0].axhline(1, color='gray', lw=0.5)
ax[0].axhline(0, color='gray', lw=0.5)
ax[0].set_ylabel('L / Lsun')
ax[0].set_title(f'flux components: solid {os.path.basename(sys.argv[2])} ({na} profiles), dashed {os.path.basename(sys.argv[5])} ({nb})',
                fontsize=9)
ax[0].legend(fontsize=8)
ax[1].plot(-A[:, 2]/1e3, A[:, 8]/100, '-', color='#1f77b4', label='rms v_r')
ax[1].plot(-A[:, 2]/1e3, A[:, 9]/100, '-', color='#2ca02c', label='rms v_h')
ax[1].plot(-B[:, 2]/1e3, B[:, 8]/100, '--', color='#1f77b4')
ax[1].plot(-B[:, 2]/1e3, B[:, 9]/100, '--', color='#2ca02c')
ax[1].set_ylabel('m/s')
ax[1].set_title('rms velocities (solid before, dashed after)')
ax[1].legend(fontsize=8)
for a in ax:
    a.set_xlabel('depth [Mm]')
    a.set_xlim(200, 25)
fig.suptitle(title)
fig.tight_layout()
fig.savefig(out, dpi=100)
print('figure:', out)
