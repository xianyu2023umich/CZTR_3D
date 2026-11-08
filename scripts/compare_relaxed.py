#!/usr/bin/env python3
"""Relaxed states of two runs: profiles averaged over a step window of each.

Usage: /usr/bin/python3 scripts/compare_relaxed.py <run A> <stepA0> <stepA1> <run B> <stepB0> <stepB1> [out.png]

A profile's step comes from its header, or from its t through the run's log.
Plots L by component, <T> (A/B - 1), rms(v_h), rms_h(s1) against depth.
"""
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def mean_profile(run, s0, s1):
    lg = np.loadtxt(os.path.join(run, 'log.txt'), comments='#')
    ps = []
    for f in sorted(glob.glob(os.path.join(run, 'profile_0*.dat'))):
        h = open(f).readline().split()
        step = int(h[h.index('step') + 1].rstrip(';')) if 'step' in h else \
            int(round(np.interp(float(h[3]), lg[:, 1], lg[:, 0])))
        if s0 <= step <= s1:
            ps.append(np.loadtxt(f)[:, :16])
    return np.mean(ps, axis=0), len(ps)


a, sa0, sa1, b, sb0, sb1 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], int(sys.argv[5]), int(sys.argv[6])
out = sys.argv[7] if len(sys.argv) > 7 else os.path.join(a, 'compare_relaxed.png')
A, na = mean_profile(a, sa0, sa1)
B, nb = mean_profile(b, sb0, sb1)
z = -A[:, 2]
comp = {'radiation': 10, 'MLT': 11, 'enthalpy (resolved)': None, 'kinetic (resolved)': 13, 'diffusion': 14, 'cooling': 15}
colors = ['#d62728', '#9467bd', '#1f77b4', '#2ca02c', '#ff7f0e', '#17becf']
fig, ax = plt.subplots(2, 2, figsize=(13, 9))
for (name, c), col in zip(comp.items(), colors):
    for P, ls in [(A, '-'), (B, '--')]:
        y = P[:, 12] - P[:, 13] if c is None else P[:, c]
        ax[0, 0].plot(-P[:, 2], y, ls, color=col, label=name if ls == '-' else None)
for P, ls, lab in [(A, '-', f'{a} ({na} profiles)'), (B, '--', f'{b} ({nb} profiles)')]:
    ax[0, 0].plot(-P[:, 2], P[:, 10] + P[:, 11] + P[:, 12] + P[:, 14] + P[:, 15], ls, color='k', lw=2, label='total: ' + lab)
    ax[1, 0].semilogy(-P[:, 2], P[:, 9], ls, color='k', label=lab)
    ax[1, 1].semilogy(-P[:, 2], P[:, 6], ls, color='k', label=lab)
ax[0, 1].plot(z, A[:, 4]/B[:, 4] - 1, color='k')
ax[0, 0].set_title('L / Lsun by component (solid A, dashed B)')
ax[0, 1].set_title('<T>_A / <T>_B - 1')
ax[1, 0].set_title('rms(v_h) [cm/s]')
ax[1, 1].set_title('rms_h(s1) [erg/g/K]')
for x in ax.flat:
    x.set_xscale('log')
    x.invert_xaxis()
    x.set_xlabel('depth [km]')
ax[0, 0].legend(fontsize=7)
ax[1, 0].legend(fontsize=7)
fig.suptitle(f'relaxed: {a} steps {sa0}-{sa1} vs {b} steps {sb0}-{sb1}')
fig.tight_layout()
fig.savefig(out, dpi=100)
print('figure:', out)
d = A[:, 4]/B[:, 4] - 1
print(f'max |dT/T| {np.abs(d).max():.2e} (at {z[np.abs(d).argmax()]/1e3:.0f} Mm); '
      f'rms(v_h) ratio A/B: {np.min(A[:,9]/B[:,9]):.2f}..{np.max(A[:,9]/B[:,9]):.2f}; '
      f'rms(v_r) ratio: {np.min(A[:,8]/B[:,8]):.2f}..{np.max(A[:,8]/B[:,8]):.2f}')
