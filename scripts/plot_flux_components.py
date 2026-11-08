#!/usr/bin/env python3
"""Luminosity by component, averaged over a run's profiles.

Usage: /usr/bin/python3 scripts/plot_flux_components.py <run dir> [first profile] [out.png]

Averages profile_NNNN.dat from `first` (default 1: the first profile after a
start or restart has no diffusion flux yet) and plots L_rad, L_conv, L_enth,
L_kin, L_diff, L_cool and their total against depth.
"""
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

run = sys.argv[1]
first = int(sys.argv[2]) if len(sys.argv) > 2 else 1
out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(run, 'flux_components.png')

fs = sorted(glob.glob(os.path.join(run, 'profile_0*.dat')))[first:]
P = np.mean([np.loadtxt(f) for f in fs], axis=0)
t0 = float(open(fs[0]).readline().split()[3])/86400.0
t1 = float(open(fs[-1]).readline().split()[3])/86400.0
z = P[:, 2]
L = {'radiation': P[:, 10], 'MLT': P[:, 11], 'enthalpy (resolved)': P[:, 12] - P[:, 13],
     'kinetic (resolved)': P[:, 13], 'numerical diffusion': P[:, 14], 'cooling layer': P[:, 15]}
total = sum(L.values())

fig, ax = plt.subplots(figsize=(9, 6))
styles = {'radiation': '#d62728', 'MLT': '#9467bd', 'enthalpy (resolved)': '#1f77b4', 'kinetic (resolved)': '#2ca02c',
          'numerical diffusion': '#ff7f0e', 'cooling layer': '#17becf'}
for name, y in L.items():
    if np.abs(y).max() < 1e-6:
        continue
    ax.plot(-z, y, color=styles[name], lw=1.8, label=name)
ax.plot(-z, total, color='black', lw=2.6, label='total')
ax.axhline(1.0, color='gray', lw=1, ls='--')
ax.text(-z[0], 1.05, 'L$_\\odot$', color='gray', ha='left', va='bottom')
ax.axhline(0.0, color='gray', lw=0.5)
ax.set_xscale('log')
ax.invert_xaxis()
ax.set_xlabel('depth [km]')
ax.set_ylabel('L / L$_\\odot$ (horizontal mean)')
ax.set_title(f'{run}: luminosity by component, mean of {len(fs)} profiles, day {t0:.2f}–{t1:.2f}')
ax.legend(loc='upper left', frameon=False)
fig.tight_layout()
fig.savefig(out, dpi=110)
print('figure:', out)
