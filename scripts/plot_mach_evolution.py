#!/usr/bin/env python3
"""Mach number evolution of one run: rms Mach against depth (true sound speed
and the reduced one c_s/xi), one curve per profile coloured by step; space-time
maps of both; the log's maxima against the step.

Usage: /usr/bin/python3 scripts/plot_mach_evolution.py <run dir> [every n-th profile] [out.png]

Needs out_dir/background.dat (c_s0, xi; written by cztr3d.x at every start).
rms |v| = sqrt(rms(v_r)^2 + rms(v_h)^2) of each profile; the background
sound speed (the state's differs by its perturbation, small here).
"""
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from matplotlib.cm import ScalarMappable

run = sys.argv[1]
every = int(sys.argv[2]) if len(sys.argv) > 2 else 1
out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(run, 'mach_evolution.png')
B = np.loadtxt(os.path.join(run, 'background.dat'))
cs, xi = B[:, 8], B[:, 9]
lg = np.loadtxt(os.path.join(run, 'log.txt'), comments='#')
P, st = [], []
for f in sorted(glob.glob(os.path.join(run, 'profile_0*.dat'))):
    h = open(f).readline().split()
    s = int(h[h.index('step') + 1].rstrip(';')) if 'step' in h else int(round(np.interp(float(h[3]), lg[:, 1], lg[:, 0])))
    if st and s == st[-1]:
        P[-1] = np.loadtxt(f)
        continue
    st.append(s)
    P.append(np.loadtxt(f))
st = np.array(st)
z = -P[0][:, 2]/1e3
V = np.array([np.sqrt(p[:, 8]**2 + p[:, 9]**2) for p in P])      # (profile, cell)
M, Mr = V/cs, V*xi/cs
keep = list(range(0, len(P), every))
if keep[-1] != len(P) - 1:
    keep.append(len(P) - 1)

fig, ax = plt.subplots(2, 3, figsize=(18, 10))
norm = Normalize(st.min(), st.max())
for a, Q, name in [(ax[0, 0], M, 'rms Mach (true c_s)'), (ax[0, 1], Mr, 'rms Mach (reduced c_s/xi)')]:
    for i in keep:
        a.semilogy(z, np.maximum(Q[i], 1e-8), color=plt.cm.turbo(norm(st[i])), lw=1.2)
    a.set_xscale('log')
    a.invert_xaxis()
    a.set_xlabel('depth [Mm]')
    a.set_title(name)
    a.set_ylim(1e-7, 1)
fig.colorbar(ScalarMappable(norm=norm, cmap=plt.cm.turbo), ax=ax[0, :2], label='step', pad=0.01)
a = ax[0, 2]
a.semilogy(lg[:, 0], lg[:, 3], label='max Mach (true c_s)')
a.semilogy(lg[:, 0], lg[:, 12], '--', label='max Mach (reduced)')
a.axhline(0.7, color='gray', lw=0.6, ls=':')
a.text(lg[0, 0], 0.75, 'RSST validity ~0.7', fontsize=8, color='gray')
a.set_xlabel('step')
a.set_title('log maxima')
a.legend(fontsize=8)
for a, Q, name in [(ax[1, 0], M, 'rms Mach (true c_s)'), (ax[1, 1], Mr, 'rms Mach (reduced)')]:
    im = a.pcolormesh(st, z, np.maximum(Q.T, 1e-7), norm=LogNorm(1e-6, 1e-1), cmap='magma', shading='nearest')
    a.set_yscale('log')
    a.set_ylim(z.max(), z.min())
    a.set_xlabel('step')
    a.set_ylabel('depth [Mm]')
    a.set_title(name + ': depth vs step')
    fig.colorbar(im, ax=a, pad=0.01)
a = ax[1, 2]
for d, c in zip([0.1, 1, 5, 20, 50, 150], plt.cm.viridis(np.linspace(0, 0.95, 6))):
    k = np.argmin(abs(z - d))
    a.semilogy(st, np.maximum(M[:, k], 1e-8), color=c, label=f'{z[k]:.1f} Mm (true)')
    a.semilogy(st, np.maximum(Mr[:, k], 1e-8), ':', color=c)
a.set_xlabel('step')
a.set_title('rms Mach at fixed depths (solid true, dotted reduced)')
a.legend(fontsize=7)
fig.suptitle(run)
fig.savefig(out, dpi=100, bbox_inches='tight')
print('figure:', out)
