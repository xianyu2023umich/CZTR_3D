#!/usr/bin/env python3
"""Mass exchange between up- and downflows in several runs (profiles with the
<rho v_r>, <rho |v_r|> columns, written by cztr3d.x since 2026-10-13): the mean
mass flux, the exchange rate, and the exchange times the entropy contrast
4 pi r^2 <rho |v_r|> T rms(s1) (the scale of the enthalpy the exchange moves,
in Lsun), against the total L of each run (its own relaxed profile file set).

Usage: /usr/bin/python3 scripts/compare_massflux.py out.png <dir>:<label> [<dir>:<label> ...]
"""
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LSUN = 3.828e33
out = sys.argv[1]
runs = [a.split(':', 1) for a in sys.argv[2:]]
fig, ax = plt.subplots(2, 2, figsize=(14, 10))
colors = ['k', '#1f77b4', '#d62728', '#2ca02c']
print(f'{"run":28s} {"depth":>6s}  <rho|vr|>  <rho vr>/<rho|vr|>  exch.scale[Lsun]  L_total')
for (d, lab), c in zip(runs, colors):
    P = np.mean([np.loadtxt(f) for f in sorted(glob.glob(os.path.join(d, 'profile_0*.dat')))], axis=0)
    z = -P[:, 2]/1e3
    r, T, rms_s = P[:, 1], P[:, 4], P[:, 6]
    m, ma = P[:, 18], P[:, 19]
    L = P[:, 10] + P[:, 11] + P[:, 12] + P[:, 14] + P[:, 15]
    ex = 4*np.pi*r**2*ma*T*rms_s/LSUN
    ax[0, 0].semilogy(z, ma, color=c, label=lab)
    ax[0, 1].plot(z, m/ma, color=c, label=lab)
    ax[1, 0].semilogy(z, ex, color=c, label=lab)
    ax[1, 1].plot(z, L, color=c, label=lab)
    for dd in [150, 100, 50, 20, 10, 8]:
        k = np.argmin(abs(z - dd))
        print(f'{lab:28s} {z[k]:6.0f}  {ma[k]:9.2e}  {m[k]/ma[k]:+9.2e}        {ex[k]:9.2f}       {L[k]:5.2f}')
ax[0, 0].set_title('<rho |v_r|> [g/cm2/s]: up/down mass exchange')
ax[0, 1].set_title('<rho v_r> / <rho |v_r|>: net mean mass flux')
ax[1, 0].set_title('4 pi r^2 <rho |v_r|> T rms(s1) [Lsun]: enthalpy moved by the exchange')
ax[1, 1].set_title('L_total [Lsun]')
ax[1, 1].axhline(1, color='gray', lw=0.5)
for a in ax.flat:
    a.set_xscale('log')
    a.invert_xaxis()
    a.set_xlabel('depth [Mm]')
    a.legend(fontsize=8)
fig.tight_layout()
fig.savefig(out, dpi=100)
print('figure:', out)
