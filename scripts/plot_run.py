#!/usr/bin/env python3
"""Overview of a cztr3d.x run: time series from log.txt and horizontally
averaged profiles (profile_NNNN.dat).

Usage: /usr/bin/python3 scripts/plot_run.py <run dir> [out.png]
"""
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

run = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(run, 'overview.png')

lg = np.loadtxt(os.path.join(run, 'log.txt'), comments='#')
# time unit: hours, or days for runs longer than 3 days
tscale, unit = (86400.0, 'd') if lg[-1, 1] > 3*86400.0 else (3600.0, 'h')
th = lg[:, 1]/tscale
profs = sorted(glob.glob(os.path.join(run, 'profile_0*.dat')))
P = [np.loadtxt(f) for f in profs]
pt = [float(open(f).readline().split()[3])/tscale for f in profs]
# at most ~12 profile curves, evenly spaced, always the last
keep = sorted(set(np.linspace(0, len(P) - 1, min(12, len(P))).round().astype(int)))
P = [P[i] for i in keep]
pt = [pt[i] for i in keep]

fig, ax = plt.subplots(2, 3, figsize=(16, 9))
a = ax[0, 0]
a.semilogy(th, lg[:, 5], label='max |v_h| [cm/s]')
a.semilogy(th, lg[:, 4], label='max |v_r| [cm/s]')
a.set_xlabel(f't [{unit}]')
a.legend()
a.set_title('velocities')
a = ax[0, 1]
a.semilogy(th, lg[:, 3], label='max Mach')
if lg.shape[1] > 12:
    a.semilogy(th, lg[:, 12], label='max Mach (reduced c_s)')
a.semilogy(th, lg[:, 9], label='max rms_h(s1) [erg/g/K]')
a.set_xlabel(f't [{unit}]')
a.legend()
a.set_title('Mach, horizontal entropy fluctuation')
a = ax[0, 2]
a.plot(th, lg[:, 6], label='L_top')
if lg.shape[1] > 15:
    a.plot(th, lg[:, 14], label='min L total')
    a.plot(th, lg[:, 15], label='max L total')
    a.plot(th, lg[:, 13], ':', label='max |L resolved|')
else:
    a.plot(th, lg[:, 7], label='min L')
    a.plot(th, lg[:, 8], label='max L')
a.set_xlabel(f't [{unit}]')
a.set_ylabel('L / Lsun (horizontal mean)')
a.legend()
cm = plt.cm.viridis(np.linspace(0, 1, len(P)))
for p, t, c in zip(P, pt, cm):
    z = p[:, 2]
    ax[1, 0].semilogy(-z, np.maximum(p[:, 9], 1e-12), color=c, label=f'{t:.1f} {unit}')
    ax[1, 1].semilogy(-z, np.maximum(p[:, 6], 1e-12), color=c)
    Lres = p[:, 12] if p.shape[1] > 12 else 0.0*z
    ax[1, 2].plot(-z, p[:, 10] + p[:, 11] + Lres, color=c)
    if p.shape[1] > 12:
        ax[1, 2].plot(-z, Lres, ':', color=c)
for a, lab in zip(ax[1], ['rms v_h [cm/s]', 'rms_h(s1) [erg/g/K]', 'L total (solid), L_enth+L_kin resolved (dotted) [Lsun]']):
    a.set_xscale('log')
    a.invert_xaxis()
    a.set_xlabel('depth [km]')
    a.set_ylabel(lab)
ax[1, 0].legend(fontsize=7, ncol=2)
fig.suptitle(run)
fig.tight_layout()
fig.savefig(out, dpi=100)
print('figure:', out)
i = -1
print(f't = {lg[i,1]:.0f} s ({th[i]:.2f} h), step {int(lg[i,0])}: max Mach {lg[i,3]:.3e}, max |v_h| {lg[i,5]:.3e} cm/s, '
      f'L_top {lg[i,6]:.4f}, max rms_h(s1) {lg[i,9]:.3e}')
