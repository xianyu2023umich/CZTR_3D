#!/usr/bin/env python3
"""Compare a cztr3d.x run's horizontally averaged profiles with a CZTR (1D) run.

Usage: /usr/bin/python3 scripts/compare_1d.py <3D run dir> <1D run dir> [out.png]

Profiles are matched by number (profile_NNNN.dat, same dt_prof / dt_out),
cells by radial index. Prints the largest differences per profile and plots
the last pair.
"""
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d3, d1 = sys.argv[1], sys.argv[2]
out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(d3, 'compare_1d.png')


def t_of(f):
    return float(open(f).readline().split()[3])


pairs = []
for f3 in sorted(glob.glob(os.path.join(d3, 'profile_0*.dat'))):
    f1 = os.path.join(d1, os.path.basename(f3))
    if os.path.exists(f1):
        pairs.append((f3, f1))
if not pairs:
    sys.exit('no matching profiles')

print(' t3[s]    t1[s]   max|drho/rho|  max|dT/T|   max|dv_r|[cm/s]  max|dL|/Lsun  max rms_h(s1)  max rms(v_h)')
for f3, f1 in pairs:
    a, b = np.loadtxt(f3), np.loadtxt(f1)
    n = a.shape[0]
    b = b[:n]
    L3 = a[:, 10] + a[:, 11]
    L1 = b[:, 10] + b[:, 11]
    print(f'{t_of(f3):7.0f}  {t_of(f1):7.0f}   {np.abs(a[:,3]/b[:,3]-1).max():.3e}    {np.abs(a[:,4]/b[:,4]-1).max():.3e}'
          f'    {np.abs(a[:,7]-b[:,6]).max():.3e}       {np.abs(L3-L1).max():.3e}     {a[:,6].max():.3e}     {a[:,9].max():.3e}')

a, b = np.loadtxt(pairs[-1][0]), np.loadtxt(pairs[-1][1])
b = b[:a.shape[0]]
z = a[:, 2]
fig, ax = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
ax[0, 0].plot(z, a[:, 3]/b[:, 3] - 1)
ax[0, 0].set_ylabel('<rho>_3D / rho_1D - 1')
ax[0, 1].plot(z, a[:, 4]/b[:, 4] - 1)
ax[0, 1].set_ylabel('<T>_3D / T_1D - 1')
ax[1, 0].plot(z, a[:, 7], label='<v_r> 3D')
ax[1, 0].plot(z, b[:, 6], '--', label='v 1D')
ax[1, 0].plot(z, a[:, 9], ':', label='rms v_h 3D')
ax[1, 0].set_ylabel('cm/s')
ax[1, 0].legend(fontsize=8)
ax[1, 1].plot(z, a[:, 10] + a[:, 11], label='3D')
ax[1, 1].plot(z, b[:, 10] + b[:, 11], '--', label='1D')
ax[1, 1].set_ylabel('L_rad + L_conv [Lsun]')
ax[1, 1].legend(fontsize=8)
for x in ax.flat:
    x.set_xscale('symlog', linthresh=100)
    x.set_xlabel('z [km]')
fig.suptitle(f'3D (horizontal mean) vs 1D at t = {t_of(pairs[-1][0]):.0f} s')
fig.tight_layout()
fig.savefig(out, dpi=110)
print('figure:', out)
