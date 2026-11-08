#!/usr/bin/env python3
"""Local time stepping (dt1) against a time-accurate run, at equal step numbers.

Usage: /usr/bin/python3 scripts/compare_lts.py <lts run> <reference run> [out.png]

Top row: log time series against the step (max |v_h|, max rms_h(s1), the
range of L_total over depth). Bottom row: profiles at a few common steps
(the reference's profile nearest in step; its steps from its log's t),
L_total and rms(v_h) against depth, and the shells' own clocks of the lts run.
"""
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

lts, ref = sys.argv[1], sys.argv[2]
out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(lts, 'compare_lts.png')


def log(run):
    return np.loadtxt(os.path.join(run, 'log.txt'), comments='#')


def profiles(run, lg):
    """[(step, t, profile)], the step from the header or from the log's t"""
    res = []
    for f in sorted(glob.glob(os.path.join(run, 'profile_0*.dat'))):
        h = open(f).readline().split()
        t = float(h[3])
        if 'step' in h:
            step = int(h[h.index('step') + 1].rstrip(';'))
        else:
            step = int(round(np.interp(t, lg[:, 1], lg[:, 0])))
        res.append((step, t, np.loadtxt(f)))
    return res


def L_total(p):
    # radiation + MLT + resolved enthalpy and kinetic + diffusion + cooling
    return p[:, 10] + p[:, 11] + p[:, 12] + p[:, 14] + p[:, 15]


La, Lb = log(lts), log(ref)
Pa, Pb = profiles(lts, La), profiles(ref, Lb)
n_max = La[-1, 0]
Lb = Lb[Lb[:, 0] <= n_max]

fig, ax = plt.subplots(2, 3, figsize=(17, 9))
for lg, lab, ls in [(La, 'local stepping', '-'), (Lb, 'time accurate', '--')]:
    ax[0, 0].semilogy(lg[:, 0], lg[:, 5], ls, color='k', label=lab)
    ax[0, 1].semilogy(lg[:, 0], lg[:, 9], ls, color='k', label=lab)
    ax[0, 2].plot(lg[:, 0], lg[:, 14], ls, color='#1f77b4', label=f'min L ({lab})')
    ax[0, 2].plot(lg[:, 0], lg[:, 15], ls, color='#d62728', label=f'max L ({lab})')
ax[0, 0].set_title('max |v_h| [cm/s]')
ax[0, 1].set_title('max rms_h(s1) [erg/g/K]')
ax[0, 2].set_title('L_total over depth [Lsun]')
ax[0, 2].axhline(1.0, color='gray', lw=0.5)
for a in ax[0]:
    a.set_xlabel('step')
    a.legend(fontsize=8)

# profiles at up to 5 common steps
steps = [s for s, _, _ in Pa if s > 0]
pick = [steps[i] for i in sorted(set(np.linspace(0, len(steps) - 1, min(5, len(steps))).round().astype(int)))]
cm = plt.cm.viridis(np.linspace(0, 0.9, len(pick)))
for s, c in zip(pick, cm):
    pa = next(p for st, _, p in Pa if st == s)
    sb, tb, pb = min(Pb, key=lambda x: abs(x[0] - s))
    if sb > n_max + 600:
        continue
    z = -pa[:, 2]
    ax[1, 0].plot(z, L_total(pa), '-', color=c, label=f'step {s}')
    ax[1, 0].plot(-pb[:, 2], L_total(pb), '--', color=c)
    ax[1, 1].semilogy(z, pa[:, 9], '-', color=c, label=f'step {s}')
    ax[1, 1].semilogy(-pb[:, 2], np.maximum(pb[:, 9], 1e-6), '--', color=c)
    if pa.shape[1] > 16:
        ax[1, 2].plot(z, pa[:, 16]/86400.0, color=c, label=f'step {s} (top clock {pa[-1, 16]/86400:.1f} d)')
ax[1, 0].set_title('L_total (solid: local stepping, dashed: time accurate)')
ax[1, 0].axhline(1.0, color='gray', lw=0.5)
ax[1, 1].set_title('rms(v_h) [cm/s]')
ax[1, 2].set_title("local stepping: each shell's own clock [days]")
for a in ax[1]:
    a.set_xscale('log')
    a.invert_xaxis()
    a.set_xlabel('depth [km]')
    a.legend(fontsize=8)
fig.suptitle(f'{lts} vs {ref}, equal steps')
fig.tight_layout()
fig.savefig(out, dpi=100)
print('figure:', out)
