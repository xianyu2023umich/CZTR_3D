#!/usr/bin/env python3
"""Horizontal flow structure on spherical shells (cztr3d.x dump_shells files).

Usage: /usr/bin/python3 scripts/analyze_shells.py <dump dir> [out.png]

Per shell: the weighted rms of v_h; the part that is a rigid rotation of the
shell (v_h = Omega x r_hat, least squares) and the part that is the l = 1
poloidal flow (v_h = a - (a.r_hat) r_hat, from pole to antipode); the angular
correlation of v_h (<v_h(x).v_h(x')> / <v_h^2> against the angle between x
and x', from a random subsample). Maps (Mollweide, global frame) of |v_h| and
its direction for two shells.
"""
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(d, 'shells.png')
A = np.vstack([np.loadtxt(f) for f in sorted(glob.glob(os.path.join(d, 'shell_dump_r*.dat')))])
B = np.loadtxt(os.path.join(d, 'background.dat'))
depth = {int(r[0]): -r[2] for r in B}
shells = sorted(set(A[:, 0].astype(int)), reverse=True)
rng = np.random.default_rng(1)
bins = np.linspace(0, np.pi, 37)
corr = {}

print(' shell  depth[km]  rms v_h   max|v_h| [m/s]  rigid rot.  l=1 pol.  corr. angle(1/2)  rms v_r  min rho1/rho0')
for i in shells:
    S = A[(A[:, 0] == i) & (A[:, 5] > 0)]
    x, w, vr, vh, dr = S[:, 2:5], S[:, 5], S[:, 6], S[:, 7:10], S[:, 10]
    W = w.sum()
    e_tot = np.sum(w*np.sum(vh**2, 1))
    # rigid rotation: v = Omega x x = -[x]_x Omega; normal equations
    M = np.zeros((3, 3))
    b = np.zeros(3)
    P = np.zeros((3, 3))
    c = np.zeros(3)
    for xi, wi, vi in zip(x, w, vh):
        X = np.array([[0, -xi[2], xi[1]], [xi[2], 0, -xi[0]], [-xi[1], xi[0], 0]])   # X Omega = x cross Omega
        R = -X                                                                          # Omega cross x
        M += wi*R.T@R
        b += wi*R.T@vi
        Q = np.eye(3) - np.outer(xi, xi)                                                # a - (a.x) x
        P += wi*Q.T@Q
        c += wi*Q.T@vi
    Om = np.linalg.solve(M, b)
    a = np.linalg.solve(P, c)
    e_rot = np.sum(w*np.sum(np.cross(Om, x)**2, 1))
    e_pol = np.sum(w*np.sum((a - (x@a)[:, None]*x)**2, 1))
    # angular correlation from a subsample
    k = rng.choice(len(x), size=min(2500, len(x)), replace=False)
    xs, vs = x[k], vh[k]
    g = np.arccos(np.clip(xs@xs.T, -1, 1)).ravel()
    pr = (vs@vs.T).ravel()
    h, _ = np.histogram(g, bins, weights=pr)
    n, _ = np.histogram(g, bins)
    cc = h/np.maximum(n, 1)/np.mean(np.sum(vs**2, 1))
    corr[i] = cc
    centers = 0.5*(bins[1:] + bins[:-1])
    half = centers[np.argmax(cc < 0.5)] if np.any(cc < 0.5) else np.pi
    print(f'{i:5d}  {depth[i]:9.1f}  {np.sqrt(e_tot/W)/100:7.1f}  {np.sqrt(np.max(np.sum(vh**2, 1)))/100:8.1f}       '
          f'{e_rot/e_tot:8.3f}  {e_pol/e_tot:8.3f}     {np.degrees(half):6.1f} deg      {np.sqrt(np.sum(w*vr**2)/W)/100:6.1f}  {dr.min():+.2e}')

fig = plt.figure(figsize=(16, 10))
for n, i in enumerate([shells[0], shells[-3]]):
    S = A[(A[:, 0] == i) & (A[:, 5] > 0)]
    x, vh = S[:, 2:5], S[:, 7:10]
    lat = np.arcsin(np.clip(x[:, 2], -1, 1))
    lon = np.arctan2(x[:, 1], x[:, 0])
    el = np.stack([-np.sin(lon), np.cos(lon), 0*lon], 1)
    en = np.stack([-np.sin(lat)*np.cos(lon), -np.sin(lat)*np.sin(lon), np.cos(lat)], 1)
    ue, un = np.sum(vh*el, 1), np.sum(vh*en, 1)
    ax = fig.add_subplot(2, 2, 1 + 2*n, projection='mollweide')
    sc = ax.scatter(lon, lat, c=np.hypot(ue, un)/100, s=2, cmap='viridis')
    q = rng.choice(len(lon), 1200, replace=False)
    ax.quiver(lon[q], lat[q], ue[q], un[q], color='w', scale=None, width=0.0015)
    fig.colorbar(sc, ax=ax, shrink=0.6, label='|v_h| [m/s]')
    ax.set_title(f'shell {i} ({depth[i]:.0f} km depth): |v_h| and direction')
    ax.grid(True, alpha=0.3)
ax = fig.add_subplot(2, 2, 2)
centers = 0.5*(bins[1:] + bins[:-1])
for i in shells:
    ax.plot(np.degrees(centers), corr[i], label=f'{depth[i]:.0f} km')
ax.axhline(0, color='gray', lw=0.5)
ax.set_xlabel('angle between points [deg]')
ax.set_title('angular correlation of v_h')
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(out, dpi=100)
print('figure:', out)
