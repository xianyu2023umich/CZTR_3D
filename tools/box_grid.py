#!/usr/bin/env python3
"""Box output (out_dir/box_rNNNN.dat, see cz3d_app_mod box_write) on a
regular grid: one patch's cells (Yin by default) as arrays (snapshot, cell i,
theta j, phi k) in an .npz, for field-line tracing and plots.

Columns of the box files: step t cell z[km] lat lon patch rho T vx vy vz
Bx By Bz (global Cartesian, cm/s and G) beta Mach Q_heat Q_loss.
Assumes the Yin-Yang theta grid uniform over [pi/4, 3pi/4] (yy_margin = 0)
with n_theta cells, and phi over [-3pi/4, 3pi/4] with 3 n_theta cells (the
patch-local angles; for Yin they are the global colatitude and longitude).
Cells outside the box are NaN.

Usage: box_grid.py out_dir out.npz [--n-theta 192] [--patch 1]
"""
import argparse, glob, os
import numpy as np
import pandas as pd

COLS = ['step', 't', 'cell', 'z', 'lat', 'lon', 'patch', 'rho', 'T', 'vx', 'vy', 'vz',
        'Bx', 'By', 'Bz', 'beta', 'Mach', 'Qh', 'Ql']
VARS = ['rho', 'T', 'vx', 'vy', 'vz', 'Bx', 'By', 'Bz', 'beta', 'Mach', 'Qh', 'Ql']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out_dir'); ap.add_argument('out')
    ap.add_argument('--n-theta', type=int, default=192)
    ap.add_argument('--patch', type=int, default=1)
    a = ap.parse_args()
    parts = []
    for fn in sorted(glob.glob(os.path.join(a.out_dir, 'box_r*.dat'))):
        if os.path.getsize(fn) == 0:
            continue
        d = pd.read_csv(fn, sep=r'\s+', comment='#', header=None, names=COLS, dtype=np.float64)
        parts.append(d[d.patch == a.patch])
        print('read', os.path.basename(fn), len(d))
    d = pd.concat(parts, ignore_index=True)
    if a.patch == 2:
        raise SystemExit('Yang: lat/lon are global; the patch angles would need the Yang rotation (not done)')
    dth = 90.0/a.n_theta
    th = 90.0 - d.lat.values                      # colatitude [deg]
    j = np.rint((th - 45.0)/dth - 0.5).astype(int)
    k = np.rint((d.lon.values + 135.0)/dth - 0.5).astype(int)
    i = d.cell.values.astype(int)
    steps = np.unique(d.step.values)
    s = np.searchsorted(steps, d.step.values)
    j0, k0, i0 = j.min(), k.min(), i.min()
    shape = (len(steps), i.max() - i0 + 1, j.max() - j0 + 1, k.max() - k0 + 1)
    out = {}
    for v in VARS:
        arr = np.full(shape, np.nan, dtype=np.float32 if v not in ('rho',) else np.float64)
        arr[s, i - i0, j - j0, k - k0] = d[v].values
        out[v] = arr
    z = np.full(shape[1], np.nan)
    z[i - i0] = d.z.values
    t = np.zeros(len(steps))
    t[s] = d.t.values
    out.update(step=steps, t=t, cell=np.arange(i0, i0 + shape[1]), z_km=z,
               theta_deg=45.0 + (np.arange(j0, j0 + shape[2]) + 0.5)*dth,
               phi_deg=-135.0 + (np.arange(k0, k0 + shape[3]) + 0.5)*dth)
    np.savez(a.out, **out)
    print('%d snapshots, %d cells x %d theta x %d phi -> %s' % (shape + (a.out,)))


if __name__ == '__main__':
    main()
