#!/usr/bin/env python3
"""A state export (export_state: out_dir/state_export_rNNNN.bin and
state_export_meta.bin) made coarser in theta and phi: every block refined to
horizontal level L becomes its level-(L - levels) ancestor, each new cell the
solid-angle weighted mean of the 2^levels x 2^levels cells it covers (rho1,
v_r, v_theta, v_phi, s1, B: means in the patch's own frame; psi set to 0).
The radial direction is kept. For local tests of a fine run's state, e.g.
CZ_3D's 0.47 deg mhd5 (refine_h = 2) at 1.875 deg (levels = 2): import it
with import_dir and refine_h = levels lower. The mass of every radial cell is
kept (rho1 weighted by solid angle); momentum and energy are not exactly
(means of v, s1). The meta file (time, step, clocks) is copied.

Assumes the Yin-Yang patch theta grid is uniform over [pi/4, 3pi/4] (yy_margin
= 0, as all our runs) and all blocks at the same horizontal level.

Usage: coarsen_export.py in_dir out_dir [--levels 2] [--files 8]
"""
import argparse, glob, os, shutil
import numpy as np

KEY_LEN = 7          # ARBOR key: patch, level r/theta/phi, index r/theta/phi
I_PSI = 8            # psi (0-based variable index with the field)


def read_export(path):
    with open(path, 'rb') as f:
        n, nb1, nb2, nb3, nv = np.frombuffer(f.read(20), dtype='<i4')
        rec = np.dtype([('key', '<i4', KEY_LEN), ('s', '<f8', nb1 * nb2 * nb3 * nv)])
        d = np.frombuffer(f.read(), dtype=rec, count=n)
    return (nb1, nb2, nb3, nv), d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('in_dir'); ap.add_argument('out_dir')
    ap.add_argument('--levels', type=int, default=2); ap.add_argument('--files', type=int, default=8)
    a = ap.parse_args()
    files = sorted(glob.glob(os.path.join(a.in_dir, 'state_export_r*.bin')))
    if not files:
        raise SystemExit('no state_export files in ' + a.in_dir)
    f2 = 2 ** a.levels
    parents = {}                         # parent key -> (sum w*state, sum w)
    shape = lev = None
    nbt = 0                              # theta blocks per patch at the fine level
    for fn in files:
        sh, d = read_export(fn)
        if shape is None:
            shape = sh
        elif sh != shape:
            raise SystemExit('the export files differ in block size or variables')
        nb1, nb2, nb3, nv = shape
        for key, s in zip(d['key'], d['s']):
            if lev is None:
                lev = int(key[2])
                if key[2] != key[3] or lev < a.levels:
                    raise SystemExit('blocks at theta/phi levels %d/%d: need equal levels >= %d' % (key[2], key[3], a.levels))
            if key[2] != lev or key[3] != lev:
                raise SystemExit('blocks at different horizontal levels')
            st = s.reshape((nb1, nb2, nb3, nv), order='F')
            jt, kp = int(key[5]), int(key[6])
            pk = (int(key[0]), int(key[1]), lev - a.levels, lev - a.levels, int(key[4]), jt // f2, kp // f2)
            if pk not in parents:
                parents[pk] = [np.zeros((nb1, nb2 * f2, nb3 * f2, nv)), 0]
            P = parents[pk]
            j0, k0 = (jt % f2) * nb2, (kp % f2) * nb3
            P[0][:, j0:j0 + nb2, k0:k0 + nb3, :] = st
            P[1] += 1
            nbt = max(nbt, jt + 1)
        print('read', os.path.basename(fn), len(d), 'blocks')
    nb1, nb2, nb3, nv = shape
    # the parents: means over f2 x f2 cells (theta weights; phi uniform)
    out = []
    ncell = nbt * nb2                    # theta cells per patch at the fine level (uniform over [pi/4, 3pi/4])
    dth = 0.5 * np.pi / ncell
    for pk, (S, nchild) in sorted(parents.items()):
        if nchild != f2 * f2:
            raise SystemExit('parent %s: %d of its %d children in the export' % (pk, nchild, f2 * f2))
        lo = 0.25 * np.pi + (pk[5] * f2 * nb2 + np.arange(nb2 * f2)) * dth
        w = (np.cos(lo) - np.cos(lo + dth)).reshape(nb2, f2)
        num = (S.reshape(nb1, nb2, f2, nb3, f2, nv) * w[None, :, :, None, None, None]).sum(axis=(2, 4))
        C = num / (w.sum(axis=1)[None, :, None, None] * f2)
        if nv > I_PSI:
            C[..., I_PSI] = 0.0
        out.append((pk, C))
    os.makedirs(a.out_dir, exist_ok=True)
    rec = np.dtype([('key', '<i4', KEY_LEN), ('s', '<f8', nb1 * nb2 * nb3 * nv)])
    chunks = np.array_split(np.arange(len(out)), a.files)
    for r, idx in enumerate(chunks):
        arr = np.zeros(len(idx), dtype=rec)
        for m, i in enumerate(idx):
            arr['key'][m] = out[i][0]
            arr['s'][m] = out[i][1].reshape(-1, order='F')
        with open(os.path.join(a.out_dir, 'state_export_r%04d.bin' % r), 'wb') as f:
            f.write(np.array([len(idx), nb1, nb2, nb3, nv], dtype='<i4').tobytes())
            f.write(arr.tobytes())
    meta = os.path.join(a.in_dir, 'state_export_meta.bin')
    if os.path.exists(meta):
        shutil.copy(meta, a.out_dir)
    print('%d blocks at level %d (%d theta cells per patch) -> %d blocks at level %d in %d files (%s)'
          % (len(parents) * f2 * f2, lev, ncell, len(out), lev - a.levels, len(chunks), a.out_dir))


if __name__ == '__main__':
    main()
