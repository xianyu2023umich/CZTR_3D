#!/usr/bin/env python3
"""Plot CZTR_3D synthetic EUV images (out_dir/euv_NNNNNNNN.bin, see euv_write).

One PNG per view: the channels in a grid, log scale, AIA colour maps (sunpy
if installed, else grey), the limb (r = 1) dotted. x: increasing longitude
at the sub-observer point; y: solar north.

Usage: euv_plot.py euv_00012345.bin [--out prefix] [--channels 171 193]
                   [--range 1e-1 1e3] (DN/s/pix; default: per channel from percentiles)
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def read(fn):
    with open(fn, 'rb') as f:
        nv, nch, npix = np.fromfile(f, '<i4', 3)
        t, fov = np.fromfile(f, '<f8', 2)
        lon = np.fromfile(f, '<f8', nv); lat = np.fromfile(f, '<f8', nv)
        ch = np.fromfile(f, '<i4', nch)
        img = np.fromfile(f, '<f4', npix*npix*nch*nv).reshape((nv, nch, npix, npix))   # (view, ch, y, x)
    return dict(t=t, fov=fov, lon=lon, lat=lat, ch=ch, img=img)


def cmap(ch):
    try:
        import sunpy.visualization.colormaps  # noqa: F401 (registers sdoaiaNNN)
        return plt.get_cmap(f'sdoaia{ch}')
    except Exception:
        return plt.get_cmap('gray')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('file'); ap.add_argument('--out', default=None)
    ap.add_argument('--channels', type=int, nargs='*', default=None)
    ap.add_argument('--range', type=float, nargs=2, default=None)
    a = ap.parse_args()
    d = read(a.file)
    chs = list(d['ch']) if not a.channels else a.channels
    prefix = a.out or a.file.rsplit('.', 1)[0]
    f = d['fov']; ext = [-f, f, -f, f]
    th = np.linspace(0, 2*np.pi, 400)
    for v in range(len(d['lon'])):
        n = len(chs); nc = min(3, n); nr = (n + nc - 1)//nc
        fig, axs = plt.subplots(nr, nc, figsize=(4.6*nc, 4.4*nr), squeeze=False)
        for m, ch in enumerate(chs):
            ax = axs[m//nc, m % nc]
            im = d['img'][v, list(d['ch']).index(ch)]
            pos = im[im > 0]
            lo, hi = a.range if a.range else ((np.percentile(pos, 1), np.percentile(pos, 99.9)) if pos.size else (1e-3, 1))
            p = ax.imshow(np.maximum(im, lo), origin='lower', extent=ext, cmap=cmap(ch), norm=LogNorm(lo, hi))
            ax.plot(np.cos(th), np.sin(th), 'w:', lw=0.5)
            ax.set_title(f'{ch} A', fontsize=10); ax.set_xticks([]); ax.set_yticks([])
            plt.colorbar(p, ax=ax, shrink=0.8, label='DN/s/pix')
        for m in range(n, nr*nc):
            axs[m//nc, m % nc].axis('off')
        fig.suptitle(f"t = {d['t']:.1f} s, view lon {d['lon'][v]:.1f}, lat {d['lat'][v]:.1f} deg "
                     f"(x: increasing lon, y: north; fov +-{f:.2f} R_sun)", fontsize=10)
        fig.tight_layout()
        out = f'{prefix}_v{v + 1}.png'
        fig.savefig(out, dpi=110); plt.close(fig)
        print('wrote', out)


if __name__ == '__main__':
    main()
