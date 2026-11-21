#!/usr/bin/env python3
"""Plot CZTR_3D synthetic continuum images (see phot_write).

out_dir/phot_NNNNNNNN.bin: monochromatic intensities [erg s^-1 cm^-2 sr^-1 A^-1]
(channels = wavelengths); out_dir/photb_NNNNNNNN.bin: AIA bands (phot_band_file)
in DN/s/pix (channels = 1600, 1700, 4500). The format is the EUV one
(tools/euv_plot.py read).

With --euv euv_NNNNNNNN.bin (made with euv_table data/aia_response_uv.dat)
the optically thin TR lines of its 1600/1700 channels (C IV ...) are added
to the band continua, and the 1600 lines alone get a panel.

One PNG per view: each channel normalised to the disk centre (the median
inside r < 0.1), the limb dotted, the limb darkening I(mu)/I(1) of each
channel (annuli), and with --zoom z the disk centre +-z R_sun of the first.

Usage: phot_plot.py photb_00012345.bin [--euv euv_00012345.bin] [--out prefix]
       [--range 0 1.3] [--zoom 0.3]
"""
import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from euv_plot import read  # noqa: E402


def cmap(ch):
    try:
        import sunpy.visualization.colormaps  # noqa: F401
        if ch in (1600, 1700, 4500):
            return plt.get_cmap(f'sdoaia{ch}')
    except Exception:
        pass
    return plt.get_cmap('gray')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('file'); ap.add_argument('--out', default=None)
    ap.add_argument('--euv', default=None)
    ap.add_argument('--range', type=float, nargs=2, default=None)
    ap.add_argument('--zoom', type=float, default=0.0)
    a = ap.parse_args()
    d = read(a.file)
    band = os.path.basename(a.file).startswith('photb_')
    unit = 'DN/s/pix' if band else 'erg/s/cm2/sr/A'
    prefix = a.out or a.file.rsplit('.', 1)[0]
    f = d['fov']; npix = d['img'].shape[-1]
    x = (np.arange(npix) + 0.5)/npix*2*f - f
    X, Y = np.meshgrid(x, x)
    R = np.hypot(X, Y)
    th = np.linspace(0, 2*np.pi, 400)
    chs = list(d['ch'])
    lines = {}
    if a.euv:
        e = read(a.euv)
        if e['img'].shape[-1] != npix or not np.allclose(e['lon'], d['lon']):
            sys.exit('--euv: different views or pixels')
        for ch in (1600, 1700):
            if ch in chs and ch in list(e['ch']):
                lines[ch] = e['img'][:, list(e['ch']).index(ch)]
    for v in range(len(d['lon'])):
        panels = []
        for m, ch in enumerate(chs):
            cont = d['img'][v, m].astype(float)
            tot = cont + (lines[ch][v] if ch in lines else 0.0)
            panels.append((f'{ch}' + (' cont + TR lines' if ch in lines else ''), ch, tot))
        if 1600 in lines:
            panels.append(('1600 TR lines only', 1600, lines[1600][v].astype(float)))
        nz = 1 if a.zoom > 0 else 0
        ncol = len(panels) + 1 + nz
        fig, axs = plt.subplots(1, ncol, figsize=(4.4*ncol, 4.2), squeeze=False)
        axs = axs[0]
        for k, (title, ch, im) in enumerate(panels):
            ax = axs[k]
            ic = np.median(im[R < 0.1])
            if 'only' in title:
                pos = im[im > 0]
                lo, hi = (np.percentile(pos, 1), np.percentile(pos, 99.5)) if pos.size else (0, 1)
                p = ax.imshow(im, origin='lower', extent=[-f, f, -f, f], cmap=cmap(ch), vmin=0, vmax=hi)
                ax.set_title(f'{title}; disk median {np.median(im[R < 0.95]):.3g} {unit}', fontsize=8)
            else:
                lo, hi = a.range if a.range else (0.0, 1.3)
                p = ax.imshow(im/ic, origin='lower', extent=[-f, f, -f, f], cmap=cmap(ch), vmin=lo, vmax=hi)
                ax.set_title(f'{title}, I/I(centre); I(centre) = {ic:.3g} {unit}', fontsize=8)
            ax.plot(np.cos(th), np.sin(th), 'w:', lw=0.5)
            ax.set_xticks([]); ax.set_yticks([])
            plt.colorbar(p, ax=ax, shrink=0.8)
        # limb darkening (annuli well inside the limb)
        ax = axs[len(panels)]
        edges = np.linspace(0, 0.995, 80)
        rc = 0.5*(edges[1:] + edges[:-1]); mu = np.sqrt(1 - rc**2)
        for title, ch, im in panels:
            if 'only' in title:
                continue
            ic = np.median(im[R < 0.1])
            prof = np.array([np.mean(im[(R >= e0) & (R < e1)]) for e0, e1 in zip(edges[:-1], edges[1:])])/ic
            ax.plot(mu, prof, label=title)
            print(f'view {v + 1} {title}: I(centre) {ic:.4g} {unit}, I(mu=0.5)/I(1) = '
                  f'{np.interp(0.5, mu[::-1], prof[::-1]):.3f}, I(mu=0.2)/I(1) = {np.interp(0.2, mu[::-1], prof[::-1]):.3f}')
        ax.set_xlabel('mu'); ax.set_ylabel('I(mu)/I(1)'); ax.set_xlim(0, 1); ax.set_ylim(0, 1.3)
        ax.grid(alpha=0.3); ax.legend(fontsize=7)
        if nz:
            z = a.zoom; title, ch, im = panels[0]
            sel = np.abs(x) < z
            sub = (im/np.median(im[R < 0.1]))[np.ix_(sel, sel)]
            ax2 = axs[-1]
            q = ax2.imshow(sub, origin='lower', extent=[-z, z, -z, z], cmap=cmap(ch),
                           vmin=np.percentile(sub, 0.5), vmax=np.percentile(sub, 99.5))
            ax2.set_title(f'{title}, centre +-{z} R_sun (rms {sub.std()/sub.mean()*100:.2f}%)', fontsize=8)
            plt.colorbar(q, ax=ax2, shrink=0.8)
        fig.suptitle(f"t = {d['t']:.1f} s, view lon {d['lon'][v]:.1f}, lat {d['lat'][v]:.1f} deg", fontsize=10)
        fig.tight_layout()
        out = f'{prefix}_v{v + 1}.png'
        fig.savefig(out, dpi=110); plt.close(fig)
        print('wrote', out)


if __name__ == '__main__':
    main()
