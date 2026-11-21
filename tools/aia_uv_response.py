#!/usr/bin/env python3
"""AIA UV/visible band data for CZTR_3D's synthetic continuum and TR lines.

1. The continuum band file (phot_band_file, read by phot_write): for each
   band (AIA 1600, 1700, 4500) a few wavelengths lambda_k [A] with weights
   w_k [DN s^-1 pix^-1 per erg s^-1 cm^-2 sr^-1 A^-1] so that the image
   DN/s/pix = sum_k w_k I(lambda_k):
       w_k = Omega_pix * int_{bin k} R(lambda) lambda/(h c) dlambda,
   R aiapy's wavelength response (cm^2 DN ph^-1, nominal), bins of equal
   width over where R > 1% of its peak (4500: 4300-4800 A only; the UV leak
   meets a solar UV intensity ~1e-3 of the visible).
   Format: '#' comments, then n_band, then per band a line "id n" and n lines
   "lambda w".
2. Optically thin TR lines in the 1600 and 1700 bands (C IV 1548/1551 etc.):
   K(T) [DN cm^5 s^-1 pix^-1] per n_e n_H as tools/aia_response.py (CHIANTI
   via fiasco, Feldman92 coronal abundances, n_e T = 1e15), from the ions
   listed in UV_IONS, lines 1300-2000 A; appended as two columns to an EUV
   table (the same log T rows): the code adds them as channels 1600/1700.
   The chromospheric lines (optically thick, non-LTE: Fe II, C I, ...) are
   not included.

Usage (the 'euv' conda env):
  python tools/aia_uv_response.py data/aia_uv_bands.dat \
      --euv-in data/aia_response.dat --euv-out data/aia_response_uv.dat
"""
import argparse, sys
import numpy as np
import astropy.units as u
import astropy.constants as const
from aiapy.response import Channel

BANDS = {1600: (24, None), 1700: (24, None), 4500: (6, (4300.0, 4800.0))}
UV_IONS = ['c_4', 'si_4', 'he_2', 'si_2', 'o_4', 'o_3', 'c_3', 'n_5', 'n_4', 'si_3', 'o_5', 'c_2']
WMIN, WMAX = 1300.0, 2000.0
PRESSURE = 1.0e15
HC = (const.h*const.c).to(u.erg*u.angstrom).value


def response(ch):
    c = Channel(ch*u.angstrom)
    R = c.wavelength_response(correction_table=None)
    return c.wavelength.to_value(u.angstrom), R.to_value(u.cm**2*u.DN/u.ph), c.plate_scale.to_value(u.sr/u.pix)


def bands(out):
    lines = ['# AIA continuum bands for phot_write: DN/s/pix = sum_k w_k I(lambda_k), I in erg s^-1 cm^-2 sr^-1 A^-1',
             '# w_k = Omega_pix int_bin R lambda/(h c) dlambda (aiapy nominal response); made by tools/aia_uv_response.py',
             f'{len(BANDS)}']
    for ch, (n, rng) in BANDS.items():
        wl, R, om = response(ch)
        if rng is None:
            s = np.where(R > 0.01*R.max())[0]
            rng = (wl[s[0]], wl[s[-1]])
        edges = np.linspace(rng[0], rng[1], n + 1)
        fine = np.linspace(rng[0], rng[1], 20*n + 1)
        g = om*np.interp(fine, wl, R)*fine/HC
        lam, w = [], []
        for k in range(n):
            m = (fine >= edges[k]) & (fine <= edges[k+1])
            w.append(np.trapezoid(g[m], fine[m])); lam.append(0.5*(edges[k] + edges[k+1]))
        tot = om*np.trapezoid(np.interp(wl, wl, R)*wl/HC, wl)
        lines.append(f'{ch} {n}')
        lines += [f'{l:9.2f} {x:12.5e}' for l, x in zip(lam, w)]
        print(f'{ch}: {rng[0]:.0f}-{rng[1]:.0f} A, sum w {sum(w):.4e} (whole response {tot:.4e})', file=sys.stderr)
    open(out, 'w').write('\n'.join(lines) + '\n')


def tr_lines(euv_in, euv_out):
    import fiasco
    head, rows = [], []
    with open(euv_in) as f:
        for line in f:
            if line.startswith('#'):
                head.append(line.rstrip('\n'))
            else:
                rows.append(line.split())
    nT, nch = int(rows[0][0]), int(rows[0][1])
    chans = rows[1]
    data = np.array(rows[2:2+nT], float)
    LGT = data[:, 0]
    T = 10**LGT*u.K
    ne = PRESSURE/T.value*u.cm**-3
    K = np.zeros((nT, 2))
    used = []
    resp = [response(1600), response(1700)]
    for name in UV_IONS:
        ion = fiasco.Ion(name, T)
        wave = ion.transitions.wavelength[ion.transitions.is_bound_bound].to_value(u.angstrom)
        sel = (wave >= WMIN) & (wave <= WMAX)
        if not sel.any():
            continue
        try:
            g = ion.contribution_function(ne, couple_density_to_temperature=True)
        except Exception as e:
            print(f'  skip {name}: {type(e).__name__}', file=sys.stderr)
            continue
        g = g.to_value(u.erg*u.cm**3/u.s)[:, 0, :][:, sel]
        w = wave[sel]
        for j in range(2):
            wl, R, om = resp[j]
            K[:, j] += om*(g*(np.interp(w, wl, R, left=0.0, right=0.0)*w/(4*np.pi*HC))[None, :]).sum(axis=1)
        used.append(name)
        print(f'  {name}: {sel.sum()} lines', file=sys.stderr)
    with open(euv_out, 'w') as f:
        for h in head:
            f.write(h + '\n')
        f.write(f'# + channels 1600, 1700: optically thin TR lines {WMIN:.0f}-{WMAX:.0f} A of {", ".join(used)} '
                f'(tools/aia_uv_response.py; no continuum, no chromospheric lines)\n')
        f.write(f'{nT} {nch + 2}\n')
        f.write(' '.join(chans) + ' 1600 1700\n')
        for i in range(nT):
            f.write(f'{LGT[i]:6.3f} ' + ' '.join(f'{k:12.5e}' for k in data[i, 1:]) +
                    f' {K[i, 0]:12.5e} {K[i, 1]:12.5e}\n')
    for j, ch in enumerate((1600, 1700)):
        i = np.argmax(K[:, j])
        print(f'  {ch}: peak {K[i, j]:.3e} at log T {LGT[i]:.2f}', file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('bands_out')
    ap.add_argument('--euv-in'); ap.add_argument('--euv-out')
    a = ap.parse_args()
    bands(a.bands_out)
    if a.euv_in and a.euv_out:
        tr_lines(a.euv_in, a.euv_out)


if __name__ == '__main__':
    main()
