#!/usr/bin/env python3
"""AIA temperature responses K_c(T) for CZTR_3D's synthetic EUV (euv_table).

K_c(T) [DN cm^5 s^-1 pix^-1] = Omega_pix * int R_c(lambda) I(lambda, T) lambda/(h c) dlambda
per unit column emission measure EM = int n_e n_H dl, so the code's image is
sum n_e n_H K_c(T) dV / A_pixel. I = lines (CHIANTI G(T)/4pi, all ions with
transitions in 50-400 A) + free-free + free-bound continuum (per sr).
R_c: aiapy's wavelength response (cm^2 DN ph^-1; nominal, no degradation,
crosstalk included); Omega_pix from the channel's plate scale.
Atomic data: CHIANTI 10.1 via fiasco, coronal abundances (Feldman 1992,
extended), CHIANTI ionization equilibrium, constant pressure n_e T = 1e15 cm^-3 K.

Usage (the 'euv' conda env: numpy, aiapy, fiasco + its CHIANTI database):
  python tools/aia_response.py data/aia_response.dat [--plot out.png]
"""
import argparse, sys, time
import numpy as np
import astropy.units as u
import astropy.constants as const
import fiasco
import aiapy
from aiapy.response import Channel

CHANNELS = [94, 131, 171, 193, 211, 335]
LGT = np.round(np.arange(4.0, 8.0001, 0.05), 4)
PRESSURE = 1.0e15            # n_e T [cm^-3 K]
WMIN, WMAX = 50.0, 400.0     # [A] the lines and continuum included


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--plot', default=None)
    a = ap.parse_args()
    T = 10**LGT*u.K
    ne = PRESSURE/T.value*u.cm**-3
    # the channels' wavelength responses and pixel solid angles
    resp = {}
    for ch in CHANNELS:
        c = Channel(ch*u.angstrom)
        R = c.wavelength_response(correction_table=None)     # (obstime None: no degradation, the table unused)
        omega = c.plate_scale.to_value(u.sr/u.pix)
        resp[ch] = (c.wavelength.to_value(u.angstrom), R.to_value(u.cm**2*u.DN/u.ph), omega)
    hc = (const.h*const.c).to(u.erg*u.angstrom).value
    K = np.zeros((len(LGT), len(CHANNELS)))
    # lines: sum over transitions of R(lambda) G(T) lambda/(4 pi h c)
    t0 = time.time()
    names = fiasco.list_ions()
    nion = 0
    for name in names:
        try:
            ion = fiasco.Ion(name, T)
            wave = ion.transitions.wavelength[ion.transitions.is_bound_bound].to_value(u.angstrom)
        except Exception:
            continue
        sel = (wave >= WMIN) & (wave <= WMAX)
        if not sel.any():
            continue
        try:
            g = ion.contribution_function(ne, couple_density_to_temperature=True)
        except Exception as e:
            print(f'  skip {name}: {type(e).__name__}', file=sys.stderr)
            continue
        g = g.to_value(u.erg*u.cm**3/u.s)[:, 0, :][:, sel]            # (nT, k)
        w = wave[sel]
        for j, ch in enumerate(CHANNELS):
            wl, R, om = resp[ch]
            Rk = np.interp(w, wl, R, left=0.0, right=0.0)
            K[:, j] += om*(g*(Rk*w/(4.0*np.pi*hc))[None, :]).sum(axis=1)
        nion += 1
        print(f'  {name}: {sel.sum()} lines ({time.time() - t0:.0f} s)', file=sys.stderr)
    # continuum (per sr; abundances and ionization fractions included)
    lam = np.arange(WMIN, WMAX + 1e-9, 0.1)*u.angstrom
    coll = fiasco.IonCollection(*[fiasco.Ion(n, T) for n in names
                                  if fiasco.Ion(n, T).charge_state > 0])
    cont = np.zeros((len(LGT), lam.size))
    for kind in ('free_free', 'free_bound'):
        try:
            cont += getattr(coll, kind)(lam).to_value(u.erg*u.cm**3/u.s/u.angstrom)
        except Exception as e:
            print(f'  continuum {kind} failed: {e}', file=sys.stderr)
    Kc = np.zeros_like(K)
    for j, ch in enumerate(CHANNELS):
        wl, R, om = resp[ch]
        Rl = np.interp(lam.value, wl, R, left=0.0, right=0.0)
        Kc[:, j] = om*np.trapezoid(cont*(Rl*lam.value/hc)[None, :], lam.value, axis=1)
    K += Kc
    with open(a.out, 'w') as f:
        f.write('# AIA temperature responses K(T) [DN cm^5 s^-1 pix^-1] per unit EM = int n_e n_H dl\n')
        f.write(f'# CHIANTI {fiasco.Ion("Fe 9", T).version if hasattr(fiasco.Ion("Fe 9", T), "version") else "10.1"} via fiasco {fiasco.__version__}; '
                f'aiapy {aiapy.__version__} (nominal response, no degradation, crosstalk on)\n')
        f.write('# abundances sun_coronal_1992_feldman_ext, CHIANTI ionization equilibrium, n_e T = 1e15 cm^-3 K\n')
        f.write(f'# lines {WMIN}-{WMAX} A from {nion} ions + free-free + free-bound; made by tools/aia_response.py\n')
        i6, i7 = int(np.argmin(abs(LGT - 6.0))), int(np.argmin(abs(LGT - 7.0)))
        f.write(f'# continuum share at log T = 6, 7: ' + ', '.join(
            f'{ch}: {Kc[i6, j]/max(K[i6, j], 1e-99):.2f} {Kc[i7, j]/max(K[i7, j], 1e-99):.2f}' for j, ch in enumerate(CHANNELS)) + '\n')
        f.write(f'{len(LGT)} {len(CHANNELS)}\n')
        f.write(' '.join(str(c) for c in CHANNELS) + '\n')
        for i, lt in enumerate(LGT):
            f.write(f'{lt:6.3f} ' + ' '.join(f'{k:12.5e}' for k in K[i]) + '\n')
    print(f'wrote {a.out}: {nion} ions, {time.time() - t0:.0f} s', file=sys.stderr)
    for j, ch in enumerate(CHANNELS):
        i = np.argmax(K[:, j])
        print(f'  {ch:4d} A: peak {K[i, j]:.3e} at log T {LGT[i]:.2f}', file=sys.stderr)
    if a.plot:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, axx = plt.subplots(figsize=(7, 5))
        for j, ch in enumerate(CHANNELS):
            axx.semilogy(LGT, K[:, j], label=f'{ch} A')
        axx.set_ylim(1e-29, 1e-23); axx.set_xlabel('log10 T [K]'); axx.set_ylabel('K [DN cm^5 s^-1 pix^-1]')
        axx.legend(); axx.grid(alpha=0.3); fig.tight_layout(); fig.savefig(a.plot, dpi=110)


if __name__ == '__main__':
    main()
