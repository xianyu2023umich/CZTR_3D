#!/usr/bin/env python3
"""Figures for a CZTR_3D run directory (e.g. the casc047_l1 cascade level).

  python3 tools/casc_plots.py <run_dir> [--avg-days 45]

Reads log.txt, profile_NNNN.dat, rotation_NNNN.dat and the shell_dump_rNNNN.dat
files (the final state) and writes PNGs to <run_dir>/figs:
  1_timeseries.png   E_mag/E_kin, max|B|, dt, Mach, L_total range, Lorentz work
  2_profiles.png     rms v, B_rms (B_r, B_h), B_eq, B/B_eq, luminosity budget
  3_B_over_Beq.png   B_rms/B_eq in (time, depth)
  4_rotation.png     Omega - Omega0 and meridional flow in the (r, lat) plane
  5_shells.png       v_r, B_r, B_phi maps on the dumped shells (final state)
  6_bphi_mean.png    longitude-mean B_phi, B_r vs latitude on the dumped shells
"""
import argparse, glob, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RSUN = 6.957e10
DAY = 86400.0


def read_t(f):
    return float(open(f).readline().split('=')[1].split('s')[0])


def load_series(pattern):
    fs = sorted(glob.glob(pattern))
    t = np.array([read_t(f) for f in fs]) / DAY
    return t, [np.loadtxt(f) for f in fs]


def fig_timeseries(run, out):
    a = np.loadtxt(os.path.join(run, 'log.txt'), comments='#')
    a = a[a[:, 0] > 0]
    d = a[:, 1] / DAY
    fig, ax = plt.subplots(3, 2, figsize=(12, 9), sharex=True)
    ax = ax.ravel()
    ax[0].semilogy(d, a[:, 21]); ax[0].set_ylabel('E_mag / E_kin')
    ax[1].semilogy(d, a[:, 20], label='max|B|'); ax[1].semilogy(d, a[:, 27] / 1e2, label='max v_A [m/s]')
    ax[1].legend(); ax[1].set_ylabel('max|B| [G]')
    ax[2].plot(d, a[:, 2]); ax[2].set_ylabel('dt [s]'); ax[2].set_ylim(0, 1.1 * a[:, 2].max())
    ax[3].plot(d, a[:, 12], label='max reduced Mach'); ax[3].plot(d, 10 * a[:, 3], label='10 x max Mach')
    ax[3].legend(); ax[3].set_ylabel('Mach')
    ax[4].plot(d, a[:, 14], label='min L_total'); ax[4].plot(d, a[:, 15], label='max L_total')
    ax[4].plot(d, a[:, 6], 'k--', lw=0.8, label='L_top'); ax[4].set_ylim(0.5, 1.5)
    ax[4].legend(); ax[4].set_ylabel('L / Lsun')
    ax[5].plot(d, -a[:, 23], label='-W_Lorentz'); ax[5].plot(d, a[:, 22] * 100, label='100 x rms(divB dx)/rms|B|')
    ax[5].legend(); ax[5].set_ylabel('[Lsun] / [-]')
    for x in ax[4:]:
        x.set_xlabel('t [days]')
    fig.suptitle(os.path.basename(os.path.normpath(run)) + ': time series')
    fig.tight_layout(); fig.savefig(os.path.join(out, '1_timeseries.png'), dpi=130); plt.close(fig)


def bvals(p):
    rho, vr, vh, br, bh = p[:, 3], p[:, 8], p[:, 9], p[:, -2], p[:, -1]
    v = np.hypot(vr, vh)
    return v, np.hypot(br, bh), v * np.sqrt(4 * np.pi * rho)


def fig_profiles(run, out, avg):
    t, P = load_series(os.path.join(run, 'profile_*.dat'))
    P = np.array(P)
    m = P[t > t[-1] - avg].mean(0)
    z = m[:, 2] / 1e3
    v, B, Beq = bvals(m)
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    ax[0].plot(z, m[:, 8] / 100, label='rms v_r'); ax[0].plot(z, m[:, 9] / 100, label='rms v_h')
    ax[0].set_ylabel('[m/s]'); ax[0].legend(); ax[0].set_title('velocity')
    ax[1].semilogy(z, B, 'k', label='B_rms'); ax[1].semilogy(z, m[:, -2], label='rms B_r')
    ax[1].semilogy(z, m[:, -1], label='rms B_h'); ax[1].semilogy(z, Beq, 'k--', label='B_eq')
    a2 = ax[1].twinx(); a2.plot(z, B / Beq, 'r:', label='B/B_eq'); a2.set_ylabel('B_rms / B_eq', color='r')
    a2.set_ylim(0, max(0.3, 1.2 * (B / Beq).max()))
    ax[1].set_ylabel('[G]'); ax[1].legend(loc='lower left'); ax[1].set_title('field')
    for c, lab in [(10, 'L_rad'), (11, 'L_conv (MLT)'), (12, 'L_enth+L_kin'), (13, 'L_kin'), (14, 'L_diff'),
                   (15, 'L_cool')]:
        if c != 11 or np.abs(m[:, 11]).max() > 1e-3:     # (L_conv only when MLT carries flux)
            ax[2].plot(z, m[:, c], label=lab)
    ax[2].plot(z, m[:, 10] + m[:, 11] + m[:, 12] + m[:, 14] + m[:, 15], 'k', label='sum')
    ax[2].set_ylabel('L / Lsun'); ax[2].legend(); ax[2].set_title('luminosity budget')
    for x in ax:
        x.set_xlabel('z [Mm]')
    fig.suptitle('profiles, mean over days %.0f-%.0f' % (t[-1] - avg, t[-1]))
    fig.tight_layout(); fig.savefig(os.path.join(out, '2_profiles.png'), dpi=130); plt.close(fig)

    # B/B_eq(t, z)
    R = []
    for p in P:
        v, B, Beq = bvals(p)
        R.append(np.where(Beq > 0, B / np.maximum(Beq, 1e-30), np.nan))
    R = np.array(R)
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    pc = ax[0].pcolormesh(t, z, R.T, shading='auto', vmin=0, vmax=0.2, cmap='viridis')
    fig.colorbar(pc, ax=ax[0], label='B_rms / B_eq'); ax[0].set_xlabel('t [days]'); ax[0].set_ylabel('z [Mm]')
    for zz in (-190, -120, -60, -35):
        i = np.argmin(abs(z - zz)); ax[1].plot(t, R[:, i], label='z = %d Mm' % zz)
    ax[1].set_xlabel('t [days]'); ax[1].set_ylabel('B_rms / B_eq'); ax[1].set_ylim(0, 0.2); ax[1].legend()
    fig.tight_layout(); fig.savefig(os.path.join(out, '3_B_over_Beq.png'), dpi=130); plt.close(fig)


def fig_rotation(run, out, avg):
    t, Rr = load_series(os.path.join(run, 'rotation_*.dat'))
    a = np.mean([r for tt, r in zip(t, Rr) if tt > t[-1] - avg], 0)
    ii = np.unique(a[:, 0]).astype(int)
    lat = np.unique(a[:, 3])
    r = np.array([a[a[:, 0] == i][0, 1] for i in ii]) / RSUN
    om = np.array([a[a[:, 0] == i][:, 5] for i in ii])
    vt = np.array([a[a[:, 0] == i][:, 6] for i in ii]) / 100
    L, Rg = np.meshgrid(np.radians(lat), r)
    X, Y = Rg * np.cos(L), Rg * np.sin(L)
    fig, ax = plt.subplots(1, 3, figsize=(15, 6), gridspec_kw={'width_ratios': [1, 1, 1.3]})
    lim = np.nanpercentile(abs(om), 98)
    pc = ax[0].pcolormesh(X, Y, om, shading='gouraud', cmap='RdBu_r', vmin=-lim, vmax=lim)
    fig.colorbar(pc, ax=ax[0], label='Omega - Omega0 [nHz]'); ax[0].set_title('differential rotation')
    lim = np.nanpercentile(abs(vt), 98)
    pc = ax[1].pcolormesh(X, Y, vt, shading='gouraud', cmap='PuOr_r', vmin=-lim, vmax=lim)
    fig.colorbar(pc, ax=ax[1], label='v_theta [m/s] (>0 southward)'); ax[1].set_title('meridional flow')
    for x in ax[:2]:
        x.set_aspect('equal'); x.set_xlabel('r cos(lat) / Rsun'); x.set_ylabel('r sin(lat) / Rsun')
    for k in np.linspace(0, len(r) - 1, 6).astype(int):
        ax[2].plot(lat, om[k], label='r = %.3f Rsun' % r[k])
    ax[2].axhline(0, color='k', lw=0.5); ax[2].set_xlabel('latitude [deg]'); ax[2].set_ylabel('Omega - Omega0 [nHz]')
    ax[2].legend(fontsize=8); ax[2].set_title('Omega0 = 414 nHz')
    fig.suptitle('rotation, mean over days %.0f-%.0f' % (t[-1] - avg, t[-1]))
    fig.tight_layout(); fig.savefig(os.path.join(out, '4_rotation.png'), dpi=130); plt.close(fig)


def load_shells(run):
    fs = sorted(glob.glob(os.path.join(run, 'shell_dump_r*.dat')))
    a = np.concatenate([np.loadtxt(f, ndmin=2) for f in fs if os.path.getsize(f) > 200])
    a = a[a[:, 5] > 0]
    out = {}
    for i in np.unique(a[:, 0]).astype(int):
        s = a[a[:, 0] == i]
        x, y, zc = s[:, 2], s[:, 3], s[:, 4]
        lat = np.degrees(np.arcsin(np.clip(zc, -1, 1))); lon = np.degrees(np.arctan2(y, x))
        rho_c = np.hypot(x, y)
        ephi = np.stack([-y, x, 0 * x], 1) / rho_c[:, None]
        er = np.stack([x, y, zc], 1)
        Bv = s[:, 12:15]
        out[i] = dict(lat=lat, lon=lon, w=s[:, 5], vr=s[:, 6] / 100,
                      br=(Bv * er).sum(1), bphi=(Bv * ephi).sum(1))
    return out


def fig_shells(run, out, bg_r):
    from scipy.interpolate import griddata
    S = load_shells(run)
    lon_g, lat_g = np.meshgrid(np.linspace(-180, 180, 361), np.linspace(-89.5, 89.5, 180))
    keys = sorted(S)
    fig, ax = plt.subplots(len(keys), 3, figsize=(16, 2.9 * len(keys)))
    ax = np.atleast_2d(ax)
    for row, i in enumerate(keys):
        s = S[i]
        pts = np.stack([s['lon'], s['lat']], 1)
        zlab = '' if bg_r is None else ' (z = %.0f Mm)' % bg_r.get(i, np.nan)
        for col, (key, lab, cm) in enumerate([('vr', 'v_r [m/s]', 'RdBu_r'), ('br', 'B_r [G]', 'RdBu_r'),
                                              ('bphi', 'B_phi [G]', 'RdBu_r')]):
            g = griddata(pts, s[key], (lon_g, lat_g), method='nearest')
            lim = np.percentile(abs(s[key]), 99)
            pc = ax[row, col].pcolormesh(lon_g, lat_g, g, cmap=cm, vmin=-lim, vmax=lim, shading='auto')
            fig.colorbar(pc, ax=ax[row, col], label=lab)
            ax[row, col].set_title('cell %d%s: %s' % (i, zlab, lab.split()[0]))
            ax[row, col].set_aspect('equal')
    fig.suptitle('final state, longitude-latitude maps (global frame)')
    fig.tight_layout(); fig.savefig(os.path.join(out, '5_shells.png'), dpi=110); plt.close(fig)

    # longitude means
    edges = np.linspace(-90, 90, 37); c = 0.5 * (edges[1:] + edges[:-1])
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for i in keys:
        s = S[i]; k = np.digitize(s['lat'], edges) - 1
        lab = 'cell %d' % i + ('' if bg_r is None else ' (%.0f Mm)' % bg_r.get(i, np.nan))
        for j, key in enumerate(('bphi', 'br')):
            mean = np.array([np.sum(s['w'][k == n] * s[key][k == n]) / max(np.sum(s['w'][k == n]), 1e-30)
                             for n in range(len(c))])
            ax[j].plot(c, mean, label=lab)
    for j, key in enumerate(('<B_phi>', '<B_r>')):
        ax[j].axhline(0, color='k', lw=0.5); ax[j].set_xlabel('latitude [deg]'); ax[j].set_ylabel(key + ' [G]')
        ax[j].legend(fontsize=8)
    fig.suptitle('longitude-mean field, final state (rms B_h ~ 1-1.5 kG for comparison)')
    fig.tight_layout(); fig.savefig(os.path.join(out, '6_bphi_mean.png'), dpi=130); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run'); ap.add_argument('--avg-days', type=float, default=45.0)
    a = ap.parse_args()
    out = os.path.join(a.run, 'figs'); os.makedirs(out, exist_ok=True)
    fig_timeseries(a.run, out)
    fig_profiles(a.run, out, a.avg_days)
    fig_rotation(a.run, out, a.avg_days)
    bg_r = None
    pf = sorted(glob.glob(os.path.join(a.run, 'profile_*.dat')))
    if pf:
        p = np.loadtxt(pf[-1]); bg_r = {int(i): zz / 1e3 for i, zz in zip(p[:, 0], p[:, 2])}
    if glob.glob(os.path.join(a.run, 'shell_dump_r*.dat')):
        fig_shells(a.run, out, bg_r)
    print('figures in', out)


if __name__ == '__main__':
    main()
