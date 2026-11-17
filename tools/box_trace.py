#!/usr/bin/env python3
"""Field lines through a gridded box (tools/box_grid.py .npz): traced from a
start point both ways along B until they reach z_stop (a footpoint) or leave
the box, with the state sampled along them.

trace(f, n, lat, lon, z_km) -> dict of arrays along the line, ordered from
one end to the other: s [Mm] (arc length), z [km], lat, lon, rho, T, p
(rho k T/(mu m_H), mu from T: 1.25 below 8e3 K, 0.6 above 2e4 K, linear in
between: a rough pressure, not the EOS one), v_par [km/s] (v along B, positive
along the direction of increasing s), Mach, beta, B [G], Qh, Ql, and 'end'
('foot' or 'out') for each end.

Usage as a script:
  box_trace.py box.npz snapshot lat lon z_km line.npz      (needs scipy)
  box_trace.py --plot line.npz out.png ["title"]           (needs matplotlib)
"""
import sys
import numpy as np

RSUN_KM = 6.957e5
KB, MH = 1.380649e-16, 1.6726e-24


def _mu(T):
    return np.interp(np.log10(T), [np.log10(8e3), np.log10(2e4)], [1.25, 0.6])


class Box:
    def __init__(self, f, n):
        from scipy.interpolate import RegularGridInterpolator
        self.f, self.n = f, n
        z = f['z_km']
        ok = np.isfinite(z)
        self.z = z[ok]
        self.th = np.radians(f['theta_deg'])
        self.ph = np.radians(f['phi_deg'])
        TH, PH = np.meshgrid(self.th, self.ph, indexing='ij')
        st, ct, sp, cp = np.sin(TH), np.cos(TH), np.sin(PH), np.cos(PH)
        er = (st*cp, st*sp, ct)
        et = (ct*cp, ct*sp, -st)
        ep = (-sp, cp, 0.0*sp)
        B = [f[c][n][ok] for c in ('Bx', 'By', 'Bz')]
        V = [f[c][n][ok] for c in ('vx', 'vy', 'vz')]
        def proj(X, e):
            return X[0]*e[0] + X[1]*e[1] + X[2]*e[2]
        self.Bs = [proj(B, e) for e in (er, et, ep)]
        self.Vs = [proj(V, e) for e in (er, et, ep)]
        grid = (self.z, self.th, self.ph)
        mk = lambda a: RegularGridInterpolator(grid, np.nan_to_num(a, nan=0.0), bounds_error=False, fill_value=np.nan)
        self.iB = [mk(b) for b in self.Bs]
        self.iV = [mk(v) for v in self.Vs]
        self.iS = {v: mk(np.log(np.abs(f[v][n][ok]) + 1e-300) if v in ('rho', 'T', 'beta', 'Qh', 'Ql') else f[v][n][ok])
                   for v in ('rho', 'T', 'Mach', 'beta', 'Qh', 'Ql')}
        self.dz = np.gradient(self.z)
        self.dh = (self.th[1] - self.th[0])*RSUN_KM

    def inside(self, x):
        z, th, ph = x
        return (self.z[0] <= z <= self.z[-1] and self.th[0] <= th <= self.th[-1] and self.ph[0] <= ph <= self.ph[-1])

    def dir(self, x):
        b = np.array([i([x])[0] for i in self.iB])
        nb = np.linalg.norm(b)
        if not np.isfinite(nb) or nb == 0.0:
            return None
        r = RSUN_KM + x[0]
        return np.array([b[0], b[1]/r, b[2]/(r*np.sin(x[1]))])/nb   # d(z, th, ph)/ds, s in km


def trace(f, n, lat, lon, z_km, z_stop=0.0, box=None, max_steps=200000):
    bx = box or Box(f, n)
    x0 = np.array([z_km, np.radians(90.0 - lat), np.radians(lon)])
    halves, ends = [], []
    for sgn in (-1.0, 1.0):
        x = x0.copy()
        pts = [x.copy()]
        end = 'out'
        for _ in range(max_steps):
            ds = sgn*0.25*min(np.interp(x[0], bx.z, bx.dz), bx.dh)
            k1 = bx.dir(x)
            if k1 is None: break
            k2 = bx.dir(x + 0.5*ds*k1)
            if k2 is None: break
            k3 = bx.dir(x + 0.5*ds*k2)
            if k3 is None: break
            k4 = bx.dir(x + ds*k3)
            if k4 is None: break
            x = x + ds*(k1 + 2*k2 + 2*k3 + k4)/6.0
            if not bx.inside(x):
                break
            pts.append(x.copy())
            if x[0] <= z_stop:
                end = 'foot'
                break
        halves.append(np.array(pts)); ends.append(end)
    P = np.vstack([halves[0][::-1], halves[1][1:]])
    r = RSUN_KM + P[:, 0]
    dxyz = np.diff(np.stack([r*np.sin(P[:, 1])*np.cos(P[:, 2]), r*np.sin(P[:, 1])*np.sin(P[:, 2]), r*np.cos(P[:, 1])], 1), axis=0)
    s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(dxyz, axis=1))])/1e3
    out = {'s': s, 'z': P[:, 0], 'lat': 90.0 - np.degrees(P[:, 1]), 'lon': np.degrees(P[:, 2]), 'end': (ends[0], ends[1])}
    for v, it in bx.iS.items():
        val = it(P)
        out[v] = np.exp(val) if v in ('rho', 'T', 'beta', 'Qh', 'Ql') else val
    b = np.stack([i(P) for i in bx.iB], 1)
    vv = np.stack([i(P) for i in bx.iV], 1)
    out['B'] = np.linalg.norm(b, axis=1)
    # v along B, signed along increasing s (the line runs along +B on its second half)
    t_hat = np.gradient(np.stack([P[:, 0], r*P[:, 1], r*np.sin(P[:, 1])*P[:, 2]], 1), axis=0)
    bhat = b/np.maximum(out['B'][:, None], 1e-30)
    sgn = np.sign(np.sum(bhat*t_hat, axis=1))
    out['v_par'] = np.sum(vv*bhat, axis=1)*sgn/1e5
    out['p'] = out['rho']*KB*out['T']/(_mu(out['T'])*MH)
    return out


def plot(L, title, png):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(4, 2, figsize=(11, 11), sharex=True)
    s = L['s']
    rows = [('z', 'z [km]', False), ('T', 'T [K]', True), ('rho', 'rho [g/cm3]', True), ('p', 'p [dyn/cm2] (rough mu)', True),
            ('v_par', 'v along s [km/s]', False), ('Mach', 'Mach', False), ('beta', 'beta', True), ('Q', 'Q_heat, Q_loss [erg/cm3/s]', True)]
    for a, (k, lab, lg) in zip(ax.T.ravel(), rows):
        if k == 'Q':
            a.plot(s, L['Qh'], label='heat'); a.plot(s, L['Ql'], label='loss'); a.legend(fontsize=8)
        else:
            a.plot(s, L[k])
        if lg: a.set_yscale('log')
        a.set_ylabel(lab, fontsize=9); a.grid(alpha=0.3); a.tick_params(labelsize=8)
    for a in ax[-1]: a.set_xlabel('s [Mm]', fontsize=9)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(png, dpi=110)


if __name__ == '__main__':
    if sys.argv[1] == '--plot':
        L = dict(np.load(sys.argv[2]))
        plot(L, sys.argv[4] if len(sys.argv) > 4 else sys.argv[2], sys.argv[3])
        sys.exit()
    f = np.load(sys.argv[1])
    n, lat, lon, z = int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
    L = trace(f, n, lat, lon, z)
    print('ends:', L['end'], 'length %.1f Mm' % L['s'][-1])
    print('end points: (%.2f, %.2f, z %.0f) and (%.2f, %.2f, z %.0f)' % (L['lat'][0], L['lon'][0], L['z'][0], L['lat'][-1], L['lon'][-1], L['z'][-1]))
    np.savez(sys.argv[6], **{k: np.asarray(v) for k, v in L.items()})
