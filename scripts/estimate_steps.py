#!/usr/bin/env python3
"""Steps to relax a cztr3d.x run from t = 0, estimated before running it, for
a series of local-stepping thresholds dt1.

Usage: /usr/bin/python3 scripts/estimate_steps.py <run.nml> [key=value ...] [dt1=a,b,...] [calibrate]

key=value overrides namelist entries (r_top_rsun, c_target, rsst_q, cfl, ...);
dt1 takes a comma list (0 = time accurate, inf = no cap). The script runs
cztr3d.x for zero steps (1 rank) to get the background and the RSST factor
(background.dat), and takes the convective flux of the 1D model from
prof1d (default: CZTR's step-1 profile behind the 3D initial state).

Model (convection starting from a small random s1, then a front of
convective flux moving down from the top to the base):
  v_mlt(z) = (F_conv/rho0)^(1/3)          F_conv from the 1D model's L_conv
  v_f      = K_F v_mlt                    front speed (in the shell's own time)
  dt_own   = cfl / [(v + c_red)/dr + 2 (v + c_red)/(r dtheta)],  v = K_V v_f,  c_red = cs0/xi
  dt       = min over depth of dt_own     (the global step: the top)
  dt_eff   = max(dt, min(dt_own, dt1))    each shell's step
  N_onset  = C_ON / max_z [v_f dt_eff / max(H_p, r dtheta)]   the seed's growth near the top
             (z within 20 Mm of the top of the convective layers)
             (growth rate ~ v_f / the size of the smallest resolved eddy)
  N_front  = integral from the top to 164 Mm depth of dz / (v_f dt_eff)
             (onset and front: convective cells only, L_conv > 1% Lsun)
  N_relax  = K_R (N_onset + N_front)
The constants are calibrated on deep_c2 (0.99 R, time accurate: onset
~5500 steps, front top -> 164 Mm in ~21.5 days, relaxed at ~36000 steps,
dt ~85 s) and checked on lts_c2 (dt1 = 300 s: onset ~4000, front at the base
~10500, relaxed ~13000); 'calibrate' recomputes them from deep_c2's numbers.
Expect +-30%, worse far from the calibration (other tops, larger alpha).
"""
import os
import re
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LSUN = 3.828e33
Z_ONSET = 20.0e8         # the seed grows within 20 Mm of the convective top (lts_c2_inf: dt1 = inf no faster than 300 s)
Z_END = 164.0e8          # the front 'reaches the base' at 164 Mm depth (2-3 cells above the bottom)
# calibration (see the docstring)
K_F, K_V, C_ON, K_R = 0.615, 0.281, 7.385, 1.316     # from 'calibrate' on runs/deep_c2.nml (2026-10-12)
CAL = dict(dt=85.0, onset=5500.0, front_days=21.5, relax=36000.0)


VAL = r'\s*=\s*(\'[^\']*\'|"[^"]*"|[^,\n/]+)'     # a namelist value (quoted strings may hold / and ,)


def read_nml(path):
    txt = open(path).read()
    body = txt[txt.index('&cz3d'):]
    return txt, body


def set_key(body, key, val):
    """replace or add `key = val` in the namelist body"""
    pat = re.compile(r'(^|[\s,])' + key + VAL, re.I | re.M)
    if pat.search(body):
        return pat.sub(lambda m: f'{m.group(1)}{key} = {val}', body, count=1)
    return body.replace('&cz3d', f'&cz3d\n  {key} = {val}', 1)


def get_key(body, key, default):
    m = re.search(r'(^|[\s,])' + key + VAL, body, re.I | re.M)
    return m.group(2).strip().strip("'\"") if m else default


def background(nml, over):
    _, body = read_nml(nml)
    for k, v in over.items():
        body = set_key(body, k, v)
    tmp = tempfile.mkdtemp(prefix='estimate_')
    for k, v in [('out_dir', f"'{tmp}'"), ('n_step_max', '0'), ('restart_from', "''"), ('dt_prof', '0.0'),
                 ('prof_every', '0'), ('plot_dt', '0.0'), ('restart_dt', '0.0'), ('dt_budget', '0.0'), ('refine_h', '0')]:
        body = set_key(body, k, v)
    f = os.path.join(tmp, 'run.nml')
    open(f, 'w').write(body)
    r = subprocess.run(['mpirun', '-np', '1', os.path.join(HERE, 'cztr3d.x'), f], cwd=HERE, capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(os.path.join(tmp, 'background.dat')):
        sys.exit('cztr3d.x failed:\n' + r.stdout[-2000:] + r.stderr[-2000:])
    head = open(os.path.join(tmp, 'background.dat')).readline().split()
    return np.loadtxt(os.path.join(tmp, 'background.dat')), body, dict(zip(head[1::2], head[2::2]))


def model(B, body, prof1d):
    cfl = float(get_key(body, 'cfl', '0.5'))
    n_theta = int(get_key(body, 'n_theta', '54'))
    margin = int(get_key(body, 'yy_margin', '-1'))
    dth = 0.5*np.pi/(n_theta - 2*(margin if margin >= 0 else 3))
    i = B[:, 0].astype(int)
    r, z, dr, rho, p, g, cs, xi = B[:, 1], -B[:, 2]*1e5, B[:, 3], B[:, 4], B[:, 5], B[:, 7], B[:, 8], B[:, 9]
    P1 = np.loadtxt(prof1d)
    Lc = np.interp(i, P1[:, 0], P1[:, 11])*LSUN       # 1D convective luminosity (same radial index)
    vm = np.maximum(Lc/(4*np.pi*r**2)/rho, 0.0)**(1/3)
    vm = np.maximum(vm, 1.0e2)                       # radiative layers: no convection (a floor for the integrals)
    vf = K_F*vm
    v = K_V*vf
    c = cs/xi
    dt_own = cfl/((v + c)/dr + 2*(v + c)/(r*dth))
    Hp = p/(rho*g)
    # the front and the onset only in the convective layers (L_conv > 1% Lsun;
    # not the radiative top of a 1.0 R domain)
    conv = Lc > 0.01*LSUN
    return dict(z=z, dr=dr, vf=vf, dt_own=dt_own, dt=dt_own.min(), L=np.maximum(Hp, r*dth), conv=conv)


def steps(M, dt1):
    dt_eff = np.maximum(M['dt'], np.minimum(M['dt_own'], dt1))
    near_top = M['conv'] & (M['z'] <= M['z'][M['conv']].min() + Z_ONSET)
    onset = C_ON/np.max((M['vf']*dt_eff/M['L'])[near_top])
    sel = (M['z'] <= Z_END) & M['conv']
    front = np.sum(M['dr'][sel]/(M['vf'][sel]*dt_eff[sel]))
    return onset, front, K_R*(onset + front), dt_eff


args = sys.argv[1:]
if not args:
    sys.exit(__doc__)
nml = args[0]
over, dt1s, do_cal = {}, [0.0, 300.0, np.inf], False
for a in args[1:]:
    if a == 'calibrate':
        do_cal = True
    elif a.startswith('dt1='):
        dt1s = [float(x) for x in a[4:].split(',')]
    elif a.startswith('prof1d='):
        pass
    else:
        k, v = a.split('=', 1)
        over[k] = v
prof1d = next((a[7:] for a in args[1:] if a.startswith('prof1d=')),
              os.path.join(HERE, '../CZTR/output/step1_cz_alone_dtcz/profile_0006.dat'))

B, body, head = background(nml, over)
if do_cal:
    # K_V: the model's global dt matches CAL['dt'] (dt_own is monotonic in K_V; bisect)
    lo, hi = 0.1, 50.0
    for _ in range(60):
        K_V = 0.5*(lo + hi)
        if model(B, body, prof1d)['dt'] > CAL['dt']:
            lo = K_V
        else:
            hi = K_V
    # K_F: time-accurate front travel (days) = CAL['front_days']; t = int dz/(K_F vm) is 1/K_F
    K_F = 1.0
    M = model(B, body, prof1d)
    sel = M['z'] <= Z_END
    t1 = np.sum(M['dr'][sel]/M['vf'][sel])/86400
    K_F = t1/CAL['front_days']
    M = model(B, body, prof1d)
    C_ON = 1.0
    C_ON = CAL['onset']/steps(M, 0.0)[0]
    on, fr, _, _ = steps(M, 0.0)
    K_R = CAL['relax']/(on + fr)
    print(f'K_F = {K_F:.3f}  K_V = {K_V:.3f}  C_ON = {C_ON:.4g}  K_R = {K_R:.3f}')

M = model(B, body, prof1d)
z = M['z']
print(f"{nml}  {' '.join(f'{k}={v}' for k, v in over.items())}")
print(f"  top {z.min()/1e8:.2f} Mm below the photosphere, {len(z)} cells; global dt {M['dt']:.1f} s; "
      f"own dt at the bottom {M['dt_own'][0]:.0f} s (alpha max {M['dt_own'].max()/M['dt']:.1f})")
print('     dt1   alpha(bottom)   onset   front    relax     (steps)   wall @0.11 s/step x cells/53')
for d in dt1s:
    on, fr, rl, dte = steps(M, d if d > 0 else 0.0)
    lab = 'inf' if np.isinf(d) else ('0' if d == 0 else f'{d:.0f}')
    wall = rl*0.11*len(z)/53/60
    print(f'  {lab:>6s}   {dte[0]/M["dt"]:9.1f}   {on:8.0f} {fr:8.0f} {rl:9.0f}            {wall:7.0f} min')
