"""Step 1b: Cp, Cv, alpha, Kt, Vp and Vs of every phase on the same (s, T) grid,
kept only where the page stores that phase (same masks as the density)."""
import os, sys, time, warnings
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, 'work'); os.makedirs(WORK, exist_ok=True)   # intermediate files, not committed
if os.environ.get('SEAFREEZE_PYTHON'):        # SeaFreeze 1.2 beta checkout (its Python folder)
    sys.path.insert(0, os.environ['SEAFREEZE_PYTHON'])
import seafreeze as sf

d = np.load(os.path.join(WORK, 'grid.npz')); st = np.load(os.path.join(WORK, 'stored.npz'))
P, T, names = d['P'], d['T'], list(d['names'])
stored = np.isfinite(st['Ls'])                       # (phase, s, T) nodes the page keeps
grid = np.array([P, T], dtype=object)
PROPS = ['Cp', 'Cv', 'alpha', 'Kt', 'Vp', 'Vs']
out = {k: np.full(stored.shape, np.nan, np.float32) for k in PROPS}
t0 = time.time()
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    for k, n in enumerate(names):
        if n in ('vap', 'liq'):
            o = sf.getProp(grid, 'water3', sf.seafreeze.defpath, 'Cp', 'Cv', 'alpha', 'Kt', 'vel',
                           branch='vapor' if n == 'vap' else 'liquid')
            vals = dict(Cp=o.Cp, Cv=o.Cv, alpha=o.alpha, Kt=o.Kt, Vp=o.vel, Vs=None)
        else:
            o = sf.getProp(grid, n, sf.seafreeze.defpath, 'Cp', 'Cv', 'alpha', 'Kt', 'Vp', 'Vs')
            vals = dict(Cp=o.Cp, Cv=o.Cv, alpha=o.alpha, Kt=o.Kt, Vp=o.Vp, Vs=o.Vs)
        m = stored[k]
        for p in PROPS:
            if vals[p] is None: continue
            a = np.asarray(vals[p], float)
            out[p][k][m] = a[m]
        print(n, 'done', round(time.time() - t0, 1), 's', {p: (float(np.nanmin(out[p][k])), float(np.nanmax(out[p][k]))) if np.isfinite(out[p][k]).any() else None for p in PROPS})
np.savez_compressed(os.path.join(WORK, 'props.npz'), **out)
print('saved', round(time.time() - t0, 1), 's')
