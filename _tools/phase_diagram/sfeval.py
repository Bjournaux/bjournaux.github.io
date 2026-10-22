"""Shared: grid, s<->P mapping and direct SeaFreeze evaluation with the grid's validity masks."""
import os, sys, json, warnings
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, 'work'); os.makedirs(WORK, exist_ok=True)   # intermediate files, not committed
if os.environ.get('SEAFREEZE_PYTHON'):        # SeaFreeze 1.2 beta checkout (its Python folder)
    sys.path.insert(0, os.environ['SEAFREEZE_PYTHON'])
import seafreeze as sf
from seafreeze.phasediagram import melt_T_dq2026, _CP_MAX, triple_points
from matplotlib.figure import Figure

d = np.load(os.path.join(WORK, 'grid.npz'))
S, P, T, G, stable, PK = d['S'], d['P'], d['T'], d['G'], d['stable'], float(d['PK'])
names = list(d['names'])
PMIN, PMAX, TMIN, TMAX = P[0], P[-1], T[0], T[-1]
s_of_P = lambda p: np.log(p) + p / PK
def P_of_s(s):
    s = np.asarray(s, float); lo = np.full(s.shape, np.log(1e-3)); hi = np.full(s.shape, np.log(1e4))
    for _ in range(200):
        m = 0.5 * (lo + hi); up = s_of_P(np.exp(m)) > s
        hi = np.where(up, m, hi); lo = np.where(up, lo, m)
    return np.exp(0.5 * (lo + hi))
SA, TB = S[-1] - S[0], TMAX - TMIN

def props(name, p, t):
    """G, S, rho of one phase at scatter points, with the same validity masks as the grid."""
    p = np.atleast_1d(p).astype(float); t = np.atleast_1d(t).astype(float)
    pts = np.empty(p.size, dtype=object)
    for q in range(p.size): pts[q] = (float(p[q]), float(t[q]))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        if name in ('vap', 'liq'):
            o = sf.getProp(pts, 'water3', sf.seafreeze.defpath, 'G', 'S', 'rho',
                           branch='vapor' if name == 'vap' else 'liquid')
            g = np.asarray(o.G, float).ravel().copy()
            g[(p > 611.657e-6) & (t < melt_T_dq2026(p) - 40.0)] = np.nan
        else:
            o = sf.getProp(pts, name, sf.seafreeze.defpath, 'G', 'S', 'rho', 'Cp', 'Kt')
            g = np.asarray(o.G, float).ravel().copy(); g[g == 0] = np.nan
            cp, kt, r = (np.asarray(getattr(o, k), float).ravel() for k in ('Cp', 'Kt', 'rho'))
            g[~((r > 0) & (kt > 0) & (cp > 0) & (cp < _CP_MAX))] = np.nan
    return g, np.asarray(o.S, float).ravel(), np.asarray(o.rho, float).ravel()

_props_raw = props
def props(name, p, t):
    try:
        return _props_raw(name, p, t)
    except AttributeError:                       # every point outside this spline
        n = np.atleast_1d(p).size
        return (np.full(n, np.nan),) * 3

def refine(a, b, p, t, fix=None, it=12):
    """Newton on G_a = G_b. fix='P' or 'T' holds that coordinate."""
    s = s_of_P(p).astype(float); t = t.astype(float).copy()
    for _ in range(it):
        p = P_of_s(s)
        ga, sa, ra = props(a, p, t); gb, sb, rb = props(b, p, t)
        dG = ga - gb
        gs = 1e6 * (1 / ra - 1 / rb) / (1 / p + 1 / PK)      # d dG / ds
        gt = -(sa - sb)                                      # d dG / dT
        gu, gw = SA * gs, TB * gt
        if fix == 'T': gw = 0 * gw
        if fix == 'P': gu = 0 * gu
        den = gu ** 2 + gw ** 2
        s = s - SA * dG * gu / den; t = t - TB * dG * gw / den
    p = P_of_s(s)
    ga, _, _ = props(a, p, t); gb, _, _ = props(b, p, t)
    return p, t, np.abs(ga - gb)

