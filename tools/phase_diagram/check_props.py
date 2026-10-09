"""Step 4b: check the page's property values against direct SeaFreeze calls.

At random points, takes the stable phase from SeaFreeze directly and compares the value the
page shows (16-bit storage + bilinear interpolation, as in template.html) with getProp.
Writes work/props_check.json, used by build.py for the accuracy note under the diagram.
"""
import json
import os
import warnings

import numpy as np

from sfeval import S, T, names, props, P_of_s, sf, WORK

pr = np.load(os.path.join(WORK, 'props.npz'))
PROPS = ['Cp', 'Cv', 'alpha', 'Kt', 'Vp', 'Vs']
nS, nT = S.size, T.size


def quantized(a):
    """What the page decodes from the 16-bit data file."""
    m = np.isfinite(a)
    out = np.full(a.shape, np.nan)
    if not m.any():
        return out
    lo, hi = float(a[m].min()), float(a[m].max())
    f = (hi - lo) / 65534 if hi > lo else 0.0
    out[m] = lo + np.round((a[m] - lo) / (f * 65534 if f else 1) * 65534) * f
    return out


def interp(A, s, t):
    fi = (s - S[0]) / (S[1] - S[0]); fj = (t - T[0]) / (T[1] - T[0])
    i = np.clip(np.floor(fi).astype(int), 0, nS - 2); j = np.clip(np.floor(fj).astype(int), 0, nT - 2)
    u, w = fi - i, fj - j
    c = [(A[i, j], (1 - u) * (1 - w)), (A[i + 1, j], u * (1 - w)), (A[i, j + 1], (1 - u) * w), (A[i + 1, j + 1], u * w)]
    num = sum(np.where(np.isfinite(v), v * wt, 0) for v, wt in c)
    den = sum(np.where(np.isfinite(v), wt, 0) for v, wt in c)
    return np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)


rng = np.random.default_rng(2)
N = 6000
s = rng.uniform(S[0], S[-1], N); t = rng.uniform(T[0], T[-1], N); p = P_of_s(s)
Gd = np.vstack([props(n, p, t)[0] for n in names]); Rd = np.vstack([props(n, p, t)[2] for n in names])
Gd[0][np.isfinite(Gd[0]) & np.isfinite(Gd[1]) & (np.abs(Rd[0] / Rd[1] - 1) < 1e-6)] = np.nan
stab = np.nanargmin(np.where(np.isnan(Gd), np.inf, Gd), axis=0)

err = {q: [] for q in PROPS}
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    for k, n in enumerate(names):
        m = stab == k
        if not m.any():
            continue
        pts = np.empty(m.sum(), dtype=object)
        for q, (pp, tt) in enumerate(zip(p[m], t[m])):
            pts[q] = (float(pp), float(tt))
        if n in ('vap', 'liq'):
            o = sf.getProp(pts, 'water3', sf.seafreeze.defpath, 'Cp', 'Cv', 'alpha', 'Kt', 'vel',
                           branch='vapor' if n == 'vap' else 'liquid')
            direct = dict(Cp=o.Cp, Cv=o.Cv, alpha=o.alpha, Kt=o.Kt, Vp=o.vel)
        else:
            o = sf.getProp(pts, n, sf.seafreeze.defpath, 'Cp', 'Cv', 'alpha', 'Kt', 'Vp', 'Vs')
            direct = dict(Cp=o.Cp, Cv=o.Cv, alpha=o.alpha, Kt=o.Kt, Vp=o.Vp, Vs=o.Vs)
        for q, dv in direct.items():
            dv = np.asarray(dv, float).ravel()
            pv = interp(quantized(pr[q][k]), s[m], t[m])
            # alpha of liquid water crosses zero: compare to the typical size of alpha instead
            scale = np.abs(dv) if q != 'alpha' else np.maximum(np.abs(dv), 1e-4)
            e = np.abs(pv - dv) / scale
            err[q].append(e[np.isfinite(e)])
            print(f'{n:4s} {q:5s} n={m.sum():5d} max rel err {np.nanmax(e):.2e}')

worst = {q: float(np.max(np.concatenate(v))) for q, v in err.items()}
print('worst relative error:', worst)


def bound(x):
    """Round up to one significant figure, as HTML: 4.4e-4 -> 5 × 10<sup>−4</sup>."""
    e = int(np.floor(np.log10(x)))
    m = int(np.ceil(x / 10 ** e))
    if m == 10:
        m, e = 1, e + 1
    return f'{m} × 10<sup>−{-e}</sup>'


others = max(v for q, v in worst.items() if q != 'alpha')
text = (f'the other properties within {bound(others)} '
        f'(thermal expansivity within {bound(worst["alpha"])}).')
json.dump(dict(worst=worst, n=N, text=text), open(os.path.join(WORK, 'props_check.json'), 'w'))
print(text)
