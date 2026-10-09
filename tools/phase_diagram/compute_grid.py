"""Step 1: SeaFreeze 1.2.0b1 Gibbs energy and density of every phase on a (s, T) grid.
s = ln P + P/PK (log-like at low P, linear at high P); fluid split into liquid and vapour branches."""
import os, sys, time, warnings
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, 'work'); os.makedirs(WORK, exist_ok=True)   # intermediate files, not committed
if os.environ.get('SEAFREEZE_PYTHON'):        # SeaFreeze 1.2 beta checkout (its Python folder)
    sys.path.insert(0, os.environ['SEAFREEZE_PYTHON'])
import seafreeze as sf
from seafreeze.phasediagram import melt_T_dq2026, _CP_MAX

PMIN, PMAX, PK = 0.1, 2000.0, 80.0
TMIN, TMAX, DT = 150.0, 400.0, 1.0
DS_TARGET = 0.06

def s_of_P(P): return np.log(P) + P / PK
def P_of_s(s):
    s = np.asarray(s, float); lo = np.full(s.shape, np.log(1e-3)); hi = np.full(s.shape, np.log(1e4))
    for _ in range(200):
        m = 0.5 * (lo + hi); up = s_of_P(np.exp(m)) > s
        hi = np.where(up, m, hi); lo = np.where(up, lo, m)
    return np.exp(0.5 * (lo + hi))

s0, s1 = s_of_P(PMIN), s_of_P(PMAX)
nS = int(np.ceil((s1 - s0) / DS_TARGET)) + 1
S = np.linspace(s0, s1, nS); P = P_of_s(S); P[0], P[-1] = PMIN, PMAX
T = np.arange(TMIN, TMAX + 1e-9, DT); nT = T.size
print('grid', nS, 'x', nT, '=', nS * nT, 'nodes; dP at 2 GPa', P[-1] - P[-2], 'ratio at 0.1', P[1] / P[0])
grid = np.array([P, T], dtype=object)

names = ['vap', 'liq', 'Ih', 'II', 'III', 'V', 'VI']
G = np.full((len(names), nS, nT), np.nan); R = np.full_like(G, np.nan)
t0 = time.time()
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    Tm_melt = melt_T_dq2026(P)[:, None]
    meltbad = (P[:, None] > 611.657e-6) & (T[None, :] < Tm_melt - 40.0)
    for k, br in ((0, 'vapor'), (1, 'liquid')):
        o = sf.getProp(grid, 'water3', sf.seafreeze.defpath, 'G', 'rho', branch=br)
        g = np.asarray(o.G, float).copy(); r = np.asarray(o.rho, float).copy()
        g[meltbad] = np.nan
        G[k], R[k] = g, r
        print(names[k], 'done', round(time.time() - t0, 1), 's; finite', np.isfinite(g).mean())
    for k, ice in enumerate(names[2:], start=2):
        o = sf.getProp(grid, ice, sf.seafreeze.defpath, 'G', 'rho', 'Cp', 'Kt')
        g = np.asarray(o.G, float).copy(); g[g == 0] = np.nan
        cp, kt, r = np.asarray(o.Cp), np.asarray(o.Kt), np.asarray(o.rho)
        g[~((r > 0) & (kt > 0) & (cp > 0) & (cp < _CP_MAX))] = np.nan
        G[k], R[k] = g, r
        print(ice, 'done', round(time.time() - t0, 1), 's; finite', np.isfinite(g).mean())

# vapour and liquid branch must be distinct roots where both exist; where the branch
# solver returns the same root (no vapour root at high P), drop the vapour branch
same = np.isfinite(G[0]) & np.isfinite(G[1]) & (np.abs(R[0] / R[1] - 1) < 1e-6)
G[0][same] = np.nan
allnan = np.all(np.isnan(G), axis=0)
stable = np.where(allnan, -1, np.nanargmin(np.where(allnan[None], 0.0, G), axis=0))
print('missing nodes', allnan.sum())
for k, n in enumerate(names):
    print(f'  {n:4s} stable nodes {(stable == k).sum()}')

# cross-check against the beta's own phase_map (fluid merged)
pm = sf.phase_map(P, T)
mine = np.where(stable <= 1, 0, stable - 1)
mine[stable < 0] = -1
print('agreement with sf.phase_map:', (mine == pm.stable).mean(), 'mismatches', (mine != pm.stable).sum())
np.savez_compressed(os.path.join(WORK, 'grid.npz'), S=S, P=P, T=T, G=G, R=R, stable=stable,
                    names=np.array(names), PK=PK)
print('saved', round(time.time() - t0, 1), 's')
