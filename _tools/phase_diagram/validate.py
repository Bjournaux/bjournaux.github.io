"""Step 3: storage masks + check the page's interpolation against direct SeaFreeze calls."""
import os, sys, json, gzip
import numpy as np
from scipy.ndimage import binary_dilation
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
WORK = os.path.join(HERE, 'work'); os.makedirs(WORK, exist_ok=True)   # intermediate files, not committed
from sfeval import S, P, T, G, stable, names, props, s_of_P, PK
R = np.load(os.path.join(WORK, 'grid.npz'))['R']
DIL = 3
nph, nS, nT = G.shape
Gs = np.full(G.shape, np.nan, np.float32); Ls = np.full(G.shape, np.nan, np.float32)
for k in range(nph):
    m = binary_dilation(stable == k, iterations=DIL, structure=np.ones((3, 3))) & np.isfinite(G[k])
    Gs[k][m] = G[k][m]; Ls[k][m] = np.log(R[k][m])
    print(f'{names[k]:4s} stored nodes {m.sum():6d}')

def interp(s, t):
    """Exactly what the page does: candidate phases = all 4 corners stored; bilinear in (s, T)."""
    fi = (s - S[0]) / (S[1] - S[0]); fj = (t - T[0]) / (T[1] - T[0])
    i = np.clip(np.floor(fi).astype(int), 0, nS - 2); j = np.clip(np.floor(fj).astype(int), 0, nT - 2)
    u = fi - i; w = fj - j
    def bil(A):
        return ((1 - u) * (1 - w) * A[i, j] + u * (1 - w) * A[i + 1, j]
                + (1 - u) * w * A[i, j + 1] + u * w * A[i + 1, j + 1])
    g = np.vstack([bil(Gs[k].astype(float)) for k in range(nph)])     # NaN if any corner missing
    cls = np.nanargmin(np.where(np.isnan(g), np.inf, g), axis=0)
    lr = np.vstack([bil(Ls[k].astype(float)) for k in range(nph)])
    return cls, np.exp(lr[cls, np.arange(s.size)]), g, lr

rng = np.random.default_rng(1)
N = 20000
s = rng.uniform(S[0], S[-1], N); t = rng.uniform(T[0], T[-1], N)
from sfeval import P_of_s
p = P_of_s(s)
cls, rho_i, gI, _ = interp(s, t)
Gd = np.vstack([props(n, p, t)[0] for n in names]); Rd = np.vstack([props(n, p, t)[2] for n in names])
Gd[0][np.isfinite(Gd[0]) & np.isfinite(Gd[1]) & (np.abs(Rd[0] / Rd[1] - 1) < 1e-6)] = np.nan   # same rule as the grid
cd = np.nanargmin(np.where(np.isnan(Gd), np.inf, Gd), axis=0)
rho_d = Rd[cd, np.arange(N)]
mis = cls != cd
print('random points', N, 'class mismatches', mis.sum())
if mis.any():
    # how far from the boundary? |dG| between the two phases / |dS| -> K
    for q in np.flatnonzero(mis)[:10]:
        a, b = cls[q], cd[q]
        _, Sa, _ = props(names[a], p[q:q+1], t[q:q+1]); _, Sb, _ = props(names[b], p[q:q+1], t[q:q+1])
        print(f'   P {p[q]:.3f} T {t[q]:.3f}: page {names[a]} vs direct {names[b]}, dG {Gd[a,q]-Gd[b,q]:.3g} J/kg, ~{abs(Gd[a,q]-Gd[b,q])/abs(Sa[0]-Sb[0]):.2e} K')
ok = ~mis
rel = np.abs(rho_i[ok] / rho_d[ok] - 1)
print(f'rho rel error: max {rel.max():.2e}  99.9% {np.quantile(rel, .999):.2e}  median {np.median(rel):.2e}')
for k, n in enumerate(names):
    m = ok & (cd == k)
    if m.any():
        r = np.abs(rho_i[m] / rho_d[m] - 1); print(f'   {n:4s} n={m.sum():5d} max {r.max():.2e}')

# boundary consistency: interpolated dG at the refined boundary vertices, expressed as a T or P offset
B = json.load(open(os.path.join(WORK, 'boundaries.json')))
worstT = worstP = 0
for o in B['boundaries']:
    a, b = names.index(o['a']), names.index(o['b'])
    pb, tb = np.array(o['P']), np.array(o['T'])
    _, _, g, _ = interp(s_of_P(pb), tb)
    dg = g[a] - g[b]
    _, Sa, ra = props(o['a'], pb, tb); _, Sb, rb = props(o['b'], pb, tb)
    offT = np.abs(dg / (Sa - Sb)); offP = np.abs(dg / (1e6 * (1 / ra - 1 / rb)))
    off = np.minimum(offT / 250, offP / 2000)        # in fractions of the full view
    print(f"{o['a']:4s}-{o['b']:4s} field edge vs line: max {np.nanmax(off)*1e4:.2f} e-4 of full view; T-off max {np.nanmax(offT):.3g} K, P-off max {np.nanmax(offP):.3g} MPa")
np.savez_compressed(os.path.join(WORK, 'stored.npz'), Gs=Gs, Ls=Ls)

def pack(a):
    b = np.ascontiguousarray(a, '<f4').view(np.uint8).reshape(-1, 4).T.copy()   # byte shuffle
    return gzip.compress(b.tobytes(), 9)
tot = sum(len(pack(Gs[k])) + len(pack(Ls[k])) for k in range(nph))
print('packed size (shuffled gzip):', tot / 1e6, 'MB; base64', tot * 4 / 3 / 1e6, 'MB')
