"""Step 2: phase boundaries (G_i = G_j) and triple points, refined against SeaFreeze directly."""
from sfeval import *
from sfeval import props

# triple points (beta's own Newton refinement), in the window
pm = sf.phase_map(P, T)
tps = []
for tp in triple_points(pm):
    if not (PMIN <= tp['P'] <= PMAX and TMIN <= tp['T'] <= TMAX): continue
    ph = ['liq' if n == 'water3' else n for n in tp['phases']]
    tps.append(dict(phases=ph, P=float(tp['P']), T=float(tp['T'])))
print('triple points:'); [print('  ', t) for t in tps]

fig = Figure(); ax = fig.subplots()
out = []
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        if not ((stable == i).any() and (stable == j).any()): continue
        pair = (stable == i) | (stable == j)
        Z = np.where(pair, G[i] - G[j], np.nan)
        cs = ax.contour(S, T, Z.T, levels=[0.0])
        for seg in cs.allsegs[0]:
            if len(seg) < 3: continue
            a, b = names[i], names[j]
            p, t, res = refine(a, b, P_of_s(seg[:, 0]), seg[:, 1])
            # drop vertices where a third phase is more stable (metastable extension)
            ga, _, _ = props(a, p, t)
            others = [n for n in names if n not in (a, b)]
            gmin = np.nanmin(np.vstack([props(n, p, t)[0] for n in others] + [np.full(p.size, np.inf)]), axis=0)
            keep = ~(gmin < ga - 1e-6)
            p, t = p[keep], t[keep]
            if p.size < 2: continue
            # snap ends to triple points of this pair, or extend to the window edge
            pts = [list(x) for x in zip(p, t)]
            for end in (0, -1):
                pe, te = pts[end]
                best = None
                for tp in tps:
                    if a in tp['phases'] and b in tp['phases']:
                        dist = np.hypot((s_of_P(tp['P']) - s_of_P(pe)) / SA, (tp['T'] - te) / TB)
                        if dist < 0.03 and (best is None or dist < best[0]): best = (dist, tp)
                if best is not None:
                    q = [best[1]['P'], best[1]['T']]
                    if end == 0: pts.insert(0, q)
                    else: pts.append(q)
                    continue
                # extend to the nearest window edge if within ~2 cells
                ds = S[1] - S[0]
                sE, tE = s_of_P(pe), te
                cand = [('s0', abs(sE - S[0]) / ds), ('s1', abs(S[-1] - sE) / ds),
                        ('t0', abs(tE - TMIN) / (T[1] - T[0])), ('t1', abs(TMAX - tE) / (T[1] - T[0]))]
                edge, dcell = min(cand, key=lambda c: c[1])
                if dcell > 2.5 or dcell < 1e-9: continue
                nb = pts[1] if end == 0 else pts[-2]
                s_nb, t_nb = s_of_P(nb[0]), nb[1]
                # linear extrapolation to the edge, then Newton along it
                if edge in ('s0', 's1'):
                    sx = S[0] if edge == 's0' else S[-1]
                    tx = tE + (t_nb - tE) * (sx - sE) / (s_nb - sE) if s_nb != sE else tE
                    pp, tt, r = refine(a, b, np.array([P_of_s(sx)]), np.array([tx]), fix='P')
                else:
                    tx = TMIN if edge == 't0' else TMAX
                    sx = sE + (s_nb - sE) * (tx - tE) / (t_nb - tE) if t_nb != tE else sE
                    pp, tt, r = refine(a, b, np.array([P_of_s(sx)]), np.array([tx]), fix='T')
                if r[0] < 1e-3 and PMIN * (1 - 1e-9) <= pp[0] <= PMAX * (1 + 1e-9):
                    q = [float(np.clip(pp[0], PMIN, PMAX)), float(np.clip(tt[0], TMIN, TMAX))]
                    if end == 0: pts.insert(0, q)
                    else: pts.append(q)
            arr = np.array(pts)
            out.append(dict(a=a, b=b, P=arr[:, 0].tolist(), T=arr[:, 1].tolist()))
            print(f'{a:4s}-{b:4s} {len(arr):4d} pts  P {arr[:,0].min():9.4f}-{arr[:,0].max():9.3f}  '
                  f'T {arr[:,1].min():7.3f}-{arr[:,1].max():7.3f}  max|dG| {res.max():.2e} J/kg  dropped {int((~keep).sum())}')

# liquid-vapour check against sf.saturation
lv = [o for o in out if {o['a'], o['b']} == {'vap', 'liq'}]
for o in lv:
    sat = sf.saturation(np.array(o['T']))
    print('L-V vs sf.saturation: max rel dP', np.nanmax(np.abs(np.array(o['P']) / sat.P - 1)))
json.dump(dict(boundaries=out, triple=tps), open(os.path.join(WORK, 'boundaries.json'), 'w'))
print('saved', len(out), 'boundary polylines')
