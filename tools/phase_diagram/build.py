"""Step 5: write the web page and its data file for the SeaFreeze page.

    python tools/phase_diagram/build.py

Reads work/grid.npz, work/stored.npz, work/props.npz and work/boundaries.json and writes
    files/seafreeze/phase-diagram.html       the interactive diagram (embedded on /seafreeze/)
    files/seafreeze/phase-diagram-data.bin   its data, loaded when the diagram is shown

Data file layout: 4-byte little-endian header length, the JSON header, then gzip blobs.
Gibbs energy and log density are float32 (exact phase stability); the other properties
are 16-bit values scaled between each phase's minimum and maximum (0 = no value).
"""
import gzip
import json
import os
import struct

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, 'work')
OUT = os.path.join(HERE, '..', '..', 'files', 'seafreeze')

d = np.load(os.path.join(WORK, 'grid.npz'))
st = np.load(os.path.join(WORK, 'stored.npz'))
pr = np.load(os.path.join(WORK, 'props.npz'))
B = json.load(open(os.path.join(WORK, 'boundaries.json')))
S, P, T, PK, stable = d['S'], d['P'], d['T'], float(d['PK']), d['stable']
names = list(d['names'])
assert names == ['vap', 'liq', 'Ih', 'II', 'III', 'V', 'VI']
PROPS = ['Cp', 'Cv', 'alpha', 'Kt', 'Vp', 'Vs']

blobs, pos = [], 0


def add(raw):
    global pos
    z = gzip.compress(raw, 9, mtime=0)
    blobs.append(z)
    ref = [pos, len(z)]
    pos += len(z)
    return ref


def f32(a):
    b = np.ascontiguousarray(a, '<f4').view(np.uint8).reshape(-1, 4).T.copy()   # byte shuffle
    return add(b.tobytes())


def u16(a):
    m = np.isfinite(a)
    if not m.any():
        return None
    lo, hi = float(a[m].min()), float(a[m].max())
    q = np.zeros(a.shape, '<u2')
    q[m] = 1 + np.round((a[m] - lo) / ((hi - lo) or 1.0) * 65534).astype('<u2')
    b = q.view(np.uint8).reshape(-1, 2).T.copy()
    return dict(at=add(b.tobytes()), lo=lo, hi=hi)


# color-scale range of each property over the stable liquid and ice nodes (vapour is off scale)
rng = {}
cond = (stable >= 1)
for name in ['rho'] + PROPS:
    vals = []
    for k in range(1, len(names)):
        m = cond & (stable == k)
        a = np.exp(st['Ls'][k]) if name == 'rho' else pr[name][k]
        v = a[m & np.isfinite(a)]
        if v.size:
            vals.append(v)
    v = np.concatenate(vals)
    rng[name] = [float(np.quantile(v, 0.002)), float(np.quantile(v, 0.998))]

header = dict(
    version='SeaFreeze 1.2.0b1',
    grid=dict(s0=float(S[0]), ds=float((S[-1] - S[0]) / (S.size - 1)), nS=int(S.size), T0=float(T[0]),
              dT=float(T[1] - T[0]), nT=int(T.size), PK=PK, Pmin=float(P[0]), Pmax=float(P[-1]),
              Tmin=float(T[0]), Tmax=float(T[-1])),
    G=[f32(st['Gs'][k]) for k in range(len(names))],
    L=[f32(st['Ls'][k]) for k in range(len(names))],
    props={name: [u16(pr[name][k]) for k in range(len(names))] for name in PROPS},
    range=rng,
    boundaries=[dict(a=o['a'], b=o['b'], P=[float(f'{x:.8g}') for x in o['P']], T=[round(x, 5) for x in o['T']])
                for o in B['boundaries']],
    triple=[dict(phases=t['phases'], P=float(f"{t['P']:.8g}"), T=round(t['T'], 5)) for t in B['triple']],
)
head = json.dumps(header, separators=(',', ':')).encode()
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, 'phase-diagram-data.bin'), 'wb') as f:
    f.write(struct.pack('<I', len(head)) + head + b''.join(blobs))

page = open(os.path.join(HERE, 'template.html')).read()
check = os.path.join(WORK, 'props_check.json')   # from check_props.py
accuracy = json.load(open(check))['text'] if os.path.exists(check) else 'the other properties were not checked.'
page = page.replace('__PROP_ACCURACY__', accuracy)
open(os.path.join(OUT, 'phase-diagram.html'), 'w').write(page)
size = os.path.getsize(os.path.join(OUT, 'phase-diagram-data.bin'))
print(f'phase-diagram-data.bin {size / 1e6:.2f} MB; phase-diagram.html {len(page) / 1e3:.0f} KB')
print('color ranges:', {k: [float(f'{x:.4g}') for x in v] for k, v in rng.items()})
