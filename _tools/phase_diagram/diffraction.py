"""Reflection list of ices Ih, II, III, V and VI for the SeaFreeze diffraction tool.

    python _tools/phase_diagram/diffraction.py

Reads the reference structures in cif/ (cell, symmetry operators, atom sites with occupancies)
and writes files/seafreeze/ice-reflections.json, loaded by tools/seafreeze-diffraction/.

For each powder peak (reflections of equal d, equivalent or not) it stores, at the reference cell:
    d         spacing in angstrom
    m         number of hkl in the peak (multiplicity)
    OO DD OD  sums over those hkl of |A_O|^2, |A_H|^2 and 2 Re(A_O conj(A_H)), where
              A_X = sum over atoms X in the cell of occupancy * exp(2 pi i h.r)
    hkl       a representative hkl of each family of equivalent reflections in the peak,
              as 'h k l', families separated by ';'
so that m |F|^2 = f_O^2 OO + f_H^2 DD + f_O f_H OD for any scattering factors f_O and f_H
(X-ray form factors, or neutron scattering lengths of H or D). The page scales every d by
(V / V_ref)^(1/3), with V from the SeaFreeze density: isotropic scaling keeps the fractional
coordinates and axial ratios of the reference structure.

If work/grid.npz exists (step 1 of the phase diagram, compute_grid.py), it also writes
files/seafreeze/ice-metastable.bin: for each ice, wherever SeaFreeze evaluates it (stable or
metastable), ln(density) and the Gibbs energy above the stable phase, G - G_stable (J/kg), on
every second node of the phase-diagram grid. Same layout as phase-diagram-data.bin: 4-byte
little-endian header length, JSON header, then gzip blobs of byte-shuffled float32.
"""
import gzip
import itertools
import struct
import json
import os
import re
import shlex

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', 'files', 'seafreeze', 'ice-reflections.json')
DMIN = 0.6           # smallest d kept at the reference cell (angstrom)
M_H2O = 18.01528     # g/mol
NA = 6.02214076e23

PHASES = [
    dict(id='Ih', cif='iceIh_KuhsLehmann.cif', setting='Hexagonal axes',
         ref='Kuhs & Lehmann (1983), J. Phys. Chem. 87, 4312; cell of Röttger et al. (1994), Acta Cryst. B 50, 644', link=''),
    dict(id='II', cif='iceII_COD1566636.cif', setting='Hexagonal axes of the rhombohedral lattice (Z = 36)',
         ref='Lobban, Finney & Kuhs (2002), J. Chem. Phys. 117, 3928; COD 1566636', link='https://www.crystallography.net/cod/1566636.html'),
    dict(id='III', cif='iceIII_COD1566638.cif', setting='Tetragonal axes',
         ref='Lobban, Finney & Kuhs (2000), J. Chem. Phys. 112, 7169; COD 1566638', link='https://www.crystallography.net/cod/1566638.html'),
    dict(id='V', cif='iceV_COD1566639.cif', setting='Monoclinic, A2/a setting (unique axis b)',
         ref='Lobban, Finney & Kuhs (2000), J. Chem. Phys. 112, 7169; COD 1566639', link='https://www.crystallography.net/cod/1566639.html'),
    dict(id='VI', cif='iceVI_COD1566641.cif', setting='Tetragonal axes, origin choice 2',
         ref='Kuhs, Ahsbahs, Londono & Finney (1989), Physica B 156–157, 684; COD 1566641', link='https://www.crystallography.net/cod/1566641.html'),
]


def num(s):
    """CIF number without its uncertainty, e.g. '0.2223(7)' -> 0.2223."""
    return float(re.sub(r'\(\d+\)$', '', s))


def read_cif(path):
    tags, loops = {}, []
    lines = [l.strip() for l in open(path) if l.strip() and not l.lstrip().startswith('#')]
    i = 0
    while i < len(lines):
        l = lines[i]
        if l == 'loop_':
            head, rows = [], []
            i += 1
            while i < len(lines) and lines[i].startswith('_'):
                head.append(lines[i].split()[0]); i += 1
            while i < len(lines) and not lines[i].startswith('_') and lines[i] != 'loop_' \
                    and not lines[i].startswith(';') and not lines[i].startswith('data_'):
                rows.append(shlex.split(lines[i])); i += 1
            loops.append((head, rows))
            continue
        if l.startswith(';'):                       # multi-line text field
            i += 1
            while i < len(lines) and not lines[i].startswith(';'):
                i += 1
            i += 1
            continue
        if l.startswith('_'):
            parts = shlex.split(l)
            if len(parts) > 1:
                tags[parts[0]] = parts[1]
        i += 1
    return tags, loops


def loop_with(loops, tag):
    for head, rows in loops:
        if tag in head:
            return [dict(zip(head, r)) for r in rows if len(r) == len(head)]
    return []


def parse_op(s):
    """'x+1/2, -y, z' -> (R, t) with x' = R x + t."""
    parts = [p.strip() for p in s.lower().split(',')]
    f = lambda x, y, z: np.array([eval(p, {}, dict(x=x, y=y, z=z)) for p in parts], float)
    t = f(0, 0, 0)
    R = np.column_stack([f(*e) - t for e in np.eye(3)])
    return np.rint(R).astype(int), t


def cell_metric(a, b, c, al, be, ga):
    al, be, ga = np.radians([al, be, ga])
    return np.array([[a * a, a * b * np.cos(ga), a * c * np.cos(be)],
                     [a * b * np.cos(ga), b * b, b * c * np.cos(al)],
                     [a * c * np.cos(be), b * c * np.cos(al), c * c]])


def expand(sites, ops):
    """All positions in the cell of each site (duplicates from special positions removed)."""
    out = []
    for s in sites:
        r0 = np.array([s['x'], s['y'], s['z']])
        pos = []
        for R, t in ops:
            p = (R @ r0 + t) % 1.0
            if not any(np.all(np.abs((p - q + 0.5) % 1.0 - 0.5) < 2e-3) for q in pos):
                pos.append(p)
        out += [(s['el'], s['occ'], p) for p in pos]
    return out


def build(ph):
    tags, loops = read_cif(os.path.join(HERE, 'cif', ph['cif']))
    cell = [num(tags['_cell_length_' + k]) for k in 'abc'] + [num(tags['_cell_angle_' + k]) for k in ('alpha', 'beta', 'gamma')]
    Z = int(num(tags['_cell_formula_units_Z']))
    oprows = loop_with(loops, '_space_group_symop_operation_xyz') or loop_with(loops, '_symmetry_equiv_pos_as_xyz')
    ops = [parse_op(r.get('_space_group_symop_operation_xyz') or r['_symmetry_equiv_pos_as_xyz']) for r in oprows]
    sites = []
    for r in loop_with(loops, '_atom_site_fract_x'):
        el = r.get('_atom_site_type_symbol') or re.match(r'[A-Za-z]+', r['_atom_site_label']).group()
        el = re.match(r'[A-Z][a-z]?', el).group()
        sites.append(dict(el='H' if el in ('D', 'H') else el, occ=num(r.get('_atom_site_occupancy', '1')),
                          x=num(r['_atom_site_fract_x']), y=num(r['_atom_site_fract_y']), z=num(r['_atom_site_fract_z'])))
    atoms = expand(sites, ops)
    nO = sum(o for e, o, _ in atoms if e == 'O')
    nH = sum(o for e, o, _ in atoms if e == 'H')
    assert abs(nO - Z) < 1e-6 and abs(nH - 2 * Z) < 0.02, (ph['id'], nO, nH, Z)

    G = cell_metric(*cell)
    V = float(np.sqrt(np.linalg.det(G)))
    Gs = np.linalg.inv(G)
    hmax = [int(np.floor(x / DMIN)) for x in cell[:3]]
    H = np.array(list(itertools.product(*[range(-m, m + 1) for m in hmax])), float)
    H = H[np.any(H != 0, axis=1)]
    d = 1 / np.sqrt(np.einsum('ij,jk,ik->i', H, Gs, H))
    H, d = H[d >= DMIN], d[d >= DMIN]

    def amp(el):
        P = np.array([p for e, _, p in atoms if e == el]); w = np.array([o for e, o, _ in atoms if e == el])
        return np.exp(2j * np.pi * H @ P.T) @ w
    AO, AH = amp('O'), amp('H')
    OO, HH, OH = np.abs(AO) ** 2, np.abs(AH) ** 2, 2 * np.real(AO * np.conj(AH))
    keep = (OO + HH) > 1e-6
    H, d, OO, HH, OH = H[keep], d[keep], OO[keep], HH[keep], OH[keep]

    # families of equivalent reflections: point-group rotations of the operators, plus Friedel pairs
    rots = {R.tobytes(): R for R, _ in ops}
    rots = list(rots.values()) + [-R for R in rots.values()]
    Hi = H.astype(int)
    index = {tuple(h): i for i, h in enumerate(Hi)}
    fam = -np.ones(len(H), int)
    reps = []
    for i, h in enumerate(Hi):
        if fam[i] >= 0:
            continue
        orbit = {tuple(h @ R) for R in rots}
        for o in orbit:
            if o in index:
                fam[index[o]] = len(reps)
        reps.append(max(orbit, key=lambda o: (sum(x >= 0 for x in o), o)))

    # powder peaks: families with the same d
    key = np.round(d, 5)
    peaks = []
    for kd in np.unique(key)[::-1]:
        m = key == kd
        fams = sorted(set(fam[m]))
        peaks.append([float(f'{d[m].mean():.6g}'), int(m.sum()), float(f'{OO[m].sum():.6g}'),
                      float(f'{HH[m].sum():.6g}'), float(f'{OH[m].sum():.6g}'),
                      ';'.join(' '.join(str(x) for x in reps[f]) for f in fams)])
    rho_ref = Z * M_H2O / (NA * V * 1e-24) * 1e3   # kg/m3 if the cell held H2O
    print(f"{ph['id']:4s} Z={Z:3d} V={V:8.3f} A3  atoms O={nO:.0f} H={nH:.2f}  hkl={len(H):6d}  peaks={len(peaks):5d}  "
          f"rho(H2O)={rho_ref:.1f}  d_max={peaks[0][0]:.4f}")
    return dict(id=ph['id'], setting=ph['setting'], ref=ph['ref'], link=ph['link'], Z=Z,
                cell=[float(f'{x:.6g}') for x in cell], V=float(f'{V:.6g}'), peaks=peaks)


data = dict(
    note='Reflections of ice reference structures for the SeaFreeze diffraction tool; see _tools/phase_diagram/diffraction.py',
    dmin=DMIN, M=M_H2O,
    phases=[build(p) for p in PHASES],
)
with open(OUT, 'w') as f:
    json.dump(data, f, separators=(',', ':'), ensure_ascii=False)
print(f'{os.path.relpath(OUT)} {os.path.getsize(OUT) / 1e3:.0f} KB')


# ---------- metastable ices: density and G above the stable phase over the whole grid ----------
GRID = os.path.join(HERE, 'work', 'grid.npz')
if os.path.exists(GRID):
    g = np.load(GRID)
    G, R, S, T, stable, PK = g['G'], g['R'], g['S'], g['T'], g['stable'], float(g['PK'])
    names = list(g['names'])
    Gst = np.take_along_axis(np.where(np.isnan(G), np.inf, G), np.maximum(stable, 0)[None], 0)[0]
    Gst[stable < 0] = np.nan
    iS = np.r_[0:S.size:2] if (S.size - 1) % 2 == 0 else np.r_[np.r_[0:S.size - 1:2], S.size - 1]
    iT = np.r_[0:T.size:2] if (T.size - 1) % 2 == 0 else np.r_[np.r_[0:T.size - 1:2], T.size - 1]
    assert np.allclose(np.diff(S[iS]), S[iS][1] - S[iS][0]) and np.allclose(np.diff(T[iT]), T[iT][1] - T[iT][0])
    blobs, pos = [], 0

    def f32(a):
        global pos
        b = np.ascontiguousarray(a, '<f4').view(np.uint8).reshape(-1, 4).T.copy()   # byte shuffle
        z = gzip.compress(b.tobytes(), 9, mtime=0)
        blobs.append(z)
        ref = [pos, len(z)]
        pos += len(z)
        return ref

    ices = [p['id'] for p in PHASES]
    L, D = {}, {}
    for pid in ices:
        k = names.index(pid)
        ok = np.isfinite(G[k]) & np.isfinite(Gst)
        lr = np.where(ok, np.log(R[k]), np.nan)[np.ix_(iS, iT)]
        dg = np.where(ok, np.maximum(G[k] - Gst, 0.0), np.nan)[np.ix_(iS, iT)]
        L[pid], D[pid] = f32(lr), f32(dg)
    head = dict(grid=dict(s0=float(S[iS][0]), ds=float(S[iS][1] - S[iS][0]), nS=int(iS.size),
                          T0=float(T[iT][0]), dT=float(T[iT][1] - T[iT][0]), nT=int(iT.size), PK=PK),
                phases=ices, L=L, dG=D)
    hb = json.dumps(head, separators=(',', ':')).encode()
    MOUT = os.path.join(HERE, '..', '..', 'files', 'seafreeze', 'ice-metastable.bin')
    with open(MOUT, 'wb') as f:
        f.write(struct.pack('<I', len(hb)) + hb + b''.join(blobs))
    print(f'{os.path.relpath(MOUT)} {os.path.getsize(MOUT) / 1e3:.0f} KB')
else:
    print('work/grid.npz not found: ice-metastable.bin not updated (run compute_grid.py first)')
