"""
Generate the data for the interactive 3D plot on the SeaFreeze page.

    pip install SeaFreeze scipy
    python tools/seafreeze_3d_data.py

Writes files/seafreeze/water-ice-properties.json: for each stable phase of H2O
(liquid water, ices Ih, II, III, V, VI) a triangle mesh in pressure-temperature
space with its density, sound speed (Vp for ices), heat capacity, thermal
expansivity and isothermal bulk modulus at every mesh point.

The mesh is adaptive: coarse inside each stability field, three times finer in a
band along the phase transitions, and dense points placed exactly on the phase
boundaries (found by bisection with SeaFreeze). Change the settings below, run
again, commit and push.
"""
import json
import warnings
from pathlib import Path

import numpy as np
from scipy.ndimage import binary_dilation
from scipy.spatial import Delaunay
from seafreeze import seafreeze as sf

warnings.filterwarnings("ignore")

P_RANGE = (0.1, 2200.0)  # MPa (ice VII only becomes stable above ~2220 MPa)
T_RANGE = (200.0, 355.0)  # K  ('water1' is recommended for 200-355 K up to 2300 MPa)
COARSE = 3                # interior points every COARSE grid steps
P = np.linspace(*P_RANGE, 36 * COARSE + 1)  # medium grid: used near the boundaries
T = np.linspace(*T_RANGE, 33 * COARSE + 1)
P_FINE = np.linspace(*P_RANGE, 2 * (len(P) - 1) + 1)  # fine grid: to find the boundaries
T_FINE = np.linspace(*T_RANGE, 2 * (len(T) - 1) + 1)
BAND = 2                  # width of the finer band along the boundaries, in grid steps
BISECTIONS = 10           # boundary located to 1/1024 of a fine grid step

PHASES = [  # (SeaFreeze phase index, SeaFreeze material code, label)
    (0, "water1", "Liquid water"),
    (1, "Ih", "Ice Ih"),
    (2, "II", "Ice II"),
    (3, "III", "Ice III"),
    (5, "V", "Ice V"),
    (6, "VI", "Ice VI"),
]
# property name in the JSON -> (function of SeaFreeze output, decimals)
PROPS = {
    "rho": (lambda o, liquid: o.rho, 1),                        # density, kg/m3
    "vel": (lambda o, liquid: o.vel if liquid else o.Vp, 0),    # sound speed (Vp for ices), m/s
    "Cp": (lambda o, liquid: o.Cp, 0),                          # heat capacity, J/kg/K
    "alpha": (lambda o, liquid: o.alpha * 1e6, 1),              # thermal expansivity, 1e-6 / K
    "Kt": (lambda o, liquid: o.Kt / 1000, 3),                   # isothermal bulk modulus, GPa
}
OUT = Path(__file__).resolve().parent.parent / "files" / "seafreeze" / "water-ice-properties.json"


def grid_input(p, t):
    arr = np.empty((2,), dtype=object)
    arr[0], arr[1] = p, t
    return arr


def scatter(points):
    """1-D SeaFreeze scatter input from an (n, 2) array of (P, T)."""
    arr = np.empty(len(points), dtype=object)
    arr[:] = [tuple(map(float, p)) for p in points]
    return arr


def phase_at(points):
    return sf.whichphase(scatter(points)) if len(points) else np.array([])


def boundary_points(idx, fine):
    """Exact boundary points of phase `idx`, one per fine-grid edge it crosses."""
    inside = fine == idx
    pairs_in, pairs_out = [], []
    for axis in (0, 1):  # edges along pressure, then along temperature
        a = inside[:-1, :] if axis == 0 else inside[:, :-1]
        b = inside[1:, :] if axis == 0 else inside[:, 1:]
        for i, j in np.argwhere(a != b):
            p1 = (P_FINE[i], T_FINE[j])
            p2 = (P_FINE[i + 1], T_FINE[j]) if axis == 0 else (P_FINE[i], T_FINE[j + 1])
            pairs_in.append(p1 if a[i, j] else p2)
            pairs_out.append(p2 if a[i, j] else p1)
    lo, hi = np.array(pairs_in, float), np.array(pairs_out, float)
    for _ in range(BISECTIONS):
        mid = (lo + hi) / 2
        is_in = phase_at(mid) == idx
        lo[is_in], hi[~is_in] = mid[is_in], mid[~is_in]
    return lo


def main():
    medium = sf.whichphase(grid_input(P, T))           # (len(P), len(T))
    fine = sf.whichphase(grid_input(P_FINE, T_FINE))

    # band of medium-grid points along any phase transition
    change = np.zeros(medium.shape, bool)
    change[:-1, :] |= medium[:-1, :] != medium[1:, :]
    change[1:, :] |= medium[:-1, :] != medium[1:, :]
    change[:, :-1] |= medium[:, :-1] != medium[:, 1:]
    change[:, 1:] |= medium[:, :-1] != medium[:, 1:]
    band = binary_dilation(change, iterations=BAND)
    ii, jj = np.meshgrid(np.arange(len(P)), np.arange(len(T)), indexing="ij")
    coarse = (ii % COARSE == 0) & (jj % COARSE == 0)
    keep = band | coarse
    norm = np.array([P_RANGE[1] - P_RANGE[0], T_RANGE[1] - T_RANGE[0]])

    out = {
        "source": "Computed with SeaFreeze (Journaux et al. 2020; Bollengier et al. 2019): "
                  "stable phase of H2O at each pressure and temperature.",
        "P_range": list(P_RANGE),
        "T_range": list(T_RANGE),
        "phases": [],
    }
    for idx, code, label in PHASES:
        grid_pts = np.array([(P[i], T[j]) for i, j in np.argwhere(keep & (medium == idx))])
        edge_pts = boundary_points(idx, fine)
        if not len(grid_pts) and not len(edge_pts):
            continue
        pts = np.vstack([p for p in (grid_pts, edge_pts) if len(p)])
        # drop near-duplicates (a grid point sitting on the boundary)
        _, unique = np.unique(np.round((pts - [P_RANGE[0], T_RANGE[0]]) / norm, 5), axis=0, return_index=True)
        pts = pts[np.sort(unique)]

        # triangulate in normalised P-T space, keep triangles whose centre is in the phase
        tri = Delaunay((pts - [P_RANGE[0], T_RANGE[0]]) / norm)
        simplices = tri.simplices
        centres = pts[simplices].mean(axis=1)
        simplices = simplices[phase_at(centres) == idx]

        props = sf.getProp(scatter(pts), code)
        phase = {
            "id": idx, "name": label,
            "P": [round(float(x), 2) for x in pts[:, 0]],
            "T": [round(float(x), 3) for x in pts[:, 1]],
            "i": simplices[:, 0].tolist(), "j": simplices[:, 1].tolist(), "k": simplices[:, 2].tolist(),
        }
        for name, (get, digits) in PROPS.items():
            phase[name] = [round(float(v), digits) for v in get(props, idx == 0)]
        out["phases"].append(phase)
        print(f"{label:13} {len(pts):5d} points ({len(edge_pts)} on the boundary), {len(simplices):5d} triangles")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, separators=(",", ":")))
    print(f"wrote {OUT} ({OUT.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
