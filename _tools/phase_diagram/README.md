# Interactive water phase diagram (SeaFreeze page)

These scripts make `files/seafreeze/phase-diagram.html` and its data file
`files/seafreeze/phase-diagram-data.bin`, shown on the SeaFreeze page by
`_includes/seafreeze-phase-diagram.html`.

To change only the page (text, layout, colors), edit `template.html`, then run step 5.
To change the physics (SeaFreeze version, pressure or temperature range), run all steps.

## Requirements

- SeaFreeze 1.2 (beta with `water3` and `seafreeze.phasediagram`). If it is not installed,
  point `SEAFREEZE_PYTHON` to the `Python` folder of a SeaFreeze checkout.
- `pip install numpy scipy matplotlib`

## Steps

Run from this folder. Intermediate files go to `work/` (not committed).

```bash
export SEAFREEZE_PYTHON=/path/to/SeaFreeze/Python   # only if SeaFreeze 1.2 is not installed
python compute_grid.py    # 1. Gibbs energy and density of every phase on the (P, T) grid
python boundaries.py      # 2. phase boundaries and triple points, solved with SeaFreeze
python validate.py        # 3. storage masks; checks phases and densities at random points
python compute_props.py   # 4. Cp, Cv, alpha, Kt, Vp, Vs on the same grid
python check_props.py     # 4b. checks those properties at random points (a few minutes)
python build.py           # 5. writes the page and its data file into files/seafreeze/
python diffraction.py     # 6. diffraction tool: ice reflections and metastable-ice densities
```

Then preview the site and commit the files in `files/seafreeze/`.

## Diffraction tool (tools/seafreeze-diffraction)

`diffraction.py` reads the reference ice structures in `cif/` and writes
`files/seafreeze/ice-reflections.json` (reflections and structure-factor terms of ices Ih, II,
III, V and VI). When `work/grid.npz` exists (step 1), it also writes
`files/seafreeze/ice-metastable.bin` (density and Gibbs energy above the stable phase of every
ice wherever SeaFreeze evaluates it). The page is `_includes/seafreeze-diffraction.html`; it
frames the phase diagram with `dx=1` (flatter layout, estimated lattice parameters in the side
panel, pinned points shown by the framing page). `lattice=1` alone adds the lattice parameters.

## Page options

The page reads options from its address, e.g. `phase-diagram.html?color=Cp`:
`color=` starts with a property map (`rho`, `Cp`, `Cv`, `alpha`, `Kt`, `Vp`, `Vs`),
`map=0` hides the "Color by" buttons, `skin=original` uses the first design's fonts and colors,
`embed=1` is set by the SeaFreeze page.
