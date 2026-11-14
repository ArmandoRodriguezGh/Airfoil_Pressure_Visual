# Airfoil_Pressure_Visual

Find the forces acting on an airfoil at its operating conditions, then generate a
discretized airfoil outline, on request, that is ready to import into CAD and FEA.

<img src="Example%20photos/Screenshot%202026-11-13%20at%2012.53.29%E2%80%AFPM.png" alt="Example output" width="500">

## The two parts

| Part | Tool | What it does |
|---|---|---|
| 1. Forces | `cp_visualizer.py` | Reads an airfoil `.dat` and an XFOIL `.cp` file (Cp vs x), turns Cp into real pressure for your airspeed and rib geometry, draws the pressure arrows, and averages the load into buckets for FEA. |
| 2. Geometry | `Airfoil_polyline_discretization/` | Reads the same `.dat`, puts the leading edge at the origin, scales to your chord, splits the outline into N splines, and adds an inner wall as one closed polyline. Writes a DXF, and an IGES file for Patran. Python and MATLAB versions give the same result. |

## Quick start

Install once (details, including a virtual environment, are in `SETUP.md`):

```bash
python3 -m pip install -r requirements.txt
```

Part 1, pressures and bucket loads (needs numpy and matplotlib only):

```bash
python3 cp_visualizer.py e423.dat e423_a12.cp --V 34 --spacing 3 --thickness 0.125 --buckets 6 --bucket-spacing adaptive --alpha 12 --save out/e423
```

Part 2, discretized outline (run inside `Airfoil_polyline_discretization/`):

```bash
python3 airfoil_discretize.py e423.dat --chord 25.5 --wall 0.3 --n 40 --spacing cosine --units in
```

or in MATLAB / Octave:

```matlab
airfoil_discretize('e423.dat', 25.5, 40, 'cosine', 'in', 'out', [], 0.3);
```

The MATLAB version can also read the Cp file and graph the pressure on the segments
(last number is the Cp-to-pressure multiplier n·K·q in Pa):

```matlab
airfoil_discretize('e423.dat', 25.5, 40, 'cosine', 'in', 'out', '../e423_a12.cp', 0.3, 200, 16993.2);
```

## What you get

| Output | From | Contents |
|---|---|---|
| `out/e423.png` | Part 1 | pressure arrows, bucket-averaged load, Cp vs x |
| `out/e423_points.csv`, `out/e423_buckets.csv` | Part 1 | pressure at every point; average pressure and break points per bucket |
| `out/e423_panels.csv` | Part 1 | one row per panel: relative position on the chord, average pressure, and the local pressure at the panel start, middle and end. By default the panels are the same 40 splines Part 2 writes to its DXF, in the same order (`--splines N` to change the count, `--no-match-splines` to use the buckets instead) |
| `Airfoil_polyline_discretization/out/outline.dxf` | Part 2 | N splines (layer `SEGMENTS`) + inner wall as one closed polyline (layer `INNER`) |
| `.../out/outline_patran.igs` | Part 2, Python only | the same geometry as IGES curves and points for Patran: 40 outer splines, the inner wall cut into pieces that pair with them, a connector across the wall at each break point, and a point at every break point for attaching internal members |
| `.../out/wall_pairs.csv` | Part 2, Python only | coordinates of the paired outer / inner break points |
| `.../out/segments.csv`, `.../out/inner_polyline.csv` | Part 2 | per-segment length, midpoint, normal; inner wall vertices |
| `.../out/outline.png` | Part 2 | check plot of the outline, break points and inner wall |
| `.../out/summary.png` | Part 2, MATLAB only | relative geometry (x/c, y/c, thickness, camber, inner wall) and, with a Cp file, the pressure on each segment |

## Where things are

| Path | What |
|---|---|
| `cp_visualizer.py` | Part 1 tool |
| `Airfoil_polyline_discretization/` | Part 2 tool (`airfoil_discretize.py`, `airfoil_discretize.m`) and its `out/` |
| `e423.dat`, `e423_a12.cp`, `e423_a12.dmp` | example airfoil and XFOIL output |
| `How to for People/command_line_guide.txt` | beginner's guide to the Part 1 command |
| `How to for AI/visualized_command_spec_AI.md` | exact Part 1 command spec |
| `xfoil_to_cp_visualizer_example.html` | walkthrough from XFOIL to the visualizer |
| `SETUP.md`, `requirements.txt` | packages needed and the install commands |
| `PRESSURE_ASSUMPTIONS.md` | every value used and assumed in turning Cp into pressure |
| `CLAUDE.md` | project context for an AI assistant |
