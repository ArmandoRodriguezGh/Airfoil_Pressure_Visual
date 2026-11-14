# Airfoil_Pressure_Visual

Find the forces acting on an airfoil at its operating conditions, then generate a
discretized airfoil outline, on request, that is ready to import into CAD and FEA.

<img src="Example%20photos/Screenshot%202026-11-13%20at%2012.53.29%E2%80%AFPM.png" alt="Example output" width="500">

## The two parts

| Part | Tool | What it does |
|---|---|---|
| 1. Forces | `cp_visualizer.py` | Reads an airfoil `.dat` and an XFOIL `.cp` file (Cp vs x), turns Cp into real pressure for your airspeed and rib geometry, draws the pressure arrows, and averages the load into buckets for FEA. |
| 2. Geometry | `Airfoil_polyline_discretization/` | Reads the same `.dat`, puts the leading edge at the origin, scales to your chord, splits the outline into N splines, and adds an inner wall as one closed polyline. Writes a DXF. Python and MATLAB versions give the same result. |

## Quick start

Install once:

```bash
python3 -m pip install numpy matplotlib scipy ezdxf
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

## What you get

| Output | From | Contents |
|---|---|---|
| `out/e423.png` | Part 1 | pressure arrows, bucket-averaged load, Cp vs x |
| `out/e423_points.csv`, `out/e423_buckets.csv` | Part 1 | pressure at every point; average pressure and break points per bucket |
| `Airfoil_polyline_discretization/out/outline.dxf` | Part 2 | N splines (layer `SEGMENTS`) + inner wall as one closed polyline (layer `INNER`) |
| `.../out/segments.csv`, `.../out/inner_polyline.csv` | Part 2 | per-segment length, midpoint, normal; inner wall vertices |
| `.../out/outline.png` | Part 2 | check plot of the outline, break points and inner wall |

## Where things are

| Path | What |
|---|---|
| `cp_visualizer.py` | Part 1 tool |
| `Airfoil_polyline_discretization/` | Part 2 tool (`airfoil_discretize.py`, `airfoil_discretize.m`) and its `out/` |
| `e423.dat`, `e423_a12.cp`, `e423_a12.dmp` | example airfoil and XFOIL output |
| `How to for People/command_line_guide.txt` | beginner's guide to the Part 1 command |
| `How to for AI/visualized_command_spec_AI.md` | exact Part 1 command spec |
| `xfoil_to_cp_visualizer_example.html` | walkthrough from XFOIL to the visualizer |
| `CLAUDE.md` | project context for an AI assistant |
