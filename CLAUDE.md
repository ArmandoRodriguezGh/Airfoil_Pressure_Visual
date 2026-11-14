# CLAUDE.md: Airfoil_Pressure_Visual handoff

Context file for an AI assistant (Claude Code or similar) picking up this project.
Read it fully before changing code. It records what the user wants, why the tool is
built the way it is, the physics assumptions baked in, the tech stack, and how to
change things safely.

Repo: `ArmandoRodriguezGh/Airfoil_Pressure_Visual` (branch `main`)

---

## 1. The user and how to work with them

- Engineering student (mechanical / aero; has taken FEM labs, fluid mechanics,
  instrumentation, controls). Comfortable with physics, newer to software tooling.
- Prefers **layered explanations** (short answer first, then depth) and
  **first-principles reasoning**: derive factors like K from force balance instead of
  asserting them.
- Machine: **MacBook Pro M3 (Apple Silicon), macOS**. Give install and command advice
  for macOS/zsh (Homebrew, `python3`, `pip`). XFOIL on Mac draws through XQuartz/X11.
- Uses **Autodesk Fusion** (student license) for CAD and static stress FEA.
- Before large changes, say what will change. Keep the tool a single
  dependency-light Python file unless the user asks otherwise.

---

## 2. The engineering problem

The user is designing an RC aircraft wing and wants to stress-test a single rib in
Fusion under realistic aerodynamic loading.

| Item | Value (as stated by user) |
|---|---|
| Airfoil | Eppler **E423** (high-lift, low-Re, strongly cambered, cusped aft lower surface) |
| Span | ~89 in |
| Rib spacing | 3 in |
| Rib thickness | 1/8 in (0.125 in), balsa |
| Covering | thin heat-shrink plastic film (membrane, carries no bending) |
| XFOIL run | viscous, **Re = 9e5**, **Mach = 0.1**, α sweep −4° … 23° in 1° steps |
| Airspeed used in examples | 34 m/s (≈ Mach 0.1 at sea level). **Not confirmed as real flight speed.** |
| Chord | **unknown**: ask the user |
| Load factor / safety factor | **not decided**: ask (typical RC: 3–5 g) |

### Workflow the tool supports
```
e423.dat ──► XFOIL (VISC, ALFA, CPWR, DUMP) ──► e423_aXX.cp (+ .dmp)
                                                   │
                         cp_visualizer.py ◄────────┘
                                │
          ┌─────────────────────┼──────────────────────────┐
     plot (.png)        points CSV             buckets CSV + panels CSV
                                                       │
            Fusion: rib profile drawn as 10–15 splines, extruded to 1/8 in,
            one uniform pressure per edge face = bucket p_avg
```

### Why buckets exist (the original pain point)
The E423 DXF/polyline came into Fusion as a discretized outline, so the extrusion had
one sliver face per segment. That made load assignment painful and drove meshing RAM
up, because element size follows the smallest face. The fix: rebuild the outline as
~10–15 fit-point splines, which gives ~10–15 smooth faces, and apply one
**arc-length-averaged** pressure per face. The tool computes those averages and
reports each bucket's start/end coordinates as the **spline break points**.

---

## 3. Physics and conventions (do not change silently)

1. **Dimensional pressure:** `p = n · K · q · Cp`, with `q = ½ ρ V²`.
   Default formula string: `n*K*q*Cp`. Freestream static pressure cancels around a
   closed body, so Cp alone carries the net load.
2. **Rib edge correction K = spacing / thickness** (= 3 / 0.125 = **24** here).
   Derivation: the film transfers pressure from a tributary strip one rib-spacing wide
   onto the rib's edge face, which is only one thickness wide. Per unit arc length:
   `p_edge · t = p · spacing`, so `p_edge = p · spacing / t`.
   - End ribs (bay on one side only): **K/2 = 12**.
   - XFOIL's Cp is the true surface pressure and is **never** pre-scaled; K applies
     only to the FEA edge load. The user once asked about this; keep it explicit.
3. **Sign convention everywhere:** p > 0 pushes INTO the surface (blue arrows pointing
   inward); p < 0 is suction pulling OUTWARD (red arrows pointing outward).
4. **Integrated coefficients** (Cn, Ca, Cm about c/4, nose-up positive) come from raw
   Cp only, never from K-scaled p. With `--alpha` they convert to Cl and
   Cd_pressure, which should be compared against XFOIL's printed CL/CM. Skin friction
   is excluded.
5. **Skin friction / shear** (Cf from XFOIL `DUMP`) is ~100–1000× smaller than normal
   pressure and is carried mostly by the skin, so it is neglected for rib stress.
   Drag matters for spar chordwise bending, not for this tool.
6. **Stall check:** trailing-edge Cp < −0.2 triggers a warning (massive separation;
   XFOIL Cp unreliable). The user's `e423_a28.cp` trips it (TE Cp ≈ −0.61). Their
   `e423_a12.cp` is attached flow (TE Cp ≈ −0.06) and is a usable design case.
7. **Design load case:** near CLmax from the polar (high-g pull-up), multiplied by
   the load factor. A 2D section result also ignores the spanwise lift distribution:
   root ribs are critical, and Schrenk or XFLR5 3D can refine this.
8. **Known physical caveats** to mention when relevant: film shrink tension loads
   ribs in-plane (can rival aero loads on thin ribs); film scallops between ribs, so
   the real section differs from the true E423; balsa is orthotropic, but Fusion
   defaults to isotropic.

---

## 4. Tech stack

| Component | Version tested | Purpose |
|---|---|---|
| Python | 3.8+ (tested 3.11) | runtime |
| numpy | tested 2.x (2.5.x) | arrays, interpolation, integration |
| matplotlib | tested 3.11.x | plotting (PolyCollection envelopes, quiver arrows) |
| XFOIL | 6.99 class | produces `.cp` (CPWR) and `.dmp` (DUMP) inputs |

Install: `python3 -m pip install numpy matplotlib`.
No other dependencies. Keep it that way unless the user agrees.

NumPy 2 compatibility: the code uses `np.trapezoid` when available and falls back to
`np.trapz`. Keep that shim.

---

## 5. Files in this project

| File | Role |
|---|---|
| `cp_visualizer.py` | the tool (single file, CLI) |
| `e423.dat` | E423 coordinates, Selig format, 72 points, CRLF line endings, header `E423` |
| `e423_a28.cp` | XFOIL CPWR, 2 columns `x Cp`, stalled case |
| `e423_a12.cp` | XFOIL CPWR, attached case |
| `e423_a12.dmp` | XFOIL DUMP: `s x y Ue/Vinf Dstar Theta Cf H H* P m K` on the surface; wake rows after the TE have only 8 columns. **Not parsed by the tool yet.** |
| `How to for People/command_line_guide.txt` | novice explanation of the CLI |
| `How to for AI/visualized_command_spec_AI.md` | precise CLI spec (parameter table, grammar, outputs) |
| `xfoil_to_cp_visualizer_example.html` | end-to-end walkthrough: XFOIL commands → visualizer |
| `README.md` | repo front page (user is adding images with relative paths) |
| `Airfoil_polyline_discretization/` | Part 2 of the mission: discretized outline generator, see §5a |
| `requirements.txt`, `SETUP.md` | packages for both Python tools and install commands (venv or global); update when a dependency changes |
| `.gitignore` | keeps `.DS_Store`, `__pycache__` and `.venv` out of the repo |
| `Example photos/…png` | screenshots for the README (file names contain a macOS U+202F narrow no-break space; renaming was recommended) |

When the CLI changes, update **all three docs** (guide .txt, AI spec .md, HTML
walkthrough) plus this file.

---

## 5a. The discretizer (`Airfoil_polyline_discretization/`)

**Project mission (stated by the user):** (1) find the forces acting on an airfoil
from its operating conditions (`cp_visualizer.py`), then (2) create an automatic
discretized airfoil on request (this folder). Both parts belong; do not drop either.

- `airfoil_discretize.py` (needs numpy, scipy, ezdxf, matplotlib) and
  `airfoil_discretize.m` (base MATLAB, also runs in Octave) do the same job and must
  be kept in step. Change one, change the other, and run both.
- Output coordinates: the **leading edge (minimum-x point of the spline) is the
  origin**; the outline is scaled so the chord along x from LE to TE equals `--chord`
  (default **25.5 in**, the user's chord). No rotation, so the E423 TE sits at
  y ≈ −0.055 in.
- Two features per run, in one `out/outline.dxf`: the outer outline as N SPLINE
  entities (layer `SEGMENTS`, break point exactly on the LE) and the **inner wall**,
  offset `--wall` (default **0.3 in**) inward, as ONE closed LWPOLYLINE (layer
  `INNER`, `--inner-n` vertices, default 200). The raw offset is trimmed wherever it
  is closer than the wall to the outline, so it ends in a corner near the TE
  (x ≈ 23.21 in for E423) and the section aft of that is solid.
- Also written: `segments.csv`, `inner_polyline.csv`, `outline.png`, and from MATLAB
  `dm_curves.txt` (DesignModeler).
- **MATLAB only:** `summary.png` (`make_summary`), a graph of the relative geometry
  (x/c, y/c: outline, break points, inner wall, camber line, max thickness) and, when
  pressure is given, the average pressure per segment as arrows plus local vs average
  pressure along the chord. Pressure comes from the 7th input: an XFOIL Cp file name
  (`read_cp`, with `pscale` = n·K·q as the 10th input, default 1 = plain Cp) or a
  function handle. Sign convention as §3. Its segment averages equal the visualizer's
  spline-matched `PREFIX_panels.csv` (checked: within 0.01 Pa, Cn = 2.0152). The Python
  discretizer does not read Cp; `cp_visualizer.py` covers that side. Arrows are drawn
  with plain lines (`draw_arrow`) because `quiver` heads differ between MATLAB and
  Octave. To test in Octave on Linux: needs gnuplot, ghostscript and FreeSans fonts,
  and `graphics_toolkit('gnuplot')`.
- Reference command: `run_command_template.txt`. After a change, read the DXF back
  (expect 40 splines + 1 closed polyline) and look at `outline.png`.
- Link to Part 1: `cp_visualizer.py` writes `PREFIX_panels.csv` on these same splines
  by default (see §6), so each DXF spline has a pressure row. The match is by
  construction, not by reading files: **if the spline fit, LE definition, upper/lower
  split or spacing changes here, change `spline_panels` in `cp_visualizer.py` too**,
  then confirm panels × 25.5 equals `segments.csv` (`s0, s1, x_mid, y_mid`, normals).
  The MATLAB version uses a not-a-knot spline, so its break points differ from the
  Python ones by about 1e-5 in. Bucket break points are still not passed the other way.

---

## 6. Architecture of `cp_visualizer.py`

Single module, top to bottom:

| Function | Responsibility |
|---|---|
| `Log` | collects messages; prints `[ OK ]`, `[NOTE]` (stdout), `[WARN]`, `[ERROR]` (stderr) |
| `InputError` | raised for fatal input problems; caught in `run()` → exit 2 |
| `read_numeric_rows(path, kind, log)` | shared reader: existence, empty, binary, UTF-8 / Latin-1, CRLF, comments (`# ! %`), header lines, comma/semicolon delimiters, Fortran `D` exponents, NaN/Inf, skipped corrupt lines |
| `parse_airfoil` | Selig vs Lednicer (with or without count line) detection → canonical **CCW Selig order** (TE → upper → LE → lower → TE); duplicate removal; TE-start check; orientation via shoelace area (reverses if CW); normalization, monotonicity, crossing-surface and TE-gap checks; returns `x, y, chord, ile` |
| `parse_cp` | 2-column (`x Cp`) or 3-column (`x y Cp`); rejects other counts (e.g. DUMP files); header metadata regex for alpha / Re / Mach; Cp range sanity |
| `match_cp_to_airfoil` | x-range compatibility; exact point-for-point match within `--tol` (forward or reversed); otherwise **interpolation fallback** that splits the Cp file at its own LE and interpolates per surface, with warnings |
| `geometry` | arc length `s`, nodal outward normals (CCW → normal = (dy, −dx)) |
| `integrate_coeffs` | panel integration of nodal Cp (trapezoidal per segment) → Ca, Cn, Cm(c/4) |
| `lift_drag` | rotate by α → Cl, Cd_p |
| `evaluate_formula` | restricted `eval` (no builtins; numpy math names only; bans `__` / `import`); names `Cp q K n rho V M x y s` |
| `build_buckets` | per-surface edges in arc length: `arc`, `cosine`, `adaptive` (70% cumulative \|ΔCp\| + 30% arc length), `x` (equal x), or custom `--upper/lower-x-edges` (x/c); split at LE or stagnation; trapezoidal arc-length average of Cp and p, including interpolated bucket end points |
| `bucket_coeffs` | Cn/Ca/Cm when each bucket carries its uniform average, for the **conservation check** (warns if ΔCn > 5% or ΔCm > 10%) |
| `offset_length`, `draw_load` | load envelope offset from the surface ∝ \|p\| (or √\|p\|); filled quads colored by sign; inward arrows for p > 0, outward for p < 0 |
| `make_figure` | 3 panels: point-by-point load, bucket-averaged load with break-point ticks and labels, Cp vs x (inverted y) with bucket averages as bars |
| `build_parser`, `validate_options`, `run` | CLI, option validation and conflicts, orchestration, CSV/PNG output, exit codes |

Exit codes: `0` success, `2` input/argument error, `3` `--strict` failure (any warning).

Outputs with `--save PREFIX`: `PREFIX.png`, `PREFIX_points.csv`
(`index, surface, x, y, s, nx_out, ny_out, Cp, p_<unit>`), `PREFIX_buckets.csv`
(`id, surface, x_start, y_start, x_end, y_end, arc_length_chords, Cp_avg, Cp_min,
Cp_max, p_avg_<unit>[, p_avg_psi]`), and `PREFIX_panels.csv` (relative position
`s_c`, `x_c`, `y_c` at start/mid/end in chords, outward normal at mid, then
`Cp_avg, p_avg` and the local, not averaged, `Cp`/`p` at start/mid/end plus
`p_min, p_max`). **By default (`--match-splines`, on) its rows are the discretizer's
splines**: `spline_panels` rebuilds the discretizer's natural cubic spline, LE and
arc-length spacing with numpy only (`NaturalSpline`), so row k = `segments.csv` i = k
= k-th SPLINE in `outline.dxf`, and position × chord = DXF coordinates. `--splines`
(default 40) and `--spline-spacing` (default cosine) must equal the discretizer's
`--n` and `--spacing`. `--no-match-splines` makes the rows the buckets instead
(`panel_rows`). Bucket ids are `U1..Un` and `L1..Ln`, numbered
from the split point toward the TE.

Full CLI table: see `cp_visualizer_command_spec_for_AI.md`. Key flags: `--V --rho --q
--K --spacing --thickness --load-factor --mach --formula --chord --alpha --buckets
--buckets-upper --buckets-lower --bucket-spacing --upper-x-edges --lower-x-edges --split
--match-splines --no-match-splines --splines --spline-spacing
--arrow-scale --scale-mode --arrow-every --save --no-show --dpi --tol --strict --quiet`.

Reference command:
```bash
python3 cp_visualizer.py e423.dat e423_a12.cp --V 34 --spacing 3 --thickness 0.125 \
        --buckets 6 --bucket-spacing adaptive --alpha 12 --save out/e423_a12
```

---

## 7. Verified behavior (regression cases to re-run after edits)

On the user's files (a28 case, V = 34, K = 24, 6 adaptive buckets): exact point match;
q = 708.05 Pa; Cn = +1.7157, Ca = −0.5816, Cm(c/4) = −0.2167; bucketed Cn within 0.2%,
Cm within 8.5%; stall warning fires.

Synthetic cases tested and expected outcomes:

| Case | Expected |
|---|---|
| Lednicer airfoil (count line, blank lines) | converted to Selig, matches Cp exactly |
| Airfoil listed clockwise, comma-delimited | reversed with warning, matches |
| Cp file with half the points, 3 columns + `# Alfa = 12.0` header | interpolation warnings, runs, alpha read from header |
| Cp file with corrupt line, < 10 rows | skipped-line warning, then error, exit 2 |
| Cp x scaled ×100 | x-range mismatch error |
| Empty airfoil file / missing file | clear error, exit 2 |
| `--split stagnation` + `--bucket-spacing x` | rejected with explanation |
| Custom formula `K*q*Cp/sqrt(1-M**2)` + custom edges | runs |
| Unknown formula name | error listing the available names |
| `--strict` with any warning | exit 3, nothing written |

Add a `tests/` folder with these as pytest cases if the user wants automation (not
done yet).

---

## 8. Known limitations / open items

- `.dmp` (boundary layer) files are rejected by the Cp parser. A natural next feature:
  a `--dump FILE` option that reads Cf (surface rows only, stopping where the wake rows
  begin with fewer columns), computes `τ = n·K·q·Cf`, and reports per-bucket tangential
  load as a comparison with the normal load.
- Polar (`.pol`) files from PACC are not parsed. A possible feature: auto-find CLmax α.
- Bucket averages assume Fusion applies **one uniform pressure per face** (true as far
  as known). The arc-length average preserves each face's normal force exactly.
  Moment is only approximately preserved, hence the conservation check.
- `adaptive` spacing is a heuristic (70/30 weighting), so it can produce very short nose
  buckets for sharp suction peaks (seen on the stalled a28 case).
- Mach is not applied as a compressibility correction unless the user writes it into
  `--formula`. XFOIL's Cp already includes its Karman–Tsien correction.
- Lengths column prints in meters only when `--chord` is given; otherwise in chords.
- Airfoil coordinates are assumed normalized chord ~1. Non-normalized files warn but
  work if the Cp file uses the same coordinates.
- No spanwise (3D) load distribution; it is a single-section tool.
- Possible future features the user hinted at: DXF/spline export of the bucket break
  points for direct import into Fusion, and a batch mode over many α files.

---

## 9. XFOIL notes the user has needed

- Polar data is saved only if `PACC` (with a filename) is on **before** `ASEQ`. Syntax:
  `ASEQ first last step` (e.g. `ASEQ -4 23 1`), not MATLAB-style `-4:1:23`.
- Per-α surface data: `ALFA 12` then `CPWR file.cp` and `DUMP file.dmp`.
- Convergence: `ITER 200`; on divergence `INIT` and re-run.
- Don't `PANE` unless you also `SAVE` the paneled airfoil and give that to the tool.
  Otherwise the points won't match (the tool interpolates and warns).
- Script mode: `xfoil < commands.txt`. No comments inside; blank lines act as "Enter"
  to leave menus; `PLOP` → `G` → blank disables graphics; delete old output files
  first or XFOIL's overwrite prompt stalls the script.
- Plot window size: top-level `PLOP` → `S` (screen fraction) / `A` (aspect);
  `WDEF` saves `xfoil.def`.

---

## 10. Fusion side (context for advice)

- Import DXF → sketch → replace the polyline with fit-point splines through the
  airfoil points (denser near the LE) → extrude to 1/8 in → Simulation, Static Stress.
- One pressure load per edge face, value = bucket `p_avg` (already ×K×n). Constrain the
  rib at the spar notch faces. Use Split Face if a sub-region load is needed.
- Meshing: global element size from rib thickness (2–3 quadratic elements through
  thickness), local refinement only at holes and notches; truncate the sharp TE to
  ~0.3–0.5 mm to avoid degenerate elements.
- Education license cloud credits: allocation not guaranteed by Autodesk. The user
  should check Preferences → Tokens & Cloud Credits and pick the Education account.

---

## 11. Rules for making changes

1. Keep the sign convention, K definition, and "coefficients from raw Cp" rule unless
   the user explicitly changes them.
2. Every new input path needs validation that prints a specific `[WARN]` or
   `[ERROR]` saying what is wrong and how to fix it. The user asked for "plenty of checks."
3. Don't add dependencies beyond numpy/matplotlib without asking.
4. After any change, re-run the reference command on `e423.dat` + `e423_a12.cp` and
   `e423_a28.cp`, confirm the numbers in §7 are unchanged (unless intentionally
   changed), and look at the PNG.
5. Update the three doc files and this CLAUDE.md when flags or outputs change.
6. Commit and push finished, verified changes to `main` automatically, without
   waiting to be asked (user's instruction, 2026-11-14; replaces the earlier "commit
   only when the user asks"). Keep commit messages descriptive. Do not force-push or
   rewrite history on `main`.
