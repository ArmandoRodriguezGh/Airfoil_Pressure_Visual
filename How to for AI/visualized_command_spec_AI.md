# cp_visualizer.py: command specification (for AI assistants)

Purpose: give an assistant enough precise context to explain, validate, or adapt a
`cp_visualizer.py` command for a user. The tool is a single-file Python 3 CLI
(dependencies: `numpy`, `matplotlib`). It plots an airfoil Cp distribution as
distributed-load arrows, dimensionalizes Cp into pressure, and averages the pressure
over surface "buckets" so each bucket can be applied as one uniform FEA pressure load.

## Reference command

```bash
python3 cp_visualizer.py e423.dat e423_a28.cp --V 34 --spacing 3 --thickness 0.125 \
        --buckets 6 --bucket-spacing adaptive --alpha 23 --save out/e423
```

## Grammar

```
python3 cp_visualizer.py AIRFOIL CP [OPTIONS]

AIRFOIL  := path, Selig or Lednicer coordinate file         (positional 1, required)
CP       := path, XFOIL CPWR output: "x Cp" or "x y Cp"     (positional 2, required)
OPTIONS  := any order; all optional
```

## Core model

```
q  = --q                         if given
   = 0.5 * rho * V**2            if --V given   (rho default 1.225 kg/m^3)
   = 1  (non-dimensional mode)   otherwise
K  = --K                         if given
   = spacing / thickness         if --spacing and --thickness given
   = 1                           otherwise
n  = --load-factor               (default 1)
p  = eval(--formula)             default "n*K*q*Cp"; evaluated per surface node
sign: p > 0 pushes INTO the surface; p < 0 is suction (pulls outward)
```

Integrated Cn, Ca, Cm(c/4) always come from raw Cp (never from K or q). With an
alpha they are converted to Cl and Cd_pressure. Skin friction is excluded.

## Parameter table

| Token | Type / unit | Default | Meaning | Constraints / interactions |
|---|---|---|---|---|
| `AIRFOIL` | file path | (required) | Airfoil (x, y) loop, TE→upper→LE→lower→TE (Selig); Lednicer auto-detected and converted | ≥10 points; must start/end at TE; CRLF, commas, `D` exponents tolerated |
| `CP` | file path | (required) | Cp per surface node, XFOIL CPWR | 2 or 3 consistent columns; ≥10 rows; x range must match airfoil within 2% chord |
| `--V` | float, m/s | none | Airspeed → q | >0; ignored if `--q` given (warns) |
| `--rho` | float, kg/m³ | 1.225 | Air density | >0 |
| `--q` | float, Pa | none | Dynamic pressure, overrides V/rho | >0 |
| `--K` | float | 1 | Load correction multiplier | >0; must equal spacing/thickness if both given, else error |
| `--spacing` | float, any length unit | none | Rib spacing (tributary width) | requires `--thickness`; same unit |
| `--thickness` | float, same unit | none | Rib thickness | requires `--spacing` |
| `--load-factor` | float | 1 | Multiplier `n` (g-load, safety factor) | >0 |
| `--mach` | float | 0 | Exposed as `M` in formula only | 0 ≤ M < 1 |
| `--formula` | Python-like expr | `n*K*q*Cp` | Pressure formula | names: `Cp q K n rho V M x y s` + `sqrt abs sin cos tan exp log log10 where minimum maximum clip sign pi`; no `__`/`import`; using `V` requires `--V` |
| `--chord` | float, m | none | Physical chord; bucket lengths printed in m | >0 |
| `--alpha` | float, deg | from Cp header if present | Converts Cn/Ca to Cl, Cd_p for comparison with XFOIL | does not affect pressures; NOT parsed from filename |
| `--buckets` | int | 0 (off) | Buckets per surface | ≥0 |
| `--buckets-upper` / `--buckets-lower` | int | = `--buckets` | Per-surface override | both surfaces need a count or edges |
| `--bucket-spacing` | `arc`/`cosine`/`adaptive`/`x` | `arc` | Edge placement: equal arc length / clustered at split / 70% Cp-variation + 30% arc-length weighting / equal x | `x` incompatible with `--split stagnation` |
| `--upper-x-edges` / `--lower-x-edges` | floats, x/c | none | Custom interior break points (ends added) | values in [0, 1]; incompatible with `--split stagnation` |
| `--split` | `le`/`stagnation` | `le` | Where upper/lower buckets begin | — |
| `--arrow-scale` | float, chords | 0.30 | Length of largest arrow | >0 |
| `--scale-mode` | `linear`/`sqrt` | `linear` | Arrow length ∝ \|p\| or √\|p\| | — |
| `--arrow-every` | int | 1 | Draw every Nth arrow | ≥1 |
| `--save` | path prefix | none | Writes `PREFIX.png`, `PREFIX_points.csv`, `PREFIX_buckets.csv` (if buckets) | parent dirs created |
| `--no-show` | flag | off | No GUI window (Agg backend) | — |
| `--dpi` | int | 160 | PNG resolution | — |
| `--tol` | float, chords | 1e-4 | Point-match tolerance between files | >0 |
| `--strict` | flag | off | Any warning → exit 3, nothing written | — |
| `--quiet` | flag | off | Suppress `[ OK ]`/`[NOTE]` lines | — |

Shell note: a trailing `\` is bash/zsh line continuation, not a tool argument
(PowerShell: `` ` ``, CMD: `^`).

## Reference command decoded

| Token | Resolves to |
|---|---|
| `e423.dat` | airfoil geometry |
| `e423_a28.cp` | Cp data; "a28" is only a filename label |
| `--V 34` | q = 0.5·1.225·34² ≈ 708.05 Pa |
| `--spacing 3 --thickness 0.125` | K = 24 (interior rib; end rib would be 12) |
| `--buckets 6 --bucket-spacing adaptive` | 6 upper + 6 lower buckets, short where dCp/ds is large |
| `--alpha 23` | Cl, Cd_p reported for α = 23° |
| `--save out/e423` | out/e423.png, out/e423_points.csv, out/e423_buckets.csv |
| implied defaults | rho 1.225, n 1, formula `n*K*q*Cp`, split `le`, arrow-scale 0.30, linear |

Result: p = 1·24·708.05·Cp ≈ 16,993·Cp Pa per node.

## Outputs

- Terminal: tagged log lines `[ OK ]`, `[NOTE]`, `[WARN]` (stderr), `[ERROR]` (stderr);
  integrated coefficients; bucket table (id U1..Un / L1..Ln, start/end x/c and y/c =
  spline break points, length, Cp avg/min/max, p avg in Pa and psi); conservation
  check (bucketed vs original Cn and Cm, warns above 5% / 10%).
- `PREFIX_points.csv`: `index, surface, x, y, s, nx_out, ny_out, Cp, p_<unit>`
- `PREFIX_buckets.csv`: `id, surface, x_start, y_start, x_end, y_end, arc_length_chords, Cp_avg, Cp_min, Cp_max, p_avg_<unit>[, p_avg_psi]`
- Exit codes: 0 success, 2 input/argument error, 3 `--strict` failure.

## Built-in validation the assistant can point to

File missing/empty/binary/non-UTF-8; non-numeric lines inside data; NaN/Inf; column
count; Selig vs Lednicer detection; duplicate points; orientation (auto-reversed);
non-normalized chord; non-monotonic surfaces; crossing surfaces; x-range mismatch
between files; point-count mismatch (falls back to interpolation, with warning);
reversed Cp order; stagnation far from LE (wrong pairing); trailing-edge Cp < −0.2
(stall/separation warning); Cp max > 1.05 or < 0.5; buckets mixing strong suction
and positive pressure; bucket conservation error.

## How to adapt a command for a user (procedure)

1. Ask for or confirm: airfoil file, matching Cp file, real α of that Cp file,
   airspeed (convert to m/s) or q (Pa), air density if not sea level, rib spacing
   and thickness in one unit, whether the rib is an end rib (K/2), load factor
   (g-load × safety factor), desired number of FEA faces per surface or existing
   spline break points (x/c).
2. Map to options with the table above; prefer `--upper-x-edges/--lower-x-edges`
   when the user already has splines drawn.
3. Check constraints: `--K` vs spacing/thickness conflict; `x` spacing or custom
   edges cannot be combined with `--split stagnation`; formula names valid.
4. Warn the user that α must come from XFOIL, not the filename, and that stalled
   cases (the tool's TE-Cp warning) give unreliable loads.
5. Template:

```bash
python3 cp_visualizer.py <airfoil.dat> <case.cp> --V <m/s> --spacing <s> --thickness <t> \
        --buckets <N> --bucket-spacing <arc|cosine|adaptive|x> --alpha <deg> \
        [--load-factor <n>] [--chord <m>] [--split <le|stagnation>] \
        [--upper-x-edges <x/c ...> --lower-x-edges <x/c ...>] --save <dir>/<prefix>
```
