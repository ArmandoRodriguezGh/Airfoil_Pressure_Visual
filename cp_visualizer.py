#!/usr/bin/env python3
"""
cp_visualizer.py - visualize a pressure-coefficient distribution on any airfoil,
dimensionalize it, and discretize it into average-pressure buckets for FEA loads.

Inputs
  airfoil file : Selig (.dat: name line, then TE -> upper -> LE -> lower -> TE) or
                 Lednicer (name, "Nupper Nlower" line, upper LE->TE, lower LE->TE)
  Cp file      : XFOIL CPWR output, 2 columns (x Cp) or 3 columns (x y Cp)

Requires: Python 3.8+, numpy, matplotlib      (pip install numpy matplotlib)

Examples
  python cp_visualizer.py e423.dat e423_a12.cp
  python cp_visualizer.py e423.dat e423_a12.cp --V 34 --rho 1.225 --K 24 --buckets 10
  python cp_visualizer.py e423.dat e423_a12.cp --V 34 --spacing 3 --thickness 0.125 \
         --buckets 6 --bucket-spacing adaptive --alpha 12 --chord 0.386 --save out/e423_a12

Sign convention for every pressure this tool prints:
  positive = pushes INTO the surface, negative = suction pulling the surface OUTWARD.
"""
import argparse
import csv
import math
import os
import re
import sys

import numpy as np

PA_PER_PSI = 6894.757293168


# ----------------------------------------------------------------------------
# Message log: everything is echoed to the terminal as it happens
# ----------------------------------------------------------------------------
class Log:
    def __init__(self, quiet=False):
        self.quiet = quiet
        self.warnings = []
        self.errors = []

    def info(self, msg):
        if not self.quiet:
            print(f"[ OK ] {msg}")

    def note(self, msg):
        if not self.quiet:
            print(f"[NOTE] {msg}")

    def warn(self, msg):
        self.warnings.append(msg)
        print(f"[WARN] {msg}", file=sys.stderr)

    def error(self, msg):
        self.errors.append(msg)
        print(f"[ERROR] {msg}", file=sys.stderr)


class InputError(Exception):
    """Raised for a problem that makes the inputs unusable."""


# ----------------------------------------------------------------------------
# Low-level file reading shared by both parsers
# ----------------------------------------------------------------------------
def read_numeric_rows(path, kind, log):
    """Return (rows, header_lines, comment_lines). rows = list of (line_no, [floats])."""
    if not os.path.exists(path):
        raise InputError(f"{kind} file not found: {path}")
    if os.path.isdir(path):
        raise InputError(f"{kind} path is a directory, not a file: {path}")
    try:
        with open(path, "rb") as fh:
            raw = fh.read()
    except OSError as exc:
        raise InputError(f"cannot read {kind} file {path}: {exc}")
    if not raw.strip():
        raise InputError(f"{kind} file is empty: {path}")
    if b"\x00" in raw:
        raise InputError(f"{kind} file looks binary (contains NUL bytes): {path}")
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = raw.decode("latin-1")
        log.warn(f"{kind} file is not UTF-8; decoded as Latin-1.")
    if b"\r\n" in raw:
        log.note(f"{kind} file has Windows (CRLF) line endings - handled.")

    rows, headers, comments, skipped = [], [], [], []
    for n, line in enumerate(text.splitlines(), start=1):
        s = line.strip()
        if not s:
            continue
        if s[0] in "#!%":
            comments.append(s)
            continue
        tokens = s.replace(",", " ").replace(";", " ").split()
        try:
            vals = [float(t.replace("D", "E").replace("d", "e")) for t in tokens]
        except ValueError:
            if rows:
                skipped.append(n)
            else:
                headers.append(s)
            continue
        if len(vals) < 2:
            skipped.append(n)
            continue
        rows.append((n, vals))
    if skipped:
        shown = ", ".join(str(n) for n in skipped[:8]) + (" ..." if len(skipped) > 8 else "")
        log.warn(f"{kind} file: {len(skipped)} non-numeric/short line(s) inside the data "
                 f"were skipped (line {shown}). Check the file for corruption.")
    if not rows:
        raise InputError(f"{kind} file has no numeric data rows: {path}")
    bad = [n for n, v in rows if not all(math.isfinite(x) for x in v)]
    if bad:
        raise InputError(f"{kind} file contains NaN/Inf on line(s) {bad[:8]}.")
    return rows, headers, comments


# ----------------------------------------------------------------------------
# Airfoil parser
# ----------------------------------------------------------------------------
def parse_airfoil(path, log):
    """Return dict(name, x, y) in canonical order: TE -> upper -> LE -> lower -> TE."""
    rows, headers, _ = read_numeric_rows(path, "Airfoil", log)
    name = headers[0] if headers else os.path.splitext(os.path.basename(path))[0]
    if len(headers) > 1:
        log.note(f"Airfoil file has {len(headers)} header lines; using the first as the name.")
    if any(len(v) > 2 for _, v in rows):
        log.warn("Airfoil file has more than 2 columns on some lines; using the first two as x, y.")
    pts = np.array([v[:2] for _, v in rows], dtype=float)

    # --- Lednicer detection: first numeric row is the two point counts ---
    fmt = "Selig"
    a, b = pts[0]
    if a > 1.5 and b > 1.5 and abs(a - round(a)) < 1e-6 and abs(b - round(b)) < 1e-6:
        nu, nl = int(round(a)), int(round(b))
        body = pts[1:]
        if len(body) != nu + nl:
            raise InputError(f"Airfoil file looks like Lednicer format (counts {nu}, {nl}) but has "
                             f"{len(body)} coordinate rows instead of {nu + nl}.")
        upper, lower = body[:nu], body[nu:]
        fmt = "Lednicer"
    else:
        # Lednicer without a count line: starts at LE and x jumps back to LE mid-file
        x = pts[:, 0]
        span = x.max() - x.min()
        if span <= 0:
            raise InputError("Airfoil has zero chord (all x values identical).")
        if abs(x[0] - x.min()) < 0.05 * span:
            jump = np.where(np.diff(x) < -0.5 * span)[0]
            if len(jump) == 1:
                upper, lower = pts[:jump[0] + 1], pts[jump[0] + 1:]
                fmt = "Lednicer (no count line)"
            else:
                raise InputError("Airfoil points start at the leading edge but the file is not in "
                                 "Lednicer layout. Expected Selig order: TE -> upper -> LE -> lower -> TE.")
    if fmt.startswith("Lednicer"):
        if upper[:, 1].mean() < lower[:, 1].mean():
            upper, lower = lower, upper
            log.warn("Lednicer file listed the lower surface first; surfaces swapped.")
        if np.allclose(upper[0], lower[0], atol=1e-9):
            lower = lower[1:]
        pts = np.vstack([upper[::-1], lower])
        log.note(f"Airfoil file is in {fmt} format; converted to Selig ordering.")

    if len(pts) < 10:
        raise InputError(f"Airfoil has only {len(pts)} points; at least 10 are needed.")

    # --- consecutive duplicates ---
    d = np.hypot(*np.diff(pts, axis=0).T)
    dup = np.where(d < 1e-12)[0]
    if len(dup):
        pts = np.delete(pts, dup + 1, axis=0)
        log.warn(f"Airfoil: removed {len(dup)} duplicated consecutive point(s).")

    x, y = pts[:, 0], pts[:, 1]
    chord = x.max() - x.min()
    if chord <= 0:
        raise InputError("Airfoil has zero chord (all x values identical).")

    # --- must start and end at the trailing edge ---
    if not (x.max() - x[0] < 0.05 * chord and x.max() - x[-1] < 0.05 * chord):
        raise InputError("Airfoil points do not start and end at the trailing edge. "
                         f"First x = {x[0]:.5f}, last x = {x[-1]:.5f}, max x = {x.max():.5f}. "
                         "Expected a closed loop TE -> LE -> TE.")
    ile = int(np.argmin(x))
    if ile < 3 or ile > len(x) - 4:
        raise InputError("Leading edge (minimum x) is at the start or end of the point list; "
                         "one surface is missing.")

    # --- orientation: make it counter-clockwise (upper surface first) ---
    area = 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)
    reversed_order = False
    if abs(area) < 1e-9 * chord ** 2:
        raise InputError("Airfoil encloses zero area; the two surfaces coincide.")
    if area < 0:
        pts = pts[::-1].copy()
        x, y = pts[:, 0], pts[:, 1]
        ile = int(np.argmin(x))
        reversed_order = True
        log.warn("Airfoil points were ordered lower-surface-first (clockwise); reversed to "
                 "upper-first. The Cp file will be matched in whichever order fits.")

    # --- coordinate sanity ---
    if abs(chord - 1.0) > 1e-3 or abs(x.min()) > 1e-3:
        log.warn(f"Airfoil is not normalized: x runs {x.min():.5f} .. {x.max():.5f} "
                 f"(chord = {chord:.5f}). Coefficients are normalized by this chord; make sure the "
                 "Cp file uses the same coordinates.")
    for label, seg in (("upper", x[:ile + 1][::-1]), ("lower", x[ile:])):
        nbad = int(np.sum(np.diff(seg) <= 0))
        if nbad:
            log.warn(f"Airfoil {label} surface: x is not strictly increasing LE -> TE at {nbad} "
                     "place(s) (kinked, noisy or self-overlapping data?).")
    te_gap = math.hypot(x[0] - x[-1], y[0] - y[-1])
    te = "closed (sharp)" if te_gap < 1e-6 * chord else f"open, gap = {te_gap / chord:.5f} c"
    tmax = np.interp(np.linspace(x.min(), x.max(), 201), x[:ile + 1][::-1], y[:ile + 1][::-1]) - \
        np.interp(np.linspace(x.min(), x.max(), 201), x[ile:], y[ile:])
    if tmax.min() < -1e-4 * chord:
        log.warn("Airfoil surfaces cross each other (negative thickness somewhere).")
    log.info(f"Airfoil '{name.strip()}': {len(x)} points, {fmt} format, chord = {chord:.5f}, "
             f"max thickness = {tmax.max() / chord * 100:.2f}% c, trailing edge {te}.")
    return dict(name=name.strip(), x=x.copy(), y=y.copy(), chord=chord, ile=ile,
                reversed=reversed_order)


# ----------------------------------------------------------------------------
# Cp parser
# ----------------------------------------------------------------------------
def parse_cp(path, log):
    """Return dict(x, y or None, cp, meta) in file order."""
    rows, headers, comments = read_numeric_rows(path, "Cp", log)
    ncols = {len(v) for _, v in rows}
    if len(ncols) > 1:
        raise InputError(f"Cp file has an inconsistent number of columns per line: {sorted(ncols)}.")
    nc = ncols.pop()
    data = np.array([v for _, v in rows], dtype=float)
    if nc == 2:
        x, ygeo, cp = data[:, 0], None, data[:, 1]
        layout = "x, Cp"
    elif nc == 3:
        x, ygeo, cp = data[:, 0], data[:, 1], data[:, 2]
        layout = "x, y, Cp"
    else:
        raise InputError(f"Cp file has {nc} columns; expected 2 (x Cp) or 3 (x y Cp). "
                         "A boundary-layer DUMP file is not a Cp file - use XFOIL's CPWR.")
    if len(x) < 10:
        raise InputError(f"Cp file has only {len(x)} rows; at least 10 are needed.")

    meta = {}
    for line in comments + headers:
        for key, pat in (("alpha", r"al(?:fa|pha)\s*[=:]\s*([-+]?\d*\.?\d+)"),
                         ("Re", r"\bRe\s*[=:]\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"),
                         ("Mach", r"\bMa(?:ch)?\s*[=:]\s*([-+]?\d*\.?\d+)")):
            m = re.search(pat, line, flags=re.IGNORECASE)
            if m and key not in meta:
                meta[key] = float(m.group(1))

    # --- value sanity ---
    if cp.max() > 1.05:
        log.warn(f"Cp maximum is {cp.max():.3f} (> 1). In incompressible flow the stagnation value "
                 "is exactly 1; check that the columns are in the right order.")
    if cp.max() < 0.5:
        log.warn(f"Cp maximum is only {cp.max():.3f}; no stagnation point (Cp ~ 1) was found. "
                 "The file may be incomplete or the columns may be swapped.")
    if np.abs(cp).max() > 25:
        log.warn(f"|Cp| reaches {np.abs(cp).max():.1f}, which is unusually large.")
    if np.ptp(cp) < 1e-9:
        raise InputError("Cp is constant over the whole surface; nothing to plot.")
    if np.ptp(x) <= 0:
        raise InputError("Cp file x column is constant; the columns may be swapped.")
    log.info(f"Cp file: {len(x)} rows, columns = ({layout}), Cp range {cp.min():.3f} .. {cp.max():.3f}"
             + (f", header metadata {meta}" if meta else "") + ".")
    return dict(x=x, y=ygeo, cp=cp, meta=meta)


# ----------------------------------------------------------------------------
# Compatibility: put Cp onto the airfoil's points
# ----------------------------------------------------------------------------
def match_cp_to_airfoil(af, cpd, tol, log):
    xa, ya, c = af["x"], af["y"], af["chord"]
    xc, cpv, yc = cpd["x"], cpd["cp"], cpd["y"]

    span_a, span_c = np.ptp(xa), np.ptp(xc)
    if abs(span_c - span_a) > 0.02 * span_a or abs(xc.min() - xa.min()) > 0.02 * span_a:
        raise InputError(f"Cp file x range ({xc.min():.5f} .. {xc.max():.5f}) does not match the "
                         f"airfoil x range ({xa.min():.5f} .. {xa.max():.5f}). The two files use "
                         "different normalization or belong to different geometry.")

    def geometry_dev(order):
        dev = np.abs(xc[order] - xa).max()
        if yc is not None:
            dev = max(dev, np.abs(yc[order] - ya).max())
        return dev

    if len(xc) == len(xa):
        fwd, rev = np.arange(len(xc)), np.arange(len(xc))[::-1]
        dev_f, dev_r = geometry_dev(fwd), geometry_dev(rev)
        if dev_f <= tol * c:
            log.info(f"Point-for-point match: Cp x-stations equal the airfoil's "
                     f"(max deviation {dev_f / c:.2e} c, tolerance {tol:.1e} c).")
            return cpv.copy(), "direct"
        if dev_r <= tol * c:
            log.warn("Cp file is ordered opposite to the airfoil points; reversed it to match.")
            return cpv[::-1].copy(), "reversed"
        log.warn(f"Same point count ({len(xa)}) but the x-stations differ by up to "
                 f"{min(dev_f, dev_r) / c:.4f} c (tolerance {tol:.1e} c). Likely a different "
                 "paneling or a different airfoil. Falling back to interpolation.")
    else:
        log.warn(f"Point counts differ: airfoil has {len(xa)}, Cp file has {len(xc)}. XFOIL probably "
                 "re-paneled the airfoil (PANE/PPAR). Interpolating Cp onto the airfoil points "
                 "surface-by-surface; for an exact match, SAVE the paneled airfoil from XFOIL.")

    # --- interpolation fallback: split the Cp sequence at its own leading edge ---
    jle = int(np.argmin(xc))
    if jle < 2 or jle > len(xc) - 3:
        raise InputError("Cannot split the Cp file into upper and lower surfaces: its minimum x is "
                         "at the start or end. Expected TE -> LE -> TE ordering.")
    first, second = (xc[:jle + 1][::-1], cpv[:jle + 1][::-1]), (xc[jle:], cpv[jle:])
    for label, (xs, _) in (("first", first), ("second", second)):
        if np.any(np.diff(xs) <= 0):
            log.warn(f"Cp file: x is not monotonic on its {label} surface; interpolation there "
                     "is approximate.")
    # decide which half of the Cp file is the upper surface
    if yc is not None:
        first_is_upper = yc[:jle + 1].mean() >= yc[jle:].mean()
    else:
        first_is_upper = True
        log.note("Cp file has no y column; assuming XFOIL ordering (upper surface first).")
    up, lo = (first, second) if first_is_upper else (second, first)
    ile = af["ile"]
    out = np.empty_like(xa)

    def interp_sorted(xq, xs, vs):
        o = np.argsort(xs, kind="stable")
        return np.interp(xq, xs[o], vs[o])

    out[:ile + 1] = interp_sorted(xa[:ile + 1], *up)
    out[ile:] = interp_sorted(xa[ile:], *lo)
    out[ile] = cpv[jle]
    log.warn("Cp values were INTERPOLATED onto the airfoil points. Near the leading edge, where Cp "
             "changes fastest, this can smear the suction peak - inspect the plot.")
    return out, "interpolated"


# ----------------------------------------------------------------------------
# Geometry and integration
# ----------------------------------------------------------------------------
def geometry(x, y):
    """Arc length and unit outward normals at the nodes (points are CCW)."""
    dx, dy = np.diff(x), np.diff(y)
    seg = np.hypot(dx, dy)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    nseg = np.column_stack([dy, -dx]) / seg[:, None]          # outward for CCW ordering
    nn = np.zeros((len(x), 2))
    nn[1:-1] = nseg[:-1] + nseg[1:]
    nn[0], nn[-1] = nseg[0], nseg[-1]
    nn /= np.hypot(nn[:, 0], nn[:, 1])[:, None]
    return s, nn


def integrate_coeffs(x, y, cp, chord, xref, yref=0.0):
    """Integrate a nodal Cp (linear between nodes) -> (Ca, Cn, Cm_nose_up)."""
    dx, dy = np.diff(x), np.diff(y)
    cpm = 0.5 * (cp[:-1] + cp[1:])
    fx, fy = -cpm * dy, cpm * dx                                # force = -Cp * n * ds
    xm, ym = 0.5 * (x[:-1] + x[1:]) - xref, 0.5 * (y[:-1] + y[1:]) - yref
    mz = xm * fy - ym * fx
    return fx.sum() / chord, fy.sum() / chord, -mz.sum() / chord ** 2


def lift_drag(ca, cn, alpha_deg):
    a = math.radians(alpha_deg)
    return cn * math.cos(a) - ca * math.sin(a), cn * math.sin(a) + ca * math.cos(a)


# ----------------------------------------------------------------------------
# Dimensionalization
# ----------------------------------------------------------------------------
SAFE_FUNCS = {k: getattr(np, k) for k in
              ("sqrt", "abs", "sin", "cos", "tan", "exp", "log", "log10", "where",
               "minimum", "maximum", "clip", "sign", "pi")}


def evaluate_formula(expr, variables, log):
    if "__" in expr or "import" in expr:
        raise InputError("Formula may not contain '__' or 'import'.")
    env = dict(SAFE_FUNCS)
    env.update(variables)
    try:
        out = eval(expr, {"__builtins__": {}}, env)             # user's own formula, local tool
    except Exception as exc:
        raise InputError(f"Could not evaluate --formula '{expr}': {type(exc).__name__}: {exc}. "
                         f"Available names: {', '.join(sorted(variables))}, "
                         f"{', '.join(sorted(SAFE_FUNCS))}.")
    try:
        out = np.broadcast_to(np.asarray(out, dtype=float), variables["Cp"].shape).copy()
    except Exception:
        raise InputError("--formula must produce one value per surface point (use Cp in it).")
    if not np.all(np.isfinite(out)):
        raise InputError("--formula produced NaN or Inf at some points.")
    return out


# ----------------------------------------------------------------------------
# Buckets
# ----------------------------------------------------------------------------
def build_buckets(x, y, cp, p, s, isplit, n_upper, n_lower, spacing, x_edges_u, x_edges_l,
                  chord, log):
    """Arc-length-average Cp and p over buckets on each surface. Returns list of dicts."""
    buckets = []
    surfaces = (("U", np.arange(isplit, -1, -1), n_upper, x_edges_u, -1.0),
                ("L", np.arange(isplit, len(x)), n_lower, x_edges_l, +1.0))
    for tag, idx, n, custom, sign in surfaces:
        xs, ys, cps, ps = x[idx], y[idx], cp[idx], p[idx]
        sl = np.abs(s[idx] - s[idx[0]])
        total = sl[-1]
        if custom is not None or spacing == "x":
            if np.any(np.diff(xs) <= 0):
                raise InputError(f"x-based bucket edges need x to increase monotonically from the "
                                 f"split point to the trailing edge on the "
                                 f"{'upper' if tag == 'U' else 'lower'} surface. Use "
                                 "--split le together with --bucket-spacing arc or cosine instead.")
            if custom is not None:
                xe = np.array(sorted(set(custom)), dtype=float) * chord + x.min()
                inside = xe[(xe > xs[0] + 1e-9) & (xe < xs[-1] - 1e-9)]
                if len(inside) < len(xe):
                    log.note(f"Surface {tag}: {len(xe) - len(inside)} custom edge(s) at or outside "
                             "the surface ends were dropped (the ends are added automatically).")
                xe = np.concatenate([[xs[0]], inside, [xs[-1]]])
            else:
                xe = np.linspace(xs[0], xs[-1], n + 1)
            se = np.interp(xe, xs, sl)
        elif spacing == "adaptive":
            # edges at equal steps of a blend of Cp variation (70%) and arc length (30%),
            # so buckets are short where Cp changes fastest
            var = np.concatenate([[0.0], np.cumsum(np.abs(np.diff(cps)))])
            w = 0.7 * var / max(var[-1], 1e-12) + 0.3 * sl / total
            se = np.interp(np.linspace(0, 1, n + 1), w, sl)
        elif spacing == "cosine":
            se = total * (1.0 - np.cos(0.5 * np.pi * np.linspace(0, 1, n + 1)))
        else:
            se = total * np.linspace(0, 1, n + 1)
        nb = len(se) - 1
        if nb < 1:
            raise InputError(f"Surface {tag}: no buckets were defined.")
        if nb > len(idx) - 1:
            log.warn(f"Surface {tag}: {nb} buckets but only {len(idx)} data points; some buckets "
                     "contain no data point and are filled purely by interpolation.")
        for k in range(nb):
            sa, sb = se[k], se[k + 1]
            if sb - sa <= 1e-12 * max(total, 1.0):
                log.warn(f"Bucket {tag}{k + 1} has zero length and was skipped.")
                continue
            inner = sl[(sl > sa) & (sl < sb)]
            ss = np.concatenate([[sa], inner, [sb]])
            xx, yy = np.interp(ss, sl, xs), np.interp(ss, sl, ys)
            cc, pp = np.interp(ss, sl, cps), np.interp(ss, sl, ps)
            length = sb - sa
            trap = getattr(np, "trapezoid", None) or np.trapz
            cp_avg, p_avg = trap(cc, ss) / length, trap(pp, ss) / length
            buckets.append(dict(
                id=f"{tag}{k + 1}", surface="upper" if tag == "U" else "lower",
                sa=sa, sb=sb, length=length, xa=xx[0], ya=yy[0], xb=xx[-1], yb=yy[-1],
                cp_avg=cp_avg, p_avg=p_avg, cp_min=cc.min(), cp_max=cc.max(),
                npts=int(len(inner)), xx=xx, yy=yy, ss=ss, sign=sign, sl=sl, idx=idx))
    return buckets


def bucket_coeffs(buckets, chord, xref, yref=0.0):
    """Force/moment coefficients when each bucket carries its uniform average Cp."""
    fx = fy = mz = 0.0
    for b in buckets:
        dx, dy = np.diff(b["xx"]) * b["sign"], np.diff(b["yy"]) * b["sign"]
        bfx, bfy = -b["cp_avg"] * dy, b["cp_avg"] * dx
        xm = 0.5 * (b["xx"][:-1] + b["xx"][1:]) - xref
        ym = 0.5 * (b["yy"][:-1] + b["yy"][1:]) - yref
        fx, fy, mz = fx + bfx.sum(), fy + bfy.sum(), mz + (xm * bfy - ym * bfx).sum()
    return fx / chord, fy / chord, -mz / chord ** 2


# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------
C_SUCTION, C_PRESSURE, C_BODY = "#c0392b", "#1f5fa8", "#2b2b2b"


def offset_length(pmag, pmax, arrow_scale, chord, mode):
    r = np.asarray(pmag, dtype=float) / pmax if pmax > 0 else np.zeros_like(pmag, dtype=float)
    if mode == "sqrt":
        r = np.sqrt(r)
    return arrow_scale * chord * r


def draw_load(ax, xs, ys, normals, pvals, lengths, every=1, lw=1.0, alpha_fill=0.16):
    """Envelope + arrows. Positive p: arrow points INTO the surface. Negative: out of it."""
    from matplotlib.collections import PolyCollection
    ex, ey = xs + normals[:, 0] * lengths, ys + normals[:, 1] * lengths
    quads, cols = [], []
    for i in range(len(xs) - 1):
        quads.append([(xs[i], ys[i]), (xs[i + 1], ys[i + 1]), (ex[i + 1], ey[i + 1]), (ex[i], ey[i])])
        cols.append(C_PRESSURE if 0.5 * (pvals[i] + pvals[i + 1]) >= 0 else C_SUCTION)
    ax.add_collection(PolyCollection(quads, facecolors=cols, edgecolors="none", alpha=alpha_fill))
    ax.plot(ex, ey, color="0.35", lw=0.8)
    sel = np.arange(0, len(xs), max(1, every))
    for mask, color, inward in ((pvals[sel] >= 0, C_PRESSURE, True), (pvals[sel] < 0, C_SUCTION, False)):
        k = sel[mask]
        k = k[lengths[k] > 0]
        if not len(k):
            continue
        if inward:      # tail on the envelope, head on the surface
            ax.quiver(ex[k], ey[k], xs[k] - ex[k], ys[k] - ey[k], color=color, angles="xy",
                      scale_units="xy", scale=1, width=0.0022, headwidth=4, headlength=5, zorder=3)
        else:           # tail on the surface, head on the envelope
            ax.quiver(xs[k], ys[k], ex[k] - xs[k], ey[k] - ys[k], color=color, angles="xy",
                      scale_units="xy", scale=1, width=0.0022, headwidth=4, headlength=5, zorder=3)
    return ex, ey


def make_figure(af, cp, p, nn, buckets, opts, unit_label, title_bits):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    x, y, chord, ile = af["x"], af["y"], af["chord"], af["ile"]
    pmax = float(np.abs(p).max())
    nrows = 3 if buckets else 2
    fig, axes = plt.subplots(nrows, 1, figsize=(12, 4.6 * nrows),
                             gridspec_kw=dict(height_ratios=[1.25] * (nrows - 1) + [1.0]))
    fig.suptitle("  |  ".join(title_bits), fontsize=11.5, y=0.995)
    legend = [Line2D([0], [0], color=C_SUCTION, lw=2, label="suction (p < 0): pulls surface outward"),
              Line2D([0], [0], color=C_PRESSURE, lw=2, label="pressure (p > 0): pushes into surface")]

    def style(ax, title):
        ax.plot(x, y, color=C_BODY, lw=1.6, zorder=4)
        ax.fill(x, y, color="0.93", zorder=1)
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_title(title, fontsize=10.5, loc="left")
        ax.set_xlabel("x"), ax.set_ylabel("y")
        ax.grid(True, lw=0.3, alpha=0.5)
        ax.legend(handles=legend, loc="lower right", fontsize=8, framealpha=0.9)
        ax.margins(0.06)

    # --- panel 1: continuous distribution ---
    ax = axes[0]
    L = offset_length(np.abs(p), pmax, opts.arrow_scale, chord, opts.scale_mode)
    draw_load(ax, x, y, nn, p, L, every=opts.arrow_every)
    ipk = int(np.argmax(np.abs(p)))
    ax.annotate(f"peak {p[ipk]:.4g} {unit_label}\n(Cp = {cp[ipk]:.3f})",
                xy=(x[ipk] + nn[ipk, 0] * L[ipk], y[ipk] + nn[ipk, 1] * L[ipk]),
                xytext=(nn[ipk, 0] * 14 - 30, nn[ipk, 1] * 22 - 6), textcoords="offset points",
                fontsize=8, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.7", lw=0.5))
    ist = int(np.argmax(cp))
    ax.plot(x[ist], y[ist], "o", color="k", ms=4, zorder=5)
    ax.annotate("stagnation", xy=(x[ist], y[ist]), xytext=(-nn[ist, 0] * 16 + 4, -nn[ist, 1] * 16),
                textcoords="offset points", fontsize=8, zorder=6)
    style(ax, f"Point-by-point load (arrow length {'~ sqrt|p|' if opts.scale_mode == 'sqrt' else '~ |p|'}; "
              f"longest = {opts.arrow_scale:g} chord)")

    # --- panel 2: bucket averages ---
    if buckets:
        ax = axes[1]
        for b in buckets:
            m = max(4, int(round(18 * b["length"] / chord)) + 2)
            sq = np.linspace(b["sa"], b["sb"], m)
            xq, yq = np.interp(sq, b["sl"], x[b["idx"]]), np.interp(sq, b["sl"], y[b["idx"]])
            nq = np.column_stack([np.interp(sq, b["sl"], nn[b["idx"], 0]),
                                  np.interp(sq, b["sl"], nn[b["idx"], 1])])
            nq /= np.hypot(nq[:, 0], nq[:, 1])[:, None]
            Lb = float(offset_length(abs(b["p_avg"]), pmax, opts.arrow_scale, chord, opts.scale_mode))
            pv = np.full(m, b["p_avg"])
            ex, ey = draw_load(ax, xq, yq, nq, pv, np.full(m, Lb), every=1, alpha_fill=0.22)
            for j in (0, -1):                               # bucket boundary ticks
                ax.plot([xq[j], ex[j]], [yq[j], ey[j]], color="0.35", lw=0.8)
                ax.plot(xq[j], yq[j], "|", color="k", ms=6, zorder=5)
            mid = m // 2
            lx = xq[mid] + nq[mid, 0] * (Lb + 0.018 * chord)
            ly = yq[mid] + nq[mid, 1] * (Lb + 0.018 * chord)
            text = b["id"] if len(buckets) > 14 else f"{b['id']}\n{b['p_avg']:,.0f}" if abs(b['p_avg']) >= 100 else f"{b['id']}\n{b['p_avg']:.3g}"
            ax.text(lx, ly, text, fontsize=7.5, ha="center", va="center",
                    color=C_PRESSURE if b["p_avg"] >= 0 else C_SUCTION)
        style(ax, f"Bucket-averaged load ({sum(b['surface'] == 'upper' for b in buckets)} upper + "
                  f"{sum(b['surface'] == 'lower' for b in buckets)} lower, values in {unit_label}); "
                  "ticks = spline break points")

    # --- last panel: Cp vs x ---
    ax = axes[-1]
    ax.plot(x[:ile + 1], cp[:ile + 1], color=C_SUCTION, lw=1.4, marker=".", ms=3.5, label="upper surface")
    ax.plot(x[ile:], cp[ile:], color=C_PRESSURE, lw=1.4, marker=".", ms=3.5, label="lower surface")
    for b in buckets:
        col = C_SUCTION if b["surface"] == "upper" else C_PRESSURE
        ax.plot([b["xa"], b["xb"]], [b["cp_avg"]] * 2, color=col, lw=2.6, alpha=0.45,
                solid_capstyle="butt")
    if buckets:
        ax.plot([], [], color="0.4", lw=2.6, alpha=0.5, label="bucket average")
    ax.axhline(0, color="0.5", lw=0.6)
    ax.invert_yaxis()
    ax.set_xlabel("x"), ax.set_ylabel("Cp  (negative up)")
    ax.grid(True, lw=0.3, alpha=0.5)
    ax.legend(fontsize=8, loc="upper right")
    ax.set_title("Pressure coefficient vs chordwise position", fontsize=10.5, loc="left")
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    return fig


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def build_parser():
    ap = argparse.ArgumentParser(
        description="Visualize, dimensionalize and bucket-average a Cp distribution on an airfoil.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("airfoil", help="airfoil coordinate file (Selig or Lednicer .dat)")
    ap.add_argument("cp", help="Cp file (XFOIL CPWR: 'x Cp' or 'x y Cp')")
    g = ap.add_argument_group("dimensionalization  (p = formula; default K*q*Cp)")
    g.add_argument("--V", type=float, help="airspeed [m/s]; q = 0.5*rho*V^2")
    g.add_argument("--rho", type=float, default=1.225, help="air density [kg/m^3] (default 1.225)")
    g.add_argument("--q", type=float, help="dynamic pressure [Pa], overrides --V/--rho")
    g.add_argument("--K", type=float, help="correction factor multiplying the pressure (default 1)")
    g.add_argument("--spacing", type=float, help="rib spacing; with --thickness sets K = spacing/thickness")
    g.add_argument("--thickness", type=float, help="rib thickness (same unit as --spacing)")
    g.add_argument("--load-factor", type=float, default=1.0, help="extra multiplier n (g-load x safety), default 1")
    g.add_argument("--mach", type=float, default=0.0, help="Mach number, available to --formula as M")
    g.add_argument("--formula", default="n*K*q*Cp",
                   help="pressure formula; names: Cp q K n rho V M x y s + numpy funcs "
                        "(default 'n*K*q*Cp')")
    g.add_argument("--chord", type=float, help="physical chord [m], for lengths and forces")
    g.add_argument("--alpha", type=float, help="angle of attack [deg] to report Cl / Cd_pressure")
    b = ap.add_argument_group("buckets (average pressure per surface patch)")
    b.add_argument("--buckets", type=int, default=0, help="buckets per surface (0 = none)")
    b.add_argument("--buckets-upper", type=int, help="override bucket count on the upper surface")
    b.add_argument("--buckets-lower", type=int, help="override bucket count on the lower surface")
    b.add_argument("--bucket-spacing", choices=("arc", "cosine", "adaptive", "x"), default="arc",
                   help="arc = equal arc length; cosine = clustered at the split point; "
                        "adaptive = short buckets where Cp changes fastest; x = equal x")
    b.add_argument("--upper-x-edges", type=float, nargs="+", metavar="X",
                   help="custom upper-surface break points as x/c (ends added automatically)")
    b.add_argument("--lower-x-edges", type=float, nargs="+", metavar="X",
                   help="custom lower-surface break points as x/c")
    b.add_argument("--split", choices=("le", "stagnation"), default="le",
                   help="where upper/lower buckets start: geometric LE or the stagnation point")
    v = ap.add_argument_group("plot / output")
    v.add_argument("--arrow-scale", type=float, default=0.30,
                   help="arrow length of the largest |p| as a fraction of chord (default 0.30)")
    v.add_argument("--scale-mode", choices=("linear", "sqrt"), default="linear",
                   help="sqrt compresses a tall suction peak so small loads stay visible")
    v.add_argument("--arrow-every", type=int, default=1, help="draw every Nth arrow")
    v.add_argument("--save", metavar="PREFIX", help="write PREFIX.png, PREFIX_points.csv, PREFIX_buckets.csv")
    v.add_argument("--no-show", action="store_true", help="do not open the plot window")
    v.add_argument("--dpi", type=int, default=160)
    c = ap.add_argument_group("checks")
    c.add_argument("--tol", type=float, default=1e-4, help="x/y match tolerance in chords (default 1e-4)")
    c.add_argument("--strict", action="store_true", help="treat any warning as a fatal error")
    c.add_argument("--quiet", action="store_true", help="print warnings/errors and tables only")
    return ap


def validate_options(o, log):
    def positive(name, val):
        if val is not None and not (val > 0 and math.isfinite(val)):
            raise InputError(f"--{name} must be a positive number (got {val}).")
    for name in ("V", "rho", "q", "K", "spacing", "thickness", "chord", "arrow_scale", "tol",
                 "load_factor"):
        positive(name.replace("_", "-"), getattr(o, name))
    if (o.spacing is None) != (o.thickness is None):
        raise InputError("--spacing and --thickness must be given together.")
    if o.K is not None and o.spacing is not None:
        k2 = o.spacing / o.thickness
        if abs(k2 - o.K) > 1e-9 * max(1, o.K):
            raise InputError(f"--K {o.K:g} conflicts with spacing/thickness = {k2:g}. Give one.")
    for name in ("buckets", "buckets_upper", "buckets_lower"):
        val = getattr(o, name)
        if val is not None and val < 0:
            raise InputError(f"--{name.replace('_', '-')} cannot be negative.")
    if o.arrow_every < 1:
        raise InputError("--arrow-every must be 1 or more.")
    if not 0 <= o.mach < 1:
        raise InputError("--mach must be in [0, 1).")
    for name in ("upper_x_edges", "lower_x_edges"):
        e = getattr(o, name)
        if e is not None and (min(e) < 0 or max(e) > 1):
            raise InputError(f"--{name.replace('_', '-')} values are x/c and must lie in [0, 1].")
    if o.split == "stagnation" and (o.bucket_spacing == "x" or o.upper_x_edges or o.lower_x_edges):
        raise InputError("--split stagnation cannot be combined with x-based bucket edges "
                         "(x is not monotonic around the nose). Use --bucket-spacing arc or cosine.")
    if o.q is not None and o.V is not None:
        log.warn("--q and --V both given; using --q and ignoring --V/--rho for the dynamic pressure.")


def run(argv=None):
    o = build_parser().parse_args(argv)
    log = Log(o.quiet)
    try:
        validate_options(o, log)
        af = parse_airfoil(o.airfoil, log)
        cpd = parse_cp(o.cp, log)
        cp, how = match_cp_to_airfoil(af, cpd, o.tol, log)
    except InputError as exc:
        log.error(str(exc))
        return 2

    x, y, chord, ile = af["x"], af["y"], af["chord"], af["ile"]
    s, nn = geometry(x, y)

    # --- physical sanity of the matched distribution ---
    ist = int(np.argmax(cp))
    if abs(x[ist] - x.min()) > 0.25 * chord:
        log.warn(f"Stagnation point (max Cp) is at x/c = {(x[ist] - x.min()) / chord:.3f}, far from "
                 "the leading edge. The Cp file may not belong to this airfoil.")
    cp_te = 0.5 * (cp[0] + cp[-1])
    if cp_te < -0.2:
        log.warn(f"Trailing-edge Cp is {cp_te:.2f} (attached flow recovers to roughly 0 .. +0.2). "
                 "This indicates massive separation / stall; XFOIL's Cp is not reliable here, so "
                 "treat loads from this case as rough.")
    if abs(cp[0] - cp[-1]) > 0.3:
        log.warn(f"Cp differs between the two trailing-edge points ({cp[0]:.3f} vs {cp[-1]:.3f}).")
    alpha = o.alpha if o.alpha is not None else cpd["meta"].get("alpha")
    if o.alpha is None and alpha is not None:
        log.note(f"Using alpha = {alpha:g} deg read from the Cp file header.")

    # --- dimensional pressure ---
    K = o.K if o.K is not None else (o.spacing / o.thickness if o.spacing else 1.0)
    dimensional = o.q is not None or o.V is not None
    q = o.q if o.q is not None else (0.5 * o.rho * o.V ** 2 if o.V is not None else 1.0)
    V = o.V if o.V is not None else float("nan")
    try:
        p = evaluate_formula(o.formula, dict(Cp=cp, q=q, K=K, n=o.load_factor, rho=o.rho, V=V,
                                             M=o.mach, x=x, y=y, s=s), log)
    except InputError as exc:
        log.error(str(exc))
        return 2
    if "V" in re.findall(r"[A-Za-z_]\w*", o.formula) and o.V is None:
        log.error("--formula uses V but --V was not given.")
        return 2
    unit = "Pa" if dimensional else "(x q)"
    if dimensional:
        log.info(f"Dynamic pressure q = {q:.3f} Pa ({q / PA_PER_PSI:.5f} psi)"
                 + (f" from rho = {o.rho:g} kg/m^3, V = {o.V:g} m/s" if o.q is None else "")
                 + f"; K = {K:g}; load factor n = {o.load_factor:g}; p = {o.formula}.")
    else:
        log.note(f"No --V or --q given: running NON-DIMENSIONAL (q = 1). Printed 'pressures' are "
                 f"{o.formula} with q = 1, K = {K:g}, n = {o.load_factor:g}; multiply by your q.")
    if K != 1 and o.spacing:
        log.note(f"K = spacing/thickness = {o.spacing:g}/{o.thickness:g} = {K:g} "
                 "(interior rib; use half for an end rib with a bay on one side only).")

    # --- integrated coefficients (always from raw Cp, never from K-scaled p) ---
    xref = x.min() + 0.25 * chord
    ca, cn, cm = integrate_coeffs(x, y, cp, chord, xref)
    print("\nIntegrated pressure coefficients (from Cp only; skin friction excluded)")
    print(f"  Cn (normal to chord) = {cn:+.4f}    Ca (along chord) = {ca:+.4f}    "
          f"Cm c/4 (nose-up +) = {cm:+.4f}")
    if alpha is not None:
        cl, cdp = lift_drag(ca, cn, alpha)
        print(f"  at alpha = {alpha:g} deg:  Cl = {cl:+.4f}    Cd_pressure = {cdp:+.4f}"
              "    <- compare Cl and Cm with XFOIL's own output")
    else:
        print("  (give --alpha to convert these to Cl and pressure drag and compare with XFOIL)")

    # --- buckets ---
    nu = o.buckets_upper if o.buckets_upper is not None else o.buckets
    nl = o.buckets_lower if o.buckets_lower is not None else o.buckets
    buckets = []
    want = nu > 0 or nl > 0 or o.upper_x_edges or o.lower_x_edges
    if want:
        if (nu == 0 and not o.upper_x_edges) or (nl == 0 and not o.lower_x_edges):
            log.error("Buckets requested on one surface only; give a count or edges for both.")
            return 2
        isplit = ile if o.split == "le" else ist
        try:
            buckets = build_buckets(x, y, cp, p, s, isplit, nu, nl, o.bucket_spacing,
                                    o.upper_x_edges, o.lower_x_edges, chord, log)
        except InputError as exc:
            log.error(str(exc))
            return 2
        scale_len = o.chord / chord if o.chord else 1.0
        len_unit = "m" if o.chord else "chords"
        print(f"\nBucket averages (arc-length weighted; start/end are the spline break points, "
              f"split at the {'leading edge' if o.split == 'le' else 'stagnation point'})")
        hdr = (f"  {'id':<5}{'x/c start':>10}{'y/c start':>10}{'x/c end':>10}{'y/c end':>10}"
               f"{'length ' + len_unit:>15}{'Cp avg':>10}{'Cp min':>9}{'Cp max':>9}"
               f"{'p avg ' + unit:>15}")
        if dimensional:
            hdr += f"{'p avg psi':>12}"
        print(hdr)
        for b in buckets:
            row = (f"  {b['id']:<5}{(b['xa'] - x.min()) / chord:>10.5f}{b['ya'] / chord:>10.5f}"
                   f"{(b['xb'] - x.min()) / chord:>10.5f}{b['yb'] / chord:>10.5f}"
                   f"{b['length'] / chord * (o.chord or 1.0):>15.5f}{b['cp_avg']:>10.4f}"
                   f"{b['cp_min']:>9.3f}{b['cp_max']:>9.3f}{b['p_avg']:>15.4f}")
            if dimensional:
                row += f"{b['p_avg'] / PA_PER_PSI:>12.5f}"
            print(row)
            if b["cp_max"] - b["cp_min"] > 2.0 and b["cp_min"] * b["cp_max"] < 0:
                log.warn(f"Bucket {b['id']} spans both strong suction and positive pressure "
                         f"(Cp {b['cp_min']:.2f} .. {b['cp_max']:.2f}); its average hides a sign "
                         "change. Try --bucket-spacing adaptive, --split stagnation, more buckets, or a custom break point.")
        bca, bcn, bcm = bucket_coeffs(buckets, chord, xref)
        print("\nConservation check: uniform bucket loads vs the original distribution")
        print(f"  {'':<14}{'Cn':>10}{'Ca':>10}{'Cm c/4':>10}")
        print(f"  {'original':<14}{cn:>+10.4f}{ca:>+10.4f}{cm:>+10.4f}")
        print(f"  {'bucketed':<14}{bcn:>+10.4f}{bca:>+10.4f}{bcm:>+10.4f}")

        def pct(a, b_):
            return abs(a - b_) / max(abs(b_), 1e-9) * 100
        ecn, ecm = pct(bcn, cn), pct(bcm, cm)
        print(f"  {'difference':<14}{ecn:>9.1f}%{'':>10}{ecm:>9.1f}%")
        if ecn > 5 or ecm > 10:
            log.warn(f"Bucketed loads differ from the original by {ecn:.1f}% in normal force and "
                     f"{ecm:.1f}% in pitching moment. Use more buckets, --bucket-spacing adaptive, or custom "
                     "edges around the suction peak.")
        else:
            log.info("Bucketed loads reproduce the original normal force within 5% and moment within 10%.")

    # --- outputs ---
    if o.strict and log.warnings:
        log.error(f"--strict: {len(log.warnings)} warning(s) treated as fatal. Nothing was written.")
        return 3
    if o.save:
        outdir = os.path.dirname(os.path.abspath(o.save))
        try:
            os.makedirs(outdir, exist_ok=True)
            with open(o.save + "_points.csv", "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["index", "surface", "x", "y", "s", "nx_out", "ny_out", "Cp", f"p_{unit}"])
                for i in range(len(x)):
                    w.writerow([i, "upper" if i <= ile else "lower", f"{x[i]:.6f}", f"{y[i]:.6f}",
                                f"{s[i]:.6f}", f"{nn[i, 0]:.6f}", f"{nn[i, 1]:.6f}",
                                f"{cp[i]:.6f}", f"{p[i]:.6f}"])
            if buckets:
                with open(o.save + "_buckets.csv", "w", newline="") as fh:
                    w = csv.writer(fh)
                    w.writerow(["id", "surface", "x_start", "y_start", "x_end", "y_end",
                                "arc_length_chords", "Cp_avg", "Cp_min", "Cp_max", f"p_avg_{unit}"]
                               + (["p_avg_psi"] if dimensional else []))
                    for b in buckets:
                        w.writerow([b["id"], b["surface"], f"{b['xa']:.6f}", f"{b['ya']:.6f}",
                                    f"{b['xb']:.6f}", f"{b['yb']:.6f}", f"{b['length'] / chord:.6f}",
                                    f"{b['cp_avg']:.6f}", f"{b['cp_min']:.6f}", f"{b['cp_max']:.6f}",
                                    f"{b['p_avg']:.6f}"]
                                   + ([f"{b['p_avg'] / PA_PER_PSI:.6f}"] if dimensional else []))
        except OSError as exc:
            log.error(f"Could not write output files for prefix '{o.save}': {exc}")
            return 2

    try:
        import matplotlib
        if o.no_show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        log.error("matplotlib is not installed. Install it with:  python3 -m pip install matplotlib")
        return 2
    bits = [af["name"], os.path.basename(o.cp)]
    if alpha is not None:
        bits.append(f"alpha = {alpha:g} deg")
    bits.append(f"q = {q:.1f} Pa, K = {K:g}, n = {o.load_factor:g}" if dimensional
                else f"non-dimensional (q = 1), K = {K:g}")
    if how != "direct":
        bits.append(f"Cp {how}")
    fig = make_figure(af, cp, p, nn, buckets, o, unit, bits)
    if o.save:
        fig.savefig(o.save + ".png", dpi=o.dpi)
        wrote = [o.save + ".png", o.save + "_points.csv"] + ([o.save + "_buckets.csv"] if buckets else [])
        log.info("Wrote " + ", ".join(wrote))
    if not o.no_show:
        if matplotlib.get_backend().lower() == "agg":
            log.warn("No interactive display is available; use --save PREFIX to get a PNG.")
        else:
            plt.show()
    print(f"\nDone: {len(log.warnings)} warning(s), {len(log.errors)} error(s).")
    return 0


if __name__ == "__main__":
    sys.exit(run())
