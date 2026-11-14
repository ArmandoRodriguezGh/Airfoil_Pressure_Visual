#!/usr/bin/env python3
"""
airfoil_discretize.py
Read an airfoil .dat (Selig or Lednicer format), fit ONE continuous spline
through the whole outline, split it into N arcs of controlled length, and
export per-segment geometry (points, tangents, outward normals, lengths)
plus a DXF with each arc as its own SPLINE entity.

Usage:
    python3 airfoil_discretize.py e423.dat --chord 25.5 --n 40 --spacing cosine
    (chord in whatever units you want the CAD in, e.g. inches; default 25.5)

Output coordinates: the leading edge (the minimum-x point of the spline) is
the origin (0, 0), and the outline is scaled so the chord, measured along x
from the leading edge to the trailing edge, equals --chord. The input file
does not have to be normalized to chord = 1 or start at x = 0.

Two features are created in the same run and written to the same DXF:
  1. the outer outline, as N SPLINE entities (layer SEGMENTS)
  2. the inner wall, offset --wall (default 0.3) inward from the outline, as
     ONE continuous closed polyline (layer INNER)
so the imported profile is a shell of constant wall thickness. Where the
airfoil is thinner than twice the wall (near the trailing edge) the upper and
lower offsets meet; the inner polyline ends in a corner there and the section
aft of it is solid. Use --wall 0 to skip the inner wall.

Outputs (in --out directory):
    segments.csv        per-segment: i, s0, s1, length, x_mid, y_mid, nx, ny, tx, ty
    inner_polyline.csv  x, y vertices of the inner wall (closed: last -> first)
    outline.dxf         N SPLINE entities + 1 closed LWPOLYLINE (inner wall)
    outline.png         quick plot to eyeball the discretization
"""
import argparse
import os
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar


# ----------------------------------------------------------------- reading
def read_dat(path):
    """Return Nx2 array ordered TE(upper) -> LE -> TE(lower) (Selig order)."""
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.replace(",", " ").split()
            try:
                vals = [float(p) for p in parts]
            except ValueError:
                continue  # header / name line
            if len(vals) == 2:
                rows.append(vals)
    pts = np.array(rows)

    # Lednicer format: first numeric row is point counts (e.g. "31. 31."),
    # then upper LE->TE, then lower LE->TE.
    if pts[0, 0] > 1.5 and pts[0, 1] > 1.5:
        nu, nl = int(round(pts[0, 0])), int(round(pts[0, 1]))
        upper = pts[1:1 + nu]            # LE -> TE
        lower = pts[1 + nu:1 + nu + nl]  # LE -> TE
        pts = np.vstack([upper[::-1], lower[1:]])  # TE(up) -> LE -> TE(low)
    return pts


# ---------------------------------------------------------------- geometry
def build_spline(pts, close_te=True):
    """
    One continuous cubic spline x(s), y(s) parametrised by cumulative
    chord length s through the data points. If the trailing edge is blunt
    (first point != last point) and close_te is True, a straight TE segment
    is appended as a separate closing line (not part of the spline).
    """
    d = np.hypot(*np.diff(pts, axis=0).T)
    s = np.concatenate([[0], np.cumsum(d)])
    # 'natural' end conditions: no spurious end slopes at the trailing edge
    cs = CubicSpline(s, pts, bc_type="natural")
    return cs, s[-1]


def arclength_table(cs, s_max, m=20000):
    """Dense table so we can invert true arc length -> parameter s."""
    s = np.linspace(0, s_max, m)
    xy = cs(s)
    ds = np.hypot(*np.diff(xy, axis=0).T)
    L = np.concatenate([[0], np.cumsum(ds)])
    return s, L


def find_le(cs, s_max, s_tab=None):
    """
    Spline parameter of the leading edge = the minimum-x point of the spline
    (not just the nearest data point): coarse search on a dense table, then
    a bounded 1-D minimisation between the neighbouring table points.
    """
    if s_tab is None:
        s_tab = np.linspace(0, s_max, 20000)
    i = int(np.argmin(cs(s_tab)[:, 0]))
    lo, hi = s_tab[max(i - 1, 0)], s_tab[min(i + 1, len(s_tab) - 1)]
    res = minimize_scalar(lambda s: cs(s)[0], bounds=(lo, hi), method="bounded",
                          options={"xatol": 1e-13 * max(s_max, 1.0)})
    return float(res.x) if cs(res.x)[0] <= cs(s_tab[i])[0] else float(s_tab[i])


def place_on_chord(pts, chord):
    """
    Move the leading edge to the origin and scale uniformly so the chord
    (x distance from the leading edge to the trailing edge) equals `chord`.
    The trailing edge is the midpoint of the first and last points, so a
    blunt trailing edge is handled too. No rotation is applied.
    """
    cs, s_max = build_spline(pts)
    le = cs(find_le(cs, s_max))
    te = 0.5 * (pts[0] + pts[-1])
    c0 = te[0] - le[0]
    if not c0 > 0:
        raise SystemExit("[ERROR] trailing edge is not downstream (+x) of the leading "
                         "edge; the file must run TE -> upper -> LE -> lower -> TE.")
    return (pts - le) * (chord / c0), c0


def breakpoints(cs, s_max, n, spacing="cosine"):
    """
    n+1 breakpoints along the single spline, returned as parameter values.
    A breakpoint is always placed exactly at the leading edge (min x), so
    no segment straddles the LE and upper/lower surfaces stay separate.
    uniform : equal true arc length on each surface
    cosine  : arc-length cosine spacing on EACH surface, so segments are
              smallest at the leading edge and trailing edge, where the
              pressure gradients are steepest.
    """
    s_tab, L_tab = arclength_table(cs, s_max)
    s_le = find_le(cs, s_max, s_tab)
    L_le, Ltot = np.interp(s_le, s_tab, L_tab), L_tab[-1]
    n_up = n // 2
    n_lo = n - n_up

    def spaced(m):
        u = np.linspace(0, 1, m + 1)
        return 0.5 * (1 - np.cos(np.pi * u)) if spacing == "cosine" else u

    # upper surface: TE -> LE (arc length 0 .. L_le); lower: LE -> TE
    L_up = spaced(n_up) * L_le
    L_lo = L_le + spaced(n_lo) * (Ltot - L_le)
    target = np.concatenate([L_up, L_lo[1:]])
    s_edges = np.interp(target, L_tab, s_tab)
    s_edges[n_up] = s_le                 # breakpoint exactly on the leading edge
    return s_edges, Ltot


def segment_data(cs, s_edges, npts=25):
    """Per-segment sample points, midpoint, unit tangent and outward normal."""
    segs = []
    for i in range(len(s_edges) - 1):
        s = np.linspace(s_edges[i], s_edges[i + 1], npts)
        xy = cs(s)
        ds = np.hypot(*np.diff(xy, axis=0).T)
        length = ds.sum()
        sm = 0.5 * (s_edges[i] + s_edges[i + 1])
        mid = cs(sm)
        t = cs(sm, 1)
        t = t / np.hypot(*t)
        # Selig order is TE-upper -> LE -> TE-lower, i.e. counter-clockwise
        # for an airfoil with chord along +x. For CCW traversal the outward
        # normal is the tangent rotated -90 deg: (ty, -tx).
        n = np.array([t[1], -t[0]])
        segs.append(dict(i=i, s0=s_edges[i], s1=s_edges[i + 1], length=length,
                         xy=xy, mid=mid, t=t, n=n))
    return segs


# -------------------------------------------------------------- inner wall
def _dist_to_polyline(Q, A, B, chunk=500):
    """Shortest distance from each point in Q to the segments A[i] -> B[i]."""
    AB = B - A
    den = np.maximum((AB ** 2).sum(1), 1e-300)
    out = np.empty(len(Q))
    for k in range(0, len(Q), chunk):
        q = Q[k:k + chunk, None, :]
        u = np.clip(((q - A) * AB).sum(2) / den, 0.0, 1.0)
        out[k:k + chunk] = np.sqrt((((A + u[..., None] * AB) - q) ** 2).sum(2)).min(1)
    return out


def _seg_intersection(p1, p2, p3, p4):
    """Intersection point of segments p1-p2 and p3-p4, or None."""
    r, q = p2 - p1, p4 - p3
    den = r[0] * q[1] - r[1] * q[0]
    if abs(den) < 1e-300:
        return None
    w = p3 - p1
    a = (w[0] * q[1] - w[1] * q[0]) / den
    b = (w[0] * r[1] - w[1] * r[0]) / den
    if -1e-9 <= a <= 1 + 1e-9 and -1e-9 <= b <= 1 + 1e-9:
        return p1 + a * r
    return None


def inner_wall(cs, s_max, pts, wall, n_inner=200, m=4000):
    """
    Inner wall = the outline offset `wall` inward, as ONE closed polyline.

    1. Sample the outer spline densely and step each point `wall` along the
       inward normal (raw offset curve).
    2. Drop every raw point that is closer than `wall` to ANY part of the
       outline (including the closing line of a blunt trailing edge). This
       removes the self-crossing loop near the trailing edge, and anywhere
       else the outline is too thin or too tightly curved for the wall.
    3. Where a stretch was dropped, join the two remaining ends at the point
       where the raw offset crosses itself (a sharp corner).
    4. Resample to n_inner vertices, clustered at the inner leading edge and
       at the trailing-edge corner like the outer segments.

    Returns (vertices Kx2, without repeating the first point; info dict).
    """
    s = np.linspace(0, s_max, m)
    P = cs(s)
    T = cs(s, 1)
    T = T / np.hypot(T[:, 0], T[:, 1])[:, None]
    Q = P - wall * np.column_stack([T[:, 1], -T[:, 0]])   # outward n = (ty, -tx)

    A, B = P[:-1], P[1:]
    if np.hypot(*(pts[0] - pts[-1])) > 1e-9:               # blunt TE closing line
        A, B = np.vstack([A, P[-1]]), np.vstack([B, P[0]])
    valid = _dist_to_polyline(Q, A, B) >= wall * (1 - 1e-3)
    idx = np.flatnonzero(valid)
    if len(idx) < 10:
        raise SystemExit(f"[ERROR] wall = {wall:g} is too thick for this airfoil at this "
                         "chord: no inner wall is left. Reduce --wall.")

    loop, corners = [], []
    kf, kl = idx[0], idx[-1]
    if kf > 0 and kl < m - 1:                              # corner at the trailing edge
        c = _seg_intersection(Q[kl], Q[kl + 1], Q[kf - 1], Q[kf])
        if c is not None:
            loop.append(c)
    for a, b in zip(idx[:-1], idx[1:]):
        loop.append(Q[a])
        if b > a + 1:                                      # dropped stretch mid-outline
            c = _seg_intersection(Q[a], Q[a + 1], Q[b - 1], Q[b])
            if c is not None:
                loop.append(c)
            corners.append(c if c is not None else Q[a])
    loop.append(Q[kl])
    loop = np.array(loop)

    # resample: start vertex (TE end) -> inner LE -> back to start
    closed = np.vstack([loop, loop[:1]])
    L = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(closed, axis=0).T))])
    L_le = L[int(np.argmin(loop[:, 0]))]
    n_up = n_inner // 2
    cosine = lambda k: 0.5 * (1 - np.cos(np.pi * np.linspace(0, 1, k + 1)))
    target = np.concatenate([cosine(n_up) * L_le,
                             (L_le + cosine(n_inner - n_up) * (L[-1] - L_le))[1:-1]])
    V = np.column_stack([np.interp(target, L, closed[:, 0]),
                         np.interp(target, L, closed[:, 1])])
    info = dict(x_te=float(loop[0, 0]), x_le=float(loop[:, 0].min()),
                length=float(L[-1]), mid_corners=corners,
                min_gap=float(_dist_to_polyline(V, A, B).min()))
    return V, info


# ---------------------------------------------------------------- pressure
def average_pressure(cs, segs, p_of_s, nint=200):
    """
    p_bar_i = (1/L_i) * integral p(s) ds over segment i, by trapezoid rule.
    p_of_s : callable giving pressure at spline parameter s (your own Cp*q
             interpolation goes here).
    """
    pbar = []
    for sg in segs:
        s = np.linspace(sg["s0"], sg["s1"], nint)
        xy = cs(s)
        ds = np.hypot(*np.diff(xy, axis=0).T)
        p = p_of_s(s)
        integ = np.sum(0.5 * (p[1:] + p[:-1]) * ds)
        pbar.append(integ / ds.sum())
    return np.array(pbar)


def resultant_force(segs, pbar):
    """
    Force on the body per unit span from pressure acting INTO the surface,
    with n the outward normal:  F_i = -p_i * INTEGRAL(n ds) over segment i.
    For an arc, INTEGRAL(n ds) = (dy, -dx) where (dx, dy) = end - start point,
    so this is exact for constant pressure on a curved segment, which is
    what Mechanical sees when you give each arc a constant 'Normal To'
    pressure. Use it to check the total against your expected lift/drag.
    """
    F = np.zeros(2)
    for sg, p in zip(segs, pbar):
        dx, dy = sg["xy"][-1] - sg["xy"][0]
        F += -p * np.array([dy, -dx])
    return F


# ----------------------------------------------------------------- output
def write_csv(segs, path):
    with open(path, "w") as f:
        f.write("i,s0,s1,length,x_mid,y_mid,nx,ny,tx,ty\n")
        for sg in segs:
            f.write(f"{sg['i']},{sg['s0']:.6f},{sg['s1']:.6f},{sg['length']:.6f},"
                    f"{sg['mid'][0]:.6f},{sg['mid'][1]:.6f},"
                    f"{sg['n'][0]:.6f},{sg['n'][1]:.6f},"
                    f"{sg['t'][0]:.6f},{sg['t'][1]:.6f}\n")


def write_inner_csv(inner, path):
    with open(path, "w") as f:
        f.write("i,x,y\n")
        for i, (x, y) in enumerate(inner):
            f.write(f"{i},{x:.6f},{y:.6f}\n")


def write_dxf(segs, pts_scaled, path, units="in", close_te=True, inner=None):
    import ezdxf
    doc = ezdxf.new("R2010")
    doc.units = {"in": 1, "mm": 4}[units]  # DXF $INSUNITS
    msp = doc.modelspace()
    doc.layers.add("SEGMENTS")
    for sg in segs:
        msp.add_spline(fit_points=[(x, y, 0) for x, y in sg["xy"]],
                       dxfattribs={"layer": "SEGMENTS"})
    if close_te and np.hypot(*(pts_scaled[0] - pts_scaled[-1])) > 1e-9:
        msp.add_line((*pts_scaled[-1], 0), (*pts_scaled[0], 0),
                     dxfattribs={"layer": "SEGMENTS"})
    if inner is not None:
        # inner wall: ONE closed polyline entity
        doc.layers.add("INNER")
        msp.add_lwpolyline([(x, y) for x, y in inner], close=True,
                           dxfattribs={"layer": "INNER"})
    doc.saveas(path)


def plot(segs, pts_scaled, path, inner=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 3.2))
    ax.plot(*pts_scaled.T, "k.", ms=2, label="data points")
    for sg in segs:
        ax.plot(*sg["xy"].T, "-", lw=1.2)
        ax.plot(*sg["xy"][0], "r|", ms=8)
        ax.arrow(*sg["mid"], *(sg["n"] * 0.03 * np.ptp(pts_scaled[:, 0])),
                 head_width=0.004 * np.ptp(pts_scaled[:, 0]), color="g")
    title = f"{len(segs)} segments (red ticks = breakpoints, green = outward normals)"
    if inner is not None:
        ax.plot(*np.vstack([inner, inner[:1]]).T, "k-", lw=0.9)
        title += f" + inner wall polyline ({len(inner)} vertices)"
    ax.set_aspect("equal")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=150)


# -------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dat")
    ap.add_argument("--chord", type=float, default=25.5,
                    help="chord length (LE to TE) in CAD units; default 25.5")
    ap.add_argument("--n", type=int, default=40, help="number of segments")
    ap.add_argument("--spacing", choices=["uniform", "cosine"], default="cosine")
    ap.add_argument("--wall", type=float, default=0.3,
                    help="wall thickness in CAD units: the inner polyline is offset this far "
                         "inward from the outline; default 0.3, 0 = no inner wall")
    ap.add_argument("--inner-n", type=int, default=200,
                    help="number of vertices in the inner wall polyline; default 200")
    ap.add_argument("--units", choices=["in", "mm"], default="in")
    ap.add_argument("--out", default="out")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    if not a.chord > 0:
        raise SystemExit("[ERROR] --chord must be a positive length.")
    if a.n < 2:
        raise SystemExit("[ERROR] --n must be at least 2 (one segment per surface).")
    if a.wall < 0:
        raise SystemExit("[ERROR] --wall must be positive (or 0 for no inner wall).")
    if a.wall > 0 and a.inner_n < 8:
        raise SystemExit("[ERROR] --inner-n must be at least 8.")

    # leading edge -> origin, chord (LE to TE along x) -> a.chord
    pts, c0 = place_on_chord(read_dat(a.dat), a.chord)
    cs, s_max = build_spline(pts)
    s_edges, Ltot = breakpoints(cs, s_max, a.n, a.spacing)
    segs = segment_data(cs, s_edges)

    inner = None
    if a.wall > 0:
        inner, info = inner_wall(cs, s_max, pts, a.wall, a.inner_n)
        write_inner_csv(inner, os.path.join(a.out, "inner_polyline.csv"))

    write_csv(segs, os.path.join(a.out, "segments.csv"))
    write_dxf(segs, pts, os.path.join(a.out, "outline.dxf"), units=a.units, inner=inner)
    plot(segs, pts, os.path.join(a.out, "outline.png"), inner=inner)

    L = np.array([s["length"] for s in segs])
    le = cs(s_edges[a.n // 2])
    te = 0.5 * (pts[0] + pts[-1])
    print(f"input chord: {c0:.5f}  -> scaled by {a.chord / c0:.5f} to chord {a.chord:g} {a.units}")
    print(f"LE at ({le[0]:.6f}, {le[1]:.6f})  | TE at ({te[0]:.6f}, {te[1]:.6f})")
    print(f"points read: {len(pts)}  | outline arc length: {Ltot:.4f}")
    print(f"segments: {a.n}  | length min/max: {L.min():.4f} / {L.max():.4f}")
    print(f"TE gap (first-last point): {np.hypot(*(pts[0]-pts[-1])):.5f}")
    if inner is not None:
        print(f"inner wall: {a.wall:g} {a.units} inside the outline, one closed polyline of "
              f"{len(inner)} vertices, length {info['length']:.4f}")
        print(f"inner wall spans x = {info['x_le']:.4f} to {info['x_te']:.4f} "
              f"(solid aft of that) | thinnest wall at a vertex: {info['min_gap']:.4f}")
        for c in info["mid_corners"]:
            print(f"[NOTE] the outline curves tighter than the {a.wall:g} wall near "
                  f"({c[0]:.3f}, {c[1]:.3f}); the inner wall has a sharp corner there instead of "
                  "a rounded offset (the wall is slightly thicker than nominal at that spot).")


if __name__ == "__main__":
    main()
