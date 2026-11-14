#!/usr/bin/env python3
"""
airfoil_discretize.py
Read an airfoil .dat (Selig or Lednicer format), fit ONE continuous spline
through the whole outline, split it into N arcs of controlled length, and
export per-segment geometry (points, tangents, outward normals, lengths)
plus a DXF with each arc as its own SPLINE entity.

Usage:
    python3 airfoil_discretize.py e423.dat --chord 6.0 --n 40 --spacing cosine
    (chord in whatever units you want the CAD in, e.g. inches)

Outputs (in --out directory):
    segments.csv   per-segment: i, s0, s1, length, x_mid, y_mid, nx, ny, tx, ty
    outline.dxf    N SPLINE entities (one per segment) + optional full outline
    outline.png    quick plot to eyeball the discretization
"""
import argparse
import os
import numpy as np
from scipy.interpolate import CubicSpline


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
    xy = cs(s_tab)
    i_le = int(np.argmin(xy[:, 0]))
    L_le, Ltot = L_tab[i_le], L_tab[-1]
    n_up = n // 2
    n_lo = n - n_up

    def spaced(m):
        u = np.linspace(0, 1, m + 1)
        return 0.5 * (1 - np.cos(np.pi * u)) if spacing == "cosine" else u

    # upper surface: TE -> LE (arc length 0 .. L_le); lower: LE -> TE
    L_up = spaced(n_up) * L_le
    L_lo = L_le + spaced(n_lo) * (Ltot - L_le)
    target = np.concatenate([L_up, L_lo[1:]])
    return np.interp(target, L_tab, s_tab), Ltot


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


def write_dxf(segs, pts_scaled, path, units="in", close_te=True):
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
    doc.saveas(path)


def plot(segs, pts_scaled, path):
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
    ax.set_aspect("equal")
    ax.set_title(f"{len(segs)} segments (red ticks = breakpoints, green = outward normals)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)


# -------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dat")
    ap.add_argument("--chord", type=float, default=1.0, help="chord length in CAD units")
    ap.add_argument("--n", type=int, default=40, help="number of segments")
    ap.add_argument("--spacing", choices=["uniform", "cosine"], default="cosine")
    ap.add_argument("--units", choices=["in", "mm"], default="in")
    ap.add_argument("--out", default="out")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    pts = read_dat(a.dat) * a.chord          # scale normalized coords by chord
    cs, s_max = build_spline(pts)
    s_edges, Ltot = breakpoints(cs, s_max, a.n, a.spacing)
    segs = segment_data(cs, s_edges)

    write_csv(segs, os.path.join(a.out, "segments.csv"))
    write_dxf(segs, pts, os.path.join(a.out, "outline.dxf"), units=a.units)
    plot(segs, pts, os.path.join(a.out, "outline.png"))

    L = np.array([s["length"] for s in segs])
    print(f"points read: {len(pts)}  | outline arc length: {Ltot:.4f}")
    print(f"segments: {a.n}  | length min/max: {L.min():.4f} / {L.max():.4f}")
    print(f"TE gap (first-last point): {np.hypot(*(pts[0]-pts[-1])):.5f}")


if __name__ == "__main__":
    main()
