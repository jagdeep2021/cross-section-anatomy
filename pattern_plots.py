"""Cortex pattern figures and profile features, called by summarise.py (and live by review.py).

Two views of cortical cell size, both on a 0-1 axis from the stele outward
(plan: 5_Projects/04_root_anatomy/01_Cortex patterning/6_cortical_radial_profile_analysis/):
  spoke   the path cells of each quarter's spoke, x = relative file position (file - 0.5) / files in the quarter
  radial  every cortical cell, x = r_norm, radial position (0 = stele edge, 1 = root edge along the cell's own
          direction), cells pooled into N_BINS equal bins
Curves are cubic smoothing splines, smoothness by generalised cross-validation
(scipy.interpolate.make_smoothing_spline): no assumed shape, nothing to tune.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import make_smoothing_spline

QCOL = {1: "#2a78d6", 2: "#eb6834", 3: "#1baf7a", 4: "#eda100"}   # validated CVD-safe; every line is labelled
GCOL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK, MUTED, GRID = "#222222", "#8a8a85", "#ebebe8"
N_BINS = 10
MIN_PTS = 5       # fewer points than this: no spline
DECAY_WIN = 0.2   # decay slope measured over this fraction of the x-range just past the peak
VIEWS = {"spoke": ("rel_position", "Along the spoke", "Relative file position (stele → epidermis)"),
         "radial": ("r_bin", "All cortical cells", "Radial position (stele edge → epidermis)")}
plt.rcParams.update({"font.size": 10, "font.family": "sans-serif", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.titleweight": "bold",
                     "axes.titlesize": 11, "axes.titlelocation": "left", "axes.titlecolor": INK})


# ---------------------------------------------------------------- curves and features
def spline(x, y, w=None, n=100, at=None):
    """GCV smoothing spline through (x, y); repeated x are averaged first (weight = count). None if too few points.
    at: evaluate there instead (NaN outside the data range: no extrapolation)."""
    d = pd.DataFrame({"x": np.asarray(x, float), "y": np.asarray(y, float),
                      "w": 1.0 if w is None else np.asarray(w, float)}).dropna()
    d = d.assign(wy=d.y * d.w).groupby("x")[["wy", "w"]].sum()
    if len(d) < MIN_PTS:
        return None
    xs = np.linspace(d.index.min(), d.index.max(), n) if at is None else np.asarray(at, float)
    try:
        ys = make_smoothing_spline(d.index.values, (d.wy / d.w).values, w=d.w.values)(xs)
    except Exception:          # degenerate input (e.g. all y equal)
        return None
    ys[(xs < d.index.min() - 1e-9) | (xs > d.index.max() + 1e-9)] = np.nan
    return xs, ys


def profile_features(x, y, w=None) -> dict:
    """Plan step 2: features read off the fitted spline (NaN where undefined).
    peak_height / peak_pos   max of the curve and where it sits
    rise_slope               (peak - first value) / (peak_pos - first x); NaN if the peak is the first point
    decay_slope              slope from the peak to DECAY_WIN of the range further out (local, so a late
                             upturn can't bend it); NaN if the peak is the last point
    outer_rise               last value - minimum after the peak; ~0 for a plain hump, > 0 for an uptick"""
    out = dict.fromkeys(["peak_height", "peak_pos", "rise_slope", "decay_slope", "outer_rise"], np.nan)
    sp = spline(x, y, w, n=201)
    if sp is None:
        return out
    xs, ys = sp
    i = int(np.argmax(ys))
    out.update(peak_height=ys[i], peak_pos=xs[i])
    if i > 0:
        out["rise_slope"] = (ys[i] - ys[0]) / (xs[i] - xs[0])
    if i < len(xs) - 1:
        j = min(np.searchsorted(xs, xs[i] + DECAY_WIN * (xs[-1] - xs[0])), len(xs) - 1)
        out["decay_slope"] = (ys[j] - ys[i]) / (xs[j] - xs[i])
        out["outer_rise"] = ys[-1] - ys[i:].min()
    return out


def with_rel(ft):
    if "rel_position" in ft:
        return ft
    return ft.assign(rel_position=(ft.file - 0.5) / ft.groupby("quadrant").file.transform("max"))


def with_bin(cells):
    """Cortical cells only, r_norm in N_BINS equal bins, plotted at the bin centre."""
    c = cells[cells.tissue == "cortex"]
    return c.assign(r_bin=(np.clip(c.r_norm, 0, 0.999) * N_BINS).astype(int) / N_BINS + 0.5 / N_BINS)


# ---------------------------------------------------------------- drawing helpers
def _end_labels(ax, items, fontsize=9):
    """Direct labels at line ends, nudged apart vertically so they never sit on each other."""
    items = sorted(items, key=lambda i: i[1])
    gap = 0.06 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    ys = []
    for _, y, _ in items:
        ys.append(max(y, ys[-1] + gap) if ys else y)
    for (x, _, text), y in zip(items, ys):
        ax.annotate(text, (x, y), xytext=(6, 0), textcoords="offset points", va="center", fontsize=fontsize,
                    color=INK)


def _axes(ax, title, xl, y_lab):
    ax.set_title(title, pad=8)
    ax.set_xlabel(xl)
    ax.set_ylabel(y_lab)
    ax.set_xlim(0, 1.18)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_ylim(bottom=0)


def _curve(ax, x, y, w, color, lw=2.2, dots=True, ms=4):
    if dots:
        ax.plot(x, y, lw=0, marker="o", ms=ms, color=color, alpha=0.75)
    sp = spline(x, y, w)
    if sp is not None:
        ax.plot(*sp, color=color, lw=lw, solid_capstyle="round")
        return sp[0][-1], sp[1][-1]
    return (x.iloc[-1], y.iloc[-1]) if len(x) else None


def quadrant_profiles(axs, ft, cells, y_col, y_lab, fontsize=9, compact=False):
    """One image: per quarter (colour, QC) and pooled over the quarters (black: the image's curve, the one
    the features come from), for the spoke and for all cortical cells. Also drawn live in review.py.
    ft: file_table rows; cells: derive() props (tissue, quadrant, r_norm, y_col)."""
    ft, cb = with_rel(ft), with_bin(cells)
    for ax, (view, d_all) in zip(axs, (("spoke", ft), ("radial", cb))):
        x, title, xl = VIEWS[view]
        ends = []
        for q, d in d_all.groupby("quadrant"):
            if view == "radial":
                d = d.groupby(x)[y_col].agg(["mean", "size"]).reset_index().rename(columns={"mean": y_col})
                e = _curve(ax, d[x], d[y_col], d["size"], QCOL[q], lw=1.1, ms=3)
            else:
                d = d.sort_values(x)
                e = _curve(ax, d[x], d[y_col], None, QCOL[q], lw=1.1, ms=3.5)
            if e:
                ends.append((*e, f"Q{q}"))
        if view == "radial":
            pooled = d_all.groupby(x)[y_col].agg(["mean", "size"])
            sp = spline(pooled.index, pooled["mean"], pooled["size"])
        else:
            sp = spline(d_all[x], d_all[y_col])
        if sp is not None:                      # the image's own curve: all four quarters pooled
            ax.plot(*sp, color=INK, lw=2.6, alpha=0.85)
            ends.append((sp[0][-1], sp[1][-1], "image"))
        _axes(ax, title, "Stele → epidermis (0-1)" if compact else xl, y_lab)
        _end_labels(ax, ends, fontsize)


# ---------------------------------------------------------------- figures
def image_figure(ft, props, y_col, y_lab, title, out):
    """One image's cortex pattern: cell size along the spoke and across all cortical cells, per quarter."""
    fig, bs = plt.subplots(1, 2, figsize=(10, 4.4))
    quadrant_profiles(bs, ft, props, y_col, y_lab)
    fig.suptitle(title, x=0.01, ha="left", fontsize=10, color=MUTED)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=200)
    plt.close(fig)


def group_figure(pf, cells, group_col, species, y_col, y_lab, out, note=""):
    """pf: per_file_path rows, cells: all_cells rows, one species. Plant = replicate: its quarters and images are
    averaged within each x bin first. Thin = one spline per plant; bold = spline through the plant means
    (weighted by plants), band = ± 1 SE across plants; only bins with >= 3 plants."""
    pf = with_rel(pf).assign(rel_bin=lambda d: (np.clip(d.rel_position, 0, 0.999) * N_BINS).astype(int) / N_BINS
                             + 0.5 / N_BINS)
    cb = with_bin(cells)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    groups = sorted(set(pf[group_col].dropna()) | set(cb[group_col].dropna()), key=str)
    for ax, (view, d_all, xb) in zip(axs, (("spoke", pf, "rel_bin"), ("radial", cb, "r_bin"))):
        _, title, xl = VIEWS[view]
        ends = []
        plant = d_all.groupby([group_col, "plant", xb])[y_col].mean().reset_index()
        for g, d in plant.groupby(group_col):
            c = GCOL[groups.index(g) % len(GCOL)]
            for _, p in d.groupby("plant"):
                sp = spline(p[xb], p[y_col])
                if sp is not None:
                    ax.plot(*sp, color=c, lw=0.6, alpha=0.22)
            m = d.groupby(xb)[y_col].agg(["mean", "std", "size"]).reset_index()
            m = m[m["size"] >= 3]
            if not len(m):
                continue
            se = m["std"] / np.sqrt(m["size"])
            lo, hi = spline(m[xb], m["mean"] - se, m["size"]), spline(m[xb], m["mean"] + se, m["size"])
            if lo is not None and hi is not None:
                ax.fill_between(lo[0], lo[1], hi[1], color=c, alpha=0.15, lw=0)
            e = _curve(ax, m[xb], m["mean"], m["size"], c, lw=2.8, ms=5)
            if e:
                ends.append((*e, f"{g}  (n = {d.plant.nunique()})"))
        _axes(ax, title, xl, y_lab if ax is axs[0] else "")
        _end_labels(ax, ends, 10)
    fig.suptitle(species, x=0.01, ha="left", fontsize=13, fontweight="bold", color=INK)
    fig.text(0.01, 0.005, "Bold: spline through plant means (dots), band ± 1 SE across plants. Thin: one plant. "
             f"Bins with ≥ 3 plants.{note}", fontsize=8, color=MUTED)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out, dpi=200)
    plt.close(fig)
