"""Cortex pattern figures and profile features, called by summarise.py (and live by review.py).

The cortex pattern (sectors, 2026-11-02): every cortical cell in the kept angles (counted by its centre), x = r_norm
(0 = stele edge, 1 = root edge along the cell's own direction), N_BINS equal bins, MEDIAN area per bin; the epidermis
is the curve's last point, one bin past the cortex (decided 2026-11-03: the pattern includes it; it is still
tagged, and samples.xlsx says whether it was intact). Replaced one spoke per quarter, which
on R1 was far noisier (split-half r 0.52 vs 0.87; 04_root_anatomy 12_sector_vs_spoke).
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

GCOL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]   # validated CVD-safe; every line is labelled
INK, MUTED, GRID = "#222222", "#8a8a85", "#ebebe8"
N_BINS = 10
MIN_PTS = 5       # fewer points than this: no spline
EPI_X = 1.05     # the epidermis row's x (1 + half a bin of 10): the last point of the curve
DECAY_WIN = 0.2   # decay slope measured over this fraction of the x-range just past the peak
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


def profile_features(x, y, w=None, peak_max_x=None) -> dict:
    """Plan step 2: features read off the fitted spline (NaN where undefined).
    peak_max_x: the peak is searched only at x <= this (1.0 when the epidermis is the curve's last point, so the
    peak is always a cortex peak and the epidermis shows up as the outer rise, the lip)
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
    zone = xs <= peak_max_x + 1e-9 if peak_max_x is not None else np.ones(len(xs), bool)
    if not zone.any():
        return out
    i = int(np.argmax(np.where(zone, ys, -np.inf)))
    out.update(peak_height=ys[i], peak_pos=xs[i])
    if i > 0:
        out["rise_slope"] = (ys[i] - ys[0]) / (xs[i] - xs[0])
    if i < len(xs) - 1:
        j = min(np.searchsorted(xs, xs[i] + DECAY_WIN * (xs[-1] - xs[0])), len(xs) - 1)
        out["decay_slope"] = (ys[j] - ys[i]) / (xs[j] - xs[i])
        out["outer_rise"] = ys[-1] - ys[i:].min()
    return out


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
    ax.set_xlim(0, 1.22)
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


def sector_pattern(ax, prof, y_col, y_lab, color=INK, label=None, compact=False):
    """One image's cortex pattern on ax: bin medians (dots; open = few cells), the epidermis (square) as the last
    point, and the spline through all of them (the curve the main features come from).
    prof: pipeline.sector_profile rows (+ y_col). Also drawn live by review.py and the web page."""
    d = prof.dropna(subset=[y_col])
    cx = d[d.bin != "epidermis"]
    ok = cx[~cx.few_cells]
    ax.plot(ok.x, ok[y_col], lw=0, marker="o", ms=5, color=color, alpha=0.8)
    few = cx[cx.few_cells]
    ax.plot(few.x, few[y_col], lw=0, marker="o", ms=5, mfc="none", mec=color, alpha=0.8)
    sp = spline(d.x, d[y_col], d.n_cells.where(d.n_cells > 0))
    if sp is not None:
        ax.plot(*sp, color=color, lw=2.4, solid_capstyle="round", label=label)
    e = d[d.bin == "epidermis"]
    if len(e):                                  # open square: epidermis on too few degrees (see epidermis_deg)
        few_e = bool(e.few_cells.iloc[0]) if "few_cells" in e else False
        ax.plot(e.x, e[y_col], lw=0, marker="s", ms=7, color=color, mfc="none" if few_e else color)
    _axes(ax, "" if compact else "Cortex, stele → epidermis", "Stele → epidermis (0-1)" if compact else
          "Radial position (0 = stele edge, 1 = root edge)", y_lab)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1, EPI_X + 0.06])
    ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1", "epi"])


# ---------------------------------------------------------------- figures
def image_figure(prof, y_lab, title, out):
    """One image's cortex pattern (prof: sector_profile with a y column)."""
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    sector_pattern(ax, prof, "y", y_lab)
    fig.suptitle(title, x=0.01, ha="left", fontsize=10, color=MUTED)
    fig.text(0.01, 0.005, "Dots: median per bin (open: fewer than 5 cells). Square: epidermis, the curve's last point.",
             fontsize=8, color=MUTED)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    fig.savefig(out, dpi=200)
    plt.close(fig)


def group_figure(pip, group_col, species, y_col, y_lab, out, note=""):
    """pip: per_image_profile rows of one species. Plant = replicate: its images are averaged per bin first.
    Thin = one spline per plant; bold = spline through the plant means (weighted by plants), band = ± 1 SE across
    plants, bins with >= 3 plants; square ± SE = epidermis, its own point."""
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    groups = sorted(pip[group_col].dropna().unique(), key=str)
    plant = pip.groupby([group_col, "plant", "bin", "x_position"])[y_col].mean().reset_index()
    ends = []
    for g, d in plant.groupby(group_col):
        c = GCOL[groups.index(g) % len(GCOL)]
        cx, ep = d, d[d.bin == "epidermis"]                      # the curve runs through the epidermis point
        for _, p in cx.groupby("plant"):
            sp = spline(p.x_position, p[y_col])
            if sp is not None:
                ax.plot(*sp, color=c, lw=0.6, alpha=0.22)
        m = cx.groupby("x_position")[y_col].agg(["mean", "std", "count"]).reset_index()
        m = m[m["count"] >= 3]
        if len(m):
            se = m["std"] / np.sqrt(m["count"])
            lo = spline(m.x_position, m["mean"] - se, m["count"])
            hi = spline(m.x_position, m["mean"] + se, m["count"])
            if lo is not None and hi is not None:
                ax.fill_between(lo[0], lo[1], hi[1], color=c, alpha=0.15, lw=0)
            e = _curve(ax, m.x_position, m["mean"], m["count"], c, lw=2.8, ms=5)
            if e:
                ends.append((*e, f"{g}  (n = {d.plant.nunique()})"))
        em = ep[y_col].dropna()
        if len(em):
            ax.errorbar(EPI_X + 0.008 * groups.index(g), em.mean(),
                        yerr=em.std() / np.sqrt(len(em)) if len(em) > 1 else None, fmt="s", ms=7, color=c, capsize=3)
    _axes(ax, "", "Radial position (0 = stele edge, 1 = root edge)", y_lab)
    ax.set_xlim(0, 1.3)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1, EPI_X + 0.06])
    ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1", "epi"])
    _end_labels(ax, ends, 10)
    fig.suptitle(species, x=0.01, ha="left", fontsize=13, fontweight="bold", color=INK)
    fig.text(0.01, 0.005, "Bold: spline through plant means (dots), band ± 1 SE across plants (bins with ≥ 3 plants). "
             f"Thin: one plant. Square: epidermis (the curve's last point).{note}", fontsize=8, color=MUTED)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(out, dpi=200)
    plt.close(fig)
