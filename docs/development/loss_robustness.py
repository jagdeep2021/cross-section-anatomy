"""Which way of placing cells survives missing tissue? (decision log 2026-11-04)

One synthetic root (tests/check_numbers.py: 8 cortex layers x 36 cells + 72 epidermal cells, stele r = 150 px),
damaged five ways, each cell's position computed by every method. The cells that survive did not move, so the right
answer is always "no shift". Methods:
  position          Rooh: r_norm, stele edge (0) to root outline (1) at the cell's own angle
  file              Rooh: file_rel, cells counted along rays from the stele (before 2026-11-04)
  file, gap-aware   Rooh: the same, a gap on the ray counted as the cells it would have held (since 2026-11-04)
  file number       the gap-aware file number itself (1, 2 ... from the stele) / 9, the intact file count
  distance          Rooh: distance from the stele edge, as a fraction of the intact cortex width (no outer anchor)
  Tyler min-max     Tyler's maize-cells: distance from the centre, min-max scaled over the cells of each quarter
                    (given the stele centre here; his crop corner would only add error)
Shift = |damaged - intact| for each surviving cell, mean over cells (a median hides damage on a third of
the root). Peak shift = the pattern's peak position
(10 bins, median area, the GCV spline: the same downstream for every method) damaged minus intact.

  .venv/bin/python docs/development/loss_robustness.py     -> loss_robustness.png + loss_robustness.csv beside it
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path[:0] = [str(ROOT), str(ROOT / "tests")]
import pattern_plots as PP          # noqa: E402
import pipeline as P                # noqa: E402
import check_numbers as T           # noqa: E402

TISSUE = list(P.PATTERN_TISSUES)
WIDTH0 = T.R_ROOT - T.R_STELE       # intact cortex + epidermis width (px)


def damage(case):
    cells, region, mx, gray = T.section()
    RR, TH, E = T.RR, T.TH, T.EDGES
    arcs = []
    if case == "aerenchyma (mid-cortex)":
        m = (RR >= E[3]) & (RR < E[5]) & ((TH // 10) % 3 == 0)          # layers 4-5 gone in every third sector
        cells[m] = 0
    elif case == "inner file missing":
        cells[(RR >= E[0]) & (RR < E[1])] = 0                          # layer 1 all round
    elif case == "outer lost all round":
        m = RR >= E[7]                                                  # layer 8 + epidermis
        cells[m], region[m] = 0, 0
    elif case.startswith("outer lost on one side"):
        m = (RR >= E[7]) & (TH < 60)
        cells[m], region[m] = 0, 0
        if "excluded" in case:
            arcs = [(0.0, 60.0, "torn")]
    return cells, region, mx, arcs


def positions(cells, region, mx, arcs):
    centre = P.centre_of(region, None)
    props, _, _ = P.derive(cells, region, mx, 0, centre)
    props = P.exclude_sectors(props, arcs)
    c = props[props.tissue.isin(TISSUE)].copy()
    out = pd.DataFrame({"x": c.x.round().values, "y": c.y.round().values, "area": c.area_px.values,
                        "position": c.r_norm.values, "distance": c.dist_from_stele_px.values / WIDTH0})
    for name, g in (("file", False), ("file, gap-aware", True)):
        fl = P.cell_files(cells, props, centre, arcs, region=region, gaps=g)
        out[name] = c.label.map(fl.file_rel).values
        out[name.replace("file", "file number")] = c.label.map(fl.file).values / 9   # 9 files intact: same axis
    q = P.quarter_of(np.degrees(c.angle_rad.values))
    r = c.r_px.values
    ty = np.full(len(c), np.nan)
    for k in (1, 2, 3, 4):
        s = q == k
        if s.sum() > 1:
            ty[s] = (r[s] - r[s].min()) / (r[s].max() - r[s].min())
    out["Tyler min-max"] = ty
    return out


def peak(x, a):
    ok = np.isfinite(x)
    x, a = x[ok], a[ok]
    b = np.clip((np.clip(x, 0, 0.999) * 10).astype(int), 0, 9)
    g = pd.Series(a).groupby(b).agg(["median", "size"])
    xs, ys = (g.index.values + 0.5) / 10, g["median"].values / np.median(a)
    return PP.profile_features(xs, ys, g["size"].values)["peak_pos"]


def main():
    methods = ["position", "file", "file, gap-aware", "file number, gap-aware", "distance", "Tyler min-max"]
    cases = ["aerenchyma (mid-cortex)", "inner file missing", "outer lost all round", "outer lost on one side",
             "outer lost on one side, arc excluded"]
    c0, r0, m0, _ = damage("intact")
    base = positions(c0, r0, m0, [])
    rows = []
    for case in cases:
        d = positions(*damage(case))
        j = d.merge(base, on=["x", "y"], suffixes=("", "_0"))
        for m in methods:
            rows.append({"case": case, "method": m, "cells": len(j),
                         "mean_shift": float(np.nanmean(np.abs(j[m] - j[m + "_0"]))),
                         "max_shift": float(np.nanmax(np.abs(j[m] - j[m + "_0"]))),
                         "peak_shift": peak(d[m].values, d.area.values) - peak(base[m].values, base.area.values)})
    t = pd.DataFrame(rows)
    t.to_csv(HERE / "loss_robustness.csv", index=False)
    print(t.pivot(index="case", columns="method", values="mean_shift").reindex(cases)[methods].round(3)
          .to_string())
    print(t.pivot(index="case", columns="method", values="peak_shift").reindex(cases)[methods].round(3).to_string())

    fig, axs = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
    for ax, val, title in ((axs[0], "mean_shift", "A  How far surviving cells move (mean, 0-1 axis)"),
                           (axs[1], "peak_shift", "B  How far the pattern's peak moves (0-1 axis)")):
        m = t.pivot(index="case", columns="method", values=val).reindex(cases)[methods]
        v = np.abs(m.values)
        ax.imshow(v, cmap="Reds", vmin=0, vmax=0.15, aspect="auto")
        for i in range(v.shape[0]):
            for k in range(v.shape[1]):
                ax.text(k, i, "-" if np.isnan(m.values[i, k]) else f"{m.values[i, k]:+.2f}" if val == "peak_shift"
                        else f"{m.values[i, k]:.2f}", ha="center", va="center", fontsize=9,
                        color="white" if v[i, k] > 0.09 else "#222222")
        ax.set_xticks(range(len(methods)), methods, rotation=25, ha="right")
        ax.set_yticks(range(len(cases)), cases)
        ax.set_title(title, loc="left", fontsize=10.5)
    fig.suptitle("Missing tissue vs the way cells are placed: 0 = the method is not fooled (synthetic root)",
                 x=0.01, ha="left", fontsize=11.5)
    fig.tight_layout()
    fig.savefig(HERE / "loss_robustness.png", dpi=170)
    print(HERE / "loss_robustness.png")


if __name__ == "__main__":
    main()
