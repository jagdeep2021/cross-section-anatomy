"""Do the numbers come out right? Synthetic cross-sections whose answers are known, run through the same code as
Make results (summarise.main), then edited with the editors' own rules (pipeline.click_edit / outline_edit /
sector_edit / save_edit) and summarised again.

  .venv/bin/python tests/check_numbers.py [out_dir]      -> table of checks (PASS / FAIL), exit 1 on any FAIL

The section: centre (500, 500), stele r < 150 (small cells, 6 metaxylem vessels r = 15 at r = 105), cortex
8 layers x 36 cells (10 degrees each, widths below), epidermis 72 cells to r = 475. Walls drawn in a grey image so
click-to-add works as on a real picture. Expected values are geometry (pi r^2 ...), not what the code gave.
"""
from __future__ import annotations

import json
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from skimage import segmentation

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
import pipeline as P      # noqa: E402
import summarise as SUM   # noqa: E402

N = 1000
C = (500.0, 500.0)
R_STELE, R_ROOT = 150, 475
WIDTHS = [20, 30, 40, 50, 50, 45, 40, 30]                  # cortex layers, stele -> out (sum 305)
EDGES = np.cumsum([R_STELE] + WIDTHS)                       # 150 ... 455
R_EPI = R_ROOT
SECT = 36                                                    # cortex cells per layer
MX_R, MX_AT, MX_ANG = 15, 105, [15, 75, 135, 195, 255, 315]
UM = 0.5

yy, xx = np.mgrid[:N, :N]
RR = np.hypot(yy - C[0], xx - C[1])
TH = np.mod(np.degrees(np.arctan2(-(yy - C[0]), xx - C[1])), 360)


def at(r, deg):
    """pixel (y, x) at radius r, angle deg (counter-clockwise from the right)."""
    t = np.radians(deg)
    return int(round(C[0] - r * np.sin(t))), int(round(C[1] + r * np.cos(t)))


def cortex_label(layer, sector):            # layer 1..8, sector 0..35 (sector k spans 10k .. 10k+10 deg)
    return 1000 + 100 * layer + sector


def section():
    cells = np.zeros((N, N), np.int32)
    # stele: centre disc + rings of small cells
    cells[RR < 50] = 1
    for i, (r0, r1, n) in enumerate([(50, 100, 12), (100, 150, 24)]):
        m = (RR >= r0) & (RR < r1)
        cells[m] = 10 + 100 * i + (TH[m] // (360 / n)).astype(int)
    for k in range(8):
        m = (RR >= EDGES[k]) & (RR < EDGES[k + 1])
        cells[m] = cortex_label(k + 1, 0) + (TH[m] // (360 / SECT)).astype(int)
    m = (RR >= EDGES[-1]) & (RR < R_EPI)
    cells[m] = 5000 + (TH[m] // 5).astype(int)
    mx = np.zeros((N, N), np.int32)
    for i, a in enumerate(MX_ANG):
        y, x = at(MX_AT, a)
        mx[(yy - y) ** 2 + (xx - x) ** 2 <= MX_R ** 2] = i + 1
    cells[mx > 0] = 0
    region = np.where(RR < R_ROOT, 1, 0).astype(np.uint8)
    region[RR < R_STELE] = 2
    # grey picture: dark lumens, bright 2-px walls, dark vessels, black outside
    walls = segmentation.find_boundaries(cells + (mx > 0) * 9000, mode="thick")
    gray = np.where(region > 0, 0.25, 0.0)
    gray[walls & (region > 0)] = 1.0
    gray[mx > 0] = 0.05
    cells, _, _ = segmentation.relabel_sequential(cells)    # ids as Cellpose gives them (keep a map for tests)
    return cells, region, mx, gray


def ann_area(r0, r1, n=1):
    return np.pi * (r1 ** 2 - r0 ** 2) / n


# ------------------------------------------------------------------ expected, from geometry
CELL_AREA = [ann_area(EDGES[k], EDGES[k + 1], SECT) for k in range(8)]
EXPECT = {
    "root_area_px": np.pi * R_ROOT ** 2, "stele_area_px": np.pi * R_STELE ** 2,
    "cortex_area_est_px": np.pi * (R_ROOT ** 2 - R_STELE ** 2), "cortex_width_px": R_ROOT - R_STELE,
    "n_cortex_cells": 8 * SECT, "mean_cortex_cell_area_px": float(np.mean(CELL_AREA)),
    "n_air_spaces": 0, "mx_n": 6, "mx_mean_area_px": np.pi * MX_R ** 2,
    "n_cortical_layers": 8, "angle_kept_deg": 360, "low_coverage": 0,
    "cortex_area_kept_px": np.pi * (R_ROOT ** 2 - R_STELE ** 2), "root_area_est_px": np.pi * R_ROOT ** 2,
    "n_cortex_cells_est": 8 * SECT,
}
MED_CORTEX = float(np.median(np.repeat(CELL_AREA, SECT)))
EPI_AREA0 = np.pi * (475 ** 2 - 455 ** 2) / 72
MED_TISSUE = float(np.median(np.r_[np.repeat(CELL_AREA, SECT), np.repeat(EPI_AREA0, 72)]))   # the reference for rel
EPI_AREA = ann_area(EDGES[-1], R_EPI, 72)
TOL = {"root_area_px": .005, "root_diameter_eq_px": .005, "stele_area_px": .01, "cortex_area_px": .005, "cortex_width_px": .01,
       "mean_cortex_cell_area_px": .02, "mx_mean_area_px": .04, "cortex_area_kept_px": .005,
       "cortex_area_est_px": .005, "root_area_est_px": .005, "stele_share_est": .01, "epidermis_rel": .03,
       "lacuna_area_px": .03, "lacuna_r_norm_mean": .03}


def write_project(d: Path, images: dict):
    """images: name -> (cells, region, mx, gray, extra npz fields). A Rooh project as run_batch leaves it."""
    (d / "anatomy/masks").mkdir(parents=True)
    (d / "raw").mkdir()
    (d / "anatomy/settings.json").write_text(json.dumps({"image_dir": "../raw", "diameter": None,
                                                        "quadrant_angle_deg": 0}))
    rows = []
    for name, (cells, region, mx, gray, extra) in images.items():
        P.save_state(d / f"anatomy/masks/{name}.npz", cells, region, mx, image_file=name + ".tif", **extra)
        rows.append({"image_file": name + ".tif", "plant": "P1", "um_per_px": UM, "include": "yes"})
    pd.DataFrame(rows).to_excel(d / "anatomy/samples.xlsx", index=False)


def results(d: Path) -> dict:
    SUM.main(d / "anatomy")
    out = d / "results"
    img = pd.read_excel(next(out.glob("anatomy_summary_*.xlsx")), sheet_name="per_image")
    prof_all = pd.read_excel(next(out.glob("cortex_pattern_*.xlsx")), sheet_name="per_image_profile", dtype={"bin": str})
    prof = prof_all[(prof_all.basis == "position") & (prof_all.scope == "section")]
    plant = pd.read_excel(next(out.glob("anatomy_summary_*.xlsx")), sheet_name="per_plant")
    cells = pd.read_excel(next(out.glob("all_cells_*.xlsx")), sheet_name="all_cells")
    img["quarters_used"] = img.quarters_used.astype(str)
    feats = pd.read_excel(next(out.glob("cortex_pattern_*.xlsx")), sheet_name="profile_features")
    xl = next(out.glob("anatomy_summary_*.xlsx"))
    lac = pd.read_excel(xl, sheet_name="lacunae")
    aer = pd.read_excel(next(out.glob("cortex_pattern_*.xlsx")), sheet_name="aerenchyma_profile")
    return {"img": img.set_index("image_file"), "prof": prof, "prof_all": prof_all, "plant": plant, "cells": cells,
            "feats": feats, "lac": lac, "aer": aer}


CHECKS: list[dict] = []


def check(case, what, got, want, tol=0.0, note=""):
    if isinstance(want, float) and np.isnan(want):
        ok = pd.isna(got)
    elif tol:
        ok = abs(got - want) <= tol * abs(want)
    else:
        ok = got == want
    CHECKS.append({"case": case, "check": what, "expected": want if not isinstance(want, float) else round(want, 2),
                   "got": got if not isinstance(got, float) else round(got, 2), "ok": "PASS" if ok else "FAIL",
                   "note": note})


def lab(cells, layer, deg):
    """the cell id at the middle of a cortex layer, at angle deg."""
    return int(cells[at((EDGES[layer - 1] + EDGES[layer]) / 2, deg)])


def run(out_dir: Path):
    cells, region, mx, gray = section()
    base = d = out_dir / "00_base"
    write_project(d, {"S": (cells, region, mx, gray, {}),
                      "S_rot": (np.rot90(cells).copy(), np.rot90(region).copy(), np.rot90(mx).copy(), gray, {})})
    r = results(d)
    i = r["img"].loc["S.tif"]
    for k, v in EXPECT.items():
        check("base", k, float(i[k]), float(v), TOL.get(k, 0))
    check("base", "µm²: root_area_um2 = px x 0.5²", float(i.root_area_um2), float(i.root_area_px) * UM ** 2, 1e-9)
    check("base", "µm: cortex_width_um = px x 0.5", float(i.cortex_width_um), float(i.cortex_width_px) * UM, 1e-9)
    check("base", "rotated copy: same cortex cells", int(r["img"].loc["S_rot.tif", "n_cortex_cells"]), 8 * SECT)
    check("base", "plant mean = mean of its 2 images", float(r["plant"].n_cortex_cells.iloc[0]),
          float(r["img"].n_cortex_cells.mean()), 1e-9)
    f = r["prof"][r["prof"].image_file == "S.tif"].set_index("bin")
    check("base", "profile: cortex + epidermis cells over the bins", int(f.n_cells.sum()), 8 * SECT + 72)
    check("base", "profile: epidermis is the outermost bin (10), alone", int(f.n_cells["10"]), 72)
    check("base", "profile: bin 10 median area = epidermis cell", float(f.median_area_px["10"]), EPI_AREA, .03)
    check("base", "epidermis_rel = epidermis / tissue median", float(i.epidermis_rel), EPI_AREA / MED_TISSUE, .03)
    for k in range(8):                                        # layer k's cells land in the bin of their r_norm
        r_mid = (EDGES[k] + EDGES[k + 1]) / 2
        b = str(int(min((r_mid - R_STELE) / (R_ROOT - R_STELE), 0.999) * SUM.N_BINS) + 1)
        check("base", f"layer {k + 1}: median area in bin {b}", float(f.median_area_px[b]), CELL_AREA[k], .03)

    # cell files (2026-11-03): every cortical cell gets the number of its layer
    ac = r["cells"][(r["cells"].image_file == "S.tif") & (r["cells"].tissue == "cortex")]
    layer_true = np.searchsorted(EDGES, ac.r_px.values, side="right")
    check("base", "file number = layer, every cortical cell", float((np.floor(ac.file + 0.5) == layer_true).mean()), 1.0)
    check("base", "every cortical cell has a file number", int(ac.file.notna().sum()), 8 * SECT)
    ae = r["cells"][(r["cells"].image_file == "S.tif") & (r["cells"].tissue == "epidermis")]
    check("base", "every epidermal cell is file 9", float((np.floor(ae.file + 0.5) == 9).mean()), 1.0)
    pf = r["prof_all"][(r["prof_all"].image_file == "S.tif") & (r["prof_all"].basis == "file")
                       & (r["prof_all"].scope == "section")].set_index("bin")
    check("base", "file pattern: 8 cortex files + the epidermis file", len(pf), 9)
    check("base", "file 9 = the epidermis", float(pf.median_area_px["9"]), EPI_AREA, .03)
    for k in range(8):
        check("base", f"file {k + 1}: median area", float(pf.median_area_px[str(k + 1)]), CELL_AREA[k], .03)
    # quarters: all four used, mean over 4
    fe = r["feats"][r["feats"].image_file == "S.tif"]
    check("base", "quarters used", r["img"].loc["S.tif", "quarters_used"], "1,2,3,4")
    qm = fe[(fe.basis == "position") & (fe.scope == "quarter_mean")]
    check("base", "quarter mean over 4 quarters", int(qm.n_quarters.iloc[0]), 4)
    sec = fe[(fe.basis == "position") & (fe.scope == "section")]
    check("base", "symmetric root: quarter-mean peak = section peak", float(qm.peak_pos_rel.iloc[0]),
          float(sec.peak_pos_rel.iloc[0]), .02)
    # diameter, lacunae
    check("base", "root_diameter_px (kept degrees)", float(i.root_diameter_px), 2 * R_ROOT, .01)
    check("base", "root_diameter_eq_px", float(i.root_diameter_eq_px), 2 * R_ROOT, .005)
    check("base", "no lacunae in an intact root", int(i.n_lacunae), 0)
    pd_ = r["prof_all"][(r["prof_all"].basis == "distance") & (r["prof_all"].scope == "section") &
                        (r["prof_all"].image_file == "S.tif")]
    check("base", f"distance basis: {P.DIST_BIN_UM:g} µm bins out to the epidermis (centre 157.5 µm -> bin 8)",
          len(pd_), 8)

    # missing tissue (2026-11-04, docs/development/loss_robustness.py): straight on the pipeline, no summarise
    def files_after(remove, region_cut=None):
        c, rg, m, _ = section()
        c[remove] = 0
        if region_cut is not None:
            c[region_cut], rg[region_cut] = 0, 0
        ctr = P.centre_of(rg, None)
        pr, _, _ = P.derive(c, rg, m, 0, ctr)
        fl = P.cell_files(c, pr, ctr, region=rg)
        pr = pr.join(fl, on="label")
        pr["true_layer"] = np.searchsorted(EDGES, np.hypot(pr.y - C[0], pr.x - C[1]), side="right")
        return pr
    pr = files_after((RR >= EDGES[0]) & (RR < EDGES[1]))
    check("loss", "inner file gone all round: cortex still found (was: all 'unassigned')",
          int((pr.tissue == "cortex").sum()), 7 * SECT)
    k = pr[pr.tissue.isin(["cortex", "epidermis"])]
    check("loss", "inner file gone: every cell keeps its file number (gap at the stele counted)",
          float((np.floor(k.file + 0.5) == k.true_layer).mean()), 1.0)
    pr = files_after((RR >= EDGES[3]) & (RR < EDGES[5]) & ((TH // 10) % 3 == 0))
    k = pr[pr.tissue.isin(["cortex", "epidermis"]) & (pr.true_layer >= 6)]
    check("loss", "aerenchyma (layers 4-5 gone, every 3rd sector): cells outside keep their file number",
          float((np.floor(k.file + 0.5) == k.true_layer).mean()), 1.0)
    pr = files_after(np.zeros_like(RR, bool), RR >= EDGES[7])
    k = pr[pr.tissue.isin(["cortex", "epidermis"])]
    check("loss", "outer layer + epidermis lost: file numbers unchanged (counted from the stele)",
          float((np.floor(k.file + 0.5) == k.true_layer).mean()), 1.0)
    docs = []                                    # output documentation: every column described, run info written
    for wb in sorted((base / "results").glob("*.xlsx")):
        rd = pd.read_excel(wb, sheet_name="README")
        rd = rd[rd.column.notna()]
        docs += [f"{wb.name.split('_20')[0]}:{c}" for c in rd.column[rd["what it is"].fillna("").str.strip() == ""]]
        check("docs", f"{wb.name.split('_20')[0]}: run_info sheet present",
              "run_info" in pd.ExcelFile(wb).sheet_names, True)
    check("docs", "every output column has a description", ", ".join(docs) or "none missing", "none missing")
    check("docs", "root_diameter_eq_um written", "root_diameter_eq_um" in r["img"].columns, True)
    check("docs", "methods text written", any((base / "results").glob("methods_*.txt")), True)

    def edit(case, fn, expect, note=""):
        """Copy the base project, edit image S with fn(npz, cells, region, mx), summarise, compare."""
        d = out_dir / case
        shutil.copytree(base, d)
        shutil.rmtree(d / "results", ignore_errors=True)
        npz = d / "anatomy/masks/S.npz"
        c, rg, m, _ = P.load_state(npz)
        fn(npz, c.copy(), rg.copy(), m.copy())
        i = results(d)["img"].loc["S.tif"]
        for k, v in expect.items():
            if k not in i.index:
                CHECKS.append({"case": case, "check": k, "expected": v, "got": "no such column", "ok": "FAIL",
                               "note": note})
                continue
            got = i[k]
            check(case, k, float(got) if pd.notna(got) else np.nan, float(v), TOL.get(k, 0), note)
        check(case, "status = reviewed", i.status, "reviewed")
        return i

    def sectors(*arcs):
        """Save excluded arcs on S (as the editors do)."""
        return lambda npz, c, rg, m: P.save_edit(npz, **P.sector_arrays(list(arcs)))

    # 1: delete a cortex cell (Cellpose saw a cell that is not there)
    def e1(npz, c, rg, m):
        y, x = at((EDGES[3] + EDGES[4]) / 2, 115)
        P.click_edit(gray, c, rg, m, "Cell", y, x, remove=True)
        P.save_edit(npz, cells=c)
    edit("01_delete_cell", e1, {"n_cortex_cells": 287, "n_cortical_layers": 8}, "one ray loses a layer; median stays 8")

    # 2: delete a cell, then click it back in (Cellpose missed a cell; the reviewer adds it)
    def e2(npz, c, rg, m):
        y, x = at((EDGES[3] + EDGES[4]) / 2, 135)
        P.click_edit(gray, c, rg, m, "Cell", y, x, remove=True)
        P.click_edit(gray, c, rg, m, "Cell", y, x)
        P.save_edit(npz, cells=c)
    edit("02_add_missed_cell", e2, {"n_cortex_cells": 288, "n_cortical_layers": 8})

    # 3: draw one cell over two (Cellpose split one cell; reviewer merges)
    def e3(npz, c, rg, m):
        P.outline_edit(c, rg, m, "Cell", (RR >= EDGES[3]) & (RR < EDGES[5]) & (TH >= 130) & (TH < 140))
        P.save_edit(npz, cells=c)
    edit("03_merge_two_cells", e3, {"n_cortex_cells": 287, "n_cortical_layers": 8})

    # 4: draw a cell over the outer half of one (Cellpose merged two; reviewer splits)
    def e4(npz, c, rg, m):
        mid = (EDGES[4] + EDGES[5]) / 2
        P.outline_edit(c, rg, m, "Cell", (RR >= mid) & (RR < EDGES[5]) & (TH >= 220) & (TH < 230))
        P.save_edit(npz, cells=c)
    edit("04_split_cell", e4, {"n_cortex_cells": 289, "n_cortical_layers": 8})

    # 5: exclude 180-270 (torn cortex): cortex from 270 degrees; stele, vessels, root outline whole; estimates scale
    ca = EXPECT["cortex_area_est_px"]
    edit("05_exclude_quarter_sector", sectors((180, 270, "torn")),
         {"n_cortex_cells": 216, "angle_kept_deg": 270, "n_cortical_layers": 8, "mx_n": 6,
          "stele_area_px": EXPECT["stele_area_px"], "root_area_px": EXPECT["root_area_px"],
          "cortex_width_px": R_ROOT - R_STELE, "cortex_area_kept_px": 0.75 * ca, "cortex_area_est_px": ca,
          "root_area_est_px": EXPECT["root_area_px"], "n_cortex_cells_est": 288,
          "stele_share_est": R_STELE ** 2 / R_ROOT ** 2, "low_coverage": 0,
          "mean_cortex_cell_area_px": EXPECT["mean_cortex_cell_area_px"], "epidermis_rel": EPI_AREA / MED_TISSUE},
         "stele + vessels whole; cortex from 3/4; _est back to the whole")

    # 6: an old project: dropped quarter 3 (before sectors) reads as the same 90 degree arc
    edit("06_old_dropped_quarter", lambda npz, c, rg, m: P.save_edit(npz, dropped_q=np.array([3], np.int16)),
         {"n_cortex_cells": 216, "angle_kept_deg": 270})

    # 6b: edge cells go by their centre: 3-23 degrees holds the centres at 5 and 15 (not 25): 2 per layer out
    edit("06b_edge_cells_by_centre", sectors((3, 23, "torn")), {"n_cortex_cells": 288 - 16, "angle_kept_deg": 340},
         "no half cells: a cell is in or out by its centre")

    # 6c: an arc across 0 degrees (350 -> 10)
    edit("06c_arc_across_zero", sectors((350, 10, "folded")), {"n_cortex_cells": 288 - 16, "angle_kept_deg": 340})

    # 6d: two arcs, overlapping: kept degrees count the overlap once
    edit("06d_overlapping_arcs", sectors((0, 60, "torn"), (40, 90, "debris")),
         {"n_cortex_cells": 216, "angle_kept_deg": 270})

    # 6e: a narrow kept sector: 40 degrees -> low_coverage; estimates still scale
    i = edit("06e_narrow_kept_sector", sectors((40, 360, "smeared / oblique")),
             {"n_cortex_cells": 32, "angle_kept_deg": 40, "low_coverage": 1, "n_cortex_cells_est": 288,
              "cortex_area_est_px": ca, "n_cortical_layers": 8}, f"below {P.MIN_KEPT_DEG} degrees: flagged")

    # 6f: exclude with a drag, then restore with a Ctrl/Cmd click (sector_edit, as the editors call it)
    def e6f(npz, c, rg, m):
        ctr = P.centre_of(rg)
        arcs = P.sector_edit([], ctr, at(300, 100), at(300, 170), "torn")
        arcs = P.sector_edit(arcs, ctr, at(300, 200), at(300, 250), "debris")
        arcs = P.sector_edit(arcs, ctr, at(300, 130), remove=True)        # the first comes back
        P.save_edit(npz, **P.sector_arrays(arcs))
    edit("06f_drag_then_restore", e6f, {"angle_kept_deg": 310, "n_cortex_cells": 288 - 5 * 8},
         "200-250: centres 205 ... 245 out")

    # 6g: excluded epidermis leaves the epidermis point; excluded cells stay in all_cells, tagged
    d = out_dir / "06_old_dropped_quarter"
    r6 = results(d)
    pe = r6["prof"][(r6["prof"].image_file == "S.tif") & (r6["prof"].bin == "10")]
    check("06g_excluded_cells", "outer bin (epidermis): kept cells only", int(pe.n_cells.iloc[0]), 54)
    ac = r6["cells"][r6["cells"].image_file == "S.tif"]
    check("06g_excluded_cells", "all_cells keeps excluded cortex, tagged", int((ac.tissue == "cortex_excluded").sum()), 72)
    check("06g_excluded_cells", "stele never excluded", int(ac.tissue.astype(str).str.startswith("stele_").sum()), 0)

    # 7-8: vessels
    def e7(npz, c, rg, m):
        y, x = at(MX_AT, 135)
        P.click_edit(gray, c, rg, m, "Metaxylem", y, x, remove=True)
        P.save_edit(npz, mx=m)
    edit("07_remove_vessel", e7, {"mx_n": 5})

    def e8(npz, c, rg, m):
        y, x = at(MX_AT, 45)                                 # draw a 7th vessel
        P.outline_edit(c, rg, m, "Metaxylem", (yy - y) ** 2 + (xx - x) ** 2 <= MX_R ** 2)
        P.save_edit(npz, cells=c, mx=m)
    edit("08_add_vessel", e8, {"mx_n": 7, "mx_mean_area_px": np.pi * MX_R ** 2})

    # 9: stele redrawn wider, taking cortex layer 1
    def e9(npz, c, rg, m):
        P.outline_edit(c, rg, m, "Stele", RR < EDGES[1])
        P.save_edit(npz, region=rg)
    edit("09_stele_takes_layer1", e9, {"stele_area_px": np.pi * EDGES[1] ** 2, "n_cortex_cells": 7 * SECT,
                                       "n_cortical_layers": 7, "cortex_width_px": R_ROOT - EDGES[1]})

    # 10: root redrawn inside the epidermis: outer cortex layer becomes the edge
    def e10(npz, c, rg, m):
        P.outline_edit(c, rg, m, "Root", RR < EDGES[-1])
        P.save_edit(npz, region=rg)
    edit("10_root_without_epidermis", e10, {"root_area_px": np.pi * EDGES[-1] ** 2, "n_cortex_cells": 7 * SECT,
                                             "n_cortical_layers": 7})

    # 11: aerenchyma: 4 neighbouring cells of layer 5 gone (10-50 deg) = one lacuna; one layer-6 cell gone (205 deg)
    # = a one-cell gap (a missed cell, review hint), not a lacuna
    def e11(npz, c, rg, m):
        gone = [lab(c, 5, d) for d in (15, 25, 35, 45)] + [lab(c, 6, 205)]
        c[np.isin(c, gone)] = 0
        P.save_edit(npz, cells=c)
    r0, r1 = EDGES[4], EDGES[5]
    r_bar = 2 / 3 * (r1 ** 3 - r0 ** 3) / (r1 ** 2 - r0 ** 2)          # area-weighted mean radius of the annulus
    i11 = edit("11_aerenchyma_lacuna", e11, {"n_cortex_cells": 283, "n_lacunae": 1, "lacuna_area_px": 4 * CELL_AREA[4],
                                             "lacuna_r_norm_mean": (r_bar - R_STELE) / (R_ROOT - R_STELE),
                                             "n_gaps_missed_cell": 1},
               "the one-cell gap is a missed cell, not a lacuna")
    r11 = results(out_dir / "11_aerenchyma_lacuna")
    a = r11["aer"][(r11["aer"].image_file == "S.tif") & (r11["aer"].lacuna_share > 0)]
    check("11_aerenchyma_lacuna", "lacuna sits in bins 5-6 only (r_norm 0.43-0.58)", ",".join(map(str, a.bin)), "5,6")
    lt = r11["lac"][r11["lac"].image_file == "S.tif"]
    check("11_aerenchyma_lacuna", "lacuna inner edge (r_norm)", float(lt.r_norm_inner.iloc[0]),
          (r0 - R_STELE) / (R_ROOT - R_STELE), .03)
    check("11_aerenchyma_lacuna", "lacuna outer edge (r_norm)", float(lt.r_norm_outer.iloc[0]),
          (r1 - R_STELE) / (R_ROOT - R_STELE), .03)
    check("11_aerenchyma_lacuna", "lacuna quarter", int(lt.quarter.iloc[0]), 1)

    # 11d: a tear against the stele (layer-1 cells gone over 60 deg) is not a lacuna
    def e11d(npz, c, rg, m):
        c[np.isin(c, [lab(c, 1, d) for d in range(185, 245, 10)] + [lab(c, 2, d) for d in range(185, 245, 10)])] = 0
        P.save_edit(npz, cells=c)
    edit("11d_tear_at_stele", e11d, {"n_lacunae": 0})

    # 11e: no stele (a root tip, or the stele cut off by the frame): stele / cortex traits blank, never 0
    def e11e(npz, c, rg, m):
        rg[rg == 2] = 1
        m[:] = 0
        P.save_edit(npz, region=rg, mx=m)
    edit("11e_no_stele", e11e, {"stele_found": 0, "n_cortical_layers": np.nan, "stele_area_um2": np.nan,
                                "n_cortex_cells": np.nan, "mx_n": np.nan, "n_lacunae": np.nan,
                                "root_diameter_eq_px": 2 * R_ROOT}, "root size still measured")

    # 11c: big epidermal cells (4 merged into 1): in the main curve (_rel, epidermis part of the tissue) the outer
    # end rises above the cortex peak, so the peak sits at the edge; without the tagged cells (_rel_cortex) it does not
    def e11c(npz, c, rg, m):
        e = (RR >= EDGES[-1]) & (RR < R_EPI) & (c > 0)
        c[e] = c.max() + 1 + (TH[e] // 20).astype(int)
        P.save_edit(npz, cells=c)
    edit("11c_big_epidermis", e11c, {"n_cortex_cells": 288})
    f11c = results(out_dir / "11c_big_epidermis")["feats"]
    f11c = f11c[(f11c.image_file == "S.tif") & (f11c.basis == "position") & (f11c.scope == "section")].iloc[0]
    check("11c_big_epidermis", "main curve: peak at the outer edge (peak_pos_rel > 0.9)", bool(f11c.peak_pos_rel > 0.9), True)
    check("11c_big_epidermis", "without the epidermis: peak inside (peak_pos_rel_cortex < 0.8)",
          bool(f11c.peak_pos_rel_cortex < 0.8), True)

    # 11b: a quarter with < 45 deg kept is not used for quarter features
    i11b = edit("11b_quarter_mostly_excluded", sectors((0, 50, "torn")), {"angle_kept_deg": 310})
    check("11b_quarter_mostly_excluded", "quarters used (Q1 has 40 deg)", i11b.quarters_used, "2,3,4")

    # 12: a lateral root breaking out at 70-80 deg: its cells replace the epidermis there and stick out past the
    # root outline (r 455-560). The root outline stays on the main root. None of its cells may count.
    def e12(npz, c, rg, m):
        lr = (RR >= EDGES[-1]) & (RR < 560) & (TH >= 70) & (TH < 80)
        rr = RR[lr]
        c[lr] = c.max() + 1 + ((rr - EDGES[-1]) // 25).astype(int)        # 5 lateral-root cells, stacked outward
        P.save_edit(npz, cells=c)
    edit("12_lateral_root_outside", e12, {"n_cortex_cells": 288, "n_cortical_layers": 8,
                                          "mean_cortex_cell_area_px": EXPECT["mean_cortex_cell_area_px"]},
         "cells outside the root outline must not count")

    # 13: the same lateral root, but the reviewer leaves the root outline round it (outline includes the LR)
    def e13(npz, c, rg, m):
        lr = (RR >= EDGES[-1]) & (RR < 560) & (TH >= 70) & (TH < 80)
        c[lr] = c.max() + 1 + ((RR[lr] - EDGES[-1]) // 25).astype(int)
        rg[lr] = 1
        P.save_edit(npz, cells=c, region=rg)
    edit("13_lateral_root_inside_outline", e13, {"n_cortex_cells": 288},
         "shows what happens if the outline takes in the lateral root")


if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(tempfile.mkdtemp(prefix="rooh_check_"))
    shutil.rmtree(out, ignore_errors=True)
    out.mkdir(parents=True)
    try:
        run(out)
    except Exception as e:                                    # report what ran
        CHECKS.append({"case": "ERROR", "check": repr(e)[:200], "ok": "FAIL"})
        import traceback
        traceback.print_exc()
    t = pd.DataFrame(CHECKS)
    pd.set_option("display.width", 220, "display.max_colwidth", 60, "display.max_rows", 200)
    print(t.to_string(index=False))
    t.to_excel(out / "checks.xlsx", index=False)
    print(f"\n{(t.ok == 'PASS').sum()} pass, {(t.ok == 'FAIL').sum()} fail  ->  {out / 'checks.xlsx'}")
    sys.exit(int((t.ok == "FAIL").any()))
