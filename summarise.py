"""Masks -> tables. Run any time; it reads whatever is saved (auto or reviewed).

  python summarise.py <anatomy_dir>

Uses <anatomy_dir>/samples.xlsx (one row per image; made by run_batch.py, or here if missing). Fill it in:
  plant      the replicate; images of the same plant are averaged
  group_*    columns from the folder names; rename the headers (species, genotype ...)
             or add columns (treatment ...). Every column except image_file, plant,
             um_per_px and include is treated as a grouping column.
  um_per_px  microns per pixel, from the microscope/camera: always given by the user (scale bars are not read).
             Fill it any time after measuring: the scale is used only here, never for outlining.
  include    "no" drops the image
Writes <anatomy_dir>/../results/anatomy_summary_YYYY-MM-DD.xlsx
"""
from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from skimage import measure

import pattern_plots as PP
import pipeline as P

FIXED = ["image_file", "plant", "um_per_px", "include"]
N_BINS = 10   # radial bins across the cortex: 1 = next to the stele ... 10 = next to the epidermis
# what every column means: written into each workbook's README sheet. Columns ending _px are the same
# measure in pixels (kept at the end of each sheet, for images without a scale). Areas run to mid-wall.
COLUMNS = {
    "image_file": "image path inside the image folder (the tag that links every row to its picture)",
    "plant": "plant (biological replicate), from the folder name; images of one plant share it",
    "group_1": "first folder level (e.g. species); rename in anatomy/samples.xlsx",
    "group_2": "second folder level (e.g. genotype); rename in anatomy/samples.xlsx",
    "tech_rep": "image number within its plant (technical replicate)",
    "status": "auto = computer outlines only; reviewed = checked by hand (Save + next)",
    "um_per_px": "scale given in samples.xlsx (micrometres per pixel); empty = no µm columns",
    "n_images": "images averaged into this plant row",
    "n_reviewed": "of those, images checked by hand",
    "n_files": "cortical cell files, stele to epidermis: median of the four spokes",
    "n_files_q1": "files on the Q1 spoke (top-right quarter); q2 top-left, q3 bottom-left, q4 bottom-right: a check",
    "n_files_quadrant": "files on this cell's spoke",
    "quadrant": "1 top-right, 2 top-left, 3 bottom-left, 4 bottom-right around the centre (where the spoke lies)",
    "file": "cell file number on the spoke: 1 = next to the stele",
    "rel_position": "(file - 0.5) / n_files_quadrant: 0 = stele side, 1 = epidermis side",
    "label": "cell number inside the image's mask",
    "tissue": "cortex | stele | epidermis | air_space | metaxylem | unassigned | cortex_dropped / air_space_dropped (in a dropped quarter)",
    "dropped_quarters": "quarters left out by hand (damaged tissue): cortex traits (files, cell sizes, air spaces, "
                        "cortex width) use the other quarters; root, stele and metaxylem areas stay whole-section",
    "ring": "adjacency steps from the stele (a gap check; runs high)",
    "x": "cell centre, pixels from the left", "y": "cell centre, pixels from the top",
    "angle_rad": "direction of the cell centre from the centre, counter-clockwise from the right (radians)",
    "angle_deg": "the same in degrees",
    "r_um": "distance of the cell centre from the centre (µm)",
    "r_norm": "radial position of the cell centre: 0 = stele edge, 1 = root edge, along its direction",
    "radial_bin": f"r_norm in {N_BINS} equal bins: 1 = next to the stele",
    "dist_from_stele_um": "distance from the stele edge to the cell centre (µm)",
    "area_um2": "cell area (µm², to mid-wall)", "major_um": "cell length, long axis (µm)",
    "minor_um": "cell width, short axis (µm)",
    "stele_radius_um": "stele radius along this direction (µm)", "root_radius_um": "root radius along this direction (µm)",
    "root_area_um2": "whole cross-section area (µm²)", "stele_area_um2": "stele area (µm²)",
    "cortex_area_um2": "root minus stele (µm²)", "cortex_width_um": "mean root radius minus stele radius (µm)",
    "n_cortex_cells": "cortical cells outlined", "mean_cortex_cell_area_um2": "mean cortical cell area, all cortical cells (µm²)",
    "n_air_spaces": "small spaces between cortical cells", "air_space_total_area_um2": "their total area (µm²)",
    "mx_n": "metaxylem vessels", "mx_mean_area_um2": "mean metaxylem vessel area (µm²)",
    "mx_total_area_um2": "total metaxylem area (µm²)",
    "n_cells": "cells averaged", "mean_area_um2": "mean cell area (µm²)", "mean_minor_um": "mean cell width (µm)",
    "mean_r_norm": "mean radial position",
    "area_unit": "unit of the cell sizes in this row: um2, or px when the image has no scale",
    "mean_cell_area_spoke": "mean area of the spoke cells", "mean_cell_area_cortex": "mean area of all cortical cells",
    "n_spoke_cells": "cells on the four spokes",
    "peak_height_file": "largest cell size on the fitted curve; _file = x is file number, _rel = relative position, "
                        "_radial = all cortical cells by radial position (same for the four columns below)",
    "peak_pos_file": "where the peak is (file number / relative position / radial position)",
    "rise_slope_file": "(peak - first value) / distance: how fast cells grow from the stele to the peak",
    "decay_slope_file": "slope just past the peak (over 20% of the range): how fast cells shrink outward",
    "outer_rise_file": "last value - lowest value after the peak: an upturn at the epidermis",
    "axis": "file, rel or radial: which x the curve is on", "fitted": "the curve's cell size at x",
}


def make_samples(files: list[Path]) -> pd.DataFrame:
    rows = []
    for f in files:
        rel = Path(str(np.load(f)["image_file"]))
        parts = list(rel.parent.parts)
        plant = (parts[-1] if len(parts) >= 3 else rel.stem).lower()   # "P1" and "p1" are one plant
        groups = parts[:-1] if len(parts) >= 3 else parts
        rows.append({"image_file": str(rel), "plant": plant, "um_per_px": None, "include": "yes",
                     **{f"group_{k + 1}": g for k, g in enumerate(groups)}})
    return pd.DataFrame(rows)


def write_samples(anatomy_dir: Path) -> Path:
    """samples.xlsx for the user to fill (never overwritten once it exists)."""
    sfile = anatomy_dir / "samples.xlsx"
    if not sfile.exists():
        make_samples(sorted((anatomy_dir / "masks").glob("*.npz"))).to_excel(sfile, index=False)
        print(f"\nWrote {sfile}\nFill in um_per_px (microns per pixel) for every image, rename group_1, group_2 ... "
              "(e.g. species, genotype), check plant, then run Summarise.")
    return sfile


def px_last(df):
    """Pixel columns at the end: the µm columns are the ones to read."""
    px = [c for c in df.columns if c.endswith("_px") and c != "um_per_px"]
    return df[[c for c in df.columns if c not in px] + px]


def describe(col):
    if col in COLUMNS:
        return COLUMNS[col]
    for sfx in ("_rel", "_radial"):                  # profile features: same meaning on another x
        if col.endswith(sfx) and col[:-len(sfx)] + "_file" in COLUMNS:
            return COLUMNS[col[:-len(sfx)] + "_file"].split(";")[0] + f" ({sfx[1:]} axis)"
    if col.endswith("_px"):
        base = col[:-3]
        for unit in ("_um2", "_um"):
            if base + unit in COLUMNS:
                return COLUMNS[base + unit].split("(µm")[0].strip() + " (pixels)"
        return "pixels"
    if col.startswith("n_files_q"):
        return COLUMNS["n_files_q1"]
    return ""


def main(anatomy_dir: Path):
    settings = json.loads((anatomy_dir / "settings.json").read_text())
    files = sorted((anatomy_dir / "masks").glob("*.npz"))
    sfile = write_samples(anatomy_dir)
    samples = pd.read_excel(sfile, dtype={"image_file": str, "plant": str})
    samples["include"] = samples["include"].fillna("yes").astype(str).str.lower()
    groups = [c for c in samples.columns if c not in FIXED]
    samples = samples.set_index("image_file")
    tech_rep = samples.groupby(["plant", *groups], dropna=False).cumcount() + 1   # image number within a plant

    angle = float(settings.get("quadrant_angle_deg", 0))
    fig_dir = anatomy_dir.parent / "results" / "figures" / "cortex_pattern"
    fig_dir.mkdir(parents=True, exist_ok=True)
    img_rows, cell_rows, file_rows, all_rows = [], [], [], []
    for f in files:
        z = np.load(f)
        rel = str(z["image_file"])
        if rel not in samples.index or samples.at[rel, "include"] == "no":
            continue
        cells, region, mx, status = P.load_state(f)
        centre = P.centre_of(region, P.load_centre(f))
        props, mxp, tr = P.derive(cells, region, mx, angle, centre)
        paths = P.load_paths(f) or P.auto_paths(cells, region, props, angle, centre)
        dropped = P.load_dropped(f)
        props = P.drop_quarters(props, dropped)
        ft = P.file_table(cells, region, props, paths[1], paths[2], dropped)
        s = samples.at[rel, "um_per_px"]
        s = float(s) if pd.notna(s) else np.nan
        meta = {"image_file": rel, "plant": samples.at[rel, "plant"], **{g: samples.at[rel, g] for g in groups}}

        nq = ft.groupby("quadrant").file.max().reindex([1, 2, 3, 4]) if len(ft) else pd.Series([np.nan] * 4, index=[1, 2, 3, 4])
        tr = {k: v for k, v in tr.items() if not k.startswith(("n_cortical_layers", "n_layers"))}
        # cortex traits from the kept quarters only; root, stele and metaxylem stay whole-section
        tr["n_cortex_cells"] = int((props.tissue == "cortex").sum())
        tr["n_air_spaces"] = int((props.tissue == "air_space").sum())
        tr["air_space_total_area_px"] = int(props.area_px[props.tissue == "air_space"].sum())
        if dropped:
            r_st, r_rt = P.radial_profile(region, centre)
            q_deg = (np.floor(np.mod(np.arange(360) - angle, 360) / 90) + 1).astype(int)
            keep = ~np.isin(q_deg, dropped)
            tr["cortex_width_px"] = float(np.mean((r_rt - r_st)[keep])) if keep.any() else np.nan
        tr["dropped_quarters"] = " ".join(f"Q{q}" for q in sorted(dropped))
        tr["mean_cortex_cell_area_px"] = float(props.area_px[props.tissue == "cortex"].mean()) if len(props) else np.nan
        img_rows.append({**meta, "tech_rep": int(tech_rep[rel]), "status": status, "um_per_px": s,
                         "n_files": float(nq.median()), **{f"n_files_q{q}": nq[q] for q in (1, 2, 3, 4)}, **tr,
                         **{k.replace("_px", "_um2"): v * s * s for k, v in tr.items() if k.endswith("area_px")},
                         **{k.replace("_px", "_um"): v * s for k, v in tr.items() if k.endswith("width_px")}})

        tech = int(tech_rep[rel])
        if len(ft):
            pr = props.set_index("label")
            ft = ft.join(pr[["angle_rad", "stele_radius_px", "root_radius_px"]], on="label")
            ft = ft.assign(n_files_quadrant=ft.quadrant.map(nq),
                           rel_position=(ft.file - 0.5) / ft.quadrant.map(nq),
                           area_um2=ft.area_px * s * s, major_um=ft.major_px * s, minor_um=ft.minor_px * s,
                           dist_from_stele_um=ft.dist_from_stele_px * s,
                           stele_radius_um=ft.stele_radius_px * s, root_radius_um=ft.root_radius_px * s)
            file_rows.append(ft.assign(**meta, tech_rep=tech, status=status))
            um = not np.isnan(s)
            PP.image_figure(P.load_gray(P.image_dir(anatomy_dir, settings) / rel), cells, region, ft,
                            props.assign(area_um2=props.area_px * s * s), centre, angle, "area_um2" if um else "area_px", "Cell area (µm²)" if um else "Cell area (px)",
                            f"{rel}  ·  {status}", fig_dir / (Path(rel).with_suffix("").as_posix().replace("/", "__") + ".png"))

        mxt = pd.DataFrame(measure.regionprops_table(mx, properties=("label", "area", "centroid", "axis_major_length",
                                                                     "axis_minor_length"))) if mx.any() else pd.DataFrame()
        if len(mxt):
            mxt = mxt.rename(columns={"centroid-0": "y", "centroid-1": "x", "area": "area_px",
                                      "axis_major_length": "major_px", "axis_minor_length": "minor_px"})
            mxt = mxt.join(P.radial_position(mxt, region, centre))
            mxt["tissue"], mxt["quadrant"] = "metaxylem", P.quadrant_of(mxt.y, mxt.x, centre, angle)
        allc = pd.concat([props, mxt], ignore_index=True)
        allc["file"] = allc.label.where(allc.tissue != "metaxylem").map(dict(zip(ft.label, ft.file))) if len(ft) else np.nan
        allc["ring"] = allc.layer if "layer" in allc else np.nan
        allc = allc.assign(area_um2=allc.area_px * s * s, major_um=allc.major_px * s, minor_um=allc.minor_px * s,
                           r_um=allc.r_px * s, dist_from_stele_um=allc.dist_from_stele_px * s,
                           angle_deg=np.degrees(allc.angle_rad))
        all_rows.append(allc.assign(**meta, tech_rep=tech, status=status)[
            [*meta, "tech_rep", "status", "label", "tissue", "quadrant", "file", "ring", "x", "y",
             "angle_rad", "angle_deg", "r_px", "r_um", "r_norm", "dist_from_stele_um",
             "area_px", "area_um2", "major_um", "minor_um"]])

        c = props[props.tissue == "cortex"].copy()
        c["area_um2"] = c.area_px * s * s
        c["major_um"], c["minor_um"] = c.major_px * s, c.minor_px * s
        c["dist_from_stele_um"] = c.dist_from_stele_px * s
        c["radial_bin"] = (np.clip(c.r_norm, 0, 0.999) * N_BINS).astype(int) + 1
        c["file"] = c.label.map(dict(zip(ft.label, ft.file))) if len(ft) else np.nan
        cell_rows.append(c.assign(**meta, status=status)[[*meta, "status", "quadrant", "file", "r_norm", "radial_bin",
                                                          "dist_from_stele_um", "label", "area_px", "area_um2",
                                                          "major_um", "minor_um"]])

    keys = ["plant", *groups]
    first = lambda df, cols: df[[*cols, *[c for c in df.columns if c not in cols]]]
    per_image = pd.DataFrame(img_rows)
    per_cell = pd.concat(cell_rows, ignore_index=True) if cell_rows else pd.DataFrame()
    per_file_path = first(pd.concat(file_rows, ignore_index=True), ["image_file", *keys, "status", "quadrant", "file"]) \
        if file_rows else pd.DataFrame()
    per_plant_file = (per_file_path.groupby([*keys, "file"], dropna=False)
                      .agg(n_cells=("area_um2", "size"), mean_area_um2=("area_um2", "mean"),
                           mean_area_px=("area_px", "mean"), mean_minor_um=("minor_um", "mean"),
                           mean_r_norm=("r_norm", "mean")).reset_index())
    per_radial_img = (per_cell.groupby(["image_file", *keys, "radial_bin"], dropna=False)
                      .agg(n_cells=("area_px", "size"), mean_area_um2=("area_um2", "mean"),
                           mean_area_px=("area_px", "mean")).reset_index())
    per_plant_radial = (per_radial_img.groupby([*keys, "radial_bin"], dropna=False)
                        [["mean_area_um2", "mean_area_px", "n_cells"]].mean().reset_index())
    num = [c for c in per_image.columns if c.endswith(("_px", "_um2", "_um")) or c.startswith(("n_", "mx_"))]
    per_plant = (per_image.groupby(keys, dropna=False)
                 .agg(n_images=("image_file", "size"),
                      n_reviewed=("status", lambda s: int((s == "reviewed").sum())),
                      **{c: (c, "mean") for c in num}).reset_index())

    per_file_path = first(per_file_path, ["image_file", *keys, "tech_rep", "status", "quadrant", "file",
                                          "n_files_quadrant", "rel_position"]) if len(per_file_path) else per_file_path
    all_cells = pd.concat(all_rows, ignore_index=True)
    # plan step 2-3: one spline per image, the four quarters pooled (6-9 points per quarter is too few for a
    # stable GCV fit); features read off it on three axes, and the curve itself on a fixed grid (for fPCA)
    feat_rows, curve_rows = [], []
    grid_rel = np.round(np.linspace(0, 1, 21), 2)
    for img, d in per_file_path.groupby("image_file") if len(per_file_path) else []:
        um = d.area_um2.notna().all()
        ycol = "area_um2" if um else "area_px"
        c = PP.with_bin(all_cells[all_cells.image_file == img])
        b = c.groupby("r_bin")[ycol].agg(["mean", "size"])
        meta = {k: d[k].iloc[0] for k in ["image_file", *keys, "tech_rep", "status"]}
        row = {**meta, "n_files": d.groupby("quadrant").file.max().median(),
               "stele_radius_um": d.stele_radius_um.iloc[0], "root_radius_um": d.root_radius_um.iloc[0],
               "area_unit": "um2" if um else "px", "mean_cell_area_spoke": d[ycol].mean(),
               "mean_cell_area_cortex": c[ycol].mean(), "n_spoke_cells": len(d), "n_cortex_cells": len(c)}
        axes = {"_file": (d.file, d[ycol], None), "_rel": (d.rel_position, d[ycol], None),
                "_radial": (b.index, b["mean"], b["size"])}
        for sfx, (x, y, w) in axes.items():
            row.update({k + sfx: v for k, v in PP.profile_features(x, y, w).items()})
            at = np.arange(1, d.file.max() + 1) if sfx == "_file" else grid_rel
            sp = PP.spline(x, y, w, at=at)
            curve_rows += [{**meta, "axis": sfx[1:], "x": xv, "fitted": yv, "area_unit": row["area_unit"]}
                           for xv, yv in zip(at, sp[1] if sp is not None else [np.nan] * len(at))]
        feat_rows.append(row)
    spline_curves = pd.DataFrame(curve_rows)
    profile_features = pd.DataFrame(feat_rows)
    out_dir = anatomy_dir.parent / "results"
    out_dir.mkdir(exist_ok=True)
    today = dt.date.today().isoformat()
    files_out = {
        f"cell_size_by_file_{today}.xlsx": {
            "per_file_path": (per_file_path, "one path cell per row: image x quadrant x file. file 1 = next to the "
                              "endodermis; one spoke per quadrant, picked where the cortex is most complete, edited "
                              "by hand. rel_position = (file - 0.5) / n_files_quadrant. plant = biological replicate, "
                              "tech_rep = image number within the plant. status auto = not reviewed yet"),
            "profile_features": (profile_features, "one row per image: smoothing spline (GCV) through cell size, the "
                                 "four quarters pooled, features read off the curve on three axes. _file = absolute file "
                                 "number; _rel = relative position (file - 0.5) / n_files_quadrant; _radial = all cortical cells of the image, "
                                 f"mean per {N_BINS} r_norm bins (0 = stele edge, 1 = epidermis). peak_height, "
                                 "peak_pos; rise_slope = (peak - first) / distance; decay_slope = local, over 20% of "
                                 "the range past the peak; outer_rise = last value - minimum after the peak. Slopes "
                                 "and heights in area_unit per file (or per unit relative position). NaN = fewer "
                                 f"than {PP.MIN_PTS} files or feature undefined (peak at an end)"),
            "spline_curves": (spline_curves, "the fitted image curves themselves (long format): axis file = at files "
                              "1..n, rel / radial = at 0, 0.05 ... 1; NaN outside the data (no extrapolation). "
                              "A cubic smoothing spline has no short parameter list; the curve on a fixed grid is "
                              "the input for fPCA (plan step 5)"),
            "per_plant_file": (per_plant_file, "plant x file: path cells pooled over quadrants and images "
                               "(display only; analyse per_file_path with plant as the unit)")},
        f"all_cells_{today}.xlsx": {
            "all_cells": (all_cells, "every segmented object: tissue = cortex | stele | epidermis | air_space | "
                          "unassigned | metaxylem. angle_rad counter-clockwise from the image's right, around the "
                          "stele centre. r_norm 0 = stele edge, 1 = root edge along the object's direction. "
                          "file = path file (path cells only); ring = adjacency hops from the stele (gap check, runs high)")},
        f"anatomy_summary_{today}.xlsx": {
            "per_plant": (per_plant, "plant means of per_image; n_reviewed = images checked by hand"),
            "per_image": (per_image, "whole-root traits per image; n_files = median of the 4 quadrant paths"),
            "per_plant_radial": (per_plant_radial, f"all cortical cells in {N_BINS} equal r_norm bins (1 = inner)")},
    }
    for fname, tabs in files_out.items():
        tabs = {name: (px_last(df), note) for name, (df, note) in tabs.items()}
        readme = pd.DataFrame([(k, "", n) for k, (_, n) in tabs.items()], columns=["sheet", "column", "what it is"])
        cols = pd.DataFrame([(name, c, describe(c)) for name, (df, _) in tabs.items() for c in df.columns],
                            columns=["sheet", "column", "what it is"])
        with pd.ExcelWriter(out_dir / fname) as w:
            pd.concat([readme, cols], ignore_index=True).to_excel(w, sheet_name="README", index=False)
            for name, (df, _) in tabs.items():
                df.to_excel(w, sheet_name=name, index=False)
    if len(per_file_path) and groups:
        sp, gt = groups[0], groups[-1]            # e.g. species, genotype
        (out_dir / "figures").mkdir(exist_ok=True)
        for spv, d in per_file_path.groupby(sp):
            um = d.area_um2.notna().any()          # µm² if any image has a scale; images without one are left out
            drop = d.loc[d.area_um2.isna(), "image_file"].nunique() if um else 0
            ac = all_cells[all_cells[sp] == spv]
            PP.group_figure(d[d.area_um2.notna()] if um else d, ac[ac.area_um2.notna()] if um else ac, gt, spv,
                            "area_um2" if um else "area_px", "Cortical cell area (µm²)" if um else "Cortical cell area (px)",
                            out_dir / "figures" / f"cortex_pattern_{spv}_{today}.png",
                            f"  {drop} image(s) without a scale left out." if drop else "")
    no_scale = per_image.loc[per_image.um_per_px.isna(), "image_file"].tolist()
    if no_scale:
        print(f"\nWARNING: {len(no_scale)} image(s) have no scale (um_per_px empty); their sizes "
              f"stay in pixels and they are left out of the um figures. Type um_per_px in {sfile.name} and run again:")
        print("  " + "\n  ".join(no_scale[:10]) + ("\n  ..." if len(no_scale) > 10 else ""))
    out = out_dir
    print(f"{len(per_image)} images ({(per_image.status == 'reviewed').sum()} reviewed), "
          f"{per_plant.shape[0]} plants -> {out}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(Path(sys.argv[1]).resolve())
