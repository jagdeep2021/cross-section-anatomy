"""Masks -> tables. Run any time; it reads whatever is saved (auto or reviewed).

  python summarise.py <anatomy_dir>

First run writes <anatomy_dir>/samples.xlsx (one row per image) for you to check:
  plant      the replicate; images of the same plant are averaged
  group_*    columns from the folder names; rename the headers (species, genotype ...)
             or add columns (treatment ...). Every column except image_file, plant,
             um_per_px and include is treated as a grouping column.
  um_per_px  leave empty to use the scale bar; fill to override (e.g. images with no bar)
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

import pipeline as P

FIXED = ["image_file", "plant", "um_per_px", "include"]
N_BINS = 10   # radial bins across the cortex: 1 = next to the stele ... 10 = next to the epidermis


def make_samples(files: list[Path]) -> pd.DataFrame:
    rows = []
    for f in files:
        rel = Path(str(np.load(f)["image_file"]))
        parts = list(rel.parent.parts)
        plant = parts[-1] if len(parts) >= 3 else rel.stem
        groups = parts[:-1] if len(parts) >= 3 else parts
        rows.append({"image_file": str(rel), "plant": plant, "um_per_px": None, "include": "yes",
                     **{f"group_{k + 1}": g for k, g in enumerate(groups)}})
    return pd.DataFrame(rows)


def main(anatomy_dir: Path):
    settings = json.loads((anatomy_dir / "settings.json").read_text())
    files = sorted((anatomy_dir / "masks").glob("*.npz"))
    sfile = anatomy_dir / "samples.xlsx"
    if not sfile.exists():
        make_samples(files).to_excel(sfile, index=False)
        print(f"Wrote {sfile}: check plant / group columns, then run again.")
    samples = pd.read_excel(sfile, dtype={"image_file": str, "plant": str})
    samples["include"] = samples["include"].fillna("yes").astype(str).str.lower()
    groups = [c for c in samples.columns if c not in FIXED]
    samples = samples.set_index("image_file")
    tech_rep = samples.groupby(["plant", *groups], dropna=False).cumcount() + 1   # image number within a plant

    angle = float(settings.get("quadrant_angle_deg", 0))
    img_rows, cell_rows, quad_rows, file_rows, all_rows = [], [], [], [], []
    for f in files:
        z = np.load(f)
        rel = str(z["image_file"])
        if rel not in samples.index or samples.at[rel, "include"] == "no":
            continue
        cells, region, mx, status = P.load_state(f)
        props, mxp, tr = P.derive(cells, region, mx, angle)
        paths = P.load_paths(f) or P.auto_paths(cells, region, props, angle)
        ft = P.file_table(cells, region, props, paths[1], paths[2])
        bar = int(z["scale_bar_px"]) if "scale_bar_px" in z else 0
        s = samples.at[rel, "um_per_px"]
        if pd.isna(s):
            s = settings.get("scale_bar_um", 0) / bar if bar else np.nan
        meta = {"image_file": rel, "plant": samples.at[rel, "plant"], **{g: samples.at[rel, g] for g in groups}}

        nq = ft.groupby("quadrant").file.max().reindex([1, 2, 3, 4]) if len(ft) else pd.Series([np.nan] * 4, index=[1, 2, 3, 4])
        tr = {k: v for k, v in tr.items() if not k.startswith(("n_cortical_layers", "n_layers"))}
        img_rows.append({**meta, "status": status, "um_per_px": s,
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

        mxt = pd.DataFrame(measure.regionprops_table(mx, properties=("label", "area", "centroid", "axis_major_length",
                                                                     "axis_minor_length"))) if mx.any() else pd.DataFrame()
        if len(mxt):
            mxt = mxt.rename(columns={"centroid-0": "y", "centroid-1": "x", "area": "area_px",
                                      "axis_major_length": "major_px", "axis_minor_length": "minor_px"})
            mxt = mxt.join(P.radial_position(mxt, region, P.centre_of(region)))
            mxt["tissue"], mxt["quadrant"] = "metaxylem", P.quadrant_of(mxt.y, mxt.x, P.centre_of(region), angle)
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

        qt = P.quadrant_traits(cells, region, mx, props, angle)
        qt = qt.drop(columns=[c for c in qt.columns if c.startswith("n_layers")])
        qt.insert(1, "n_files", qt.quadrant.map(nq).values)
        for col in [c for c in qt.columns if c.endswith("area_px")]:
            qt[col.replace("_px", "_um2")] = qt[col] * s * s
        quad_rows.append(qt.assign(**meta, status=status))

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
    per_quadrant = first(pd.concat(quad_rows, ignore_index=True), ["image_file", *keys, "status", "quadrant"])
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
    out_dir = anatomy_dir.parent / "results"
    out_dir.mkdir(exist_ok=True)
    today = dt.date.today().isoformat()
    files_out = {
        f"cell_size_by_file_{today}.xlsx": {
            "per_file_path": (per_file_path, "one path cell per row: image x quadrant x file. file 1 = next to the "
                              "endodermis; one spoke per quadrant, picked where the cortex is most complete, edited "
                              "by hand. rel_position = (file - 0.5) / n_files_quadrant. plant = biological replicate, "
                              "tech_rep = image number within the plant. status auto = not reviewed yet"),
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
            "per_quadrant": (per_quadrant, f"Q1 top-right, Q2 top-left, Q3 bottom-left, Q4 bottom-right around the "
                             f"stele centre (rotated {angle:g} deg); sub-samples within an image"),
            "per_plant_radial": (per_plant_radial, f"all cortical cells in {N_BINS} equal r_norm bins (1 = inner)")},
    }
    for fname, tabs in files_out.items():
        readme = pd.DataFrame([(k, n) for k, (_, n) in tabs.items()], columns=["sheet", "what it is"])
        with pd.ExcelWriter(out_dir / fname) as w:
            readme.to_excel(w, sheet_name="README", index=False)
            for name, (df, _) in tabs.items():
                df.to_excel(w, sheet_name=name, index=False)
    out = out_dir
    print(f"{len(per_image)} images ({(per_image.status == 'reviewed').sum()} reviewed), "
          f"{per_plant.shape[0]} plants -> {out}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(Path(sys.argv[1]).resolve())
