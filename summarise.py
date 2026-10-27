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

import pipeline as P

FIXED = ["image_file", "plant", "um_per_px", "include"]


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

    angle = float(settings.get("quadrant_angle_deg", 0))
    img_rows, cell_rows, quad_rows = [], [], []
    for f in files:
        z = np.load(f)
        rel = str(z["image_file"])
        if rel not in samples.index or samples.at[rel, "include"] == "no":
            continue
        cells, region, mx, status = P.load_state(f)
        props, mxp, tr = P.derive(cells, region, mx, angle)
        bar = int(z["scale_bar_px"]) if "scale_bar_px" in z else 0
        s = samples.at[rel, "um_per_px"]
        if pd.isna(s):
            s = settings.get("scale_bar_um", 0) / bar if bar else np.nan
        meta = {"image_file": rel, "plant": samples.at[rel, "plant"], **{g: samples.at[rel, g] for g in groups}}
        img_rows.append({**meta, "status": status, "um_per_px": s, **tr,
                         **{k.replace("_px", "_um2"): v * s * s for k, v in tr.items() if k.endswith("area_px")}})
        qt = P.quadrant_traits(cells, region, mx, props, angle)
        for col in [c for c in qt.columns if c.endswith("area_px")]:
            qt[col.replace("_px", "_um2")] = qt[col] * s * s
        quad_rows.append(qt.assign(**meta, status=status))
        c = props[props.tissue == "cortex"].copy()
        c["area_um2"] = c.area_px * s * s
        c["major_um"], c["minor_um"] = c.major_px * s, c.minor_px * s
        cell_rows.append(c.assign(**meta, status=status)[[*meta, "status", "quadrant", "layer", "label", "area_px", "area_um2", "major_um", "minor_um"]])

    per_image = pd.DataFrame(img_rows)
    per_cell = pd.concat(cell_rows, ignore_index=True) if cell_rows else pd.DataFrame()
    keys = ["plant", *groups]
    per_layer_img = (per_cell.groupby(["image_file", *keys, "layer"], dropna=False)
                     .agg(n_cells=("area_px", "size"), mean_area_um2=("area_um2", "mean"),
                          mean_area_px=("area_px", "mean")).reset_index())
    num = [c for c in per_image.columns if c.endswith(("_px", "_um2")) or c.startswith(("n_", "mx_"))]
    per_quadrant = pd.concat(quad_rows, ignore_index=True)
    per_quadrant = per_quadrant[["image_file", *keys, "status", *[c for c in per_quadrant.columns
                                                                  if c not in ("image_file", "status", *keys)]]]
    per_layer_quadrant = (per_cell.groupby(["image_file", *keys, "quadrant", "layer"], dropna=False)
                          .agg(n_cells=("area_px", "size"), mean_area_um2=("area_um2", "mean"),
                               mean_area_px=("area_px", "mean")).reset_index())
    per_plant = (per_image.groupby(keys, dropna=False)
                 .agg(n_images=("image_file", "size"),
                      n_reviewed=("status", lambda s: int((s == "reviewed").sum())),
                      **{c: (c, "mean") for c in num}).reset_index())
    per_plant_layer = (per_layer_img.groupby([*keys, "layer"], dropna=False)
                       [["mean_area_um2", "mean_area_px", "n_cells"]].mean().reset_index())

    readme = pd.DataFrame({"sheet": ["per_plant", "per_plant_layer", "per_image", "per_layer_image",
                                     "per_quadrant", "per_layer_quadrant", "per_cell"],
                           "one row per": ["plant (mean of its images)", "plant x cortical layer", "image",
                                           "image x cortical layer", "image x quadrant", "image x quadrant x layer",
                                           "cortical cell"],
                           "note": ["n_reviewed = images checked by hand; use only fully reviewed plants for results",
                                    "layer 1 = touching the endodermis",
                                    "status auto = not yet checked in review.py; air spaces = small objects at cell junctions",
                                    "",
                                    f"Q1 top-right, Q2 top-left, Q3 bottom-left, Q4 bottom-right around the stele centre "
                                    f"(rotated {angle:g} deg; settings.json quadrant_angle_deg). Quadrants are sub-samples "
                                    f"within an image, not comparable across images",
                                    "", "um values empty where the image has no scale"]})
    out_dir = anatomy_dir.parent / "results"
    out_dir.mkdir(exist_ok=True)
    out = out_dir / f"anatomy_summary_{dt.date.today().isoformat()}.xlsx"
    with pd.ExcelWriter(out) as w:
        for name, df in [("README", readme), ("per_plant", per_plant), ("per_plant_layer", per_plant_layer),
                         ("per_image", per_image), ("per_layer_image", per_layer_img),
                         ("per_quadrant", per_quadrant), ("per_layer_quadrant", per_layer_quadrant), ("per_cell", per_cell)]:
            df.to_excel(w, sheet_name=name, index=False)
    print(f"{len(per_image)} images ({(per_image.status == 'reviewed').sum()} reviewed), "
          f"{per_plant.shape[0]} plants -> {out}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(Path(sys.argv[1]).resolve())
