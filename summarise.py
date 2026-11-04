"""Masks -> tables. Run any time; it reads whatever is saved (auto or reviewed).

  python summarise.py <anatomy_dir> [--image-figures]     (also one cortex pattern figure per image)

Uses <anatomy_dir>/samples.xlsx (one row per image; made by run_batch.py, or here if missing). Fill it in:
  plant      the replicate; images of the same plant are averaged
  group_*    columns from the folder names; rename the headers (species, genotype ...)
             or add columns (treatment ...). Every column except image_file, plant,
             um_per_px and include is treated as a grouping column.
  um_per_px  microns per pixel, from the microscope/camera: always given by the user (scale bars are not read).
             Fill it any time after measuring: the scale is used only here, never for outlining.
  include    "no" drops the image
  notes      any column whose name starts with note / check / comment / remark: free text, copied to per_image
  epidermis  yes = intact, no = lost in sectioning (the outermost layer is then outer cortex, tagged epidermis by
             position); empty = not checked. The pattern counts every cell outside the endodermis either way.
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

FIXED = ["image_file", "plant", "um_per_px", "include", "epidermis"]
NOTE_COLS = ("note", "check", "comment", "remark")   # samples.xlsx columns starting so = free text: copied to per_image,
                                                    # never a grouping column
N_BINS = 10   # radial bins across the cortex: 1 = next to the stele ... 10 = next to the epidermis; epidermis its own point
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
    "angle_kept_deg": "degrees of the section measured: 360 minus the excluded sectors (damaged tissue taken out by hand)",
    "excluded_sectors": "the excluded arcs, 'start-end reason' in degrees counter-clockwise from the image's right",
    "stele_found": "FALSE = no stele outlined (root tip, or stele cut off by the frame): stele, metaxylem, cortex "
                   "and aerenchyma columns are blank for this image (root size still measured)",
    "low_coverage": f"TRUE when fewer than {P.MIN_KEPT_DEG} degrees are kept: cortex numbers rest on few cells",
    "label": "cell number inside the image's mask",
    "tissue": "cortex | stele | epidermis | air_space | metaxylem | unassigned; '..._excluded' = centre in an excluded "
              "sector (kept in this table, out of every cortex number)",
    "sector": "kept or excluded (by the cell's centre)",
    "ring": "adjacency steps from the stele (a gap check; runs high)",
    "x": "cell centre, pixels from the left", "y": "cell centre, pixels from the top",
    "angle_rad": "direction of the cell centre from the centre, counter-clockwise from the right (radians)",
    "angle_deg": "the same in degrees",
    "r_um": "distance of the cell centre from the centre (µm)",
    "r_norm": "radial position of the cell centre: 0 = stele edge, 1 = root edge, along its own direction",
    "radial_bin": f"r_norm in {N_BINS} equal bins: 1 = next to the stele",
    "dist_from_stele_um": "distance from the stele edge to the cell centre (µm)",
    "area_um2": "cell area (µm², to mid-wall)", "major_um": "cell length, long axis (µm)",
    "minor_um": "cell width, short axis (µm)",
    "root_area_um2": "whole cross-section area as outlined (µm²); with excluded sectors use root_area_est_um2",
    "stele_area_um2": "stele area, whole section (µm²; never sector dependent)",
    "cortex_area_kept_um2": "root outside the stele (cortex + epidermis) inside the kept angles (µm²)",
    "cortex_area_est_um2": "the same scaled to the whole circle (x 360 / angle_kept_deg): an ESTIMATE that assumes "
                           "the excluded arcs looked like the kept ones",
    "root_area_est_um2": "stele area + cortex_area_est_um2 (estimate)",
    "stele_share_est": "stele area / root_area_est (estimate; whole-section stele over the scaled root)",
    "cortex_width_um": "mean root radius minus stele radius over the kept degrees (µm)",
    "n_cortical_layers": "cortical cells crossed by a ray from the stele out, median over rays every 5 degrees in the "
                         "kept angles (_min / _max: the range: one side can carry an extra layer)",
    "n_cortical_layers_min": "fewest cortical layers on any 5-degree ray (kept angles)",
    "n_cortical_layers_max": "most cortical layers on any 5-degree ray; a wide min-max range = a wrong stele outline "
                             "or an extra layer on one side: check the image",
    "n_low_coverage": f"images of this plant with fewer than {P.MIN_KEPT_DEG} degrees kept",
    "n_cortex_cells": "cortical cells counted in the kept angles", "n_cortex_cells_est": "scaled to 360 degrees (estimate)",
    "mean_cortex_cell_area_um2": "mean cortical cell area, kept cells (µm²)",
    "median_cortex_cell_area_um2": "median cortical cell area, kept cells (µm²): the reference for rel",
    "n_air_spaces": "small spaces between kept cortical cells", "n_air_spaces_est": "scaled to 360 degrees (estimate)",
    "air_space_total_area_um2": "their total area (µm²)",
    "mx_n": "metaxylem vessels (whole section)", "mx_mean_area_um2": "mean metaxylem vessel area (µm²)",
    "mx_total_area_um2": "total metaxylem area (µm²)",
    "bin": f"position: 1..{N_BINS} = r_norm bins (1 next to the stele); file: the cell file number (the epidermis is "
           "the outermost file)",
    "x_position": "where the bin sits: position = bin centre (0-1); file = median file_rel of its cells (0-1); distance = bin centre, µm from the stele edge",
    "n_cells": "cells in this bin (kept angles, cortex + epidermis); few_cells = TRUE below " + str(P.MIN_CELLS_PER_BIN),
    "few_cells": f"TRUE when the bin has fewer than {P.MIN_CELLS_PER_BIN} cells: read it with care",
    "median_area_um2": "median cell area in this bin (µm²)",
    "rel": "median area in the bin / the image's median cell area (cortex + epidermis): shape without overall size",
    "epidermis_rel": "median area of the cells tagged epidermis (within 6 px of the outline) / the image's median",
    "min_cells_in_bin": "fewest cortical cells in any cortex bin of this image",
    "area_unit": "unit of the absolute sizes in this row: um2, or px when the image has no scale",
    "peak_height_rel": "largest value on the fitted curve; _rel = relative size (shape), all cells outside the "
                       "endodermis; _abs = the same in area_unit; _rel_cortex = relative, without the cells tagged "
                       "epidermis (a check) (same for the columns below)",
    "peak_pos_rel": "where the peak is on the 0-1 axis (stele side 0, epidermis side 1)",
    "rise_slope_rel": "(peak - first value) / distance: how fast cells grow from the stele to the peak",
    "decay_slope_rel": "slope just past the peak (over 20% of the range): how fast cells shrink outward",
    "outer_rise_rel": "last value - lowest value after the peak: an upturn at the epidermis",
    "basis": "position = cells binned by r_norm (0 stele edge - 1 root edge; the headline); file = cells grouped by "
             "cell file number (1 = next to the stele; gaps such as lacunae counted as the cells they would hold); "
             f"distance = {P.DIST_BIN_UM:g} µm bins of distance from the stele edge (images with a scale only): survives "
             "lost outer tissue but carries root size, a cross-check. For distance, x_position and peak_pos are µm "
             "from the stele edge and slopes are per µm",
    "scope": "section = all kept cells of the image; Q1..Q4 = one quarter (Q1 0-90 deg top right, counter-clockwise); "
             "quarter_mean / quarter_sd = mean and SD of the quarters' features",
    "n_quarters": "quarters averaged (a quarter needs >= 45 degrees kept)",
    "kept_deg": "degrees kept in this scope",
    "quarter": "quarter of the cell / lacuna centre: 1 = 0-90 deg (top right), counter-clockwise",
    "quarters_used": "quarters with >= 45 degrees kept (used for quarter features)",
    "file": "cell file number: median of the numbers given by rays every 0.5 deg (1 = next to the stele); a gap on a "
            "ray (lacuna, missed cell, gap at the stele) counts as the cells it would have held",
    "file_rel": "relative file position: median over rays of (k - 0.5) / cortical cells on that ray",
    "n_rays": "rays that crossed the cell (file number rests on them)",
    "epidermis": "from samples.xlsx: yes = intact, no = lost (outermost layer is outer cortex), empty = not checked",
    "root_diameter_um": "2 x mean root radius over the kept degrees, from the stele centre (µm)",
    "root_diameter_eq_um": "diameter of a circle with the root's area (µm; whole outline)",
    "n_giant_cortex_cells": "kept cortical cells > 5 x the median: possibly a lacuna Cellpose took for a cell, check",
    "n_lacunae": "aerenchyma lacunae: gaps in the cortex that are not cells, >= 2 x the median cortical cell, "
                 "touching neither the root edge nor the stele, centre in the kept angles",
    "epidermis_deg": "kept degrees (5-deg steps) holding a cell tagged epidermis: Cellpose often misses the epidermis, "
                     "so the outer bins then rest on hypodermis; below 90 compare _rel with _rel_cortex",
    "n_cells": "cells in the pattern (cortex + epidermis, kept angles)",
    "n_gaps_missed_cell": "gaps of about one cell (0.5-2 x the median) inside the cortex: usually a cell Cellpose "
                          "missed. Review hint: draw the cell",
    "lacuna_area_um2": "total lacuna area (µm²)", "aerenchyma_share": "lacuna area / cortex area kept",
    "lacuna_r_norm_mean": "area-weighted mean radial position of lacuna pixels (0 stele edge, 1 root edge)",
    "lacuna": "lacuna number in the image",
    "r_norm_inner": "innermost radial position of the lacuna", "r_norm_outer": "outermost radial position",
    "cortex_px": "cortex pixels in this bin (kept angles)", "lacuna_px": "lacuna pixels in this bin",
    "lacuna_share": "lacuna pixels / cortex pixels in this bin",
    "axis": "rel (the curve of relative size)", "fitted": "the curve's value at x",
}
FEATURE_SETS = ("_rel", "_abs", "_rel_cortex")


# Topics: the README's topic column and the column order of per_image / per_plant (first rule that matches wins).
TOPICS = (
    ("1 identity", lambda c, g: c in ("image_file", "plant", "tech_rep", "status", "um_per_px", "n_images",
                                      "n_reviewed", "epidermis", "label", "lacuna", "basis", "scope", "area_unit",
                                      "bin", "axis") or c in g or c.startswith("group_")),
    ("2 check first", lambda c, g: c in ("stele_found", "angle_kept_deg", "excluded_sectors", "low_coverage", "n_low_coverage",
                                         "quarters_used", "n_quarters", "kept_deg", "epidermis_deg",
                                         "min_cells_in_bin", "few_cells", "n_gaps_missed_cell",
                                         "n_giant_cortex_cells", "sector", "n_rays", "ring")
                                or c.startswith("n_cortical_layers_m")),
    ("3 whole root", lambda c, g: c.startswith("root_")),
    ("4 stele + metaxylem", lambda c, g: c.startswith(("stele_", "mx_"))),
    ("5 cortex", lambda c, g: "cortex" in c and "lacuna" not in c and not c.endswith("_rel_cortex")
                              or c.startswith(("n_cortical_layers", "n_air_spaces", "air_space"))),
    ("7 aerenchyma", lambda c, g: "lacuna" in c or "aerenchyma" in c),
    ("6 cortex pattern", lambda c, g: c.startswith(("peak_", "rise_", "decay_", "outer_rise", "epidermis_rel"))
                                      or c in ("rel", "x_position", "n_cells", "fitted", "radial_bin")
                                      or c.startswith("median_area")),
    ("8 cell position + size", lambda c, g: True),
)


def topic_of(col, groups=()):
    return next(t for t, rule in TOPICS if rule(col, groups))


def unit_of(col):
    c = col
    if c == "um_per_px":
        return "µm per pixel"
    if c.endswith(("_abs",)) or "_abs" in c:
        return "area_unit per unit position (per µm for basis distance)" if "slope" in c else \
            ("0-1 position (µm for basis distance)" if "pos" in c else "area_unit")
    if c.endswith("_um2"):
        return "µm²"
    if c.endswith("_um"):
        return "µm"
    if c.endswith("_px"):
        return "pixels²" if "area" in c else "pixels"
    if c.endswith("_deg") or c == "angle_deg":
        return "degrees"
    if c == "angle_rad":
        return "radians"
    if c in ("low_coverage", "few_cells", "stele_found"):
        return "TRUE / FALSE"
    if c in ("x", "y"):
        return "pixels"
    if c == "file_rel":
        return "0-1 (file position, stele side 0)"
    if c == "quarter":
        return "1-4"
    if c == "x_position" or c.startswith("peak_pos"):
        return "0-1 position (0 stele edge, 1 root edge); µm from the stele edge for basis distance"
    if "pos" in c or c.startswith(("r_norm", "lacuna_r_norm")):
        return "0-1 position (0 stele edge, 1 root edge)"
    if "slope" in c:
        return "relative size per unit position (per µm for basis distance)"
    if c in ("rel", "fitted", "epidermis_rel") or c.startswith(("peak_height", "outer_rise")):
        return "x the image's median cell (cortex + epidermis)"
    if "share" in c:
        return "fraction (0-1)"
    if c.startswith(("n_", "min_cells")) or c.endswith(("_n", "_est")) or c in ("tech_rep", "file", "bin",
                                                                                 "radial_bin", "lacuna", "label", "ring"):
        return "count / number"
    return "text"


def order_cols(df, groups=()):
    """per_image / per_plant: columns grouped by topic (identity, checks, root, stele, cortex, pattern, aerenchyma);
    within a topic µm before pixels."""
    cols = list(df.columns)
    return df[sorted(cols, key=lambda c: (topic_of(c, groups), c.endswith("_px") and c != "um_per_px",
                                          cols.index(c)))]


def rooh_version():
    """Commit of the Rooh code that made the numbers ('+ local changes' when the checkout is edited)."""
    import os
    import subprocess
    if os.environ.get("ROOH_VERSION"):
        return os.environ["ROOH_VERSION"]
    here = Path(__file__).resolve().parent
    try:
        h = subprocess.run(["git", "-C", str(here), "log", "-1", "--format=%h %cs"], capture_output=True,
                           text=True, timeout=10).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(here), "status", "--porcelain", "--untracked-files=no"],
                               capture_output=True, text=True, timeout=10).stdout.strip()
        return (h + (" + local changes" if dirty else "")) if h else "unknown"
    except Exception:
        return "unknown"


def run_info(anatomy_dir, settings):
    """Provenance: written into every workbook (sheet run_info) so a table can always be traced to its rules."""
    return pd.DataFrame([
        ("made", dt.datetime.now().isoformat(timespec="minutes"), "date and time these tables were made"),
        ("rooh_version", rooh_version(), "Rooh code commit (github.com/jagdeep2021/cross-section-anatomy)"),
        ("anatomy_folder", str(anatomy_dir), "masks + samples.xlsx these tables were made from"),
        ("cell_model", "Cellpose cpsam" + (f" (custom: {settings['model_path']})" if settings.get("model_path")
                                           else ""), "proposes the cells; corrected by hand where status = reviewed"),
        ("cell_diameter", settings.get("diameter") or "auto", "Cellpose diameter setting (pixels)"),
        ("N_BINS", N_BINS, "radial bins across the tissue outside the endodermis"),
        ("MIN_KEPT_DEG", P.MIN_KEPT_DEG, "below this many degrees kept the image is flagged low_coverage"),
        ("MIN_QUARTER_DEG", P.MIN_QUARTER_DEG, "a quarter needs this many degrees kept to give quarter features"),
        ("MIN_CELLS_PER_BIN", P.MIN_CELLS_PER_BIN, "a bin with fewer cells is flagged few_cells"),
        ("MIN_EPI_DEG", P.MIN_EPI_DEG, "epidermis_deg below this: compare _rel with _rel_cortex"),
        ("FILE_STEP_DEG", P.FILE_STEP_DEG, "ray spacing (degrees) for cell file numbers"),
        ("LACUNA_MIN_CELLS", P.LACUNA_MIN_CELLS, "a lacuna is at least this many median cortical cells"),
        ("GAP_MIN_CELLS", P.GAP_MIN_CELLS, "gaps from this size up to LACUNA_MIN_CELLS = probably a missed cell"),
        ("MIN_PTS", PP.MIN_PTS, "fewest bins with cells for a fitted curve"),
    ], columns=["setting", "value", "what it is"])


def start_here(today, n_img, n_rev, n_plant, no_scale):
    """The first sheet of anatomy_summary: what the files are, the unit, what to check, where each answer is."""
    rows = [
        ("FILES", f"anatomy_summary_{today}.xlsx", "root, stele, metaxylem, cortex, aerenchyma: per_plant, per_image, "
                                                   "per_cell, lacunae"),
        ("", f"cortex_pattern_{today}.xlsx", "cell size from the stele out: profiles, curve features, aerenchyma "
                                             "by position, fitted curves"),
        ("", f"all_cells_{today}.xlsx", "every segmented object with its tissue, position and size (raw material)"),
        ("", f"methods_{today}.txt", "the Methods paragraph for these numbers, with the thresholds used"),
        ("", "figures/", "cortex pattern per group"),
        ("", "README sheet (every workbook)", "every column: topic, unit, meaning. run_info sheet: code version "
                                              "and thresholds"),
        ("THIS RUN", f"{n_img} images, {n_plant} plants", f"{n_rev} of {n_img} images checked by hand (status = "
                                                          "reviewed); the rest are computer outlines only"),
        ("UNIT", "the plant", "images of one plant are technical replicates: analyse plant means (per_plant) or "
                              "use plant as a random effect. Never use cells or images as n."),
        ("CHECK FIRST", "status", "auto rows have not been checked by a person"),
        ("", "um_per_px", (f"{no_scale} image(s) have no scale: their sizes stay in pixels" if no_scale else
                           "every image has a scale")),
        ("", "stele_found", "FALSE = no stele (root tip or cut off): stele / cortex / aerenchyma columns blank"),
        ("", "low_coverage", f"TRUE = fewer than {P.MIN_KEPT_DEG} degrees kept: cortex numbers rest on few cells"),
        ("", "n_cortical_layers_min / _max", "a wide range = stele outline wrong or one side thicker: look"),
        ("", "epidermis_deg", f"below {P.MIN_EPI_DEG}: Cellpose missed most epidermis; compare _rel with "
                              "_rel_cortex"),
        ("", "n_gaps_missed_cell / n_giant_cortex_cells", "likely segmentation misses: review those images"),
        ("", "min_cells_in_bin, few_cells", f"bins under {P.MIN_CELLS_PER_BIN} cells are noisy"),
        ("WHERE IS", "root size", "per_image: root_area_um2, root_diameter_um (kept angles), root_diameter_eq_um"),
        ("", "stele + metaxylem", "per_image: stele_area_um2, stele_share_est, mx_n, mx_mean_area_um2 (whole "
                                  "section, never by sector)"),
        ("", "cortical cell size", "per_image: median_cortex_cell_area_um2 (kept angles)"),
        ("", "cortical cell file number (CCFN)", "per_image: n_cortical_layers"),
        ("", "cortex width / area", "per_image: cortex_width_um, cortex_area_est_um2 (_est = scaled to 360 deg)"),
        ("", "cell size pattern (shape)", "cortex_pattern > profile_features, basis = position, scope = section: "
                                          "peak_pos_rel, peak_height_rel, rise / decay slopes, outer_rise"),
        ("", "pattern by cell file", "the same with basis = file"),
        ("", "quarters", "scope = Q1..Q4, quarter_mean, quarter_sd (spread between quarters)"),
        ("", "aerenchyma amount", "per_image: aerenchyma_share, n_lacunae, lacuna_area_um2"),
        ("", "where aerenchyma forms", "cortex_pattern > aerenchyma_profile (share per radial bin); "
                                       "per_image: lacuna_r_norm_mean; anatomy_summary > lacunae (each one)"),
        ("", "curves for fPCA", "cortex_pattern > spline_curves"),
        ("", "your own analysis", "all_cells (filter sector = kept for cortex numbers)"),
    ]
    return pd.DataFrame(rows, columns=["", "item", "what / where"])


def methods_text(per_image, n_plant, today):
    n, rev = len(per_image), int((per_image.status == "reviewed").sum())
    return f"""Methods: root cross-section anatomy (Rooh {rooh_version()}; tables made {today})

Text to adapt for a Methods section. The thresholds are this run's settings (sheet run_info); edit the species-specific wording.

Segmentation. Cells were segmented with Cellpose (cpsam model; Stringer et al. 2021, Pachitariu et al. 2025). The
root outline was derived from the cell walls and cells, the stele as the cluster of small cells near the centre, and
metaxylem vessels as large round lumens in the stele (area > 5 x the median stele cell, solidity > 0.85). All
outlines were proposed automatically and checked and corrected by hand in Rooh's editor ({rev} of {n} images in
this run; the rest are automatic only). Each pixel belongs to one object. The stele is everything inside the outer edge of the endodermis; cells
touching a 6-pixel band inside the root outline were tagged epidermis; the cells in between are cortex (including
the hypodermis). Objects smaller than 25% of their neighbours' median area with at most 5 neighbours were classed as
intercellular air spaces, not cells.

Damaged tissue. Torn, folded, oblique, out-of-focus or debris-covered parts of the cortex were excluded by hand as
arcs of angle around the stele centre, each with its reason; cells were kept or excluded by their centre. Biology
(aerenchyma, collapsed cells) was never a reason to exclude. Stele and metaxylem traits use the whole section.
Cortex counts and areas from the kept angles are also given scaled to 360 degrees (_est). Images with fewer than
{P.MIN_KEPT_DEG} degrees kept were flagged.

Scale. Pixel sizes (micrometres per pixel) were entered per image from the microscope calibration; areas are to the
mid-wall.

Cortical cell file number. Rays from the stele centre every 5 degrees (kept angles only) counted the cortical cells
crossed (a cell counts when the ray stays in it for at least 2 pixels); CCFN is the median over rays. Every cell
outside the endodermis was also given a file number (1 = next to the stele) as the median of its numbers along
rays every {P.FILE_STEP_DEG} degrees.

Radial pattern of cell size. Each cell's radial position was normalised along its own direction from the stele
centre: 0 at the stele edge, 1 at the root edge (r_norm). All cells outside the endodermis (cortex and epidermis,
treated as one tissue) in the kept angles were binned into {N_BINS} equal r_norm bins; cell area per bin was the
median, divided by the image's median cell area (relative size). A smoothing spline with the smoothing chosen by
generalised cross-validation (scipy make_smoothing_spline), weighted by cells per bin, was fitted through the bin
medians (at least {PP.MIN_PTS} bins with cells), and the peak height and position, the rise slope from the stele to
the peak, the decay slope over 20% of the range past the peak, and the outer rise were read off the curve. The same
was done by cell file instead of position, without the tagged epidermis (a check), and for each quarter around the
stele centre with at least {P.MIN_QUARTER_DEG} degrees kept (quarter mean and SD). Bins with fewer than
{P.MIN_CELLS_PER_BIN} cells were flagged.

Aerenchyma. Lacunae were non-cell gaps in the cortex of at least {P.LACUNA_MIN_CELLS:g} x the median cortical cell
area that touched neither the root edge nor the stele, with their centre in the kept angles. Aerenchyma was
expressed as lacuna area / cortex area kept, and by radial position as lacuna pixels / cortex pixels per r_norm bin.

Statistics. The plant was the experimental unit: per-cell values were summarised per image (medians), images
averaged per plant ({n} images from {n_plant} plants here); cells were never pooled across images.
"""


def make_samples(files: list[Path]) -> pd.DataFrame:
    rows = []
    for f in files:
        rel = Path(str(np.load(f)["image_file"]))
        parts = list(rel.parent.parts)
        plant = (parts[-1] if len(parts) >= 3 else rel.stem).lower()   # "P1" and "p1" are one plant
        groups = parts[:-1] if len(parts) >= 3 else parts
        rows.append({"image_file": str(rel), "plant": plant, "um_per_px": None, "include": "yes", "epidermis": None,
                     **{f"group_{k + 1}": g for k, g in enumerate(groups)}})
    return pd.DataFrame(rows)


def write_samples(anatomy_dir: Path) -> Path:
    """samples.xlsx for the user to fill (never overwritten once it exists)."""
    sfile = anatomy_dir / "samples.xlsx"
    if not sfile.exists():
        make_samples(sorted((anatomy_dir / "masks").glob("*.npz"))).to_excel(sfile, index=False)
        print(f"\nWrote {sfile}\nFill in um_per_px (microns per pixel) for every image, rename group_1, group_2 ... "
              "(e.g. species, genotype), check plant, then click Make results.")
    return sfile


def px_last(df):
    """Pixel columns at the end: the µm columns are the ones to read."""
    px = [c for c in df.columns if c.endswith("_px") and c != "um_per_px"]
    return df[[c for c in df.columns if c not in px] + px]


def describe(col):
    if col in COLUMNS:
        return COLUMNS[col]
    if col.endswith("_px"):
        base = col[:-3]
        for unit in ("_um2", "_um"):
            if base + unit in COLUMNS:
                return COLUMNS[base + unit].split("(µm")[0].strip() + " (pixels)"
        return ""
    for sfx in FEATURE_SETS:                           # profile features: same meaning on another scale
        if col.endswith(sfx) and col[:-len(sfx)] + "_rel" in COLUMNS:
            return COLUMNS[col[:-len(sfx)] + "_rel"].split(";")[0] + f" ({sfx[1:]})"
    return ""


def to_um(d: dict, s: float) -> dict:
    """Pixel traits -> µm: *area*_px -> *_um2, *width_px -> *_um (NaN without a scale)."""
    out = {}
    for k, v in d.items():
        if k.endswith("_px") and "area" in k:
            out[k[:-3] + "_um2"] = v * s * s
        elif k.endswith("_px") and ("width" in k or "diameter" in k):
            out[k[:-3] + "_um"] = v * s
    return out


def main(anatomy_dir: Path, image_figures: bool = False):
    """The project's tables (xlsx) and group figures. image_figures: also one cortex pattern figure per image
    (slow for big sets; the editors draw it live for the image on screen).
    Cortex numbers come from the kept angles only (excluded sectors out, cells by their centre); stele and metaxylem
    are whole-section. Unit of analysis: the plant. Per image first (medians per bin), then plant means; cells are
    never pooled across images."""
    settings = json.loads((anatomy_dir / "settings.json").read_text())
    files = sorted((anatomy_dir / "masks").glob("*.npz"))
    sfile = write_samples(anatomy_dir)
    samples = pd.read_excel(sfile, dtype={"image_file": str, "plant": str})
    samples["include"] = samples["include"].fillna("yes").astype(str).str.lower()
    notes = [c for c in samples.columns if str(c).lower().startswith(NOTE_COLS)]
    groups = [c for c in samples.columns if c not in FIXED and c not in notes]
    samples = samples.set_index("image_file")
    tech_rep = samples.groupby(["plant", *groups], dropna=False).cumcount() + 1   # image number within a plant

    angle = float(settings.get("quadrant_angle_deg", 0))       # only to read old dropped quarters as sectors
    fig_dir = anatomy_dir.parent / "results" / "figures" / "cortex_pattern"
    if image_figures:
        fig_dir.mkdir(parents=True, exist_ok=True)
    img_rows, cell_rows, prof_rows, all_rows, feat_rows, curve_rows, lac_rows, aer_rows = [], [], [], [], [], [], [], []
    grid = np.round(np.linspace(0, 1, 21), 2)
    for f in files:
        z = np.load(f)
        rel = str(z["image_file"])
        if rel not in samples.index or samples.at[rel, "include"] == "no":
            continue
        cells, region, mx, status = P.load_state(f)
        centre = P.centre_of(region, P.load_centre(f))
        props, mxp, tr = P.derive(cells, region, mx, angle, centre)
        arcs = P.load_sectors(f, angle)
        props = P.exclude_sectors(props, arcs)
        st = P.sector_traits(cells, region, props, centre, arcs)
        s = samples.at[rel, "um_per_px"]
        s = float(s) if pd.notna(s) else np.nan
        tech = int(tech_rep[rel])
        meta = {"image_file": rel, "plant": samples.at[rel, "plant"], **{g: samples.at[rel, g] for g in groups}}
        whole = {k: tr[k] for k in ("root_area_px", "root_diameter_eq_px", "stele_area_px", "mx_n", "mx_mean_area_px",
                                    "mx_total_area_px")}

        fl = P.cell_files(cells, props, centre, arcs, region=region)
        props = props.join(fl, on="label")
        props["quarter"] = P.quarter_of(np.degrees(props.angle_rad))
        c_all = props[props.tissue.isin(P.PATTERN_TISSUES)]          # rel reference: the tissue median
        ref = c_all.area_px.median() if len(c_all) else np.nan
        qs = [int(q) for q in st["quarters_used"].split(",") if q]
        scopes = {"section": (props, st["angle_kept_deg"])}
        scopes.update({f"Q{q}": (props[props.quarter == q], P.quarter_kept_deg(arcs, q)) for q in qs})
        um = not np.isnan(s)
        unit = "um2" if um else "px"
        base = {**meta, "tech_rep": tech, "status": status}
        prof = None
        for basis in ("position", "file", "distance"):
            if basis == "distance" and not um:                 # µm bins: needs the scale
                continue
            qfeat = []
            for scope, (sub, kdeg) in scopes.items():
                pr = P.sector_profile(sub, N_BINS, ref, basis, st["n_cortical_layers"], um_per_px=s)
                pr["median_area_um2"] = pr.median_area_px * s * s
                if basis == "position" and scope == "section":
                    prof = pr
                prof_rows.append(pr.rename(columns={"x": "x_position"}).assign(**base, basis=basis, scope=scope))
                allp = pr.dropna(subset=["rel"])
                w_all = allp.n_cells.where(allp.n_cells > 0)
                pc = P.sector_profile(sub, N_BINS, ref, basis, tissues=("cortex",), um_per_px=s).dropna(subset=["rel"])
                ycol = "median_area_um2" if um else "median_area_px"
                e_sub = sub[sub.tissue == "epidermis"]
                f = {**base, "basis": basis, "scope": scope, "area_unit": unit, "kept_deg": kdeg,
                     "low_coverage": st["low_coverage"] if scope == "section" else kdeg < P.MIN_QUARTER_DEG,
                     "n_cells": int(pr.n_cells.sum()), "min_cells_in_bin": int(pr.n_cells.min()) if len(pr) else 0,
                     "epidermis_deg": P.epidermis_deg(sub, arcs),
                     "epidermis_rel": e_sub.area_px.median() / ref if len(e_sub) else np.nan,
                     "epidermis": samples.at[rel, "epidermis"] if "epidermis" in samples.columns else None}
                for sfx, (x, y, ww) in {"_rel": (allp.x, allp.rel, w_all), "_abs": (allp.x, allp[ycol], w_all),
                                        "_rel_cortex": (pc.x, pc.rel, pc.n_cells.where(pc.n_cells > 0))}.items():
                    f.update({k + sfx: v for k, v in PP.profile_features(x, y, ww).items()})
                feat_rows.append(f)
                if scope != "section":
                    qfeat.append(f)
                sp = PP.spline(allp.x, allp.rel, w_all, at=grid)
                curve_rows += [{**base, "basis": basis, "scope": scope, "axis": "rel", "x_position": xv, "fitted": yv}
                               for xv, yv in zip(grid, sp[1] if sp is not None else [np.nan] * len(grid))]
            if qfeat:
                q = pd.DataFrame(qfeat)
                cols = [c for c in q.columns if c.endswith(FEATURE_SETS)]
                for scope, agg in (("quarter_mean", "mean"), ("quarter_sd", "std")):
                    feat_rows.append({**base, "basis": basis, "scope": scope, "area_unit": unit,
                                      "n_quarters": len(q), **q[cols].agg(agg).to_dict()})
        e_all = props[props.tissue == "epidermis"]
        lac, aprof, n_gaps = P.lacunae(cells, region, mx, props, centre, arcs, N_BINS)
        lac_area = int(lac.area_px.sum()) if len(lac) else 0
        lt = {"n_lacunae": len(lac), "lacuna_area_px": lac_area, "n_gaps_missed_cell": n_gaps,
              "aerenchyma_share": lac_area / st["cortex_area_kept_px"] if st["cortex_area_kept_px"] else np.nan,
              "lacuna_r_norm_mean": float(np.average(lac.r_norm, weights=lac.area_px)) if len(lac) else np.nan}
        lac_rows.append(lac.assign(**base, area_um2=lac.area_px * s * s))
        aer_rows.append(aprof.assign(**base))
        if not st["stele_found"]:                       # no stele: whole-section stele traits blank too, not 0
            whole.update({k: np.nan for k in P.NEEDS_STELE if k in whole})
            lt.update({k: np.nan for k in P.NEEDS_STELE if k in lt})
        img_rows.append({**meta, "tech_rep": tech, "status": status, "um_per_px": s,
                         "epidermis": samples.at[rel, "epidermis"] if "epidermis" in samples.columns else None,
                         **{n: samples.at[rel, n] for n in notes},
                         **whole, **st, **lt,
                         "epidermis_rel": e_all.area_px.median() / ref if len(e_all) else np.nan,
                         "min_cells_in_bin": int(prof.n_cells.min()),
                         **to_um({**whole, **st, **lt}, s)})
        ycol = "median_area_um2" if um else "median_area_px"
        if image_figures:
            PP.image_figure(prof.assign(y=prof[ycol]), "Cell area (µm²)" if um else "Cell area (px)",
                            f"{rel}  ·  {status}  ·  {st['angle_kept_deg']:.0f}° kept",
                            fig_dir / (Path(rel).with_suffix("").as_posix().replace("/", "__") + ".png"))

        mxt = pd.DataFrame(measure.regionprops_table(mx, properties=("label", "area", "centroid", "axis_major_length",
                                                                     "axis_minor_length"))) if mx.any() else pd.DataFrame()
        if len(mxt):
            mxt = mxt.rename(columns={"centroid-0": "y", "centroid-1": "x", "area": "area_px",
                                      "axis_major_length": "major_px", "axis_minor_length": "minor_px"})
            mxt = mxt.join(P.radial_position(mxt, region, centre))
            mxt["tissue"] = "metaxylem"
        allc = pd.concat([props, mxt], ignore_index=True)
        allc["sector"] = np.where(allc.tissue.astype(str).str.endswith("_excluded"), "excluded", "kept")
        allc["ring"] = allc.layer if "layer" in allc else np.nan
        allc["quarter"] = P.quarter_of(np.degrees(allc.angle_rad))
        for k in ("file", "file_rel", "n_rays"):
            if k not in allc:
                allc[k] = np.nan
        allc = allc.assign(area_um2=allc.area_px * s * s, major_um=allc.major_px * s, minor_um=allc.minor_px * s,
                           r_um=allc.r_px * s, dist_from_stele_um=allc.dist_from_stele_px * s,
                           angle_deg=np.degrees(allc.angle_rad))
        all_rows.append(allc.assign(**meta, tech_rep=tech, status=status)[
            [*meta, "tech_rep", "status", "label", "tissue", "sector", "ring", "x", "y",
             "angle_rad", "angle_deg", "quarter", "r_px", "r_um", "r_norm", "dist_from_stele_um", "file", "file_rel",
             "n_rays", "area_px", "area_um2", "major_um", "minor_um"]])

        c = props[props.tissue == "cortex"].copy()
        c["area_um2"] = c.area_px * s * s
        c["major_um"], c["minor_um"] = c.major_px * s, c.minor_px * s
        c["dist_from_stele_um"] = c.dist_from_stele_px * s
        c["radial_bin"] = (np.clip(c.r_norm, 0, 0.999) * N_BINS).astype(int) + 1
        cell_rows.append(c.assign(**meta, status=status)[[*meta, "status", "quarter", "r_norm", "radial_bin",
                                                          "file", "file_rel", "dist_from_stele_um",
                                                          "label", "area_px", "area_um2", "major_um", "minor_um"]])

    keys = ["plant", *groups]
    first = lambda df, cols: df[[*cols, *[c for c in df.columns if c not in cols]]]
    per_image = pd.DataFrame(img_rows)
    per_cell = pd.concat(cell_rows, ignore_index=True) if cell_rows else pd.DataFrame()
    per_image_profile = first(pd.concat(prof_rows, ignore_index=True),
                              ["image_file", *keys, "tech_rep", "status", "basis", "scope", "bin"]) \
        if prof_rows else pd.DataFrame()
    pos_section = per_image_profile[(per_image_profile.basis == "position") & (per_image_profile.scope == "section")] \
        if len(per_image_profile) else per_image_profile
    per_plant_profile = (pos_section.groupby([*keys, "bin", "x_position"], dropna=False)
                         .agg(n_images=("image_file", "size"), median_area_um2=("median_area_um2", "mean"),
                              rel=("rel", "mean"), n_cells=("n_cells", "sum")).reset_index()
                         .sort_values([*keys, "x_position"])) if len(per_image_profile) else pd.DataFrame()
    num = [c for c in per_image.columns if c.endswith(("_px", "_um2", "_um", "_deg", "_est", "_rel"))
           or c.startswith(("n_", "mx_", "min_"))]
    per_plant = (per_image.groupby(keys, dropna=False)
                 .agg(n_images=("image_file", "size"),
                      n_reviewed=("status", lambda s: int((s == "reviewed").sum())),
                      n_low_coverage=("low_coverage", "sum"),
                      **{c: (c, "mean") for c in num}).reset_index())
    all_cells = pd.concat(all_rows, ignore_index=True)
    profile_features = first(pd.DataFrame(feat_rows), ["image_file", *keys, "tech_rep", "status", "basis", "scope"])
    lacunae_t = first(pd.concat(lac_rows, ignore_index=True), ["image_file", *keys, "tech_rep", "status"]) \
        if lac_rows else pd.DataFrame()
    aerenchyma_profile = first(pd.concat(aer_rows, ignore_index=True), ["image_file", *keys, "tech_rep", "status"]) \
        if aer_rows else pd.DataFrame()
    spline_curves = pd.DataFrame(curve_rows)
    out_dir = anatomy_dir.parent / "results"
    out_dir.mkdir(exist_ok=True)
    today = dt.date.today().isoformat()
    files_out = {
        f"cortex_pattern_{today}.xlsx": {
            "per_image_profile": (per_image_profile, f"one row per image x basis x scope x bin: every cell outside the "
                                  f"endodermis in the kept angles (centre rule), position: r_norm in {N_BINS} bins (1 = "
                                  "next to the stele); file: one row per cell file. rel = bin median / the image's median "
                                  f"cell area. few_cells = fewer than {P.MIN_CELLS_PER_BIN} cells in the bin"),
            "profile_features": (profile_features, "one row per image x basis (position | file) x scope (section, "
                                 "Q1..Q4, quarter_mean, quarter_sd): smoothing spline (GCV) through the bin medians "
                                 "(weights = cells per bin); features read off the curve. The pattern = every cell "
                                 "outside the endodermis (cortex + epidermis, one tissue). Basis distance = µm from the stele edge (cross-check for lost outer tissue). _rel = relative size (the "
                                 "main pattern), _abs = the same in area_unit, _rel_cortex = without the cells tagged "
                                 "epidermis (a check). "
                                 "peak_height, peak_pos; rise_slope = (peak - first) / distance; decay_slope = local, "
                                 "over 20% of the range past the peak; outer_rise = last value - minimum after the peak. "
                                 f"NaN = fewer than {PP.MIN_PTS} bins with cells or feature undefined (peak at an end). "
                                 "Check low_coverage and min_cells_in_bin before using a row"),
            "aerenchyma_profile": (aerenchyma_profile, f"one row per image x r_norm bin ({N_BINS} bins, 1 = next to the "
                                   "stele): lacuna pixels / cortex pixels in the kept angles. Tests where aerenchyma "
                                   "forms; plant = unit"),
            "spline_curves": (spline_curves, "the fitted relative-size curves at x_position 0, 0.05 ... 1 (long format; NaN outside "
                              "the data): the input for fPCA"),
            "per_plant_profile": (per_plant_profile, "plant x bin: mean of the image medians (display; analyse "
                                  "per_image_profile with plant as the unit)")},
        f"all_cells_{today}.xlsx": {
            "all_cells": (all_cells, "every segmented object: tissue = cortex | stele | epidermis | air_space | "
                          "unassigned | metaxylem, '..._excluded' when its centre is in an excluded sector. angle_rad "
                          "counter-clockwise from the image's right, around the stele centre. r_norm 0 = stele edge, "
                          "1 = root edge along the object's direction. ring = adjacency hops from the stele (gap check, runs high)")},
        f"anatomy_summary_{today}.xlsx": {
            "per_plant": (per_plant, "plant means of per_image; n_reviewed = images checked by hand; "
                          "n_low_coverage = images with too few degrees kept"),
            "per_image": (per_image, "per image. Stele, metaxylem, root_area: whole section. Cortex traits: kept angles "
                          "only; _est = scaled to 360 degrees (estimate). angle_kept_deg and excluded_sectors say what "
                          "was taken out and why"),
            "per_cell": (per_cell, "kept cortical cells (centre rule), one per row"),
            "lacunae": (lacunae_t, "one row per aerenchyma lacuna (centre in the kept angles): area, angle, quarter, "
                        "radial position (mean, inner, outer)")},
    }
    info = run_info(anatomy_dir, settings)
    n_plant = per_image.groupby(keys, dropna=False).ngroups
    for fname, tabs in files_out.items():
        tabs = {name: (order_cols(df, [*groups, *notes]) if name in ("per_image", "per_plant") else px_last(df), note)
                for name, (df, note) in tabs.items()}
        readme = pd.DataFrame([(k, "", "", "", n) for k, (_, n) in tabs.items()],
                              columns=["sheet", "column", "topic", "unit", "what it is"])
        cols = pd.DataFrame([(name, c, topic_of(c, [*groups, *notes])[2:], unit_of(c), describe(c) or
                              ("grouping column from samples.xlsx" if c in groups else
                               "your note from samples.xlsx" if c in notes else ""))
                             for name, (df, _) in tabs.items() for c in df.columns],
                            columns=["sheet", "column", "topic", "unit", "what it is"])
        with pd.ExcelWriter(out_dir / fname) as w:
            if fname.startswith("anatomy_summary"):
                start_here(today, len(per_image), int((per_image.status == "reviewed").sum()), n_plant,
                           int(per_image.um_per_px.isna().sum())).to_excel(w, sheet_name="START_HERE", index=False)
            pd.concat([readme, cols], ignore_index=True).to_excel(w, sheet_name="README", index=False)
            for name, (df, _) in tabs.items():
                df.to_excel(w, sheet_name=name, index=False)
            info.to_excel(w, sheet_name="run_info", index=False)
    paras = methods_text(per_image, n_plant, today).split("\n\n")         # one line per paragraph: pastes into Word
    (out_dir / f"methods_{today}.txt").write_text("\n\n".join(" ".join(x.split("\n")) for x in paras) + "\n")
    if len(per_image_profile) and groups:
        sp, gt = groups[0], groups[-1]            # e.g. species, genotype
        (out_dir / "figures").mkdir(exist_ok=True)
        for spv, d in pos_section.groupby(sp):
            um = d.median_area_um2.notna().any()   # µm² if any image has a scale; images without one are left out
            drop = d.loc[d.median_area_um2.isna(), "image_file"].nunique() if um else 0
            PP.group_figure(d[d.median_area_um2.notna()] if um else d, gt, spv,
                            "median_area_um2" if um else "median_area_px",
                            "Cortical cell area (µm²)" if um else "Cortical cell area (px)",
                            out_dir / "figures" / f"cortex_pattern_{spv}_{today}.png",
                            f"  {drop} image(s) without a scale left out." if drop else "")
    no_scale = per_image.loc[per_image.um_per_px.isna(), "image_file"].tolist()
    if no_scale:
        print(f"\nWARNING: {len(no_scale)} image(s) have no scale (um_per_px empty); their sizes "
              f"stay in pixels and they are left out of the um figures. Type um_per_px in {sfile.name} and run again:")
        print("  " + "\n  ".join(no_scale[:10]) + ("\n  ..." if len(no_scale) > 10 else ""))
    low = per_image.loc[per_image.low_coverage, "image_file"].tolist()
    if low:
        print(f"\nNOTE: {len(low)} image(s) keep fewer than {P.MIN_KEPT_DEG} degrees: low_coverage = TRUE.")
    print(f"{len(per_image)} images ({(per_image.status == 'reviewed').sum()} reviewed), "
          f"{per_plant.shape[0]} plants -> {out_dir}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(Path(sys.argv[1]).resolve(), image_figures="--image-figures" in sys.argv)
