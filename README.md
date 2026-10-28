# Cross-section anatomy: propose, correct, summarise

Root cross-sections with bright walls (fluorescence / autofluorescence) → per-cell,
per-layer, per-image and per-plant tables. The computer proposes, a person corrects with
clicks, the computer recomputes. Built 2026-10-27 for Hallie's SL-2026-004 (Ortiz-Ramírez
et al. 2021 maize + Setaria SHR images); written for any new image set.

```
cross_section_anatomy/
├── pipeline.py     propose() auto pass, derive() tissues/layers/areas, click edits, storage
├── run_batch.py    folder of images -> anatomy/masks/*.npz + anatomy/qc/*.png (resumable)
├── review.py       napari app: fix cells / stele / metaxylem; edit the file path per quadrant
├── summarise.py    masks -> results/{cell_size_by_file, all_cells, anatomy_summary}_YYYY-MM-DD.xlsx
│                   + results/figures/ (cortex pattern per species, per image)
├── pattern_plots.py  splines, profile features, the figures (also the live plot in review.py)
├── Install (…), Open app (…)   double-click launchers (Mac .command, Windows .bat)
└── requirements.txt
```

## Run

```bash
cd ~/lab-code/09_phenotyping_pipelines/scripts/cross_section_anatomy
.venv/bin/python run_batch.py "<experiment>/raw"               # -> <experiment>/anatomy/ (+ samples.xlsx to fill)
.venv/bin/python review.py    "<experiment>/anatomy"
.venv/bin/python summarise.py "<experiment>/anatomy"           # -> <experiment>/results/
```
Setup once: double-click `Install (Mac).command` or `Install (Windows).bat`. It uses uv to install
Python 3.14 and requirements.txt into this folder (`.uv/`, `.venv/`; no admin rights; ~1 min, 2.8 GB;
delete the folder to uninstall). Then `Open app (…)`: the top bar has Open folder… (measures new images), Make results, the µm/px scale,
the image list and the editor (1-5: cell, metaxylem, stele, file path, quarters; in quarters mode Ctrl/⌘-click drops a
damaged quarter, click brings it back). The quarter centre is always the traced stele's centre; `run_batch.py` / `summarise.py` still work from the command line.
On a Mac, a downloaded .command opens the first time only with right-click → Open.
Updates: `git pull`, then Open app as usual; it reinstalls packages by itself when requirements.txt changed.
Speed: ~10 s/image with a GPU (Apple silicon, NVIDIA); on a laptop without one (e.g. i3, Intel graphics, 8 GB) it is
tens of minutes per image. Measure there on a fast computer instead: `anatomy/` stores paths relative to itself, so the
experiment folder (OneDrive) opens anywhere with nothing left to measure.

## Editable state

One `.npz` per image: `cells` (uint16 ids), `region` (0 out / 1 root / 2 stele),
`mx` (uint16 vessel ids), `status` (`auto` | `reviewed`), `image_file`, `scale_bar_px`,
`spokes` (4 angles), `path_q` + `path_cells` (quadrant, cell id of each path cell).
Everything else is derived by `derive()`, so the editor and the summary always agree.
`run_batch.py --redo` re-proposes only `auto` images; reviewed work is never overwritten.

## Rules baked in (change here, not downstream)

- **Stele** = inside the outer edge of the endodermis (endodermis counts as stele).
- **Epidermis** = cells within 6 px of the root outline. Everything else between the stele and
  the epidermis is cortex, so the maize hypodermis/exodermis is counted as the outer cortical
  layer(s). Decide before reporting layer number.
- **Cell size by file (primary, decided with Jagdeep 2026-10-27)** = one spoke per quadrant, picked
  where the cortex is most complete (share of the stele-to-edge line inside cortical cells). Cells the
  spoke crosses (>= 30% of their minor axis) form the path; the user adds/removes cells (zig-zag OK) and
  can move the spoke. Path cells sorted by distance from the stele = files 1..n, renumbered on every
  edit. Number of files = path length; image value = median of the 4 quadrants.
- **Rings** (adjacency hops from the stele) are kept only as a gap check in the app: a missing cell
  pushes everything outside it one ring out, so ring counts run high (maize sibs ~13 auto).
- **Radial position** r_norm = 0 at the stele edge, 1 at the root edge along each cell's direction;
  10 bins in per_plant_radial. Continuous, so genotypes with different file numbers stay comparable.
- **Metaxylem** = large round lumens inside the stele (> 5× median stele cell, solidity > 0.85),
  found as dark lumens or as big round Cellpose "cells"; touching pieces merged into one vessel.
- **Air spaces** = objects < 25% of their neighbours' median area with ≤ 5 neighbours (the gaps at
  four-cell junctions), outside the stele. Kept and counted (`n_air_spaces`, area) but not cells:
  excluded from cell sizes and from layer rings. Tune `SPACE_RATIO` / `SPACE_MAX_NB` in pipeline.py.
- **Quadrants** = Q1 top-right, Q2 top-left, Q3 bottom-left, Q4 bottom-right around the stele
  centroid, image axes rotated by `settings.json:quadrant_angle_deg`. Sub-samples within an image;
  Q1 of one image is not Q1 of another.
- **Editing (one model for everything)**: pick the target (keys 1 Cell, 2 Metaxylem, 3 Stele, 4 File path;
  the name follows the cursor). click = add, ⌘-click = remove (Ctrl on Windows/Linux; napari's
  "Control" is ⌘ on Mac), Shift-drag = draw an outline and add, ⌘-Shift-drag = outline and remove,
  drag = pan, Alt-click (file path) = move the spoke. A path outline over an unsegmented gap draws the
  missing cell and adds it to the path. The dock shows the live cortex pattern for the image.
- **Drawing rules**: one owner per pixel, so the drawn
  object takes any overlap and a neighbour left with < 25% of its area is removed (no slivers,
  nothing counted twice). Metaxylem differs: a drawn vessel merges with any vessel it touches (vessels
  never touch, so overlap = same vessel; one count, union area), and clears cell pixels under it.
- **Scale** is always given by the user: `um_per_px` per image in `anatomy/samples.xlsx` (written by
  run_batch.py at the end of the auto pass). Scale bars are not read (detection misread tissue as bars and
  many images have none). Used only by the summary/app for µm values; images left blank stay in px and are
  listed in a warning. Ortiz images: maize 0.7194, Setaria 0.3584 um/px (Jagdeep 2026-10-27).
- **Multi-page TIFF** (stack): page 1 is used (`pipeline.read_plane`; the earlier PNG copies matched page 1).
- **Cortex pattern** (plan: `5_Projects/04_root_anatomy/01_Cortex patterning/6_cortical_radial_profile_analysis/`):
  two views on a 0-1 axis, stele → epidermis. *Along the spoke*: path cells, relative file position
  (file - 0.5) / files in the quarter. *All cortical cells*: r_norm of each cell's centroid, 10 bins; uses
  every cell, so it is the sanity check on the hand-edited spokes. Curves = cubic smoothing splines,
  smoothness by GCV (`scipy.interpolate.make_smoothing_spline`). Group figures: plant means per bin, bins
  with >= 3 plants, band ± 1 SE across plants.
- **Profile features** (`cell_size_by_file → profile_features`, one row per image × quarter): peak height,
  peak position, rise slope, decay slope (local: 20% of the range past the peak), outer rise (last value -
  minimum after the peak), on three axes: `_file` (absolute file), `_rel` (relative file), `_radial`.
- **Replicate** = `plant` in samples.xlsx; per_plant averages that plant's images.

## Status (2026-10-27)

- Auto pass run on all 71 SL-2026-004 images (~15–40 s each on M-series MPS).
  Maize: cells, root, stele good; MX mostly found. Setaria: stele loose, cortex cells missed → needs review.
- Not validated yet. Plan: Hallie reviews all; then compare auto vs reviewed per image
  (layer number, mean cell area per layer) to report the agreement for the auto pass alone.
- Next (v2): fine-tune Cellpose on reviewed masks (`cellpose` train on `cells` arrays) once ~10 per
  species are reviewed; pass the model with `propose(..., model_path=)`.
- Tested: the review app opened on this Mac (napari 0.9.2, PyQt6); clicks/save/next run in code;
  a hand session hasn't been tested.
