# Rooh (ਰੂਹ): root cross-section anatomy — propose, correct, summarise

Root cross-sections with bright walls (fluorescence / autofluorescence) → per-cell,
per-layer, per-image and per-plant tables. The computer proposes, a person corrects with
clicks, the computer recomputes. Built 2026-10-27 for Hallie's SL-2026-004 (Ortiz-Ramírez
et al. 2021 maize + Setaria SHR images); written for any new image set.

```
cross_section_anatomy/
├── pipeline.py     propose() auto pass, derive() tissues/layers/areas, click edits, storage
├── run_batch.py    folder of images -> anatomy/masks/*.npz + anatomy/qc/*.png (resumable)
├── review.py       napari app: fix cells / stele / metaxylem / root; exclude damaged sectors
├── summarise.py    masks -> results/{cortex_pattern, all_cells, anatomy_summary}_YYYY-MM-DD.xlsx
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
delete the folder to uninstall). Then `Open app (…)`: the top bar has Open folder… (measures new images; or drag a folder or an image onto the window; anything inside raw/ opens the whole experiment), Make results, the µm/px scale,
the image list and the editor (1-5: cell, metaxylem, stele, sectors, root (Shift-drag = this outline is the root, Ctrl/⌘-Shift-drag cuts an area out); in sectors mode
Shift-drag across damaged tissue excludes that arc, Ctrl/⌘-click an excluded arc brings it back). Sectors are measured around the traced
stele's centre; `run_batch.py` / `summarise.py` still work from the command line.
On a Mac, a downloaded .command opens the first time only with right-click → Open.
Updates: `git pull`, then Open app as usual; it reinstalls packages by itself when requirements.txt changed.
Speed: ~10 s/image with a GPU (Apple silicon, NVIDIA); on a laptop without one (e.g. i3, Intel graphics, 8 GB) it is
tens of minutes per image. Measure there on a fast computer instead: `anatomy/` stores paths relative to itself, so the
experiment folder (OneDrive) opens anywhere with nothing left to measure.

## Try it

Two test experiments, 10 images each, to try the app without OneDrive:
- `examples/hallie_SL-2026-004/`: maize sibling and *shr* double, Setaria wild type and *SvSHR1/2* knockout, from
  Ortiz-Ramirez et al. (2021) *Science* 374:1247, doi:10.1126/science.abj2327.
- `examples/amool_root_tip_to_mature/`: grayscale sections from root tip to mature zone, colour JPG and 16-bit DAPI
  TIFF (DAPI copies saved with lossless compression; pixels identical to the originals).

Both come measured (`anatomy/`), like a student's real data: open app, drag one of these folders onto the window
(or Open folder…) and it opens straight into editing, nothing to measure. Speed test: delete that folder's
`anatomy/` and open it again; it measures the 10 images and the title bar shows the time per image. Restore with
`git checkout examples`.

## Measure in the cloud (Modal)

`cloud.py` runs the same `pipeline.propose()` on a rented GPU (Modal, L4, 4 CPU; model baked into the image).
`run_batch.py <experiment>/raw --cloud`; saving and QC stay local, so outputs match a local run.
Measured 2026-10-31: 18-30 s per fluorescence/DAPI image and ~65 s per 12 MP brightfield JPG on the server, plus
~65 s once after a few idle minutes (machine start). Results travel compressed (~0.2 MB). ~$0.004-0.015 per image.
Setup (account owner): `.venv/bin/pip install modal`, `.venv/bin/modal token new`, `.venv/bin/modal deploy cloud.py`
(again after any change to pipeline.py; if old machines still answer with the old code, `modal app stop
cross-section-anatomy` first). Stopping a local run does not stop calls already on the server.

## Compare two sets of results

`python compare.py <reference.xlsx> <other.xlsx>`: two `anatomy_summary_*.xlsx` of the same images (auto pass vs a
review, or reviewer vs reviewer). Matches on image_file, compares pixel columns (scale can't fake a difference),
and writes `compare_<other>_vs_<reference>.xlsx` next to `other`: per trait n changed, mean difference with 95%
limits of agreement (Bland & Altman 1986), mean |% difference| and Lin's CCC (Lin 1989); per image the values.
For the examples, the reference is `summarise.py examples/<set>/anatomy` on the untouched (auto) copy from git.

## Web version (live on Modal; browser editing in progress)

The same pipeline behind a browser page, for laptops without a GPU and users outside the lab. `web/server.py`
(FastAPI) keeps each project in the desktop layout (`raw/`, `anatomy/`, `results/`) and calls `run_batch.batch`,
`derive()` and `summarise.main` unchanged, so web and desktop numbers are identical (checked on 6 SL-2026-004 images,
2026-10-29: every measurement column equal). One GPU worker measures projects in turn.

```bash
.venv/bin/pip install -r web/requirements-web.txt
.venv/bin/python web/server.py invite someone@example.com [--admin]   # prints their login link
.venv/bin/python web/server.py serve --port 8000                       # data in $CSA_DATA (default ./web_data)
```
Done: invite-only accounts (only the token hash is stored), projects, folder upload (folders become groups), queued
GPU measuring with progress, image + outline layers with pan/zoom, per-image numbers, scale, Make results, download zip.
Live: https://jagdeep2021--rooh.modal.run (deployed from `cloud.py` with
`modal deploy cloud.py`; a CPU machine serves the page, measuring goes to the GPU function, projects live on the
Modal volume `csa-web-data`, saved after every change).
Access: share the site address; a new person types name + email and asks, an admin approves them in the People
panel (that browser then works). Or send a personal login link: `python cloud.py invite <email> [--admin]`.
Anyone can type any email, so approve only people you know.
Tested 2026-10-31: login, project, upload, measure, Make results, download.
Layout as the desktop: one top bar in work order (project, add images, Outline, µm/px, image ▾ ‹ › Save + next,
Make results), the image, and the right panel. Cells and masks fade separately.
Editing in the browser (right panel): pick Cell / Vessel / Stele / Sector / Root (keys 1-5); drag draws, click adds,
hold Ctrl/⌘ to remove (Sector: drag around the stele across the damage, with a reason; Ctrl/⌘-click a grey arc = back); Space-drag or right-drag moves the view, scroll zooms, F fits. Below: this image's cortex
pattern (live) and its traits as small bars against the project median (orange = under 2/3 or over 1.5x).
Make results = the project's Excel files + group graphs; per-image graphs only if ticked (`summarise.py --image-figures`).
🗑 deletes a project: its dialog first offers ⬇ the whole project (images, masks, settings, results) as one zip;
the project is moved to `<data>/trash/<user>/<name>_<date-time>/`, not erased (an admin can move it back).
Same rules as the desktop (`pipeline.click_edit` / `outline_edit` / `sector_edit`, shared). Next: per-user quotas;
a container for deployment (MU server or a GPU service, whichever IT approves).

## Editable state

Saving: Save + next (Shift-S) saves and marks the image `reviewed`; leaving an edited image any other way
(< > , . the image list, another folder, Make results, closing the window) saves it too. An image you only look at stays `auto`.

One `.npz` per image: `cells` (uint16 ids), `region` (0 out / 1 root / 2 stele),
`mx` (uint16 vessel ids), `status` (`auto` | `reviewed`), `image_file`, `scale_bar_px`,
`excl_start` / `excl_end` / `excl_reason` (excluded arcs, degrees counter-clockwise from the image's right; none =
the whole section counts), `centre` (optional, else the stele centroid). Older files may carry `dropped_q` (read as
90° arcs) and `spokes` / `path_q` / `path_cells` (ignored since 2026-11-02).
Everything else is derived by `derive()`, so the editor and the summary always agree.
`run_batch.py --redo` re-proposes only `auto` images; reviewed work is never overwritten.

## Rules baked in (change here, not downstream)

- **Stele** = inside the outer edge of the endodermis (endodermis counts as stele).
- **Epidermis** = cells within 6 px of the root outline. Everything else between the stele and
  the epidermis is cortex, so the maize hypodermis/exodermis is counted as the outer cortical
  layer(s). Decide before reporting layer number.
- **Sectors (decided with Jagdeep 2026-11-02; replaced one spoke per quarter)**. The whole section counts by
  default; damaged tissue is taken out as an arc of angles around the stele centre, each with a technical reason
  (torn, folded, smeared / oblique, out of focus, debris, other). Any angle, any number of arcs (overlaps count once).
  - A cell is in or out by its **centre** (no half cells: cells cut by the edge are neither trimmed nor dropped
    as a size class, so the arc edge biases nothing).
  - Sectors touch only cortex, epidermis, air spaces and unassigned cells. **Stele and metaxylem are always
    whole-section**; an unusable stele means excluding the image, not a sector.
  - **Biology is never a reason**: aerenchyma, collapsed or irregular cells are measured, not removed.
  - Excluded cells stay in `all_cells` as `<tissue>_excluded` (`sector` = excluded).
  - Why: on R1 (04_root_anatomy `12_sector_vs_spoke`) the pattern from all cortical cells was far more
    reliable than from one file: split-half r 0.87 vs 0.52, sectors better in 16 of 19 sections (Wilcoxon
    *P* = 0.002), ~340 vs ~31 cells per section.
- **Quarters (2026-11-03, alongside sectors)**: Q1 = 0-90° (top right), counter-clockwise, around the stele centre.
  Excluded arcs take cortex out inside a quarter; a quarter with < 45° kept is not used. Every feature is given for
  the section (all kept cells), each quarter, and the quarter mean / SD (Tyler's maize-cells unit). Stele and
  metaxylem: whole section only. Editors draw quarter lines and excluded arcs on the cortex only (stele edge -> root edge).
- **Cell file number (2026-11-03)**: rays every 0.5° in the kept angles number the cortical cells 1, 2 ... from the
  stele out; a cell's `file` = median of its numbers, `file_rel` = median of (k - 0.5) / n. Gives the second pattern
  (`basis = file`: median area per file) next to the position pattern (`basis = position`).
- **Epidermis is part of the tissue (decided with Jagdeep 2026-11-03)**: the pattern = every cell outside the
  endodermis, binned by position like any cell (the epidermis lands in the outermost bin / is the outermost file),
  peak searched over the whole curve; `rel` = bin median / the median of those cells. As in Tyler's maize-cells, and
  "epidermis" is a position tag (within 6 px of the outline), not a cell type we can trust. The tag stays on each cell;
  `_rel_cortex` = the same features without the tagged cells (a check for reviewers). Cellpose often misses epidermal
  cells (auto masks covered 25-300 of 360°): `epidermis_deg` says how much; below 90 compare `_rel` with `_rel_cortex`.
  CCFN (`n_cortical_layers`) stays cortex only. `samples.xlsx` column `epidermis` (yes / no / empty) = intact or not.
- **Aerenchyma lacunae (2026-11-03)**: a non-cell gap in the cortex >= 2 x the median cortical cell, touching neither
  the root edge nor the stele, centre in the kept angles. Derived from the masks (fix a false one by drawing the missing
  cell). Per lacuna: area, angle, quarter, r_norm (mean, inner, outer) (`lacunae` sheet); per image x r_norm bin: lacuna
  share of the cortex (`aerenchyma_profile`), to test where aerenchyma forms. One-cell gaps (0.5-2 x) =
  `n_gaps_missed_cell`, a review hint. Cortex cells > 5 x median = `n_giant_cortex_cells` (a lacuna taken for a cell?).
- **Root diameter (2026-11-03)**: `root_diameter` = 2 x mean root radius over the kept degrees (from the stele centre);
  `root_diameter_eq` = equal-area circle. The outline stays walls + cells: a cells-only outline was 1.4-11.6% smaller on
  Hallie's sections because Cellpose misses epidermal cells (`root_from_cells` kept, unused).
- **Cortical layers** = cortical cells crossed by a ray from the centre out (a cell counts when the ray stays in
  it >= 2 px), rays every 5° (offset 2.5°) in the kept angles; `n_cortical_layers` = median, `_min` / `_max` the
  range (one side can carry an extra layer; a wide range on a real image is worth a look).
- **Rings** (adjacency hops from the stele) are kept only as a gap check in the app: a missing cell
  pushes everything outside it one ring out, so ring counts run high (maize sibs ~13 auto).
- **Radial position** r_norm = 0 at the stele edge, 1 at the root edge along each cell's direction (so an
  uneven or off-centre cortex is not mixed across layers). Continuous, so genotypes with different layer numbers
  stay comparable.
- **Metaxylem** = large round lumens inside the stele (> 5× median stele cell, solidity > 0.85),
  found as dark lumens or as big round Cellpose "cells"; touching pieces merged into one vessel.
- **Air spaces** = objects < 25% of their neighbours' median area with ≤ 5 neighbours (the gaps at
  four-cell junctions), outside the stele. Kept and counted (`n_air_spaces`, area) but not cells:
  excluded from cell sizes and from layer rings. Tune `SPACE_RATIO` / `SPACE_MAX_NB` in pipeline.py.
- **Editing (one model for everything)**: pick the target (keys 1 Cell, 2 Metaxylem, 3 Stele, 4 Sectors, 5 Root;
  the name follows the cursor). click = add, ⌘-click = remove (Ctrl on Windows/Linux; napari's
  "Control" is ⌘ on Mac), Shift-drag = draw an outline and add, ⌘-Shift-drag = outline and remove,
  drag = pan. Sectors: Shift-drag from one side of the damage to the other (the smaller arc; a wider one takes two
  drags), reason from the bar; ⌘-click inside an excluded arc = back. The dock shows the live cortex pattern.
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
  every kept cortical cell, r_norm in 10 equal bins (1 = next to the stele), **median** area per bin (robust to
  a few mis-segmented cells), n per bin; `rel` = bin median / the image's median cortical cell area (shape without
  overall size). The **epidermis is its own point** past 1 (a different tissue: never joined to the cortex curve;
  `epidermis_rel`). Curves = cubic smoothing splines through the bin medians weighted by cells per bin, smoothness
  by GCV (`scipy.interpolate.make_smoothing_spline`). Group figures: images -> plant means per bin, bins with
  >= 3 plants, band ± 1 SE across plants; epidermis as a square ± SE.
- **Profile features** (`cortex_pattern → profile_features`, one row per image): peak height, peak position, rise
  slope, decay slope (local: 20% of the range past the peak), outer rise (last value - minimum after the peak), on
  `_rel` (relative size: the shape), `_abs` (µm², or px without a scale) and `_rel_epi` (relative, the epidermis as
  the last point). `spline_curves` = the fitted `rel` curves on a 0-1 grid, for fPCA.
- **Scaled estimates** (`_est`): kept value × 360 / kept degrees (cortex area, root area, cortex cells, air spaces);
  `stele_share_est` = stele / estimated root area. They assume the excluded arcs looked like the kept ones.
- **Flags**: `low_coverage` = fewer than 45° kept; `few_cells` = a bin with < 5 cortical cells;
  `min_cells_in_bin` per image. Make results prints the low-coverage images.

### Reading the results (caveats to honour)

1. **Unit = the image, then the plant.** Never pool cells across images or plants: bin medians per image, then
   plant means. Cells are pseudo-replicates; tests use plant (or section) as the unit.
2. **Excluded arcs are not random.** Damage tends to hit the outer cortex and the thin side of a section; a
   pattern from 200° is a sample of the root, not the root. Report `angle_kept_deg` and check that conclusions
   hold on images with nothing excluded.
3. **`_est` columns are estimates**, not measurements: an asymmetric root (lateral root side, flattened section)
   breaks the "excluded looked like kept" assumption. Use kept values (`cortex_area_kept`, cell sizes) for tests;
   `_est` for whole-root context only.
4. **Low coverage and few-cell bins**: below 45° the curve rests on a handful of files; a bin with < 5 cells is a
   point, not a median. Drop or down-weight; say which in Methods.
5. **Oblique sections** stretch cells along one direction: a smeared / oblique arc should be excluded, not kept
   and averaged away. A whole oblique section is an image to exclude (`include` = no in samples.xlsx).
6. **r_norm is relative**: the same bin is a different absolute depth in a thick and a thin cortex. Compare shapes
   on `_rel`, sizes on `_abs`, and report layer number alongside.
7. **The epidermis point is not part of the cortex curve**: compare it as its own trait (`epidermis_rel`), and
   check it exists (epidermis lost in sectioning leaves it NaN).
8. **Layer range**: `n_cortical_layers` is the median over rays; a wide `_min`–`_max` means an uneven cortex or a
   segmentation gap. Look before using the number.
- **Replicate** = `plant` in samples.xlsx; per_plant averages that plant's images.

## Development protocols

Why the measurements are defined the way they are, and how they differ from earlier lab pipelines: [docs/development/](docs/development/README.md).

## Status (2026-11-02)

- Auto pass run on all 71 SL-2026-004 images (~15–40 s each on M-series MPS).
  Maize: cells, root, stele good; MX mostly found. Setaria: stele loose, cortex cells missed → needs review.
- 2026-11-02: spokes / quarters replaced by sectors (desktop, web, summary, checks). `tests/check_numbers.py`:
  105 pass; the 2 fails are case 11, the aerenchyma tag (not built: aerenchyma is measured as air spaces,
  there is no hand tag yet). Make results on the BW root-tip set (88 images) runs; old dropped quarters read as arcs.
- Not validated yet. Plan: Hallie reviews all; then compare auto vs reviewed per image
  (layer number, mean cell area per layer) to report the agreement for the auto pass alone.
- Next (v2): fine-tune Cellpose on reviewed masks (`cellpose` train on `cells` arrays) once ~10 per
  species are reviewed; pass the model with `propose(..., model_path=)`.
- Tested: the review app opened on this Mac (napari 0.9.2, PyQt6); clicks/save/next run in code;
  a hand session hasn't been tested.
