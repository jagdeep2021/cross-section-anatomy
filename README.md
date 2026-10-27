# Cross-section anatomy: propose, correct, summarise

Root cross-sections with bright walls (fluorescence / autofluorescence) → per-cell,
per-layer, per-image and per-plant tables. The computer proposes, a person corrects with
clicks, the computer recomputes. Built 2026-10-27 for Hallie's SL-2026-004 (Ortiz-Ramírez
et al. 2021 maize + Setaria SHR images); written for any new image set.

```
cross_section_anatomy/
├── pipeline.py     propose() auto pass, derive() tissues/layers/areas, click edits, storage
├── run_batch.py    folder of images -> anatomy/masks/*.npz + anatomy/qc/*.png (resumable)
├── review.py       napari app: fix cells / stele / metaxylem by clicks, save as reviewed
├── summarise.py    masks -> results/anatomy_summary_YYYY-MM-DD.xlsx
├── Review.command  double-click launcher (Mac): pick the anatomy folder, opens review.py
└── requirements.txt
```

## Run

```bash
cd ~/lab-code/09_phenotyping_pipelines/scripts/cross_section_anatomy
.venv/bin/python run_batch.py "<experiment>/images_png"        # -> <experiment>/anatomy/
.venv/bin/python review.py    "<experiment>/anatomy"
.venv/bin/python summarise.py "<experiment>/anatomy"           # -> <experiment>/results/
```
Setup once: `python3 -m venv .venv && .venv/bin/pip install -r requirements.txt` (~2 GB, torch).

## Editable state

One `.npz` per image: `cells` (uint16 ids), `region` (0 out / 1 root / 2 stele),
`mx` (uint16 vessel ids), `status` (`auto` | `reviewed`), `image_file`, `scale_bar_px`.
Everything else is derived by `derive()`, so the editor and the summary always agree.
`run_batch.py --redo` re-proposes only `auto` images; reviewed work is never overwritten.

## Rules baked in (change here, not downstream)

- **Stele** = inside the outer edge of the endodermis (endodermis counts as stele).
- **Epidermis** = cells within 6 px of the root outline. Everything else between the stele and
  the epidermis is cortex, so the maize hypodermis/exodermis is counted as the outer cortical
  layer(s). Decide before reporting layer number.
- **Cortical layer** = rings of touching cells counted outward from the stele (C1 = touching the
  endodermis). A missing cell makes the rings step outward too fast, so **automatic layer counts run high**
  (maize sibs ~13 auto). Only reviewed images count.
- **Metaxylem** = large round lumens inside the stele (> 5× median stele cell, solidity > 0.85),
  found as dark lumens or as big round Cellpose "cells"; touching pieces merged into one vessel.
- **Air spaces** = objects < 25% of their neighbours' median area with ≤ 5 neighbours (the gaps at
  four-cell junctions), outside the stele. Kept and counted (`n_air_spaces`, area) but not cells:
  excluded from cell sizes and from layer rings. Tune `SPACE_RATIO` / `SPACE_MAX_NB` in pipeline.py.
- **Quadrants** = Q1 top-right, Q2 top-left, Q3 bottom-left, Q4 bottom-right around the stele
  centroid, image axes rotated by `settings.json:quadrant_angle_deg`. Sub-samples within an image;
  Q1 of one image is not Q1 of another.
- **Scale**: `settings.json:scale_bar_um` ÷ detected bar length (bottom-left bright run).
  Maize Ortiz images: 100 um = 139 px (0.719 um/px). Setaria: no bar → `um_per_px` in samples.xlsx.
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
