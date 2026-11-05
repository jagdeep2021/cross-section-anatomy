# Rooh: rules for anyone editing this repo

Rooh measures root cross-sections, and its numbers go into theses and papers. So every rule here protects one
thing: **a number must always be traceable to the rule that made it, and the rule to the evidence for it.**
Read this at the start of every session. Claude Code loads it automatically.

## What goes where

| What | Where | Never |
|---|---|---|
| Code | this repo (github.com/jagdeep2021/rooh) | on OneDrive, or copied into an experiment folder |
| Raw images | `<experiment>/raw/` on OneDrive | edited, renamed or moved; work on a copy |
| Masks, `samples.xlsx`, `settings.json` | `<experiment>/anatomy/` | in git (except the small shipped sets in `examples/`) |
| Results (xlsx, figures, methods text) | `<experiment>/results/`, made by Make results | edited by hand; re-run instead |
| How to run, rules in short | `README.md` | |
| **How each number is made now** | `docs/traits.md` | out of date: the checks fail if a `per_image` column is missing |
| **Why** (decisions, evidence, dates) | `docs/development/README.md`, newest at the bottom | rewritten: add a new entry, never edit an old one |
| Figures behind a decision | `docs/development/`, beside the script that draws them | without their script |
| Known-answer tests | `tests/check_numbers.py` (synthetic sections, exact geometry) | skipped |
| Thresholds | named constants at the top of their section in `pipeline.py` | bare numbers inside functions |
| **Anything an editor shows or does** | `web/` only: the desktop runs the same page (`server.py local`, 2026-11-04). `review.py` (napari) is the frozen classic fallback; drawing shared in `pattern_plots.py` | a feature in `review.py` only |
| Rules shared by desktop and web | `pipeline.py` (`derive`, `click_edit`, `outline_edit`, `sector_edit` …) | copied into `review.py` or `web/server.py` |

## Before a session

1. `git switch main && git pull`, then make a branch for the work: `git switch -c <name>/<topic>` (e.g. `amool/brightfield-stele`).
2. Read the part of `docs/traits.md` you will touch and the latest entries of the decision log.
3. Run `.venv/bin/python tests/check_numbers.py` and confirm everything passes **before** you change anything. A failure you started with is not yours to explain later.
4. Testing on real data: copy the experiment (or the images you need) to a scratch folder and work there. Never run on someone's live `anatomy/` folder.

## During

- **If the change alters how a number is defined, write the test first.** Make a synthetic case with a known answer in `check_numbers.py` (or a figure script in `docs/development/` when the question is "which method is better"). Then change the code until it passes.
- One change at a time, so each one can be explained in a single log entry.
- A new output column needs a description in `COLUMNS` (`summarise.py`) and a row in `docs/traits.md`. The checks enforce both.
- Check a change on at least one real project. Look at the images, not just the numbers, and compare before and after.
- **Biology is never a reason to exclude** (aerenchyma, collapsed cells): measure it. Only technical damage is excluded.
- **The plant is the unit.** Never pool cells across images, and never count cells or images as n.

## After a session (the drill)

1. `tests/check_numbers.py` passes in full.
2. If a number's definition changed, update all three in **the same commit**:
   - `docs/traits.md`: the new rule;
   - `COLUMNS` in `summarise.py`;
   - a dated entry in `docs/development/README.md` with **Decision**, **Why**, **Evidence** (numbers or a figure), and which results need a re-run.
3. If how to install or run changed, update `README.md`.
4. Commit with a message that says what changed in the numbers, not just in the code. Push your branch and open a pull request into `main`. Jagdeep reviews and merges. **No force-push; no direct pushes to `main`.**
5. The website (Modal, Jagdeep's account) is deployed by Jagdeep after the merge: `modal app stop cross-section-anatomy`, then `modal deploy cloud.py`.
6. Delete scratch copies of data when you are done.

## Things that look harmless and are not

- `python3` on a lab Mac is not the Rooh Python: always use `.venv/bin/python`.
- Renaming an output column breaks anyone's analysis that reads it. Add the new one and leave the old one for a release, or say so in the log.
- Changing a threshold changes every past result. That needs a log entry with before / after numbers on real data.
- `examples/` masks are shipped data: replace them only on purpose, and say so in the commit.

## Map (where a change goes)

- `pipeline.py`: every measurement and edit rule (`derive`, sectors, `cell_files`, `lacunae`, `lacuna_edit`, …).
- `run_batch.py`: outlining (Cellpose); `locate` decides which folder is the experiment.
- `summarise.py`: Make results; `COLUMNS` describes every output column.
- `pattern_plots.py`: all graphs, including the editor's live plot and file colours.
- `web/server.py` + `web/static/index.html`: **the only editor**. The website runs it on Modal (`cloud.py`); the
  desktop runs the same page on the user's computer (`web/server.py local`, started by the Open app launchers: no
  accounts, a key cookie, projects = folders opened in place). `review.py` (napari) is a frozen fallback.
- Web display caches: `<project>/.cache/<code fingerprint>/`, refreshed when the measurement code changes.

## Environment

- Hellbender (Linux, glibc 2.28): a new compiled dependency must resolve for `manylinux_2_28`
  (`uv pip compile requirements.txt --python-version 3.14 --python-platform x86_64-manylinux_2_28`). PyQt6 is pinned
  6.9.1 on Linux for this; the Linux launcher fetches libxcb-cursor into `.venv`.
- `cloud.py` builds the GPU image from `requirements.txt`: a new package there rebuilds it (a failed rebuild took the
  site down on 2026-11-04). Server-only packages go in `SERVER_PKGS`.
- The page has no build step: check it with `node --check` on its extracted `<script>`.
- Never run `modal deploy`: Jagdeep deploys after merging.

## Open problems (2026-11-04, most important first)

1. The automatic stele misses on brightfield JPGs (the root outline is already fixed by inversion). Every position
   and file number starts at the stele.
2. Gap-aware CCFN: `n_cortical_layers` counts only cells present, so aerenchyma lowers it (B73 R1: 8–9 with lacunae vs
   12). Proposal: the median file number of the outermost cortical file over kept rays, added beside it.
3. Cortex pulled away from the stele: the gap-aware file rule counts the gap as cells (phantom files). Needs a
   synthetic case and a flag (a gap at the stele edge on more than X% of rays).
4. The lacuna rule cannot tell aerenchyma from patches of cells Cellpose missed; reviewers override by hand
   (`n_lacunae_overridden`).
5. Provenance stops at Make results: `run_info` has the Rooh commit, but masks (`.npz`) don't record the version
   that outlined them or last edited them, `run_info` lacks the Cellpose version/model, the page doesn't show the
   live version, and `run_info` still links the old repo name. Proposal: stamp `outlined_with` / `edited_with` in
   each npz (per-image columns), Cellpose version in `run_info`, version in the page footer.

## Working with this lab

Build and run things rather than describing them; judge a change on rendered images. Name a stopping rule before a
sweep. When there are several variants, recommend one and say why. Argue the statistics: deferring is the failure.
Report failures plainly, with the output.
