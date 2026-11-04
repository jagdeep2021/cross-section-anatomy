"""Root cross-section anatomy from bright-wall (fluorescence / autofluorescence) images.

The editable state of an image is three label arrays, saved as one .npz per image:
  cells   uint16  one id per cell (Cellpose; corrected by hand in review.py)
  region  uint8   0 outside, 1 root, 2 stele (stele = inside the endodermis' outer edge)
  mx      uint16  one id per metaxylem vessel
Everything else (tissue of each cell, cortical layer, areas) is derived by derive(),
so a hand correction only ever touches those three arrays.

Auto pass (propose):
  1. gray   : percentile stretch + CLAHE, max over colour channels; brightfield (dark walls on a bright field)
              is inverted to bright walls and its black frame blanked
  2. cells  : Cellpose (cpsam)
  3. root   : walls (Li threshold) + cells, closed, filled, largest object
  4. stele  : small-cell cluster near the centre (Otsu on log cell area), convex hull
  5. mx     : large round lumens (non-wall) inside the stele
"""
from __future__ import annotations

from pathlib import Path

import os
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage import exposure, filters, io, measure, morphology, segmentation

_MODEL = None
IMAGE_EXT = (".png", ".tif", ".tiff", ".jpg", ".jpeg")


# ---------------------------------------------------------------- loading

def read_plane(path) -> np.ndarray:
    """One 2-D grey plane from any image: colour -> brightest channel; a multi-page TIFF (stack) -> page 1."""
    im = io.imread(path)
    if im.ndim == 4 or (im.ndim == 3 and im.shape[-1] not in (3, 4)):
        im = im[0]                                  # stack (pages first): first page
    if im.ndim == 3:
        im = im[..., :3].max(axis=2)
    return im


def is_brightfield(im: np.ndarray) -> bool:
    """Dark walls on a bright field (brightfield, e.g. colour JPGs from a microscope camera) rather than bright walls
    on a dark background (fluorescence, DAPI). The background is the most common brightness of the blurred image,
    a black frame around a circular field left out: ~0.99 brightfield, <= 0.55 in every fluorescence image so far."""
    sm = ndi.uniform_filter(im, 31)
    h, e = np.histogram(sm[sm > 0.05], 50, (0, 1))
    return h.sum() > 0 and e[np.argmax(h)] > 0.8


def load_gray(path) -> np.ndarray:
    """0-1 grey image with bright walls on a dark background, whatever the microscope gave: the auto pass and the
    click tools look for bright walls, so a brightfield image is inverted and its black frame (outside the circular
    field) blanked, else the bright field itself would be taken for the root."""
    im = read_plane(path).astype(float)
    lo, hi = np.percentile(im, (1, 99.8))
    im = np.clip((im - lo) / max(hi - lo, 1e-9), 0, 1)
    if is_brightfield(im):
        dark = ndi.uniform_filter(im, 15) < 0.05
        lab, _ = ndi.label(dark)
        edge = np.unique(np.r_[lab[0], lab[-1], lab[:, 0], lab[:, -1]])
        frame = np.isin(lab, edge[edge > 0])            # black areas touching the image edge
        im = 1 - im
        im[frame] = 0
    return exposure.equalize_adapthist(im, clip_limit=0.02)


def wall_mask(gray: np.ndarray) -> np.ndarray:
    sm = filters.gaussian(gray, 1)
    return sm > filters.threshold_local(sm, 51, offset=-0.02)


# ---------------------------------------------------------------- auto pass

def segment_cells(gray: np.ndarray, diameter: float | None = None, model_path: str | None = None) -> np.ndarray:
    global _MODEL
    import logging
    for name in ("cellpose.vit", "cellpose.models"):
        logging.getLogger(name).setLevel(logging.ERROR)   # hides 'Could not import CPDINO' (a model we don't use)
    from cellpose import models
    if _MODEL is None:
        _MODEL = models.CellposeModel(gpu=True, pretrained_model=model_path) if model_path else models.CellposeModel(gpu=True)
    masks, *_ = _MODEL.eval(gray, diameter=diameter)
    return masks.astype(np.uint16)


def root_from_cells(cells: np.ndarray, close_px: int = 15) -> np.ndarray:
    """The root outline from the cells alone: every cell merged, gaps closed, holes (lacunae, the stele's
    unsegmented bits) filled, largest piece. The user edits cells and the outline follows."""
    fg = morphology.closing(cells > 0, morphology.disk(close_px))
    fg = ndi.binary_fill_holes(fg)
    lab = measure.label(fg)
    if lab.max() == 0:
        return fg
    biggest = np.argmax(np.bincount(lab.ravel())[1:]) + 1
    return morphology.opening(lab == biggest, morphology.disk(5))


def find_root(gray: np.ndarray, cells: np.ndarray, use_walls: bool = True) -> np.ndarray:
    if not use_walls:
        return root_from_cells(cells)
    sm = filters.gaussian(gray, 2)
    fg = (sm > filters.threshold_li(sm)) | (cells > 0)
    fg = morphology.closing(fg, morphology.disk(15))
    fg = ndi.binary_fill_holes(fg)
    lab = measure.label(fg)
    if lab.max() == 0:
        return fg
    biggest = np.argmax(np.bincount(lab.ravel())[1:]) + 1
    return morphology.opening(lab == biggest, morphology.disk(5))


def find_stele(cells: np.ndarray, root: np.ndarray) -> np.ndarray:
    cy, cx = ndi.center_of_mass(root)
    R = np.sqrt(root.sum() / np.pi)
    p = pd.DataFrame(measure.regionprops_table(cells, properties=("label", "area", "centroid")))
    if p.empty:
        return np.zeros_like(root)
    p["r"] = np.hypot(p["centroid-0"] - cy, p["centroid-1"] - cx) / R
    la = np.log(p.area)
    inner = p.r < 0.6
    cut = filters.threshold_otsu(la[inner].values) if inner.sum() > 5 else la.median()
    small = p.label[(la < cut) & (p.r < 0.7)].values
    m = morphology.closing(np.isin(cells, small), morphology.disk(7))
    m = ndi.binary_fill_holes(m)
    lab = measure.label(m)
    if lab.max() == 0:
        return np.zeros_like(root)
    regs = measure.regionprops(lab)
    best = min(regs, key=lambda g: np.hypot(g.centroid[0] - cy, g.centroid[1] - cx) - 0.01 * g.area)
    return morphology.convex_hull_image(lab == best.label) & root


def find_mx(gray: np.ndarray, cells: np.ndarray, stele: np.ndarray) -> np.ndarray:
    if not stele.any():
        return np.zeros(cells.shape, np.uint16)
    lumen = ~wall_mask(gray) & morphology.erosion(stele, morphology.disk(3))
    lumen = morphology.opening(lumen, morphology.disk(2))
    lab = measure.label(lumen)
    st_cells = np.unique(cells[stele & (cells > 0)])
    st_areas = np.bincount(cells[stele].ravel())[st_cells[st_cells > 0]] if len(st_cells) else np.array([100])
    ref = np.median(st_areas) if len(st_areas) else 100
    out = np.zeros(cells.shape, np.uint16)
    k = 0
    for g in measure.regionprops(lab):
        if g.area > 5 * ref and g.solidity > 0.85 and g.eccentricity < 0.85:
            k += 1
            out[lab == g.label] = k
    return out


def propose(path, diameter: float | None = None, model_path: str | None = None):
    """Auto pass for one image -> (gray, cells, region, mx)."""
    gray = load_gray(path)
    cells = segment_cells(gray, diameter, model_path)
    root = find_root(gray, cells)
    cells[~root] = 0
    stele = find_stele(cells, root)
    mx = find_mx(gray, cells, stele)
    # big round "cells" inside the stele are metaxylem lumens Cellpose kept as cells
    st = pd.DataFrame(measure.regionprops_table(cells * stele, properties=("label", "area", "solidity", "eccentricity")))
    if len(st) > 5:
        big = st[(st.area > 5 * st.area.median()) & (st.solidity > 0.85) & (st.eccentricity < 0.85)]
        for lbl in big.label:
            m = (cells == lbl) & (mx == 0)
            if m.any():
                mx[m] = mx.max() + 1
    if mx.any():  # a vessel is not a cell
        hit = np.unique(cells[mx > 0])
        frac = {c: (mx[cells == c] > 0).mean() for c in hit if c}
        cells[np.isin(cells, [c for c, f in frac.items() if f > 0.5])] = 0
    cells = segmentation.relabel_sequential(cells)[0].astype(np.uint16)
    mx = measure.label(morphology.closing(mx > 0, morphology.disk(2))).astype(np.uint16)  # one id per vessel
    region = root.astype(np.uint8) + stele.astype(np.uint8)
    return gray, cells, region, mx


# ---------------------------------------------------------------- click edits (review.py)

def lumen_at(gray: np.ndarray, y: int, x: int, max_area: int) -> np.ndarray:
    """Dark region enclosed by walls around a click; grown by 2 px to reach mid-wall.
    Falls back to a small disk if the walls leak (region bigger than max_area)."""
    walls = wall_mask(gray)
    if walls[y, x]:  # clicked on a wall: nudge to the darkest pixel nearby
        r = 4
        win = gray[max(y - r, 0):y + r + 1, max(x - r, 0):x + r + 1]
        dy, dx = np.unravel_index(np.argmin(win), win.shape)
        y, x = max(y - r, 0) + dy, max(x - r, 0) + dx
    reg = segmentation.flood(~walls, (y, x))
    if reg.sum() > max_area or reg.sum() < 10:
        yy, xx = np.ogrid[:gray.shape[0], :gray.shape[1]]
        return (yy - y) ** 2 + (xx - x) ** 2 <= 10 ** 2
    return morphology.dilation(reg, morphology.disk(2))


def typical_cell_area(cells: np.ndarray) -> int:
    a = np.bincount(cells.ravel())[1:]
    a = a[a > 0]
    return int(np.median(a)) if len(a) else 2000


def add_cell(gray, cells, y, x) -> int:
    """Add a cell at the click (only into unlabelled pixels). Returns its id (0 if none)."""
    reg = lumen_at(gray, y, x, 8 * typical_cell_area(cells)) & (cells == 0)
    if reg.sum() < 10:
        return 0
    new = int(cells.max()) + 1
    cells[reg] = new
    return new


def vessel_at(gray: np.ndarray, y: int, x: int, max_area: int):
    """A vessel's dark lumen around a click, or None. Vessel walls often have gaps, so a plain flood (lumen_at) runs
    out into the neighbours and fell back to a 10 px dot: here the dark pixels near the click (Otsu in a window) are
    cut at narrow gaps by an opening of growing size until the click's piece is closed and small enough, then tails
    are trimmed and holes filled."""
    r = int(np.sqrt(max_area / np.pi) * 1.5)
    y0, x0 = max(y - r, 0), max(x - r, 0)
    win = filters.gaussian(gray[y0:y + r + 1, x0:x + r + 1], 1.5)
    dark = win < filters.threshold_otsu(win)
    cy, cx = y - y0, x - x0
    for k in (0, 2, 4, 6, 9, 12):
        m = ndi.binary_opening(dark, morphology.disk(k)) if k else dark
        if not m[cy, cx]:
            return None
        lab, _ = ndi.label(m)
        reg = lab == lab[cy, cx]
        if reg.sum() <= max_area and not (reg[0].any() or reg[-1].any() or reg[:, 0].any() or reg[:, -1].any()):
            trim = ndi.binary_opening(reg, morphology.disk(max(k, 4)))      # tails along wall gaps
            if trim[cy, cx]:
                lab, _ = ndi.label(trim)
                reg = lab == lab[cy, cx]
            reg = ndi.binary_dilation(ndi.binary_fill_holes(reg), morphology.disk(2))   # out to mid-wall
            out = np.zeros(gray.shape, bool)
            out[y0:y0 + reg.shape[0], x0:x0 + reg.shape[1]] = reg
            return out
    return None


def add_mx(gray, mx, cells, y, x) -> int:
    reg = vessel_at(gray, y, x, 60 * typical_cell_area(cells))
    if reg is None:                                  # not dark under the click: the old fill (or a small dot)
        reg = lumen_at(gray, y, x, 60 * typical_cell_area(cells))
    new = int(mx.max()) + 1
    mx[reg & (mx == 0)] = new
    cells[reg] = 0
    return new


# ---------------------------------------------------------------- derive (used by auto pass and review)

def cell_adjacency(cells: np.ndarray) -> dict:
    adj: dict[int, set] = {}
    grown = ndi.grey_dilation(cells, footprint=morphology.disk(3))
    m = (grown != cells) & (cells > 0)
    code = np.unique(cells[m].astype(np.int64) * 65536 + grown[m])
    for x, y in zip((code // 65536).tolist(), (code % 65536).tolist()):
        adj.setdefault(x, set()).add(y)
        adj.setdefault(y, set()).add(x)
    return adj


SPACE_RATIO = 0.25   # object < 25% of its neighbours' median area ...
SPACE_MAX_NB = 5     # ... with at most 5 neighbours = air space (intercellular space at a cell junction)


def air_spaces(props: pd.DataFrame, adj: dict) -> set:
    area = dict(zip(props.label, props.area_px))
    out = set()
    for lbl, a in area.items():
        nb = [area[n] for n in adj.get(lbl, ()) if n in area]
        if nb and len(nb) <= SPACE_MAX_NB and a < SPACE_RATIO * np.median(nb):
            out.add(lbl)
    return out


def centre_of(region: np.ndarray, centre=None) -> tuple[float, float]:
    """The origin for sectors and radial position: the user's centre if set, else the stele centroid (root centroid if no stele)."""
    if centre is not None:
        return float(centre[0]), float(centre[1])
    m = region == 2 if (region == 2).any() else region >= 1
    return ndi.center_of_mass(m) if m.any() else (region.shape[0] / 2, region.shape[1] / 2)


def derive(cells: np.ndarray, region: np.ndarray, mx: np.ndarray, angle_deg: float = 0, centre=None):
    """Per-cell table (tissue, cortical layer, angle, radial position) and per-image traits, in pixels, whole section.
    angle_deg is unused (kept for old callers); sectors are applied afterwards (exclude_sectors, sector_traits)."""
    root = region >= 1
    stele = region == 2
    props = pd.DataFrame(measure.regionprops_table(
        cells, properties=("label", "area", "perimeter", "centroid",
                           "axis_major_length", "axis_minor_length")))
    props = props.rename(columns={"centroid-0": "y", "centroid-1": "x",
                                  "axis_major_length": "major_px", "axis_minor_length": "minor_px",
                                  "area": "area_px", "perimeter": "perimeter_px"})
    if props.empty:
        props["tissue"] = []
        props["layer"] = []
    else:
        yi = props.y.round().astype(int).clip(0, cells.shape[0] - 1)
        xi = props.x.round().astype(int).clip(0, cells.shape[1] - 1)
        in_stele = set(props.label[stele[yi, xi]])
        edge_band = root & ~morphology.erosion(root, morphology.disk(6))
        edge = set(np.unique(cells[edge_band])) - {0}
        adj = cell_adjacency(cells)
        air = air_spaces(props, adj) - in_stele
        cortex = set(props.label) - in_stele - edge - air
        boundary = morphology.dilation(stele, morphology.disk(6)) & ~stele
        seed = (set(np.unique(cells[boundary])) - {0}) & cortex
        if stele.any():
            # also the first cell inside the root on each ray from the centre (every 1 deg): a gap at the stele
            # (cortex pulled away, inner file lost in sectioning) would otherwise leave no seed on that side, and
            # with a gap all round every cortical cell came out 'unassigned' (found 2026-11-04, loss_robustness.py)
            inside = cortex & set(props.label[root[yi, xi]])
            c0 = centre_of(region, centre)
            for t in np.arange(0.5, 360, 1.0):
                seed |= set(ray_cells(cells, c0, t, inside)[:1])
        layer = {c: 1 for c in seed}
        frontier, k = set(seed), 1
        while frontier:
            k += 1
            frontier = {n for c in frontier for n in adj.get(c, ()) if n in cortex and n not in layer}
            layer.update({n: k for n in frontier})
        props["tissue"] = np.where(props.label.isin(air), "air_space",
                          np.where(props.label.isin(in_stele), "stele",
                                   np.where(props.label.isin(edge), "epidermis",
                                            np.where(props.label.isin(layer.keys()), "cortex", "unassigned"))))
        props["layer"] = props.label.map(layer)
        props = props.join(radial_position(props, region, centre_of(region, centre)))
    mxp = pd.DataFrame(measure.regionprops_table(mx, properties=("label", "area"))) if mx.any() else pd.DataFrame(columns=["label", "area"])
    traits = {
        "root_area_px": int(root.sum()),
        "root_diameter_eq_px": float(2 * np.sqrt(root.sum() / np.pi)),
        "stele_area_px": int(stele.sum()),
        "cortex_area_px": int(root.sum() - stele.sum()),
        "n_cortical_layers": float(np.median([ray_layers(cells, region, props, centre_of(region, centre), t)
                                              for t in RAY_ANGLES])) if len(props) else 0,
        "n_layers_ring": int(props.layer.max()) if props.layer.notna().any() else 0,
        "n_cortex_cells": int((props.tissue == "cortex").sum()),
        "n_air_spaces": int((props.tissue == "air_space").sum()),
        "air_space_total_area_px": int(props.area_px[props.tissue == "air_space"].sum()),
        "cortex_width_px": float(np.mean(np.subtract(*radial_profile(region, centre_of(region, centre))[::-1]))),
        "mx_n": int(len(mxp)),
        "mx_mean_area_px": float(mxp.area.mean()) if len(mxp) else np.nan,
        "mx_total_area_px": int(mxp.area.sum()) if len(mxp) else 0,
    }
    return props, mxp.rename(columns={"area": "area_px"}), traits


def radial_profile(region: np.ndarray, centre, n_angles: int = 360):
    """Stele and root radius (px) at each whole degree, measured from the centre."""
    t = np.radians(np.arange(n_angles) * 360 / n_angles)
    r = np.arange(0, max(region.shape))
    y = np.round(centre[0] - np.outer(np.sin(t), r)).astype(int)
    x = np.round(centre[1] + np.outer(np.cos(t), r)).astype(int)
    ok = (y >= 0) & (y < region.shape[0]) & (x >= 0) & (x < region.shape[1])
    v = np.where(ok, region[y.clip(0, region.shape[0] - 1), x.clip(0, region.shape[1] - 1)], 0)
    last = lambda m: np.where(m.any(1), m.shape[1] - 1 - np.argmax(m[:, ::-1], 1), 0)
    return last(v == 2).astype(float), last(v >= 1).astype(float)


def radial_position(props: pd.DataFrame, region: np.ndarray, centre) -> pd.DataFrame:
    """r_norm: 0 at the stele edge, 1 at the root edge, along each cell's own direction."""
    r_st, r_rt = radial_profile(region, centre)
    dy, dx = -(props.y - centre[0]), props.x - centre[1]
    deg = np.mod(np.round(np.degrees(np.arctan2(dy, dx))), 360).astype(int)
    r = np.hypot(dy, dx)
    width = np.maximum(r_rt[deg] - r_st[deg], 1)
    return pd.DataFrame({"angle_rad": np.mod(np.arctan2(dy, dx), 2 * np.pi), "r_px": r,
                         "r_norm": (r - r_st[deg]) / width, "dist_from_stele_px": r - r_st[deg],
                         "stele_radius_px": r_st[deg], "root_radius_px": r_rt[deg],
                         "cortex_width_px": width}, index=props.index)


RAY_ANGLES = np.arange(2.5, 360, 5)  # layer counts: one ray every 5 degrees (only the kept ones count), off the
                                     # round angles so a ray does not run along a radial wall drawn at one


def ray_runs(cells, centre, theta_deg: float, labels) -> list:
    """(cell, start px from the centre, length px) for the cells in `labels` crossed by one ray from the centre
    outwards, in order (first crossing). A cell counts when the ray stays in it >= 2 px (a ray grazing a wall corner
    is not a crossing)."""
    t = np.radians(theta_deg)
    r = np.arange(0, max(cells.shape))
    y = np.round(centre[0] - r * np.sin(t)).astype(int)
    x = np.round(centre[1] + r * np.cos(t)).astype(int)
    ok = (y >= 0) & (y < cells.shape[0]) & (x >= 0) & (x < cells.shape[1])
    lab = cells[y[ok], x[ok]]
    if not len(lab):
        return []
    start = np.r_[0, np.nonzero(lab[1:] != lab[:-1])[0] + 1]     # runs of the same label along the ray
    length = np.diff(np.r_[start, len(lab)])
    out, seen = [], set()
    for l, s0, n in zip(lab[start], start, length):
        l = int(l)
        if n >= 2 and l in labels and l not in seen:
            seen.add(l)
            out.append((l, int(s0), int(n)))
    return out


def ray_cells(cells, centre, theta_deg: float, labels) -> list:
    """Cells in `labels` crossed by one ray from the centre outwards, in order (see ray_runs)."""
    return [l for l, _, _ in ray_runs(cells, centre, theta_deg, labels)]


def ray_layers(cells, region, props, centre, theta_deg: float) -> int:
    """Distinct cortical cells crossed by one ray from the centre to the root edge (CCFN on that ray)."""
    return len(ray_cells(cells, centre, theta_deg, set(props.label[props.tissue == "cortex"])))


# ---------------------------------------------------------------- cell files (decided with Jagdeep 2026-11-03)
# Cell file number per cell, from the same kind of rays that count CCFN: along each ray (every
# FILE_STEP_DEG in the kept angles) the cortical and epidermal cells are numbered 1, 2 ... from the stele out; a cell's file is the
# median of the numbers its rays gave it (file_rel: median of (k - 0.5) / n, n = files on that ray). Since 2026-11-04
# a gap on the ray (lacuna, missed cell, gap at the stele edge) counts as the cells it would have held, so cells
# beyond it keep their number (docs/development/loss_robustness.py; on R1, lacuna sections went from epidermis =
# file 9-11 to 14, matching the lacuna-free mature section S29 at 13-14).
FILE_STEP_DEG = 0.5
# The pattern = every cell outside the endodermis (decided with Jagdeep 2026-11-03): the epidermis is the outer part of
# the same tissue, binned by position like any cell, peak anywhere (as in Tyler's maize-cells, which never split it
# out; and "epidermis" here is a position tag, not a cell type). _rel_cortex = the same without the tagged cells.
PATTERN_TISSUES = ("cortex", "epidermis")


def cell_files(cells, props, centre, arcs=(), step: float = FILE_STEP_DEG, tissues=PATTERN_TISSUES,
               region=None, gaps: bool = True) -> pd.DataFrame:
    """index label -> file (median number, 1 = next to the stele; the epidermis is the outermost file), file_rel
    (0-1), n_rays. props: after exclude_sectors (only kept cells are numbered). CCFN (n_cortical_layers) stays
    cortex only.
    gaps (2026-11-04): a gap on a ray (aerenchyma lacuna, a missed cell; with `region`, also a gap between the stele
    edge and the first cell) counts as the cells it would have held, round(gap / the neighbouring cells' length on
    the ray), so the cells beyond it keep their file number. Small gaps (walls, air spaces at junctions: < half a
    cell) add nothing. Without it a lacuna two cells deep numbered every cell outside it two files too low."""
    cortex = set(props.label[props.tissue.isin(tissues)].astype(int))
    r_st = radial_profile(region, centre)[0] if (gaps and region is not None and (region == 2).any()) else None
    nums, rels = {}, {}
    for t in np.arange(step / 2, 360, step):
        if arcs and in_arcs([t], arcs)[0]:
            continue
        runs = ray_runs(cells, centre, t, cortex)
        if not runs:
            continue
        ks, k = [], 0
        for i, (l, s0, n) in enumerate(runs):
            if gaps:
                if i == 0:
                    gap = s0 - r_st[int(t) % 360] if r_st is not None else 0
                    local = n
                else:
                    gap = s0 - (runs[i - 1][1] + runs[i - 1][2])
                    local = (n + runs[i - 1][2]) / 2
                k += int(round(max(gap, 0) / local)) if local > 0 else 0
            k += 1
            ks.append(k)
        for (l, _, _), kk in zip(runs, ks):
            nums.setdefault(l, []).append(kk)
            rels.setdefault(l, []).append((kk - 0.5) / k)
    return pd.DataFrame({"file": {l: float(np.median(v)) for l, v in nums.items()},
                         "file_rel": {l: float(np.median(v)) for l, v in rels.items()},
                         "n_rays": {l: len(v) for l, v in nums.items()}}, dtype=float)


# ---------------------------------------------------------------- sectors (decided with Jagdeep 2026-11-02)
# The whole section counts by default; the user takes out damaged arcs ("excluded sectors"), each with a technical
# reason. Biology (aerenchyma, collapsed or irregular cells) is never a reason: it is measured, not removed.
# A sector is an arc of angles around the stele centre, stored as start/end degrees (counter-clockwise from the
# image's right, the same convention as angle_rad). A cell belongs to the arc its CENTRE lies in (no size bias from
# dropping cells cut by the edge). Sectors touch only cortex, epidermis, air spaces and unassigned cells; the stele
# and metaxylem are always whole-section (an unusable stele is excluded with the image, not with a sector).
# Nothing refers to cell labels, so cell edits never leave a sector stale.

SECTOR_REASONS = ("torn", "folded", "smeared / oblique", "out of focus", "debris", "other")
SECTOR_TISSUES = ("cortex", "epidermis", "air_space", "unassigned")
MIN_KEPT_DEG = 45          # below this the image is flagged low_coverage
MIN_CELLS_PER_BIN = 5      # a profile bin with fewer cortical cells is flagged
MIN_EPI_DEG = 90           # epidermis_deg below this = the outer bins rest on hypodermis more than epidermis: Cellpose
                           # often misses the epidermis (Hallie's auto masks 2026-11-03: 25-300 of 360 deg covered)


def epidermis_deg(props, arcs=()) -> float:
    """Kept degrees (in 5-deg steps) that hold at least one epidermal cell centre."""
    e = props[props.tissue == "epidermis"]
    if not len(e):
        return 0.0
    steps = set((np.floor(np.mod(np.degrees(e.angle_rad), 360) / 5)).astype(int))
    return float(5 * sum(1 for k in steps if not (arcs and in_arcs([5 * k + 2.5], arcs)[0])))


def angle_of(y, x, centre):
    """Degrees 0-360, counter-clockwise from the image's right, around centre (y, x)."""
    return np.mod(np.degrees(np.arctan2(-(np.asarray(y, float) - centre[0]), np.asarray(x, float) - centre[1])), 360)


def in_arcs(deg, arcs) -> np.ndarray:
    """True where the angle lies in any arc (start, end): counter-clockwise from start to end."""
    deg = np.mod(np.asarray(deg, float), 360)
    out = np.zeros(deg.shape, bool)
    for a, b, *_ in arcs:
        a, b = a % 360, b % 360
        out |= ((deg >= a) & (deg < b)) if a <= b else ((deg >= a) | (deg < b))
    return out


def kept_deg(arcs) -> float:
    """Degrees of the circle not covered by any excluded arc (overlaps counted once)."""
    return float(360 - in_arcs(np.arange(0, 360, 0.25) + 0.125, arcs).sum() * 0.25)


def load_sectors(npz_path, angle_deg: float = 0) -> list:
    """Excluded arcs [(start, end, reason)]. Old projects: each dropped quarter becomes a 90 degree arc."""
    z = np.load(npz_path, allow_pickle=False)
    if "excl_start" in z:
        return [(float(a), float(b), str(r)) for a, b, r in zip(z["excl_start"], z["excl_end"], z["excl_reason"])]
    if "dropped_q" in z:
        return [(angle_deg + 90 * (int(q) - 1), angle_deg + 90 * int(q), "dropped quarter (before sectors)")
                for q in z["dropped_q"]]
    return []


def sector_arrays(arcs) -> dict:
    """The arcs as npz fields (save_edit(npz, **sector_arrays(arcs)))."""
    return {"excl_start": np.array([a for a, _, _ in arcs], float), "excl_end": np.array([b for _, b, _ in arcs], float),
            "excl_reason": np.array([r for _, _, r in arcs], dtype="U40")}


def sector_edit(arcs, centre, p0, p1=None, reason="torn", remove=False) -> list:
    """Shared by the desktop editor and the web page. Drag p0 -> p1 (image (y, x)): the smaller arc between the two
    directions is excluded (a wider one takes two drags). remove: the excluded arc(s) holding p0's direction come
    back. A click without a drag changes nothing."""
    arcs = list(arcs)
    t0 = float(angle_of(*p0, centre))
    if remove:
        return [a for a in arcs if not in_arcs([t0], [a])[0]]
    if p1 is None:
        return arcs
    t1 = float(angle_of(*p1, centre))
    d = (t1 - t0) % 360
    a, b = (t0, t1) if d <= 180 else (t1, t0)
    if min(d, 360 - d) < 2:                     # not a sector, a click
        return arcs
    return arcs + [(round(a, 1) % 360, round(b, 1) % 360, reason if reason in SECTOR_REASONS else "other")]


def exclude_sectors(props, arcs):
    """Cells whose centre is in an excluded arc get tissue '<tissue>_excluded' (kept in all_cells, out of every
    cortex number). Stele cells are never touched."""
    if not arcs or not len(props):
        return props
    props = props.copy()
    out = in_arcs(np.degrees(props.angle_rad), arcs) & props.tissue.isin(SECTOR_TISSUES)
    props.loc[out, "tissue"] = props.tissue[out] + "_excluded"
    return props


def sector_traits(cells, region, props, centre, arcs) -> dict:
    """Cortex traits from the kept angles only, plus whole-section estimates scaled by 360 / kept degrees (they
    assume the excluded arcs looked like the kept ones: marked _est). props: after exclude_sectors."""
    centre = centre_of(region, centre)
    kept = kept_deg(arcs)
    yy, xx = np.nonzero(region >= 1)
    outer = region[yy, xx] == 1                                 # root outside the stele (cortex + epidermis)
    keep_px = ~in_arcs(angle_of(yy, xx, centre), arcs)
    scale = 360 / kept if kept > 0 else np.nan
    r_st, r_rt = radial_profile(region, centre)
    keep_deg_mask = ~in_arcs(np.arange(360) + 0.5, arcs)
    rays = [ray_layers(cells, region, props, centre, t) for t in RAY_ANGLES if not in_arcs([t], arcs)[0]]
    c = props[props.tissue == "cortex"]
    a = props[props.tissue == "air_space"]
    stele_px = int((region == 2).sum())
    ctx_kept = int((outer & keep_px).sum())
    out = {"angle_kept_deg": kept,
            "excluded_sectors": "; ".join(f"{a_:.0f}-{b_:.0f} {r}" for a_, b_, r in arcs),
            "low_coverage": bool(kept < MIN_KEPT_DEG),
            "n_cortical_layers": float(np.median(rays)) if rays else np.nan,
            "n_cortical_layers_min": int(min(rays)) if rays else np.nan,
            "n_cortical_layers_max": int(max(rays)) if rays else np.nan,
            "n_cortex_cells": len(c), "n_air_spaces": len(a),
            "air_space_total_area_px": int(a.area_px.sum()),
            "mean_cortex_cell_area_px": float(c.area_px.mean()) if len(c) else np.nan,
            "median_cortex_cell_area_px": float(c.area_px.median()) if len(c) else np.nan,
            "cortex_width_px": float(np.mean((r_rt - r_st)[keep_deg_mask])) if keep_deg_mask.any() else np.nan,
            "cortex_area_kept_px": ctx_kept,
            "cortex_area_est_px": ctx_kept * scale,
            "root_area_est_px": stele_px + ctx_kept * scale,
            "stele_share_est": stele_px / (stele_px + ctx_kept * scale) if kept > 0 else np.nan,
            "n_cortex_cells_est": len(c) * scale, "n_air_spaces_est": len(a) * scale,
            "root_diameter_px": float(2 * np.mean(r_rt[keep_deg_mask])) if keep_deg_mask.any() else np.nan,
            "quarters_used": ",".join(str(q) for q, *_ in QUARTERS if quarter_kept_deg(arcs, q) >= MIN_QUARTER_DEG),
            "n_giant_cortex_cells": int((c.area_px > 5 * c.area_px.median()).sum()) if len(c) else 0,
            "epidermis_deg": epidermis_deg(props, arcs)}
    out["stele_found"] = stele_px > 0
    if not stele_px:            # root tip / stele cut off: no cortex without a stele -> blank, never 0
        out.update({k: np.nan for k in NEEDS_STELE if k in out})
    return out


# traits that only exist once the stele is outlined (blank, not 0, without one: a 0 would drag means down)
NEEDS_STELE = ("n_cortical_layers", "n_cortical_layers_min", "n_cortical_layers_max", "n_cortex_cells",
               "n_air_spaces", "air_space_total_area_px", "mean_cortex_cell_area_px", "median_cortex_cell_area_px",
               "cortex_width_px", "cortex_area_kept_px", "cortex_area_est_px", "stele_share_est",
               "n_cortex_cells_est", "n_air_spaces_est", "n_giant_cortex_cells", "stele_area_px", "mx_n",
               "mx_mean_area_px", "mx_total_area_px", "n_lacunae", "lacuna_area_px", "n_gaps_missed_cell",
               "aerenchyma_share", "lacuna_r_norm_mean")


DIST_BIN_UM = 20.0   # basis "distance": bin width in µm from the stele edge (decided 2026-11-04)


def sector_profile(props, n_bins: int = 10, ref=None, by: str = "position", n_files=None,
                   tissues=PATTERN_TISSUES, um_per_px: float = np.nan) -> pd.DataFrame:
    """The pattern from every kept cell of `tissues` (default: cortex + epidermis, one tissue).
    by = "position": r_norm in n_bins equal bins (1 = next to the stele; r_norm is normalised along each cell's own
         direction, so an uneven cortex is not mixed across layers); x = bin centre.
    by = "file": one row per cell file (props.file rounded; 1 = next to the stele; the epidermis is the outermost
         file); x = the median file_rel of its cells (0-1).
    by = "distance": distance from the stele edge in DIST_BIN_UM bins (needs um_per_px); x = bin centre in µm. No
         outer anchor, so it survives lost outer tissue (loss_robustness.py), but it carries root size: a cross-check.
    Median cell area per row, n per row; rel = median / ref (default: the median of these cells; pass the section's
    for a quarter). props: after exclude_sectors (and with file columns for by="file"). n_files: unused (kept for
    old callers)."""
    c = props[props.tissue.isin(tissues)]
    med = (c.area_px.median() if len(c) else np.nan) if ref is None else ref
    if by == "file":
        c = c.dropna(subset=["file"])
        k = np.floor(c.file + 0.5).astype(int)
        top = int(k.max()) if len(k) else 0
        g = c.groupby(k).agg(median=("area_px", "median"), size=("area_px", "size"), x=("file_rel", "median"))
        g = g.reindex(range(1, top + 1))
        rows = [{"bin": str(i), "x": r["x"], "n_cells": int(0 if pd.isna(r["size"]) else r["size"]),
                 "median_area_px": r["median"]} for i, r in g.iterrows()]
    elif by == "distance":
        w = DIST_BIN_UM / um_per_px if um_per_px and np.isfinite(um_per_px) else np.nan
        c = c[np.isfinite(c.dist_from_stele_px)] if np.isfinite(w) else c.iloc[:0]
        b = (np.clip(c.dist_from_stele_px, 0, None) // w).astype(int) + 1 if len(c) else pd.Series(dtype=int)
        top = int(b.max()) if len(b) else 0
        g = c.area_px.groupby(b).agg(["median", "size"]).reindex(range(1, top + 1))
        rows = [{"bin": str(k), "x": (k - 0.5) * DIST_BIN_UM, "n_cells": int(0 if pd.isna(r["size"]) else r["size"]),
                 "median_area_px": r["median"]} for k, r in g.iterrows()]
    else:
        b = (np.clip(c.r_norm, 0, 0.999) * n_bins).astype(int) + 1
        g = c.area_px.groupby(b).agg(["median", "size"]).reindex(range(1, n_bins + 1))
        rows = [{"bin": str(k), "x": (k - 0.5) / n_bins, "n_cells": int(0 if pd.isna(r["size"]) else r["size"]),
                 "median_area_px": r["median"]} for k, r in g.iterrows()]
    d = pd.DataFrame(rows, columns=["bin", "x", "n_cells", "median_area_px"])
    d["rel"] = d.median_area_px / med
    d["few_cells"] = d.n_cells < MIN_CELLS_PER_BIN
    return d


# ---------------------------------------------------------------- quarters (decided with Jagdeep 2026-11-03)
# Four fixed quarters around the stele centre, on the image axes (Q1 = 0-90 deg, top right, then counter-clockwise),
# as in Tyler's maize-cells quarter images. Excluded arcs take cortex out inside a quarter; a quarter with fewer than
# MIN_QUARTER_DEG kept is not used. Stele and metaxylem are never per quarter (whole section).
QUARTERS = ((1, 0, 90), (2, 90, 180), (3, 180, 270), (4, 270, 360))
MIN_QUARTER_DEG = 45


def quarter_of(deg) -> np.ndarray:
    return (np.mod(np.asarray(deg, float), 360) // 90).astype(int) + 1


def quarter_kept_deg(arcs, q: int) -> float:
    a = np.arange(90 * (q - 1), 90 * q, 0.25) + 0.125
    return float((~in_arcs(a, arcs)).sum() * 0.25)


# ---------------------------------------------------------------- aerenchyma lacunae (decided 2026-11-03)
# A lacuna = a gap in the cortex that is not a cell, >= LACUNA_MIN_CELLS x the median cortical cell, touching
# neither the root edge (torn margins, the band outside the epidermis) nor the stele (tears and stele outline errors;
# lysigenous lacunae form in the mid cortex, the inner cortex stays intact). Derived from the masks (nothing stored):
# the user fixes a false one by drawing the missing cell(s). Gaps of about one cell (GAP_MIN_CELLS .. LACUNA_MIN_CELLS
# x the median) are counted as n_gaps_missed_cell: usually a cell Cellpose missed, a review hint. The thin wall lines
# between cells are broken by an opening first; air spaces (small objects at cell junctions) stay their own class.
# Checked 2026-11-03 on Hallie's maize (no aerenchyma): the 1x-median rule found 6 false lacunae (missed cells, a
# tear at the stele); with these two rules none.
LACUNA_MIN_CELLS = 2.0
GAP_MIN_CELLS = 0.5
def pixel_position(region: np.ndarray, centre):
    """r_norm (0 stele edge, 1 root edge, along each pixel's direction) and angle (deg) for every pixel."""
    r_st, r_rt = radial_profile(region, centre)
    yy, xx = np.indices(region.shape)
    deg = angle_of(yy, xx, centre)
    d = np.mod(np.round(deg), 360).astype(int)
    r = np.hypot(yy - centre[0], xx - centre[1])
    return (r - r_st[d]) / np.maximum(r_rt[d] - r_st[d], 1), deg


def lacunae(cells, region, mx, props, centre, arcs=(), n_bins: int = 10):
    """(per-lacuna table, aerenchyma profile per r_norm bin, n_gaps_missed_cell). props: after exclude_sectors. Only
    lacunae whose centre is in the kept angles count; the profile = lacuna pixels / cortex pixels per bin, kept angles."""
    root, stele = region >= 1, region == 2
    c = props[props.tissue == "cortex"]
    med = float(c.area_px.median()) if len(c) else np.nan
    rn, deg = pixel_position(region, centre)
    ctx = root & ~stele & ~in_arcs(deg, arcs)
    bins = np.clip((np.clip(rn, 0, 0.999) * n_bins).astype(int), 0, n_bins - 1)
    empty = pd.DataFrame(columns=["lacuna", "area_px", "y", "x", "angle_deg", "quarter", "r_norm", "r_norm_inner",
                                  "r_norm_outer"])
    lab = np.zeros(region.shape, np.int32)
    n_gaps = 0
    if np.isfinite(med):
        gap = root & ~stele & (cells == 0) & (mx == 0)
        rad = max(2, int(round(0.1 * np.sqrt(med))))
        core = measure.label(morphology.opening(gap, morphology.disk(rad)))
        core = segmentation.expand_labels(core, rad) * gap
        edge = root & ~morphology.erosion(root, morphology.disk(6))
        near_stele = morphology.dilation(stele, morphology.disk(3))
        rows = []
        for g in measure.regionprops(core):
            m = core == g.label
            cy, cx = g.centroid
            a = float(angle_of(cy, cx, centre))
            if (m & edge).any() or (m & near_stele).any() or (arcs and in_arcs([a], arcs)[0]):
                continue
            if g.area < LACUNA_MIN_CELLS * med:
                n_gaps += int(g.area >= GAP_MIN_CELLS * med)
                continue
            lab[m] = len(rows) + 1
            rows.append({"lacuna": len(rows) + 1, "area_px": int(g.area), "y": cy, "x": cx, "angle_deg": a,
                         "quarter": int(quarter_of(a)), "r_norm": float(rn[m].mean()),
                         "r_norm_inner": float(rn[m].min()), "r_norm_outer": float(rn[m].max())})
        table = pd.DataFrame(rows) if rows else empty
    else:
        table = empty
    tot = np.bincount(bins[ctx], minlength=n_bins)
    hole = np.bincount(bins[ctx & (lab > 0)], minlength=n_bins)
    prof = pd.DataFrame({"bin": np.arange(1, n_bins + 1), "x": (np.arange(n_bins) + 0.5) / n_bins,
                         "cortex_px": tot, "lacuna_px": hole,
                         "lacuna_share": np.where(tot > 0, hole / np.maximum(tot, 1), np.nan)})
    return table, prof, n_gaps


def layer_image(cells: np.ndarray, props: pd.DataFrame) -> np.ndarray:
    lut = np.zeros(int(cells.max()) + 1, np.uint16)
    c = props[props.tissue == "cortex"]
    lut[c.label.values] = c.layer.values.astype(np.uint16)
    return lut[cells]


# ---------------------------------------------------------------- storage

def save_state(npz_path, cells, region, mx, status="auto", **meta):
    """Write to a temporary file, then rename: a crash mid-write never leaves a half file that looks done."""
    npz_path = Path(npz_path)
    tmp = npz_path.with_name(npz_path.name + ".part")
    with open(tmp, "wb") as fh:
        np.savez_compressed(fh, cells=cells, region=region, mx=mx, status=np.array(status),
                            **{k: np.array(v) for k, v in meta.items()})
    os.replace(tmp, npz_path)


def load_centre(npz_path):
    """The centre the user set (y, x), or None: then the stele centroid is used."""
    z = np.load(npz_path, allow_pickle=False)
    return tuple(z["centre"]) if "centre" in z else None


def outline_mask(pts, shape) -> np.ndarray:
    """A drawn outline (points as (y, x) in image pixels) -> the filled area inside it."""
    from skimage.draw import polygon
    pts = np.asarray(pts, float)
    m = np.zeros(shape, bool)
    rr, cc = polygon(pts[:, 0], pts[:, 1], shape)
    m[rr, cc] = True
    return m


def edit_region(region: np.ndarray, m: np.ndarray, target: str, remove: bool = False) -> None:
    """Apply a drawn outline m to region (in place). The one rule set for the desktop editor and the web page.
    Stele: this outline is the stele, inside the root (the old stele goes) (remove: the area leaves the stele).
    Root:  this outline is the root, the stele inside it kept (remove: the area is cut out of the root).
    Both are outlines, not painted areas: drawing one replaces it."""
    if target == "Stele":
        if remove:
            region[m & (region == 2)] = 1
        else:
            region[region == 2] = 1
            region[m & (region >= 1)] = 2
    elif target == "Root":
        if remove:
            region[m] = 0
        else:
            region[~m] = 0
            region[m & (region == 0)] = 1
    else:
        raise ValueError(f"edit_region: no rule for {target}")


def take(arr: np.ndarray, m: np.ndarray, merge: bool) -> int:
    """Give pixels m to a new object in arr (in place); one owner per pixel. merge: objects it touches join it
    (vessels never touch, so overlap = same vessel); else a cell left with < 25% of itself goes (no slivers)."""
    new = int(arr.max()) + 1
    hit = {int(l): int((arr == l).sum()) for l in np.unique(arr[m]) if l}
    arr[m] = new
    for l, before in hit.items():
        if merge:
            arr[arr == l] = new
        elif (arr == l).sum() < 0.25 * before:
            arr[arr == l] = 0
    return new


def click_edit(gray, cells, region, mx, target: str, y: int, x: int, remove: bool = False) -> str:
    """A click at (y, x) on Cell / Metaxylem / Stele, in place. The desktop editor and the web page share it.
    Cell: fill the lumen out to its walls (remove: the cell under the click goes). Metaxylem: same for a vessel.
    Stele / Root: outlines only, so a click just says how to draw one. Returns a note."""
    lbl = int(cells[y, x])
    if target == "Cell":
        if remove and lbl:
            cells[cells == lbl] = 0
        elif not remove:
            add_cell(gray, cells, y, x)
    elif target == "Metaxylem":
        if remove and mx[y, x]:
            mx[mx == mx[y, x]] = 0
        elif not remove and not mx[y, x]:
            add_mx(gray, mx, cells, y, x)
    elif target in ("Stele", "Root"):
        return f"{target}: Shift-drag around the {target.lower()} to draw its outline (it replaces the old one)."
    else:
        raise ValueError(f"click_edit: no rule for {target}")
    return ""


def outline_edit(cells, region, mx, target: str, m: np.ndarray, remove: bool = False) -> None:
    """A drawn outline m on Cell / Metaxylem / Stele / Root, in place (shared by desktop and web).
    Cell: a new cell (remove: clear the area). Metaxylem: a vessel, merged with any it overlaps, no longer a cell
    (remove: clear). Stele / Root: edit_region."""
    if target == "Cell":
        if remove:
            cells[m] = 0
        else:
            take(cells, m, merge=False)
    elif target == "Metaxylem":
        if remove:
            mx[m] = 0
        else:
            take(mx, m, merge=True)
            cells[m] = 0
    elif target in ("Stele", "Root"):
        edit_region(region, m, target, remove)
    else:
        raise ValueError(f"outline_edit: no rule for {target}")


def save_edit(npz_path, **changes) -> None:
    """Save changed arrays into an image's npz, keeping everything else in it, and mark it reviewed."""
    z = np.load(npz_path, allow_pickle=False)
    f = {k: z[k] for k in z.files}
    f.update(changes)
    save_state(npz_path, f.pop("cells"), f.pop("region"), f.pop("mx"), status="reviewed",
               **{k: v for k, v in f.items() if k != "status"})


def load_state(npz_path):
    z = np.load(npz_path, allow_pickle=False)
    return z["cells"], z["region"], z["mx"], str(z["status"])


def qc_png(gray, cells, region, mx, props, traits, out_png, title=""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    lay = layer_image(cells, props).astype(float)
    lay[lay == 0] = np.nan
    fig, ax = plt.subplots(1, 2, figsize=(12, 6.2))
    for a in ax:
        a.imshow(gray, cmap="gray")
        a.axis("off")
    ax[0].contour(region >= 1, [0.5], colors="yellow", linewidths=1)
    ax[0].contour(region == 2, [0.5], colors="cyan", linewidths=1)
    if mx.any():
        ax[0].contour(mx > 0, [0.5], colors="red", linewidths=1.2)
    air = np.isin(cells, props.label[props.tissue == "air_space"].values)
    if air.any():
        ax[1].contour(air, [0.5], colors="magenta", linewidths=0.8)
    ax[0].set_title("root (yellow), stele (cyan), metaxylem (red)")
    ax[1].imshow(lay, cmap="turbo", alpha=0.55, interpolation="nearest")
    ax[1].set_title(f"cortical layers C1..C{traits['n_cortical_layers']}, air spaces (magenta)")
    fig.suptitle(f"{title}   layers={traits['n_cortical_layers']}  MX={traits['mx_n']}  cortex cells={traits['n_cortex_cells']}")
    fig.tight_layout()
    fig.savefig(out_png, dpi=100)
    plt.close(fig)


def um_per_px(anatomy_dir, rel: str) -> float:
    """Scale for one image, always from the user: um_per_px in samples.xlsx; NaN if not filled yet."""
    import pandas as pd
    f = Path(anatomy_dir) / "samples.xlsx"
    if f.exists():
        s = pd.read_excel(f, dtype={"image_file": str}).set_index("image_file")
        if rel in s.index and pd.notna(s.at[rel, "um_per_px"]):
            return float(s.at[rel, "um_per_px"])
    return np.nan


def image_dir(anatomy_dir, settings: dict) -> Path:
    """settings image_dir is stored relative to the anatomy folder, so the experiment works on any computer
    (Mac or Windows, any OneDrive location). An absolute path still works if it exists."""
    d = Path(settings["image_dir"])
    return d if d.is_absolute() and d.exists() else (Path(anatomy_dir) / d).resolve()
