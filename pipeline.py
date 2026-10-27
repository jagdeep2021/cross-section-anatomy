"""Root cross-section anatomy from bright-wall (fluorescence / autofluorescence) images.

The editable state of an image is three label arrays, saved as one .npz per image:
  cells   uint16  one id per cell (Cellpose; corrected by hand in review.py)
  region  uint8   0 outside, 1 root, 2 stele (stele = inside the endodermis' outer edge)
  mx      uint16  one id per metaxylem vessel
Everything else (tissue of each cell, cortical layer, areas) is derived by derive(),
so a hand correction only ever touches those three arrays.

Auto pass (propose):
  1. gray   : percentile stretch + CLAHE, max over colour channels
  2. cells  : Cellpose (cpsam)
  3. root   : walls (Li threshold) + cells, closed, filled, largest object
  4. stele  : small-cell cluster near the centre (Otsu on log cell area), convex hull
  5. mx     : large round lumens (non-wall) inside the stele
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage import exposure, filters, io, measure, morphology, segmentation

_MODEL = None
IMAGE_EXT = (".png", ".tif", ".tiff", ".jpg", ".jpeg")


# ---------------------------------------------------------------- loading

def load_gray(path) -> np.ndarray:
    im = io.imread(path)
    if im.ndim == 3:
        im = im[..., :3].max(axis=2)
    im = im.astype(float)
    lo, hi = np.percentile(im, (1, 99.8))
    im = np.clip((im - lo) / max(hi - lo, 1e-9), 0, 1)
    return exposure.equalize_adapthist(im, clip_limit=0.02)


def wall_mask(gray: np.ndarray) -> np.ndarray:
    sm = filters.gaussian(gray, 1)
    return sm > filters.threshold_local(sm, 51, offset=-0.02)


# ---------------------------------------------------------------- auto pass

def segment_cells(gray: np.ndarray, diameter: float | None = None, model_path: str | None = None) -> np.ndarray:
    global _MODEL
    from cellpose import models
    if _MODEL is None:
        _MODEL = models.CellposeModel(gpu=True, pretrained_model=model_path) if model_path else models.CellposeModel(gpu=True)
    masks, *_ = _MODEL.eval(gray, diameter=diameter)
    return masks.astype(np.uint16)


def find_root(gray: np.ndarray, cells: np.ndarray) -> np.ndarray:
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


def add_mx(gray, mx, cells, y, x) -> int:
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


def centre_of(region: np.ndarray) -> tuple[float, float]:
    """Stele centroid (root centroid if no stele): the origin for quadrants."""
    m = region == 2 if (region == 2).any() else region >= 1
    return ndi.center_of_mass(m) if m.any() else (region.shape[0] / 2, region.shape[1] / 2)


def quadrant_of(y, x, centre, angle_deg: float = 0):
    """Q1 top-right, Q2 top-left, Q3 bottom-left, Q4 bottom-right (image axes, rotated by angle_deg)."""
    theta = np.degrees(np.arctan2(-(np.asarray(y) - centre[0]), np.asarray(x) - centre[1])) - angle_deg
    return (np.floor(np.mod(theta, 360) / 90) + 1).astype(int)


def quadrant_map(shape, centre, angle_deg: float = 0) -> np.ndarray:
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    return quadrant_of(yy, xx, centre, angle_deg)


def derive(cells: np.ndarray, region: np.ndarray, mx: np.ndarray, angle_deg: float = 0):
    """Per-cell table (tissue, cortical layer, quadrant) and per-image traits, in pixels."""
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
        props["quadrant"] = quadrant_of(props.y, props.x, centre_of(region), angle_deg)
    mxp = pd.DataFrame(measure.regionprops_table(mx, properties=("label", "area"))) if mx.any() else pd.DataFrame(columns=["label", "area"])
    traits = {
        "root_area_px": int(root.sum()),
        "stele_area_px": int(stele.sum()),
        "cortex_area_px": int(root.sum() - stele.sum()),
        "n_cortical_layers": int(props.layer.max()) if props.layer.notna().any() else 0,
        "n_cortex_cells": int((props.tissue == "cortex").sum()),
        "n_air_spaces": int((props.tissue == "air_space").sum()),
        "air_space_total_area_px": int(props.area_px[props.tissue == "air_space"].sum()),
        "mx_n": int(len(mxp)),
        "mx_mean_area_px": float(mxp.area.mean()) if len(mxp) else np.nan,
        "mx_total_area_px": int(mxp.area.sum()) if len(mxp) else 0,
    }
    return props, mxp.rename(columns={"area": "area_px"}), traits


def quadrant_traits(cells, region, mx, props, angle_deg: float = 0) -> pd.DataFrame:
    """One row per quadrant: areas, layers, cortex cells, air spaces, metaxylem (pixels)."""
    centre = centre_of(region)
    q = quadrant_map(region.shape, centre, angle_deg)
    mxp = pd.DataFrame(measure.regionprops_table(mx, properties=("label", "area", "centroid"))) if mx.any() else None
    rows = []
    for k in (1, 2, 3, 4):
        inq = q == k
        c = props[(props.get("quadrant") == k) & (props.tissue == "cortex")] if len(props) else props
        a = props[(props.get("quadrant") == k) & (props.tissue == "air_space")] if len(props) else props
        mq = mxp[quadrant_of(mxp["centroid-0"], mxp["centroid-1"], centre, angle_deg) == k] if mxp is not None else []
        rows.append({"quadrant": k,
                     "root_area_px": int((inq & (region >= 1)).sum()),
                     "stele_area_px": int((inq & (region == 2)).sum()),
                     "cortex_area_px": int((inq & (region == 1)).sum()),
                     "n_cortical_layers": int(c.layer.max()) if len(c) else 0,
                     "n_cortex_cells": len(c),
                     "mean_cortex_cell_area_px": float(c.area_px.mean()) if len(c) else np.nan,
                     "n_air_spaces": len(a),
                     "air_space_total_area_px": int(a.area_px.sum()) if len(a) else 0,
                     "mx_n": len(mq)})
    return pd.DataFrame(rows)


def layer_image(cells: np.ndarray, props: pd.DataFrame) -> np.ndarray:
    lut = np.zeros(int(cells.max()) + 1, np.uint16)
    c = props[props.tissue == "cortex"]
    lut[c.label.values] = c.layer.values.astype(np.uint16)
    return lut[cells]


# ---------------------------------------------------------------- storage

def save_state(npz_path, cells, region, mx, status="auto", **meta):
    np.savez_compressed(npz_path, cells=cells, region=region, mx=mx, status=np.array(status), **{k: np.array(v) for k, v in meta.items()})


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
    cy, cx = centre_of(region)
    for a in ax:
        a.axhline(cy, color="white", lw=0.6, ls="--")
        a.axvline(cx, color="white", lw=0.6, ls="--")
    for k, (dy, dx) in zip((1, 2, 3, 4), ((-1, 1), (-1, -1), (1, -1), (1, 1))):
        ax[0].text(cx + dx * 0.42 * gray.shape[1], cy + dy * 0.42 * gray.shape[0], f"Q{k}", color="white", ha="center", va="center")
    ax[0].set_title("root (yellow), stele (cyan), metaxylem (red)")
    ax[1].imshow(lay, cmap="turbo", alpha=0.55, interpolation="nearest")
    ax[1].set_title(f"cortical layers C1..C{traits['n_cortical_layers']}, air spaces (magenta)")
    fig.suptitle(f"{title}   layers={traits['n_cortical_layers']}  MX={traits['mx_n']}  cortex cells={traits['n_cortex_cells']}")
    fig.tight_layout()
    fig.savefig(out_png, dpi=100)
    plt.close(fig)
