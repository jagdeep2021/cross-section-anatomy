"""Auto pass over an image set: propose cells / root / stele / metaxylem for every image.

  python run_batch.py <image_dir> [<out_dir>] [--scale-um 100] [--diameter D] [--redo]

out_dir defaults to <image_dir>/../anatomy and gets
  settings.json       image_dir, scale bar length; edit once per image set
  masks/<key>.npz     editable state per image (status "auto" until saved in review.py)
  qc/<key>.png        overlay to eyeball the proposal
Images that already have a mask are skipped (resume), unless --redo.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from skimage import io

import pipeline as P


def image_key(rel: Path) -> str:
    return str(rel.with_suffix("")).replace("/", "__")


def list_images(image_dir: Path) -> list[Path]:
    return sorted(p.relative_to(image_dir) for p in image_dir.rglob("*")
                  if p.suffix.lower() in P.IMAGE_EXT and not p.name.startswith("."))


def scale_bar_px(path, bottom=0.06, left=0.35) -> int:
    """Longest bright horizontal run in the bottom-left corner (0 if no bar)."""
    im = io.imread(path)
    if im.ndim == 3:
        im = im[..., :3].max(axis=2)
    h, w = im.shape
    sub = im[int(h * (1 - bottom)):, :int(w * left)]
    best = 0
    for row in sub > 0.6 * max(sub.max(), 1):
        d = np.diff(np.concatenate([[0], row.astype(int), [0]]))
        starts, ends = np.where(d == 1)[0], np.where(d == -1)[0]
        if len(starts):
            best = max(best, int((ends - starts).max()))
    return best if best >= 20 else 0


def load_settings(out_dir: Path, image_dir: Path | None = None, **override) -> dict:
    f = out_dir / "settings.json"
    s = json.loads(f.read_text()) if f.exists() else {
        "image_dir": str(image_dir), "scale_bar_um": 100, "diameter": None, "quadrant_angle_deg": 0,
        "note": "scale_bar_um = length of the burned-in scale bar; 0 if the images have none "
                "(then set um_per_px per image in samples.xlsx)."}
    s.update({k: v for k, v in override.items() if v is not None})
    out_dir.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(s, indent=2))
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("image_dir", type=Path)
    ap.add_argument("out_dir", type=Path, nargs="?")
    ap.add_argument("--scale-um", type=float)
    ap.add_argument("--diameter", type=float)
    ap.add_argument("--redo", action="store_true", help="re-propose images that are still 'auto' (never overwrites reviewed ones)")
    a = ap.parse_args()
    image_dir = a.image_dir.resolve()
    out = (a.out_dir or image_dir.parent / "anatomy").resolve()
    s = load_settings(out, image_dir, scale_bar_um=a.scale_um, diameter=a.diameter)
    (out / "masks").mkdir(exist_ok=True)
    (out / "qc").mkdir(exist_ok=True)

    imgs = list_images(image_dir)
    print(f"{len(imgs)} images -> {out}")
    for i, rel in enumerate(imgs, 1):
        key = image_key(rel)
        npz = out / "masks" / f"{key}.npz"
        if npz.exists():
            if not a.redo or P.load_state(npz)[3] == "reviewed":
                continue
        t = time.time()
        gray, cells, region, mx = P.propose(image_dir / rel, s.get("diameter"))
        bar = scale_bar_px(image_dir / rel) if s.get("scale_bar_um") else 0
        P.save_state(npz, cells, region, mx, status="auto", image_file=str(rel), scale_bar_px=bar)
        props, _, tr = P.derive(cells, region, mx)
        P.qc_png(gray, cells, region, mx, props, tr, out / "qc" / f"{key}.png", str(rel))
        print(f"[{i}/{len(imgs)}] {rel}  layers={tr['n_cortical_layers']} mx={tr['mx_n']} bar={bar}px  {time.time() - t:.0f}s", flush=True)


if __name__ == "__main__":
    main()
