"""Auto pass over an image set: propose cells / root / stele / metaxylem for every image.

  python run_batch.py <image_dir> [<out_dir>] [--diameter D] [--redo]

out_dir defaults to <image_dir>/../anatomy and gets
  settings.json       image_dir (relative), cell diameter, quadrant angle
  samples.xlsx        one row per image: fill um_per_px (the scale; always given by the user), plant, groups
  masks/<key>.npz     editable state per image (status "auto" until saved in review.py)
  qc/<key>.png        overlay to eyeball the proposal
Images that already have a mask are skipped (resume), unless --redo.
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np

import pipeline as P
import summarise as S


def image_key(rel: Path) -> str:
    return str(rel.with_suffix("")).replace("/", "__")


def list_images(image_dir: Path) -> list[Path]:
    return sorted(p.relative_to(image_dir) for p in image_dir.rglob("*")
                  if p.suffix.lower() in P.IMAGE_EXT and not p.name.startswith("."))


def load_settings(out_dir: Path, image_dir: Path | None = None, **override) -> dict:
    f = out_dir / "settings.json"
    s = json.loads(f.read_text()) if f.exists() else {
        "image_dir": Path(os.path.relpath(image_dir, out_dir)).as_posix(), "diameter": None, "quadrant_angle_deg": 0}
    s.update({k: v for k, v in override.items() if v is not None})
    out_dir.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(s, indent=2))
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("image_dir", type=Path)
    ap.add_argument("out_dir", type=Path, nargs="?")
    ap.add_argument("--diameter", type=float)
    ap.add_argument("--redo", action="store_true", help="re-propose images that are still 'auto' (never overwrites reviewed ones)")
    a = ap.parse_args()
    image_dir = a.image_dir.resolve()
    out = (a.out_dir or image_dir.parent / "anatomy").resolve()
    s = load_settings(out, image_dir, diameter=a.diameter)
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
        props, _, tr = P.derive(cells, region, mx, s.get("quadrant_angle_deg", 0))
        spokes, pq, pc = P.auto_paths(cells, region, props, s.get("quadrant_angle_deg", 0))
        P.save_state(npz, cells, region, mx, status="auto", image_file=rel.as_posix(),
                     spokes=spokes, path_q=pq, path_cells=pc)
        P.qc_png(gray, cells, region, mx, props, tr, out / "qc" / f"{key}.png", str(rel))
        print(f"[{i}/{len(imgs)}] {rel.as_posix()}  files={'/'.join(str(int((pq == q).sum())) for q in (1, 2, 3, 4))} mx={tr['mx_n']}  {time.time() - t:.0f}s", flush=True)

    S.write_samples(out)


if __name__ == "__main__":
    main()
