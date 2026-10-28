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
import sys
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


def locate(folder):
    """Any folder the user picks -> (anatomy_dir or None, image_dir or None).
    Accepts the anatomy folder, the experiment folder (with anatomy/ and/or raw/), or a folder of images."""
    d = Path(folder).resolve()
    for a in (d, d / "anatomy"):
        if (a / "settings.json").exists():
            imgs = P.image_dir(a, json.loads((a / "settings.json").read_text()))
            return a, (imgs if imgs.exists() else None)
    for imgs in (d / "raw", d):
        if imgs.is_dir() and list_images(imgs):
            return imgs.parent / "anatomy", imgs
    return None, None


def n_to_measure(image_dir, anatomy_dir) -> int:
    return sum(not (Path(anatomy_dir) / "masks" / f"{image_key(r)}.npz").exists() for r in list_images(Path(image_dir)))


def batch(image_dir, out_dir=None, diameter=None, redo=False, on_start=None):
    """Outline every image; yields (i, n, line) per image so a caller can show progress.
    Writes samples.xlsx at the end. Returns (StopIteration value) the anatomy folder."""
    image_dir = Path(image_dir).resolve()
    out = Path(out_dir or image_dir.parent / "anatomy").resolve()
    s = load_settings(out, image_dir, diameter=diameter)
    (out / "masks").mkdir(exist_ok=True)
    (out / "qc").mkdir(exist_ok=True)
    imgs = list_images(image_dir)
    if not imgs:
        raise FileNotFoundError(f"No images ({', '.join(P.IMAGE_EXT)}) in {image_dir}")
    for i, rel in enumerate(imgs, 1):
        key = image_key(rel)
        npz = out / "masks" / f"{key}.npz"
        if npz.exists() and (not redo or P.load_state(npz)[3] == "reviewed"):
            yield i, len(imgs), f"{rel.as_posix()}  (done before)"
            continue
        if on_start:
            on_start(i, len(imgs), rel)
        t = time.time()
        gray, cells, region, mx = P.propose(image_dir / rel, s.get("diameter"))
        props, _, tr = P.derive(cells, region, mx, s.get("quadrant_angle_deg", 0))
        spokes, pq, pc = P.auto_paths(cells, region, props, s.get("quadrant_angle_deg", 0))
        P.save_state(npz, cells, region, mx, status="auto", image_file=rel.as_posix(),
                     spokes=spokes, path_q=pq, path_cells=pc)
        P.qc_png(gray, cells, region, mx, props, tr, out / "qc" / f"{key}.png", str(rel))
        yield i, len(imgs), (f"{rel.as_posix()}  files={'/'.join(str(int((pq == q).sum())) for q in (1, 2, 3, 4))} "
                             f"mx={tr['mx_n']}  {time.time() - t:.0f}s")
    S.write_samples(out)
    return out


def be_polite():
    """Run below normal priority with one core left free, so the app's window stays responsive."""
    try:
        import torch
        torch.set_num_threads(max(1, (os.cpu_count() or 2) - 1))
    except Exception:
        pass
    try:
        if sys.platform == "win32":
            import ctypes
            k = ctypes.windll.kernel32
            k.SetPriorityClass(k.GetCurrentProcess(), 0x4000)   # BELOW_NORMAL_PRIORITY_CLASS
        else:
            os.nice(10)
    except Exception:
        pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("image_dir", type=Path)
    ap.add_argument("out_dir", type=Path, nargs="?")
    ap.add_argument("--diameter", type=float)
    ap.add_argument("--redo", action="store_true", help="re-propose images that are still 'auto' (never overwrites reviewed ones)")
    ap.add_argument("--background", action="store_true", help="low priority, one core free (used by the app)")
    a = ap.parse_args()
    if a.background:
        be_polite()
    try:                                         # the app shows this: no GPU means minutes per image, not seconds
        import torch
        dev = "gpu" if torch.cuda.is_available() or torch.backends.mps.is_available() else "cpu"
    except Exception:
        dev = "cpu"
    print(f"[device] {dev}", flush=True)
    start = lambda i, n, rel: print(f"[start {i}/{n}] {rel.as_posix()}", flush=True)
    for i, n, line in batch(a.image_dir, a.out_dir, a.diameter, a.redo, on_start=start):
        print(f"[{i}/{n}] {line}", flush=True)


if __name__ == "__main__":
    main()
