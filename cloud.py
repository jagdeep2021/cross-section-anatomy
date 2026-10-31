"""Measure on a rented GPU (Modal, modal.com) instead of this computer: same pipeline.propose(), so the same numbers.

Once, by the account owner:
  .venv/bin/pip install modal && .venv/bin/modal token new     # links this computer to the Modal account
  .venv/bin/modal deploy cloud.py                              # again after any change to pipeline.py

Then:  .venv/bin/python run_batch.py "<experiment>/raw" --cloud

Only the outlining runs remotely (one image at a time, the file's bytes in, three label arrays out); saving, QC
images and file paths stay local, exactly as in a local run. Billing is per second of GPU while an image is
measured; the model (~1.2 GB) is built into the server image, so it never downloads per run.
"""
from __future__ import annotations

from pathlib import Path

import modal

APP, GPU, CPU = "cross-section-anatomy", "L4", 4   # L4: ~3x a T4 for ~1.4x the price; the root/stele steps
                                                    # after Cellpose run on the CPU (default 1/8 core: minutes per image)
HERE = Path(__file__).parent
PKGS = ([l.strip() for l in (HERE / "requirements.txt").read_text().splitlines()
         if l.strip() and l.split("==")[0] not in ("napari", "magicgui", "PyQt6")]  # no window on the server
        if modal.is_local() else [])            # read when building on this computer; the server has no requirements.txt

app = modal.App(APP)


def _fetch_model():
    from cellpose import models
    models.CellposeModel(gpu=False)              # downloads cpsam (~1.2 GB) into the image itself


# the model is baked into the image: a cold start (first image after a few idle minutes) then only boots the
# machine, instead of also fetching 1.2 GB (that made the first image wait ~5 min)
image = (modal.Image.debian_slim(python_version="3.13")
         .apt_install("libgl1", "libglib2.0-0")
         .pip_install(*PKGS)
         .run_function(_fetch_model)
         .add_local_python_source("pipeline"))


@app.function(gpu=GPU, cpu=CPU, image=image, timeout=900, scaledown_window=120)
def propose_remote(name: str, data: bytes, diameter: float | None = None):
    """One image's bytes in, its three label arrays out as one compressed npz (~0.2 MB instead of ~60 MB raw:
    sent raw, the download took minutes per image), plus timings."""
    import io, tempfile, time
    import numpy as np
    import pipeline as P
    p = Path(tempfile.mkdtemp()) / name
    p.write_bytes(data)
    t = time.time()
    _, cells, region, mx = P.propose(p, diameter)
    buf = io.BytesIO()
    np.savez_compressed(buf, cells=cells, region=region, mx=mx)
    return buf.getvalue(), round(time.time() - t, 1)


def propose(path, diameter=None):
    """Drop-in for pipeline.propose on the deployed app: (gray, cells, region, mx). gray is made locally (cheap,
    deterministic) so only the label arrays travel back."""
    import io, time
    import numpy as np
    import pipeline as P
    f = modal.Function.from_name(APP, "propose_remote")
    t = time.time()
    packed, gpu_s = f.remote(Path(path).name, Path(path).read_bytes(), diameter)
    z = np.load(io.BytesIO(packed))
    print(f"[cloud] measured in {gpu_s} s on the server, {time.time() - t:.0f} s round trip", flush=True)
    return P.load_gray(path), z["cells"], z["region"], z["mx"]
