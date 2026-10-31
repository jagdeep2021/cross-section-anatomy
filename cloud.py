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

APP, CPU = "cross-section-anatomy", 4   # the root/stele steps after Cellpose run on the CPU (default 1/8 core: minutes)
GPU = ["L4", "A10G", "T4"]   # first free one wins: L4 is ~3x a T4 for ~1.4x the price, but on 2026-10-31 no L4 was
                             # free for 5+ min and measuring sat at "Waiting for the GPU"
HERE = Path(__file__).parent
PKGS = ([l.strip() for l in (HERE / "requirements.txt").read_text().splitlines()
         if l.strip() and l.split("==")[0] not in ("napari", "magicgui", "PyQt6")]  # no window on the server
        if modal.is_local() else [])            # read when building on this computer; the server has no requirements.txt

app = modal.App(APP)
if not modal.is_local():
    # load the pipeline (numpy, scipy, skimage) when a machine starts, not inside a request: a request cancelled
    # mid-import left numpy half-loaded and every later call on that machine failed (RecursionError)
    import pipeline  # noqa: F401


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


# ---------------------------------------------------------------- the web version (web/server.py) on Modal
# One small CPU machine serves the page (no torch: starts in seconds, ~free while idle); measuring calls
# propose_remote above. Projects live on a Modal volume in the desktop layout, so any project downloads as a zip
# and opens in the desktop app. max_containers=1: one server owns the files and the job queue.
#   modal deploy cloud.py                                   -> prints the site's URL
#   python cloud.py invite someone@missouri.edu [--admin]   -> prints that person's login link (keep it private)

WEB_PKGS = [p for p in PKGS if p.split("==")[0] not in ("torch", "cellpose")] + \
    ([l.strip() for l in (HERE / "web/requirements-web.txt").read_text().splitlines()
      if l.strip() and not l.startswith("-r")] if modal.is_local() else [])
web_data = modal.Volume.from_name("csa-web-data", create_if_missing=True)
web_image = (modal.Image.debian_slim(python_version="3.13")
             .apt_install("libgl1", "libglib2.0-0")
             .pip_install(*WEB_PKGS)
             .env({"CSA_DATA": "/data"})
             .add_local_python_source("pipeline", "run_batch", "summarise", "pattern_plots")
             .add_local_dir(HERE / "web", "/root/web")
             .add_local_file(HERE / "icon.png", "/root/icon.png"))


def _propose_on_gpu(path, diameter=None):
    """The web server's outlining step: the image goes to propose_remote; gray is made here for the QC image."""
    import io
    import numpy as np
    import pipeline as P
    print(f"[gpu] sending {Path(path).name}", flush=True)
    packed, secs = propose_remote.remote(Path(path).name, Path(path).read_bytes(), diameter)
    print(f"[gpu] {Path(path).name} done in {secs} s", flush=True)
    z = np.load(io.BytesIO(packed))
    return P.load_gray(path), z["cells"], z["region"], z["mx"]


@app.function(image=web_image, volumes={"/data": web_data}, cpu=2, memory=4096, max_containers=1,
              scaledown_window=900, timeout=24 * 3600,
              secrets=[modal.Secret.from_name("csa-mail")])   # Gmail for "X asks for access" emails (server._notify)
@modal.concurrent(max_inputs=50)
@modal.asgi_app(label="rooh")   # the site address: https://jagdeep2021--rooh.modal.run
def web():
    import sys
    sys.path.insert(0, "/root/web")
    import server
    server.start_worker(propose=_propose_on_gpu, commit=web_data.commit, reload=web_data.reload)
    return server.app


@app.function(image=web_image, volumes={"/data": web_data})
def invite(email: str, admin: bool = False) -> str:
    import sys
    sys.path.insert(0, "/root/web")
    import server
    token = server.invite(email, admin)
    web_data.commit()
    return token


if __name__ == "__main__":                       # python cloud.py invite <email> [--admin]: on the deployed app
    import sys
    if len(sys.argv) < 3 or sys.argv[1] != "invite":
        sys.exit(__doc__)
    token = modal.Function.from_name(APP, "invite").remote(sys.argv[2], "--admin" in sys.argv)
    print(f"{modal.Function.from_name(APP, 'web').get_web_url()}/login?token={token}")
