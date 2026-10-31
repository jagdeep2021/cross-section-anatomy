"""Web version: the same pipeline behind a browser page. Users upload images, the server outlines them on its GPU,
users look (and, next, correct), set the scale and download the result tables.

  python web/server.py invite <email> [--admin]     -> prints a login link (one per person; keep it private)
  python web/server.py serve [--port 8000]          -> http://localhost:8000

Each project is a folder in the desktop layout (raw/, anatomy/, results/), so run_batch.py and summarise.py run
unchanged and a project can be copied out and opened in the desktop app. Data lives in $CSA_DATA (default ./web_data):
  users.json                          token hash -> {email, admin, [name, approved, requested]}
                                      (asked on the page: approved False until an admin approves it)
  projects/<user>/<project>/raw/      uploaded images (sub-folders kept: they become group_1, group_2 ...)
                            anatomy/  masks/*.npz, qc/, samples.xlsx, settings.json (written by run_batch)
                            results/  written by summarise
                            .cache/   display PNGs (safe to delete)
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import queue
import re
import secrets
import shutil
import sys
import threading
import time
import zipfile
from pathlib import Path, PurePosixPath

import numpy as np
from fastapi import Cookie, Depends, FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, RedirectResponse, Response, StreamingResponse
from fastapi.staticfiles import StaticFiles
from PIL import Image
from skimage.segmentation import find_boundaries
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import pipeline as P      # noqa: E402
import run_batch as RB    # noqa: E402
import summarise as SUM   # noqa: E402

DATA = Path(os.environ.get("CSA_DATA", Path(__file__).resolve().parent.parent / "web_data")).resolve()
STATIC = Path(__file__).resolve().parent / "static"
app = FastAPI(title="Root cross-section anatomy")
COMMIT = [lambda: None]                 # persists the data after a change (set on Modal: the volume's commit)
RELOAD = [lambda: None]                 # sees changes made elsewhere, e.g. an invite (set on Modal: volume reload)


@app.middleware("http")
async def _persist(request, call_next):
    response = await call_next(request)
    if request.method != "GET":         # uploads, scale, results ... saved before the next request can need them
        import asyncio
        await asyncio.to_thread(COMMIT[0])   # off the event loop: a blocking call there can stall the GPU calls
    return response


# ---------------------------------------------------------------- accounts (invite-only)

def _users() -> dict:
    f = DATA / "users.json"
    return json.loads(f.read_text()) if f.exists() else {}


def _hash(token: str) -> str:
    return hashlib.sha256(token.encode()).hexdigest()


def _save_users(users: dict):
    DATA.mkdir(parents=True, exist_ok=True)
    (DATA / "users.json").write_text(json.dumps(users, indent=2))


def invite(email: str, admin: bool = False) -> str:
    """New login token for a person. Only its hash is stored, so a lost link is replaced, never recovered."""
    users = _users()
    users = {h: u for h, u in users.items() if u["email"] != email}     # re-inviting replaces the old link
    token = secrets.token_urlsafe(24)
    users[_hash(token)] = {"email": email, "admin": admin}
    _save_users(users)
    return token


def _signed_in(csa: str | None = Cookie(default=None)) -> dict:
    """The person behind this browser's cookie, approved or still waiting."""
    u = _users().get(_hash(csa)) if csa else None
    if u is None and csa:                # a new invite was written by another machine: look once more
        try:
            RELOAD[0]()
        except Exception:
            pass
        u = _users().get(_hash(csa))
    if u is None:
        raise HTTPException(401, "Not signed in: open your login link again.")
    return u


def user(u: dict = Depends(_signed_in)) -> dict:
    if not u.get("approved", True):      # invited people have no flag: approved
        raise HTTPException(403, "Your request is waiting for Jagdeep to approve it.")
    return u


def admin(u: dict = Depends(user)) -> dict:
    if not u.get("admin"):
        raise HTTPException(403, "Only an admin can do this.")
    return u


def _uid(u: dict) -> str:
    return re.sub(r"[^a-z0-9._-]", "_", u["email"].lower())


# ---------------------------------------------------------------- projects

def _slug(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9._ -]", "", name).strip().replace(" ", "_")
    if not s or s.startswith("."):
        raise HTTPException(400, "Use letters, numbers, spaces, - or _ in the project name.")
    return s[:80]


def project_dir(u: dict, name: str) -> Path:
    return DATA / "projects" / _uid(u) / _slug(name)


def _project(u: dict, name: str) -> Path:
    d = project_dir(u, name)
    if not d.is_dir():
        raise HTTPException(404, f"No project {name}")
    return d


def _safe_rel(rel: str) -> PurePosixPath:
    """An uploaded file's path inside raw/ (folders kept), refusing anything that could escape it."""
    p = PurePosixPath(rel.replace("\\", "/"))
    parts = [x for x in p.parts if x not in ("", ".")]
    if not parts or any(x == ".." or x.startswith(".") for x in parts) or p.is_absolute():
        raise HTTPException(400, f"Bad file name: {rel}")
    return PurePosixPath(*parts)


def _images(d: Path) -> list[dict]:
    """Every uploaded image, measured or not, with its status."""
    out = []
    for rel in RB.list_images(d / "raw"):
        key = RB.image_key(rel)
        npz = d / "anatomy" / "masks" / f"{key}.npz"
        out.append({"image": rel.as_posix(), "key": key,
                    "status": P.load_state(npz)[3] if npz.exists() else "not measured"})
    return out


# ---------------------------------------------------------------- the GPU worker: one job at a time, in order

JOBS: dict[str, dict] = {}            # project path -> progress the page polls
WORK: queue.Queue = queue.Queue()


def _worker(propose=None, commit=lambda: None):
    """propose: the outlining step (None = this machine's GPU; on Modal the GPU function). commit: persist the data
    after each image (Modal volume), so a restart loses at most the image in progress; batch resumes from there."""
    while True:
        d = WORK.get()
        job = JOBS[str(d)]
        job.update(state="running", started=time.time())
        print(f"[measure] start {d}", flush=True)
        try:
            for i, n, line in RB.batch(d / "raw", d / "anatomy", propose=propose):
                job.update(i=i, n=n, line=line)
                print(f"[measure] {i}/{n} {line}", flush=True)
                commit()
            job.update(state="done", line="All images measured.")
        except Exception as e:                       # shown on the page; the server keeps running
            import traceback
            traceback.print_exc()
            job.update(state="error", line=f"{type(e).__name__}: {e}")
        finally:
            WORK.task_done()


def start_worker(propose=None, commit=lambda: None, reload=lambda: None):
    COMMIT[0], RELOAD[0] = commit, reload
    threading.Thread(target=_worker, args=(propose, commit), daemon=True).start()


def queue_measure(d: Path) -> dict:
    job = JOBS.get(str(d))
    if job and job["state"] in ("queued", "running"):
        return job
    JOBS[str(d)] = {"state": "queued", "i": 0, "n": len(RB.list_images(d / "raw")), "line": "Waiting for the GPU...",
                    "ahead": WORK.qsize()}
    WORK.put(d)
    return JOBS[str(d)]


# ---------------------------------------------------------------- display layers

def _cache(d: Path, key: str, kind: str) -> Path:
    c = d / ".cache"
    c.mkdir(exist_ok=True)
    return c / f"{key}.{kind}.png"


def gray_png(d: Path, rel: str) -> Path:
    f = _cache(d, RB.image_key(Path(rel)), "gray")
    src = d / "raw" / rel
    if not f.exists() or f.stat().st_mtime < src.stat().st_mtime:
        g = P.load_gray(src)                          # what the model saw: stretched, CLAHE
        Image.fromarray((g * 255).astype(np.uint8)).save(f, optimize=True)
    return f


OVERLAY = {"cell": (80, 255, 120, 200), "root": (255, 220, 0, 255), "stele": (0, 220, 255, 255),
           "mx": (255, 60, 60, 110), "path": (255, 150, 0, 120)}


def overlay_png(d: Path, rel: str) -> Path:
    """Outlines as a transparent layer over the image: cells, root edge, stele edge, metaxylem and the spoke cells."""
    key = RB.image_key(Path(rel))
    npz = d / "anatomy" / "masks" / f"{key}.npz"
    f = _cache(d, key, "overlay")
    if f.exists() and f.stat().st_mtime >= npz.stat().st_mtime:
        return f
    cells, region, mx, _ = P.load_state(npz)
    _, path_q, path_cells = P.load_paths(npz) or ((), np.array([]), np.array([]))
    rgba = np.zeros((*cells.shape, 4), np.uint8)
    rgba[np.isin(cells, path_cells) & (cells > 0)] = OVERLAY["path"]
    rgba[mx > 0] = OVERLAY["mx"]
    rgba[find_boundaries(cells, mode="inner") & (cells > 0)] = OVERLAY["cell"]
    w = max(2, round(max(cells.shape) / 600))      # root and stele edges thick enough to see on a whole-image view
    edge = lambda m: ndi.binary_dilation(find_boundaries(m, mode="inner") & m, iterations=w) & m
    rgba[edge(region >= 1)] = OVERLAY["root"]
    rgba[edge(region == 2)] = OVERLAY["stele"]
    Image.fromarray(rgba, "RGBA").save(f, optimize=True)
    return f


# ---------------------------------------------------------------- routes

@app.get("/login")
def login(token: str):
    if _hash(token) not in _users():
        try:
            RELOAD[0]()                  # the invite may have been written by another machine just now
        except Exception:
            pass
    if _hash(token) not in _users():
        raise HTTPException(401, "This login link is not valid (it may have been replaced). Ask Jagdeep for a new one.")
    r = RedirectResponse("/")
    r.set_cookie("csa", token, httponly=True, samesite="lax", max_age=180 * 24 * 3600)
    return r


@app.get("/")
def index():
    return FileResponse(STATIC / "index.html")


@app.get("/icon.png")
def icon():
    return FileResponse(STATIC.parent.parent / "icon.png")


@app.get("/api/me")
def me(u: dict = Depends(_signed_in)):
    return {"email": u["email"], "name": u.get("name", ""), "admin": u.get("admin", False),
            "approved": u.get("approved", True)}


@app.post("/api/join")
def join(name: str = Form(""), email: str = Form("")):
    """Ask for access from the page itself: this browser gets a cookie that works once an admin approves it.
    Anyone can ask, so an admin approves only people they know (the email is typed by the asker, not checked)."""
    name, email = name.strip()[:80], email.strip().lower()[:120]
    if not name or not re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", email):
        raise HTTPException(400, "Enter your name and email.")
    users = _users()
    if sum(not u.get("approved", True) for u in users.values()) >= 50:
        raise HTTPException(429, "Too many requests waiting. Ask Jagdeep directly.")
    token = secrets.token_urlsafe(24)
    users[_hash(token)] = {"email": email, "name": name, "admin": False, "approved": False,
                           "requested": time.strftime("%Y-%m-%d %H:%M")}
    _save_users(users)
    r = Response(status_code=204)
    r.set_cookie("csa", token, httponly=True, samesite="lax", max_age=180 * 24 * 3600)
    return r


@app.get("/api/people")
def people(_: dict = Depends(admin)):
    return [{"id": h[:16], "email": u["email"], "name": u.get("name", ""), "admin": u.get("admin", False),
             "approved": u.get("approved", True), "requested": u.get("requested", "")}
            for h, u in _users().items()]


@app.post("/api/people/{pid}/{action}")
def person(pid: str, action: str, me: dict = Depends(admin)):
    users = _users()
    h = next((h for h in users if h.startswith(pid)), None)
    if h is None or len(pid) < 16 or action not in ("approve", "remove"):
        raise HTTPException(404, "No such person.")
    if action == "approve":
        users[h]["approved"] = True
    elif users[h]["email"] == me["email"] and users[h].get("admin"):
        raise HTTPException(400, "You can't remove your own admin login.")
    else:
        del users[h]                     # their projects stay on the server; only the login goes
    _save_users(users)
    return {"ok": True}


@app.get("/api/projects")
def projects(u: dict = Depends(user)):
    root = DATA / "projects" / _uid(u)
    return sorted(p.name for p in root.iterdir() if p.is_dir()) if root.exists() else []


@app.post("/api/projects")
def new_project(name: str = Form(...), u: dict = Depends(user)):
    d = project_dir(u, name)
    (d / "raw").mkdir(parents=True, exist_ok=True)
    return {"name": d.name}


@app.post("/api/projects/{name}/upload")
async def upload(name: str, files: list[UploadFile] = File(...), u: dict = Depends(user)):
    """Files keep their folder path (a dropped folder Maize/sibs/P1/x.tif stays that way): the folders become groups."""
    d = _project(u, name)
    saved = []
    for f in files:
        rel = _safe_rel(f.filename)
        if rel.suffix.lower() not in P.IMAGE_EXT:
            continue
        dest = d / "raw" / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        with open(dest, "wb") as out:
            shutil.copyfileobj(f.file, out)
        saved.append(rel.as_posix())
    return {"saved": len(saved), "skipped": len(files) - len(saved)}


@app.get("/api/projects/{name}")
def project(name: str, u: dict = Depends(user)):
    d = _project(u, name)
    scale = {}
    if (d / "anatomy" / "samples.xlsx").exists():
        import pandas as pd
        s = pd.read_excel(d / "anatomy" / "samples.xlsx", dtype={"image_file": str})
        scale = {r.image_file: (None if pd.isna(r.um_per_px) else float(r.um_per_px)) for r in s.itertuples()}
    imgs = _images(d)
    for im in imgs:
        im["um_per_px"] = scale.get(im["image"])
    results = sorted(p.name for p in (d / "results").glob("*.xlsx")) if (d / "results").exists() else []
    return {"name": d.name, "images": imgs, "job": JOBS.get(str(d)), "results": results}


@app.post("/api/projects/{name}/measure")
def measure(name: str, u: dict = Depends(user)):
    return queue_measure(_project(u, name))


@app.get("/api/projects/{name}/image/{kind}")
def image(name: str, kind: str, rel: str, u: dict = Depends(user)):
    d = _project(u, name)
    rel = _safe_rel(rel).as_posix()
    if kind == "gray":
        return FileResponse(gray_png(d, rel), headers={"Cache-Control": "no-cache"})
    if kind == "overlay":
        if not (d / "anatomy" / "masks" / f"{RB.image_key(Path(rel))}.npz").exists():
            raise HTTPException(404, "Not measured yet")
        return FileResponse(overlay_png(d, rel), headers={"Cache-Control": "no-cache"})
    raise HTTPException(404)


@app.get("/api/projects/{name}/traits")
def traits(name: str, rel: str, u: dict = Depends(user)):
    """The numbers for one image, from derive(): the same function the results use."""
    d = _project(u, name)
    npz = d / "anatomy" / "masks" / f"{RB.image_key(Path(_safe_rel(rel)))}.npz"
    cells, region, mx, status = P.load_state(npz)
    settings = json.loads((d / "anatomy" / "settings.json").read_text())
    props, _, tr = P.derive(cells, region, mx, float(settings.get("quadrant_angle_deg", 0)),
                            P.centre_of(region, P.load_centre(npz)))
    _, path_q, _ = P.load_paths(npz) or ((), np.array([]), None)
    tr["files_per_quarter"] = [int((np.asarray(path_q) == q).sum()) for q in (1, 2, 3, 4)]
    tr["status"] = status
    return {k: (None if isinstance(v, float) and np.isnan(v) else v) for k, v in tr.items()}


@app.post("/api/projects/{name}/edit")
def edit(name: str, rel: str = Form(...), target: str = Form(...), points: str = Form(...),
         remove: bool = Form(False), u: dict = Depends(user)):
    """A drawn outline (points [[y, x], ...] in image pixels) applied with the desktop editor's own rules
    (pipeline.edit_region); the image becomes 'reviewed'. Returns the image's new numbers."""
    d = _project(u, name)
    rel = _safe_rel(rel).as_posix()
    npz = d / "anatomy" / "masks" / f"{RB.image_key(Path(rel))}.npz"
    if not npz.exists():
        raise HTTPException(404, "Not measured yet")
    if target not in ("Stele", "Root"):
        raise HTTPException(400, f"Editing {target} in the browser is not available yet: use the desktop app.")
    pts = np.asarray(json.loads(points), float)
    if pts.ndim != 2 or pts.shape[1] != 2 or len(pts) < 3:
        raise HTTPException(400, "Draw an outline: hold Shift and drag around the area.")
    _, region, _, _ = P.load_state(npz)
    region = region.copy()
    m = P.outline_mask(pts, region.shape)
    if m.sum() < 10:
        raise HTTPException(400, "That outline is too small.")
    P.edit_region(region, m, target, remove)
    P.save_edit(npz, region=region)
    return traits(name, rel, u)


@app.post("/api/projects/{name}/scale")
def set_scale(name: str, um_per_px: float = Form(...), group: str = Form(""), rel: str = Form(""),
              u: dict = Depends(user)):
    """Scale for one image (rel) or every image whose first folder is group ('' = all images)."""
    import pandas as pd
    d = _project(u, name)
    if not any((d / "anatomy" / "masks").glob("*.npz")):
        raise HTTPException(400, "Measure the images first (Outline new images), then set the scale.")
    f = SUM.write_samples(d / "anatomy")
    s = pd.read_excel(f, dtype={"image_file": str, "plant": str})
    if "image_file" not in s:               # samples.xlsx written before any image was measured: rebuild it
        f.unlink()
        f = SUM.write_samples(d / "anatomy")
        s = pd.read_excel(f, dtype={"image_file": str, "plant": str})
    rows = (s.image_file == rel) if rel else (s.image_file.str.split("/").str[0] == group) if group else s.image_file.notna()
    s.loc[rows, "um_per_px"] = um_per_px if um_per_px > 0 else np.nan
    s.to_excel(f, index=False)
    return {"set": int(rows.sum())}


@app.post("/api/projects/{name}/results")
def make_results(name: str, u: dict = Depends(user)):
    d = _project(u, name)
    if not any((d / "anatomy" / "masks").glob("*.npz")):
        raise HTTPException(400, "Measure the images first.")
    SUM.main(d / "anatomy")
    return {"results": sorted(p.name for p in (d / "results").glob("*.xlsx"))}


@app.get("/api/projects/{name}/download")
def download(name: str, u: dict = Depends(user)):
    """results/ plus anatomy/samples.xlsx as one zip."""
    d = _project(u, name)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for p in sorted((d / "results").rglob("*")) if (d / "results").exists() else []:
            if p.is_file():
                z.write(p, f"{d.name}/results/{p.relative_to(d / 'results').as_posix()}")
        if (d / "anatomy" / "samples.xlsx").exists():
            z.write(d / "anatomy" / "samples.xlsx", f"{d.name}/anatomy/samples.xlsx")
    buf.seek(0)
    return StreamingResponse(buf, media_type="application/zip",
                             headers={"Content-Disposition": f'attachment; filename="{d.name}_results.zip"'})


app.mount("/static", StaticFiles(directory=STATIC), name="static")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    inv = sub.add_parser("invite")
    inv.add_argument("email")
    inv.add_argument("--admin", action="store_true")
    inv.add_argument("--base", default="http://localhost:8000")
    srv = sub.add_parser("serve")
    srv.add_argument("--port", type=int, default=8000)
    srv.add_argument("--host", default="127.0.0.1")
    a = ap.parse_args()
    if a.cmd == "invite":
        print(f"{a.base}/login?token={invite(a.email, a.admin)}")
        return
    import uvicorn
    start_worker()
    uvicorn.run(app, host=a.host, port=a.port)


if __name__ == "__main__":
    main()
