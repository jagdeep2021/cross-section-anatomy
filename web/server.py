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
from fastapi import Cookie, Depends, FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, RedirectResponse, Response, StreamingResponse
from fastapi.staticfiles import StaticFiles
from PIL import Image
from skimage.segmentation import find_boundaries
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import pipeline as P      # noqa: E402
import run_batch as RB    # noqa: E402
import summarise as SUM   # noqa: E402
import pattern_plots as PP   # noqa: E402

DATA = Path(os.environ.get("CSA_DATA", Path(__file__).resolve().parent.parent / "web_data")).resolve()
STATIC = Path(__file__).resolve().parent / "static"
app = FastAPI(title="Rooh")
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


# ---------------------------------------------------------------- local mode (the desktop app, `server.py local`)
# The same page and rules, on this computer: no accounts; a project is any experiment folder the user opens (it
# stays where it is; the list of opened folders is kept in ~/.rooh/local_projects.json). The browser gets a private
# key in its cookie when the launcher opens it, so other web pages cannot drive this server.
LOCAL = {"on": False, "key": ""}
LOCAL_USER = {"email": "this computer", "name": "", "admin": False, "approved": True}
LOCAL_REG = Path.home() / ".rooh" / "local_projects.json"


def _reg() -> dict:
    try:
        return json.loads(LOCAL_REG.read_text())
    except (OSError, ValueError):
        return {}


def _reg_save(r: dict) -> None:
    LOCAL_REG.parent.mkdir(parents=True, exist_ok=True)
    LOCAL_REG.write_text(json.dumps(r, indent=1))


def _signed_in(csa: str | None = Cookie(default=None)) -> dict:
    """The person behind this browser's cookie, approved or still waiting."""
    if LOCAL["on"]:
        if csa != LOCAL["key"]:
            raise HTTPException(401, "Open Rooh with its Open app launcher (this page has no key).")
        return LOCAL_USER
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
    if LOCAL["on"]:
        p = _reg().get(name)
        if p is None:
            raise HTTPException(404, f"No project {name}: open its folder again.")
        return Path(p)
    return DATA / "projects" / _uid(u) / _slug(name)


def _raw(d: Path) -> Path:
    """The project's images: raw/ on the web; on this computer wherever anatomy/settings.json points (raw/, or the
    image folder itself), as the desktop app has always read it."""
    s = d / "anatomy" / "settings.json"
    return P.image_dir(d / "anatomy", json.loads(s.read_text())) if s.exists() else d / "raw"


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
    for rel in RB.list_images(_raw(d)):
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
            for i, n, line in RB.batch(_raw(d), d / "anatomy", propose=propose):
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
    JOBS[str(d)] = {"state": "queued", "i": 0, "n": len(RB.list_images(_raw(d))), "line": "Waiting for the GPU...",
                    "ahead": WORK.qsize()}
    WORK.put(d)
    return JOBS[str(d)]


# ---------------------------------------------------------------- display layers

# Display files (pictures, plots, the trait cache) depend on the code as well as the masks: a fingerprint of the
# measurement code names their folder, so after a deploy that changes a rule every project redraws once, on first view.
CODE = hashlib.sha1(b"".join((Path(__file__).resolve().parent.parent / f).read_bytes() for f in (
    "pipeline.py", "pattern_plots.py", "summarise.py", "run_batch.py", "web/server.py"))).hexdigest()[:12]


def _cache_dir(d: Path) -> Path:
    root, c = d / ".cache", d / ".cache" / CODE
    if not c.exists():
        c.mkdir(parents=True, exist_ok=True)
        for old in root.iterdir():               # drawn by older code: never shown again
            if old != c:
                shutil.rmtree(old, ignore_errors=True) if old.is_dir() else old.unlink(missing_ok=True)
    return c


def _cache(d: Path, key: str, kind: str) -> Path:
    return _cache_dir(d) / f"{key}.{kind}.png"


def gray_png(d: Path, rel: str) -> Path:
    f = _cache(d, RB.image_key(Path(rel)), "gray")
    src = _raw(d) / rel
    if not f.exists() or f.stat().st_mtime < src.stat().st_mtime:
        g = P.load_gray(src)                          # what the model saw: stretched, CLAHE
        Image.fromarray((g * 255).astype(np.uint8)).save(f, optimize=True)
    return f


OVERLAY = {"cell": (80, 255, 120, 200), "root": (255, 220, 0, 255), "stele": (0, 220, 255, 255),
           "mx": (255, 60, 60, 110), "lacuna": (255, 140, 0, 150), "lacuna_edge": (255, 140, 0, 255), "lacuna_off": (150, 150, 150, 140), "stele_fill": (0, 220, 255, 85)}


def _live(d: Path, npz: Path):
    """What every live view starts from: masks, centre, excluded arcs, and the kept cells (props after the arcs)."""
    cells, region, mx, status = P.load_state(npz)
    settings = json.loads((d / "anatomy" / "settings.json").read_text())
    centre = P.centre_of(region, P.load_centre(npz))
    arcs = P.load_sectors(npz, float(settings.get("quadrant_angle_deg", 0)))
    props, _, tr = P.derive(cells, region, mx, 0, centre)
    return cells, region, mx, status, centre, arcs, P.exclude_sectors(props, arcs), tr


def overlay_png(d: Path, rel: str, part: str) -> Path:
    """One transparent layer over the image, so the page fades them separately (as the desktop's opacities):
    part "cells" = the cell outlines; "masks" = root edge, stele (shaded + edge), metaxylem. Excluded sectors are
    drawn by the page (from traits), so changing one needs no new picture."""
    key = RB.image_key(Path(rel))
    npz = d / "anatomy" / "masks" / f"{key}.npz"
    f = _cache(d, key, part)
    if f.exists() and f.stat().st_mtime >= npz.stat().st_mtime:
        return f
    cells, region, mx, _ = P.load_state(npz)
    rgba = np.zeros((*cells.shape, 4), np.uint8)
    if part == "cells":
        rgba[find_boundaries(cells, mode="inner") & (cells > 0)] = OVERLAY["cell"]
    elif part == "lacunae":                            # aerenchyma lacunae as Make results finds them (kept angles)
        if (region == 2).any():
            *_, centre, arcs, props, _ = _live(d, npz)
            lab = P.lacunae(cells, region, mx, props, centre, arcs, mask=True, not_lacuna=P.load_not_lacuna(npz))[3]
            rgba[lab < 0] = OVERLAY["lacuna_off"]        # a reviewer said: not a lacuna
            lac = lab > 0
            rgba[lac] = OVERLAY["lacuna"]
            w = max(2, round(max(cells.shape) / 600))
            rgba[ndi.binary_dilation(find_boundaries(lac, mode="inner") & lac, iterations=w) & lac] = OVERLAY["lacuna_edge"]
    elif part == "files":                              # cells coloured by cell file number (pattern_plots.file_rgba)
        if (region == 2).any():
            *_, centre, arcs, props, _ = _live(d, npz)
            rgba = PP.file_rgba(cells, P.cell_files(cells, props, centre, arcs, region=region).file)
    else:
        rgba[region == 2] = OVERLAY["stele_fill"]      # the stele is an area: lightly shaded, its cells still visible
        rgba[mx > 0] = OVERLAY["mx"]
        w = max(2, round(max(cells.shape) / 600))      # root and stele edges thick enough to see on a whole-image view
        edge = lambda m: ndi.binary_dilation(find_boundaries(m, mode="inner") & m, iterations=w) & m
        rgba[edge(region >= 1)] = OVERLAY["root"]
        rgba[edge(region == 2)] = OVERLAY["stele"]
    Image.fromarray(rgba, "RGBA").save(f, optimize=True)
    return f


# ---------------------------------------------------------------- routes

@app.get("/local")
def local_login(key: str):
    """The launcher opens this once: the key becomes this browser's cookie (SameSite strict: other sites can't use it)."""
    if not LOCAL["on"] or not secrets.compare_digest(key, LOCAL["key"]):
        raise HTTPException(403, "Wrong key: close this tab and use the Open app launcher.")
    r = RedirectResponse("/")
    r.set_cookie("csa", key, httponly=True, samesite="strict")
    return r


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
    return FileResponse(STATIC / "index.html", headers={"Cache-Control": "no-cache"})   # a reload gets a new deploy


@app.get("/icon.png")
def icon():
    return FileResponse(STATIC.parent.parent / "icon.png")


@app.get("/api/me")
def me(u: dict = Depends(_signed_in)):
    return {"email": u["email"], "name": u.get("name", ""), "admin": u.get("admin", False),
            "approved": u.get("approved", True), "local": LOCAL["on"]}


def _notify(name: str, email: str, site: str):
    """Email the admin about a new access request (Gmail + app password from the Modal secret csa-mail:
    GMAIL_USER, GMAIL_APP_PASSWORD, optional MAIL_TO). No secret = no email; the request still waits in People."""
    user, pw = os.environ.get("GMAIL_USER"), os.environ.get("GMAIL_APP_PASSWORD")
    if not (user and pw):
        return
    import smtplib
    from email.message import EmailMessage
    m = EmailMessage()
    m["From"], m["To"] = user, os.environ.get("MAIL_TO", user)
    m["Subject"] = f"Rooh: {' '.join(name.split())} asks for access"
    m.set_content(f"{name} ({email}) asked to use Rooh.\n\n"
                  f"Approve or remove them here (opens the People list; you must be signed in as admin in that browser): {site}#people\n\n"
                  "Approve only people you know: the name and email are typed by the asker, not checked.")
    try:
        with smtplib.SMTP_SSL("smtp.gmail.com", 465, timeout=20) as smtp:
            smtp.login(user, pw)
            smtp.send_message(m)
    except Exception as e:                       # a failed email never fails the request
        print(f"[mail] not sent: {e}", flush=True)


@app.post("/api/join")
def join(request: Request, name: str = Form(""), email: str = Form(""), csa: str | None = Cookie(default=None)):
    """Ask for access from the page itself: this browser gets a cookie that works once an admin approves it.
    Anyone can ask, so an admin approves only people they know (the email is typed by the asker, not checked)."""
    name, email = name.strip()[:80], email.strip().lower()[:120]
    if not name or not re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", email):
        raise HTTPException(400, "Enter your name and email.")
    users = _users()
    if csa and _hash(csa) in users and users[_hash(csa)].get("approved", True):   # never swap an approved
        raise HTTPException(400, "This browser is already signed in: reload the page.")  # login for a waiting one
    if sum(not u.get("approved", True) for u in users.values()) >= 50:
        raise HTTPException(429, "Too many requests waiting. Ask Jagdeep directly.")
    token = secrets.token_urlsafe(24)
    users[_hash(token)] = {"email": email, "name": name, "admin": False, "approved": False,
                           "requested": time.strftime("%Y-%m-%d %H:%M")}
    _save_users(users)
    threading.Thread(target=_notify, args=(name, email, str(request.base_url)), daemon=True).start()
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
    if LOCAL["on"]:
        return sorted(n for n, p in _reg().items() if Path(p).is_dir())
    root = DATA / "projects" / _uid(u)
    return sorted(p.name for p in root.iterdir() if p.is_dir()) if root.exists() else []


@app.post("/api/projects")
def new_project(name: str = Form(...), u: dict = Depends(user)):
    if LOCAL["on"]:
        raise HTTPException(400, "On this computer, use Open folder…")
    d = project_dir(u, name)
    (d / "raw").mkdir(parents=True, exist_ok=True)
    return {"name": d.name}


@app.post("/api/projects/{name}/upload")
def upload(name: str, files: list[UploadFile] = File(...), u: dict = Depends(user)):
    """Files keep their folder path (a dropped folder Maize/sibs/P1/x.tif stays that way): the folders become groups.
    A plain def, not async: the copies block, so FastAPI runs them on a worker thread and the site keeps answering
    while a big batch is written to the volume (as async it froze every other request)."""
    d = _project(u, name)
    saved = []
    for f in files:
        rel = _safe_rel(f.filename)
        if rel.suffix.lower() not in P.IMAGE_EXT:
            continue
        dest = _raw(d) / rel
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
    figs = d / "results" / "figures"
    figures = sorted(p.relative_to(figs).as_posix() for p in figs.rglob("*.png")) if figs.exists() else []
    return {"name": d.name, "images": imgs, "job": JOBS.get(str(d)), "results": results, "figures": figures}


@app.get("/api/projects/{name}/figure")
def figure(name: str, path: str, u: dict = Depends(user)):
    """One figure from Make results (results/figures/...): group figures, and cortex_pattern/<image>.png per image."""
    f = _project(u, name) / "results" / "figures" / _safe_rel(path)
    if f.suffix != ".png" or not f.is_file():
        raise HTTPException(404, "No such figure: click Make results first.")
    return FileResponse(f, headers={"Cache-Control": "no-store"})


@app.post("/api/projects/{name}/measure")
def measure(name: str, u: dict = Depends(user)):
    return queue_measure(_project(u, name))


@app.get("/api/projects/{name}/image/{kind}")
def image(name: str, kind: str, rel: str, u: dict = Depends(user)):
    d = _project(u, name)
    rel = _safe_rel(rel).as_posix()
    if kind == "gray":
        return FileResponse(gray_png(d, rel), headers={"Cache-Control": "no-cache"})
    if kind in ("cells", "masks", "lacunae", "files"):
        if not (d / "anatomy" / "masks" / f"{RB.image_key(Path(rel))}.npz").exists():
            raise HTTPException(404, "Not measured yet")
        return FileResponse(overlay_png(d, rel, kind), headers={"Cache-Control": "no-cache"})
    raise HTTPException(404)


@app.get("/api/projects/{name}/traits")
def traits(name: str, rel: str, u: dict = Depends(user)):
    """The numbers for one image, from derive(): the same function the results use."""
    d = _project(u, name)
    npz = d / "anatomy" / "masks" / f"{RB.image_key(Path(_safe_rel(rel)))}.npz"
    if not npz.exists():                         # not measured yet, or the project was just deleted
        raise HTTPException(404, "Not measured yet")
    cells, region, mx, status, centre, arcs, props, tr = _live(d, npz)
    tr.update(P.sector_traits(cells, region, props, centre, arcs))   # cortex: kept angles; stele, vessels: whole
    lac, _, n_gaps = P.lacunae(cells, region, mx, props, centre, arcs, not_lacuna=P.load_not_lacuna(npz))
    tr.update(P.lacuna_traits(lac, n_gaps, tr["cortex_area_kept_px"]))
    if not tr.get("stele_found", (region == 2).any()):   # as in the results: no stele -> blank, not 0
        tr.update({k: None for k in P.NEEDS_STELE if k in tr})
    tr.pop("cortex_area_px", None)
    tr["status"] = status
    # for drawing the excluded arcs on the page (as review.py does)
    r_st, r_rt = P.radial_profile(region, centre)
    tr["sectors"] = {"centre": [float(c) for c in centre], "arcs": [[float(a), float(b), r] for a, b, r in arcs],
                     "rst": [float(v) for v in r_st], "rrt": [float(v) for v in r_rt]}   # stele / root radius per degree
    clean = lambda v: None if isinstance(v, float) and np.isnan(v) else v.item() if isinstance(v, np.generic) else v
    return {k: clean(v) for k, v in tr.items()}


@app.get("/api/projects/{name}/pattern")
def pattern(name: str, rel: str, u: dict = Depends(user)):
    """This image's cortex pattern, live from its current masks (pattern_plots.live_pattern, as review.py): the section
    and each quarter kept. Redrawn when the mask, its sectors, the scale or the code change."""
    import matplotlib.pyplot as plt
    import pattern_plots as PP
    d = _project(u, name)
    rel = _safe_rel(rel).as_posix()
    key = RB.image_key(Path(rel))
    npz = d / "anatomy" / "masks" / f"{key}.npz"
    if not npz.exists():
        raise HTTPException(404, "Not measured yet")
    f = _cache(d, key, "pattern")
    samples = d / "anatomy" / "samples.xlsx"
    newest = max(npz.stat().st_mtime, samples.stat().st_mtime if samples.exists() else 0)
    if not f.exists() or f.stat().st_mtime < newest:
        *_, arcs, props, _ = _live(d, npz)
        fig, ax = plt.subplots(figsize=(3.4, 2.8), dpi=150, layout="constrained")
        PP.live_pattern(ax, props, arcs, P.um_per_px(d / "anatomy", rel), SUM.N_BINS)
        fig.savefig(f)
        plt.close(fig)
    return FileResponse(f, headers={"Cache-Control": "no-cache"})


@app.get("/api/projects/{name}/all_traits")
def all_traits(name: str, u: dict = Depends(user)):
    """Every measured image's numbers (for the per-image bars that show outliers). Cached per image until its
    mask or the measurement code changes, so only edited images are recomputed."""
    d = _project(u, name)
    f = _cache_dir(d) / "traits.json"
    cache = json.loads(f.read_text()) if f.exists() else {}
    out = []
    for im in _images(d):
        if im["status"] == "not measured":
            continue
        npz = d / "anatomy" / "masks" / f"{RB.image_key(Path(im['image']))}.npz"
        t = npz.stat().st_mtime
        if cache.get(im["image"], {}).get("mtime") != t:
            cache[im["image"]] = {"mtime": t, "traits": traits(name, im["image"], u)}
        out.append({"image": im["image"], **cache[im["image"]]["traits"]})
    f.write_text(json.dumps(cache))
    return out


@app.post("/api/projects/{name}/edit")
def edit(name: str, rel: str = Form(...), target: str = Form(...), points: str = Form(...),
         remove: bool = Form(False), reason: str = Form("torn"), u: dict = Depends(user)):
    """A click (one point) or a drawn outline (points [[y, x], ...] in image pixels) on Cell / Metaxylem / Stele /
    Root, applied with the desktop editor's own rules (pipeline.click_edit / outline_edit); Sectors: a drag from
    one side of the damage to the other excludes that arc (first and last point, with a reason), Ctrl/Cmd-click
    inside an excluded arc brings it back (pipeline.sector_edit). The image becomes 'reviewed'. Returns the image's
    new numbers, plus a note when the edit needs a follow-up."""
    d = _project(u, name)
    rel = _safe_rel(rel).as_posix()
    npz = d / "anatomy" / "masks" / f"{RB.image_key(Path(rel))}.npz"
    if not npz.exists():
        raise HTTPException(404, "Not measured yet")
    if target not in ("Cell", "Metaxylem", "Stele", "Root", "Sectors", "Aerenchyma"):
        raise HTTPException(400, f"Editing {target} in the browser is not available yet: use the desktop app.")
    pts = np.asarray(json.loads(points), float)
    if pts.ndim != 2 or pts.shape[1] != 2 or len(pts) == 0 or (len(pts) == 2 and target != "Sectors"):
        raise HTTPException(400, "Click, or hold Shift and drag around the area.")
    cells, region, mx, _ = P.load_state(npz)
    if target == "Sectors":
        settings = json.loads((d / "anatomy" / "settings.json").read_text())
        centre = P.centre_of(region, P.load_centre(npz))
        arcs = P.load_sectors(npz, float(settings.get("quadrant_angle_deg", 0)))
        if len(pts) == 1 and not remove:
            return {**traits(name, rel, u), "note": "Sectors: drag around the stele across the damage to exclude it; "
                                                    "Ctrl/Cmd-click an excluded arc to bring it back."}
        new = P.sector_edit(arcs, centre, tuple(pts[0]), None if len(pts) == 1 else tuple(pts[-1]),
                            reason, remove=len(pts) == 1)
        if new == arcs:
            return {**traits(name, rel, u), "note": "Nothing changed (drag further round the stele)." if len(pts) > 1
                    else "No excluded arc there."}
        P.save_edit(npz, **P.sector_arrays(new))
        return {**traits(name, rel, u), "note": ""}
    if target == "Aerenchyma":                     # pipeline.lacuna_edit: draw = a lacuna; Cmd/Ctrl-click = not one
        cells, region, mx, _, centre, arcs, props, _ = _live(d, npz)
        cells = cells.copy()
        old = P.load_not_lacuna(npz)
        marks, note = P.lacuna_edit(cells, region, mx, props, centre, arcs, old, pts, remove)
        if marks != old or len(pts) > 2:
            P.save_edit(npz, cells=cells, **P.not_lacuna_arrays(marks))
        return {**traits(name, rel, u), "note": note}
    cells, region, mx = cells.copy(), region.copy(), mx.copy()
    note = ""
    if len(pts) == 1:
        y, x = (int(round(v)) for v in pts[0])
        if not (0 <= y < cells.shape[0] and 0 <= x < cells.shape[1]):
            raise HTTPException(400, "Click on the image.")
        gray = P.load_gray(_raw(d) / rel) if target != "Root" else None
        note = P.click_edit(gray, cells, region, mx, target, y, x, remove)
    else:
        m = P.outline_mask(pts, region.shape)
        if m.sum() < 10:
            raise HTTPException(400, "That outline is too small.")
        P.outline_edit(cells, region, mx, target, m, remove)
    P.save_edit(npz, cells=cells, region=region, mx=mx)
    return {**traits(name, rel, u), "note": note}


@app.post("/api/projects/{name}/reviewed")
def reviewed(name: str, rel: str = Form(...), u: dict = Depends(user)):
    """Save + next: the image is checked (status reviewed), edited or not, as review.py's Save + next."""
    d = _project(u, name)
    npz = d / "anatomy" / "masks" / f"{RB.image_key(Path(_safe_rel(rel)))}.npz"
    if not npz.exists():
        raise HTTPException(404, "Not measured yet")
    P.save_edit(npz)
    return {"ok": True}


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
def make_results(name: str, image_figures: bool = Form(False), u: dict = Depends(user)):
    d = _project(u, name)
    if not any((d / "anatomy" / "masks").glob("*.npz")):
        raise HTTPException(400, "Measure the images first.")
    SUM.main(d / "anatomy", image_figures)
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


@app.get("/api/projects/{name}/download_all")
def download_all(name: str, u: dict = Depends(user)):
    """The whole project as one zip in the desktop layout (raw/ images, anatomy/ masks + settings + samples,
    results/), so it opens in the desktop app. Written to a temporary file: raw images can be large."""
    import tempfile
    from starlette.background import BackgroundTask
    d = _project(u, name)
    tmp = tempfile.NamedTemporaryFile(suffix=".zip", delete=False)
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        for p in sorted(d.rglob("*")):
            rel = p.relative_to(d)
            if p.is_file() and rel.parts[0] != ".cache" and "qc" not in rel.parts[:2]:
                z.write(p, f"{d.name}/{rel.as_posix()}")
    tmp.close()
    return FileResponse(tmp.name, media_type="application/zip", filename=f"{d.name}_project.zip",
                        background=BackgroundTask(os.unlink, tmp.name))


@app.delete("/api/projects/{name}")
def delete_project(name: str, u: dict = Depends(user)):
    """Moves the project to trash/ (out of the user's list, not erased: an admin can still restore it). On this
    computer: only takes the folder off the list; nothing on disk is touched."""
    if LOCAL["on"]:
        r = _reg()
        r.pop(name, None)
        _reg_save(r)
        return {"ok": True}
    d = _project(u, name)
    if JOBS.get(str(d), {}).get("state") in ("queued", "running"):
        raise HTTPException(400, "This project is being measured: wait until it finishes.")
    dest = DATA / "trash" / _uid(u) / f"{d.name}_{time.strftime('%Y-%m-%d_%H%M%S')}"
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(d), dest)
    return {"ok": True}


def _pick_folder() -> str | None:
    """The computer's own folder picker (Mac: Finder's; Windows / Linux: Tk's). None = cancelled."""
    import subprocess
    if sys.platform == "darwin":
        r = subprocess.run(["osascript", "-e", 'POSIX path of (choose folder with prompt "Open an experiment folder '
                            '(or its raw/, anatomy/, or a folder of images)")'], capture_output=True, text=True)
    else:
        r = subprocess.run([sys.executable, "-c", "import tkinter as t, tkinter.filedialog as f; w = t.Tk(); "
                            "w.withdraw(); w.attributes('-topmost', True); "
                            "print(f.askdirectory(title='Open an experiment folder') or '')"], capture_output=True, text=True)
    p = r.stdout.strip()
    return p or None


@app.post("/api/local/open")
def local_open(path: str = Form(""), u: dict = Depends(user)):
    """This computer only: pick a folder (or pass its path) -> it joins the project list, where it is. Same rule as
    the desktop app (run_batch.locate): the experiment folder, its raw/ or anatomy/, or a folder of images."""
    if not LOCAL["on"]:
        raise HTTPException(404)
    path = path or _pick_folder()
    if not path:
        return {"name": None}
    anatomy, imgs = RB.locate(path)
    if anatomy is None:
        raise HTTPException(400, f"No images (.tif .png .jpg), raw/ or anatomy/ in {path}.")
    if imgs is not None:
        try:
            RB.load_settings(anatomy, imgs)          # anatomy/ beside the images: also tells now if it can be written
        except PermissionError:
            raise HTTPException(400, f"Rooh writes its results in an anatomy/ folder next to the images, and {Path(path).name} "
                                     "is read-only (a shared folder?). Copy the images into a folder of your own "
                                     "(e.g. <experiment>/raw) and open that.")
    d = anatomy.parent
    r = _reg()
    name = next(n for n in (d.name, *(f"{d.name} ({k})" for k in range(2, 99))) if r.get(n) in (None, str(d)))
    r[name] = str(d)
    _reg_save(r)
    return {"name": name}


app.mount("/static", StaticFiles(directory=STATIC), name="static")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    inv = sub.add_parser("invite")
    inv.add_argument("email")
    inv.add_argument("--admin", action="store_true")
    inv.add_argument("--base", default="http://localhost:8000")
    loc = sub.add_parser("local", help="the desktop app: this computer, your own folders, opens the browser")
    loc.add_argument("--port", type=int, default=0)
    loc.add_argument("--no-browser", action="store_true")
    srv = sub.add_parser("serve")
    srv.add_argument("--port", type=int, default=8000)
    srv.add_argument("--host", default="127.0.0.1")
    a = ap.parse_args()
    if a.cmd == "invite":
        print(f"{a.base}/login?token={invite(a.email, a.admin)}")
        return
    import uvicorn
    if a.cmd == "local":
        import socket
        import webbrowser
        LOCAL.update(on=True, key=secrets.token_urlsafe(24))
        port = a.port
        if not port:
            with socket.socket() as s_:
                s_.bind(("127.0.0.1", 0))
                port = s_.getsockname()[1]
        url = f"http://127.0.0.1:{port}/local?key={LOCAL['key']}"
        print(f"Rooh is running on this computer. If no browser opened, open:\n  {url}\nClose this window to stop Rooh.",
              flush=True)
        if not a.no_browser:
            threading.Timer(1.5, lambda: webbrowser.open(url)).start()
        start_worker()                              # measuring on this computer (its GPU if it has one)
        uvicorn.run(app, host="127.0.0.1", port=port, log_level="warning")
        return
    start_worker()
    uvicorn.run(app, host=a.host, port=a.port)


if __name__ == "__main__":
    main()
