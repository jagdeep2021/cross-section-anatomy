#!/bin/bash
# Hellbender (or any Linux): from a terminal inside an Open OnDemand desktop:
#   bash "Open app (Linux).sh"
# After "git pull" nothing else is needed: if the tool needs new packages, this installs them first.
cd "$(dirname "$0")"
if ! command -v nvidia-smi >/dev/null || ! nvidia-smi >/dev/null 2>&1; then
  echo "No GPU in this session: measuring new images will be very slow."
  echo "Close this desktop and start a new one on the GPU partition (see README, 'On Hellbender')."
  echo "Press Enter to open the app anyway (fine for editing images that are already measured)."
  read -r
fi
if [ -x .uv/uv ] && ! cmp -s requirements.txt .venv/installed_requirements.txt; then
  echo "Updating the tool's packages (only after an update; a few minutes)..."
  export UV_PYTHON_INSTALL_DIR="$PWD/.uv/python" UV_CACHE_DIR="$PWD/.uv/cache"
  .uv/uv pip install --python .venv/bin/python -r requirements.txt || { echo "Updating failed. Send Jagdeep a screenshot."; read -r; exit 1; }
  cp requirements.txt .venv/installed_requirements.txt
fi
if ! /sbin/ldconfig -p 2>/dev/null | grep -q "libxcb-cursor.so.0" && [ ! -e .venv/xcb-cursor/libxcb-cursor.so.0 ]; then
  echo "Getting one library this Linux is missing (xcb-cursor, needed for the window)..."
  .venv/bin/python - <<'PY' || { echo "Could not get it. Send Jagdeep a screenshot."; read -r; exit 1; }
# Hellbender's Linux (glibc 2.28) lacks libxcb-cursor, which Qt 6.5+ needs to open a window.
# Fetch the conda-forge build (glibc >= 2.17) into .venv/xcb-cursor/ once; no admin rights.
import hashlib, io, pathlib, sys, tarfile, urllib.request, zipfile
from compression import zstd
URL = "https://conda.anaconda.org/conda-forge/linux-64/xcb-util-cursor-0.1.6-hb03c661_0.conda"
SHA = "c2be9cae786fdb2df7c2387d2db31b285cf90ab3bfabda8fa75a596c3d20fc67"
out = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else ".venv/xcb-cursor")
data = urllib.request.urlopen(URL, timeout=60).read()
if hashlib.sha256(data).hexdigest() != SHA:
    sys.exit("xcb-cursor download failed its checksum; try again or send Jagdeep a screenshot.")
pkg = zipfile.ZipFile(io.BytesIO(data))
inner = next(n for n in pkg.namelist() if n.startswith("pkg-") and n.endswith(".tar.zst"))
with tarfile.open(fileobj=io.BytesIO(zstd.decompress(pkg.read(inner)))) as t:
    real = next(m for m in t.getmembers() if m.isfile() and m.name.startswith("lib/libxcb-cursor.so.0."))
    out.mkdir(parents=True, exist_ok=True)
    (out / "libxcb-cursor.so.0").write_bytes(t.extractfile(real).read())
print(f"xcb-cursor ready ({real.name})")
PY
fi
export LD_LIBRARY_PATH="$PWD/.venv/xcb-cursor${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
.venv/bin/python -W ignore web/server.py local   # the web page on this computer; classic editor: .venv/bin/python review.py
