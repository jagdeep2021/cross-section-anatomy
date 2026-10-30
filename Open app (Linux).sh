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
.venv/bin/python -W ignore review.py
