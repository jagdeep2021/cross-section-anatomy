#!/bin/bash
# Double-click: opens the app. After "git pull" nothing else is needed: if the tool needs new
# packages, this installs them first.
cd "$(dirname "$0")"
if [ -x .uv/uv ] && ! cmp -s requirements.txt .venv/installed_requirements.txt; then
  echo "Updating the tool's packages (only after an update; a few minutes)..."
  export UV_PYTHON_INSTALL_DIR="$PWD/.uv/python" UV_CACHE_DIR="$PWD/.uv/cache"
  .uv/uv pip install --python .venv/bin/python -r requirements.txt || { echo "Updating failed. Send Jagdeep a screenshot."; read -r; exit 1; }
  cp requirements.txt .venv/installed_requirements.txt
fi
.venv/bin/python -W ignore review.py
