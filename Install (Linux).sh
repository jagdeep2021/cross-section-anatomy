#!/bin/bash
# Hellbender (or any Linux): run once from a terminal inside an Open OnDemand desktop:
#   bash "Install (Linux).sh"
# Installs everything into this folder (.uv/ and .venv/): no admin rights, nothing elsewhere.
# Delete the folder to uninstall. Takes ~10 min (downloads ~4 GB: the Linux torch carries CUDA).
cd "$(dirname "$0")"
set -e
if [ ! -x .uv/uv ]; then
  echo "Getting the installer (uv)..."
  curl -LsSf https://astral.sh/uv/install.sh | env UV_INSTALL_DIR="$PWD/.uv" UV_NO_MODIFY_PATH=1 sh
fi
export UV_PYTHON_INSTALL_DIR="$PWD/.uv/python" UV_CACHE_DIR="$PWD/.uv/cache"
echo "Installing Python and the packages..."
.uv/uv venv --python 3.14 --allow-existing .venv
.uv/uv pip install --python .venv/bin/python -r requirements.txt
.venv/bin/python -c "import napari, cellpose, pipeline, review; print('\nAll set. Start the app with:  bash \"Open app (Linux).sh\"')"
cp requirements.txt .venv/installed_requirements.txt
