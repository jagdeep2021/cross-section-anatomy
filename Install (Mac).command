#!/bin/bash
# Double-click once. Installs everything into this folder (.uv/ and .venv/): no admin rights,
# nothing elsewhere on the computer. Delete the folder to uninstall. Takes ~5-10 min (downloads ~1.5 GB).
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
chmod +x *.command
.venv/bin/python -c "import napari, cellpose, pipeline, review; print('\nAll set. Double-click \"Review (Mac).command\" to start.')"
