#!/bin/bash
# Double-click: choose the folder with the images (e.g. the experiment's raw/ folder). The computer outlines
# every image and writes anatomy/ next to that folder. Stop any time; running it again continues where it left off.
cd "$(dirname "$0")"
DIR=$(osascript -e 'POSIX path of (choose folder with prompt "Choose the folder with the images (e.g. raw)")') || exit 0
.venv/bin/python -W ignore run_batch.py "$DIR"
read -n 1 -s -r -p "Done. Press any key to close"
