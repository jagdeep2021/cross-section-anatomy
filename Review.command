#!/bin/bash
# Double-click: choose an experiment's "anatomy" folder, then the review app opens.
cd "$(dirname "$0")"
DIR=$(osascript -e 'POSIX path of (choose folder with prompt "Choose the anatomy folder of the experiment")') || exit 0
.venv/bin/python -W ignore review.py "$DIR"
