#!/bin/bash
# Double-click: choose an experiment's "anatomy" folder; the tables and figures are written to its results/ folder.
cd "$(dirname "$0")"
DIR=$(osascript -e 'POSIX path of (choose folder with prompt "Choose the anatomy folder of the experiment")') || exit 0
.venv/bin/python -W ignore summarise.py "$DIR" && open "$DIR/../results"
read -n 1 -s -r -p "Press any key to close"
