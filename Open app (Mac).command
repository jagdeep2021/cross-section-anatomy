#!/bin/bash
# Double-click: opens the app (Measure images..., scale, Summarise, checking outlines).
cd "$(dirname "$0")"
.venv/bin/python -W ignore review.py
