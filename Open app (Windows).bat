@echo off
rem Double-click: opens the app (Measure images..., scale, Summarise, checking outlines).
rem Keep this black window open while you work; if something goes wrong, the message shows here.
cd /d "%~dp0"
.venv\Scripts\python.exe -W ignore review.py
if errorlevel 1 pause
