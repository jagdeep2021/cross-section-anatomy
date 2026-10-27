@echo off
rem Double-click: choose an experiment's "anatomy" folder, then the review app opens.
cd /d "%~dp0"
for /f "usebackq delims=" %%D in (`powershell -NoProfile -Command "Add-Type -AssemblyName System.Windows.Forms; $f=New-Object System.Windows.Forms.FolderBrowserDialog; $f.Description='Choose the anatomy folder of the experiment'; if($f.ShowDialog() -eq 'OK'){$f.SelectedPath}"`) do set DIR=%%D
if not defined DIR exit /b 0
.venv\Scripts\python.exe -W ignore review.py "%DIR%"
