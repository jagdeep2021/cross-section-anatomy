@echo off
rem Double-click: choose the folder with the images (e.g. the experiment's raw folder). The computer outlines
rem every image and writes anatomy\ next to that folder. Stop any time; running it again continues where it left off.
cd /d "%~dp0"
for /f "usebackq delims=" %%D in (`powershell -NoProfile -Command "Add-Type -AssemblyName System.Windows.Forms; $f=New-Object System.Windows.Forms.FolderBrowserDialog; $f.Description='Choose the folder with the images (e.g. raw)'; if($f.ShowDialog() -eq 'OK'){$f.SelectedPath}"`) do set DIR=%%D
if not defined DIR exit /b 0
.venv\Scripts\python.exe -W ignore run_batch.py "%DIR%"
echo Done.
pause
