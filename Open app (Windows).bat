@echo off
rem Double-click: opens the app. After "git pull" nothing else is needed: if the tool needs new
rem packages, this installs them first. Keep this black window open; messages show here.
cd /d "%~dp0"
if not exist .uv\uv.exe goto :run
fc /b requirements.txt .venv\installed_requirements.txt >nul 2>&1
if errorlevel 1 (
  echo Updating the tool's packages ^(only after an update; a few minutes^)...
  set UV_PYTHON_INSTALL_DIR=%CD%\.uv\python
  set UV_CACHE_DIR=%CD%\.uv\cache
  .uv\uv.exe pip install --python .venv\Scripts\python.exe -r requirements.txt || goto :fail
  copy /y requirements.txt .venv\installed_requirements.txt >nul
)
:run
.venv\Scripts\python.exe -W ignore review.py
if errorlevel 1 pause
exit /b 0
:fail
echo.
echo Updating failed. Send Jagdeep a photo of this window.
pause
