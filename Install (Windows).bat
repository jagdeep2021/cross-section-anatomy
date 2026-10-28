@echo off
rem Double-click once. Installs everything into this folder (.uv\ and .venv\): no admin rights,
rem nothing elsewhere on the computer. Delete the folder to uninstall. Takes ~5-10 min (downloads ~1.5 GB).
cd /d "%~dp0"
if not exist ".uv\uv.exe" (
  echo Getting the installer ^(uv^)...
  powershell -NoProfile -ExecutionPolicy Bypass -Command "$env:UV_INSTALL_DIR='%CD%\.uv'; $env:UV_NO_MODIFY_PATH='1'; irm https://astral.sh/uv/install.ps1 | iex"
)
set UV_PYTHON_INSTALL_DIR=%CD%\.uv\python
set UV_CACHE_DIR=%CD%\.uv\cache
echo Installing Python and the packages...
.uv\uv.exe venv --python 3.14 --allow-existing .venv || goto :fail
.uv\uv.exe pip install --python .venv\Scripts\python.exe -r requirements.txt || goto :fail
.venv\Scripts\python.exe -c "import napari, cellpose, pipeline, review; print('All set. Double-click Open app (Windows).bat to start.')" || goto :fail
copy /y requirements.txt .venv\installed_requirements.txt >nul
pause
exit /b 0
:fail
echo.
echo Something went wrong. Send Jagdeep a photo of this window.
pause
