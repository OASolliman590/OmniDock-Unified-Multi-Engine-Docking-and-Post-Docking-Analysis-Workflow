@echo off
REM Launch the OmniDock web UI using the project's virtual environment.
REM Double-click this file, or run it from any shell -- no PATH setup needed.

setlocal
cd /d "%~dp0"

set "VENV_PY=%~dp0.venv\Scripts\python.exe"

if not exist "%VENV_PY%" (
    echo.
    echo   The project virtual environment was not found at:
    echo     %VENV_PY%
    echo.
    echo   Create it once with:
    echo     py -3.12 -m venv .venv
    echo     .venv\Scripts\python.exe -m pip install -r requirements.txt
    echo.
    pause
    exit /b 1
)

echo Starting the OmniDock web UI...
echo Open http://127.0.0.1:8770 in your browser. Press Ctrl+C here to stop.
echo.

"%VENV_PY%" main.py webui %*

REM Keep the window open if it exits unexpectedly, so the error stays readable.
if errorlevel 1 pause
endlocal
