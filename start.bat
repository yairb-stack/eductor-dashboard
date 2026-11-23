@echo off
cd /d "%~dp0"
where py >nul 2>nul && set PY=py
if not defined PY where python >nul 2>nul && set PY=python
if not defined PY where python3 >nul 2>nul && set PY=python3
if not defined PY (
  echo Python is required. Install it from https://www.python.org/downloads/ and check "Add python.exe to PATH", then run this file again.
  pause
  exit /b 1
)
start "" "http://127.0.0.1:8000/dashboard/analytical-exposure.html"
%PY% -m http.server 8000
