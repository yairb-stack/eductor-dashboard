#!/bin/bash
cd "$(dirname "$0")"
if command -v xdg-open >/dev/null 2>&1; then
  xdg-open "http://127.0.0.1:8000/dashboard/analytical-exposure.html" >/dev/null 2>&1 &
fi
python3 -m http.server 8000
