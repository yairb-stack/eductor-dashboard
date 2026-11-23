#!/bin/bash
cd "$(dirname "$0")"
open "http://127.0.0.1:8000/dashboard/analytical-exposure.html"
python3 -m http.server 8000
