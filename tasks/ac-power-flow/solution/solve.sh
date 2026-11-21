#!/bin/bash
# Reference solution: installs the Newton-Raphson power flow solver at /app/powerflow.py.
set -euo pipefail
APP_DIR="${APP_DIR:-/app}"
cp "$(dirname "$0")/powerflow.py" "$APP_DIR/powerflow.py"
