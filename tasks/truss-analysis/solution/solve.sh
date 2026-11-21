#!/bin/bash
# Reference solution: installs the stiffness-method truss solver at /app/truss.py.
set -euo pipefail
APP_DIR="${APP_DIR:-/app}"
cp "$(dirname "$0")/truss.py" "$APP_DIR/truss.py"
