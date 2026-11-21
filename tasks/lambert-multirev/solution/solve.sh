#!/bin/bash
# Reference solution: installs the multi-revolution Lambert solver at /app/lambert.py.
set -euo pipefail
APP_DIR="${APP_DIR:-/app}"
cp "$(dirname "$0")/lambert.py" "$APP_DIR/lambert.py"
