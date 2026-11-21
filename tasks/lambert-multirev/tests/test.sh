#!/bin/bash
# Runs the verifier and writes the reward: 1 if every test passes, 0 otherwise.
TESTS_DIR="${TESTS_DIR:-/tests}"
LOGS_DIR="${LOGS_DIR:-/logs/verifier}"
mkdir -p "$LOGS_DIR"

if python3 -m pytest -q -rA "$TESTS_DIR/test_outputs.py" -p no:cacheprovider > "$LOGS_DIR/pytest.log" 2>&1; then
    echo 1 > "$LOGS_DIR/reward.txt"
else
    echo 0 > "$LOGS_DIR/reward.txt"
fi
tail -n 30 "$LOGS_DIR/pytest.log"
