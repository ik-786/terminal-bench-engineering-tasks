"""Runs a task's verifier on this machine, without Docker.

    python scripts/run_local.py <task> --oracle          # reference solution, expect reward 1
    python scripts/run_local.py <task> --null            # no solution, expect reward 0
    python scripts/run_local.py <task> --solution FILE   # any candidate file, installed as the task's target
    python scripts/run_local.py --all                    # every task: oracle, null and every file in validation/

It mirrors the container layout in a scratch directory (APP_DIR, TESTS_DIR,
LOGS_DIR) and prints the reward. Docker is still the reference way to run a
task; this is for fast iteration and for CI.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TARGETS = {"ac-power-flow": "powerflow.py", "truss-analysis": "truss.py", "lambert-multirev": "lambert.py"}


def run(task: str, mode: str, candidate: Path | None = None, quiet: bool = False) -> int:
    tdir = ROOT / "tasks" / task
    with tempfile.TemporaryDirectory(prefix=f"tb_{task}_") as tmp:
        tmp = Path(tmp)
        app, tests, logs = tmp / "app", tmp / "tests", tmp / "logs"
        app.mkdir()
        data = tdir / "environment" / "data"
        if data.exists():
            shutil.copytree(data, app / "data")
        shutil.copytree(tdir / "tests", tests)
        env = {**os.environ, "APP_DIR": str(app), "TESTS_DIR": str(tests), "LOGS_DIR": str(logs)}
        if mode == "oracle":
            subprocess.run(["bash", str(tdir / "solution" / "solve.sh")], env=env, check=True)
        elif mode == "candidate":
            shutil.copy(candidate, app / TARGETS[task])
        proc = subprocess.run(["bash", str(tests / "test.sh")], env=env, capture_output=True, text=True)
        reward = int((logs / "reward.txt").read_text().strip())
        if not quiet:
            print(proc.stdout[-3000:])
        return reward


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("task", nargs="?", choices=sorted(TARGETS))
    g = p.add_mutually_exclusive_group()
    g.add_argument("--oracle", action="store_true")
    g.add_argument("--null", action="store_true")
    g.add_argument("--solution", type=Path)
    g.add_argument("--all", action="store_true")
    a = p.parse_args()

    if a.all:
        failures = 0
        for task in sorted(TARGETS):
            checks = [("oracle", None, 1), ("null", None, 0)]
            checks += [("candidate", f, 0) for f in sorted((ROOT / "validation" / task).glob("*.py"))]
            for mode, f, want in checks:
                got = run(task, mode, f, quiet=True)
                label = f.name if f else mode
                status = "ok " if got == want else "BAD"
                failures += got != want
                print(f"{status} {task:18s} {label:34s} reward={got} (expected {want})")
        return 1 if failures else 0

    mode = "oracle" if a.oracle else "null" if a.null else "candidate"
    reward = run(a.task, mode, a.solution)
    print(f"reward: {reward}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
