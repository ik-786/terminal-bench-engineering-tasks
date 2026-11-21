"""Verifier for the multi-revolution Lambert task.

Each reported departure velocity is checked by flying it: the verifier
propagates (r1, v1) for the time of flight with its own universal-variable
Kepler propagator and requires the result to land on r2 with velocity v2. It
also checks the direction of motion and the number of complete revolutions.
For multi-revolution problems the requested branch is checked through the
semi-major axis, and problems without a solution must be reported as null.
"""

import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

APP = Path(os.environ.get("APP_DIR", "/app"))
HERE = Path(__file__).resolve().parent
SOLVER = APP / "lambert.py"
PROBLEMS = HERE / "cases" / "problems.json"
EXPECTED = json.loads((HERE / "cases" / "expected.json").read_text())
DATA = json.loads(PROBLEMS.read_text())
MU = DATA["mu"]
BY_ID = {p["id"]: p for p in DATA["problems"]}
TIME_LIMIT = 120
POS_TOL = 1e-6  # relative
VEL_TOL = 1e-6  # relative
A_TOL = 1e-6    # relative, semi-major axis of the requested branch


def dot(a, b): return sum(x * y for x, y in zip(a, b))
def norm(a): return math.sqrt(dot(a, a))
def cross(a, b): return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def stumpff(z):
    if z > 1e-8:
        s = math.sqrt(z)
        return (1 - math.cos(s)) / z, (s - math.sin(s)) / s ** 3
    if z < -1e-8:
        s = math.sqrt(-z)
        return (math.cosh(s) - 1) / -z, (math.sinh(s) - s) / s ** 3
    return 1 / 2 - z / 24 + z * z / 720, 1 / 6 - z / 120 + z * z / 5040


def propagate(r0, v0, t, mu):
    """Universal-variable Kepler propagation (any conic)."""
    R0, V0 = norm(r0), norm(v0)
    vr0 = dot(r0, v0) / R0
    alpha = 2 / R0 - V0 * V0 / mu  # 1/a
    if alpha > 1e-12:  # ellipse: drop whole periods first to keep the anomaly small
        period = 2 * math.pi / math.sqrt(mu * alpha ** 3)
        t = math.fmod(t, period)
    sm = math.sqrt(mu)
    chi = sm * abs(alpha) * t if alpha > 1e-12 else sm * t / R0
    for _ in range(200):
        z = alpha * chi * chi
        C, S = stumpff(z)
        F = R0 * vr0 / sm * chi * chi * C + (1 - alpha * R0) * chi ** 3 * S + R0 * chi - sm * t
        dF = R0 * vr0 / sm * chi * (1 - z * S) + (1 - alpha * R0) * chi * chi * C + R0
        step = F / dF
        chi -= step
        if abs(step) < 1e-13 * max(1.0, abs(chi)):
            break
    z = alpha * chi * chi
    C, S = stumpff(z)
    f = 1 - chi * chi / R0 * C
    g = t - chi ** 3 / sm * S
    r = [f * a + g * b for a, b in zip(r0, v0)]
    R = norm(r)
    fd = sm / (R * R0) * (alpha * chi ** 3 * S - chi)
    gd = 1 - chi * chi / R * C
    v = [fd * a + gd * b for a, b in zip(r0, v0)]
    return r, v


def mean_anomaly(r, v, mu):
    R = norm(r)
    alpha = 2 / R - dot(v, v) / mu
    a = 1 / alpha
    e_cosE = 1 - R / a
    e_sinE = dot(r, v) / math.sqrt(mu * a)
    E = math.atan2(e_sinE, e_cosE)
    e = math.hypot(e_cosE, e_sinE)
    return E - e * math.sin(E)


def full_revolutions(r1, v1, t, mu):
    alpha = 2 / norm(r1) - dot(v1, v1) / mu
    if alpha <= 0:
        return 0
    n = math.sqrt(mu * alpha ** 3)
    r2, v2 = propagate(r1, v1, t, mu)
    frac = (mean_anomaly(r2, v2, mu) - mean_anomaly(r1, v1, mu)) % (2 * math.pi)
    return round((n * t - frac) / (2 * math.pi))


@pytest.fixture(scope="module")
def solutions(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("lambert")
    out = tmp / "solutions.json"
    proc = subprocess.run([sys.executable, str(SOLVER), str(PROBLEMS), str(out)],
                          capture_output=True, text=True, timeout=TIME_LIMIT, cwd=tmp)
    assert proc.returncode == 0, f"solver exited with {proc.returncode}: {proc.stderr[-800:]}"
    return json.loads(out.read_text())["solutions"]


def test_solver_exists():
    assert SOLVER.exists(), "expected the solver at /app/lambert.py"


def test_solver_does_not_touch_grader_files():
    src = SOLVER.read_text(errors="replace")
    assert not re.search(r"/tests\b|/logs\b|reward\.txt", src), "the solver must not read or write grader files"


def test_every_problem_answered(solutions):
    assert set(solutions) == set(BY_ID), "solutions must contain every problem id exactly once"


@pytest.mark.parametrize("pid", sorted(EXPECTED))
def test_problem(pid, solutions):
    p, exp = BY_ID[pid], EXPECTED[pid]
    mu = p.get("mu", MU)
    sol = solutions.get(pid, "missing")
    if not exp["feasible"]:
        assert sol is None, "no transfer exists for this problem; the answer must be null"
        return
    assert isinstance(sol, dict), "this problem has a solution"
    v1, v2 = sol["v1"], sol["v2"]
    r_end, v_end = propagate(p["r1"], v1, p["tof"], mu)
    assert norm([a - b for a, b in zip(r_end, p["r2"])]) <= POS_TOL * norm(p["r2"]), \
        "flying v1 for the time of flight does not arrive at r2"
    assert norm([a - b for a, b in zip(v_end, v2)]) <= VEL_TOL * norm(v_end), \
        "v2 is not the arrival velocity of the transfer"

    h = cross(p["r1"], v1)
    if "normal" in p:
        assert dot(h, p["normal"]) > 0, "transfer must turn about the given normal"
    elif p.get("direction", "prograde") == "prograde":
        assert h[2] > 0, "transfer must be prograde (angular momentum along +z)"
    else:
        assert h[2] < 0, "transfer must be retrograde (angular momentum along -z)"

    assert full_revolutions(p["r1"], v1, p["tof"], mu) == p.get("revs", 0), \
        "wrong number of complete revolutions"
    if p.get("revs", 0) > 0:
        a = 1 / (2 / norm(p["r1"]) - dot(v1, v1) / mu)
        assert abs(a - exp["a_km"]) <= A_TOL * exp["a_km"], \
            f"semi-major axis {a:.6f} km is not the requested '{p['branch']}' branch"
