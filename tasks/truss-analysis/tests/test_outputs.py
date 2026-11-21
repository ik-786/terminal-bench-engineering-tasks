"""Verifier for the truss analysis task.

For a stable truss the linear-elastic solution is unique, and it is the only
set of displacements, member forces and reactions that satisfies all three of:
member constitutive law, nodal equilibrium and the support conditions. This file
checks those three directly, so no stored answers are needed. For a mechanism
the solver must report that the truss is unstable.
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
CASES = Path(__file__).resolve().parent / "cases"
SOLVER = APP / "truss.py"
EXPECTED = json.loads((CASES / "expected.json").read_text())
TIME_LIMIT = 120
REL = 1e-6  # relative tolerance on forces, scaled by the largest force in the problem


def run_solver(name, tmp_path):
    out = tmp_path / f"{name}_result.json"
    proc = subprocess.run([sys.executable, str(SOLVER), str(CASES / f"{name}.json"), str(out)],
                          capture_output=True, text=True, timeout=TIME_LIMIT, cwd=tmp_path)
    assert proc.returncode == 0, f"solver exited with {proc.returncode}: {proc.stderr[-800:]}"
    assert out.exists(), "solver did not write the result file"
    return json.loads(out.read_text())


def check_solution(case, res):
    xy = {n["id"]: (n["x"], n["y"]) for n in case["nodes"]}
    u = {k: tuple(v) for k, v in res["displacements"].items()}
    N = res["member_forces"]
    R = {k: tuple(v) for k, v in res["reactions"].items()}
    supports = {s["node"]: s for s in case["supports"]}

    assert set(u) == set(xy), "displacements must be given for every node"
    assert set(N) == {m["id"] for m in case["members"]}, "member_forces must be given for every member"
    assert set(R) == set(supports), "reactions must be given for every supported node, and only those"

    loads = {k: [0.0, 0.0] for k in xy}
    for ld in case.get("loads", []):
        loads[ld["node"]][0] += ld.get("fx", 0.0)
        loads[ld["node"]][1] += ld.get("fy", 0.0)

    fscale = max([abs(v) for v in N.values()] + [math.hypot(*l) for l in loads.values()] +
                 [m["E"] * m["A"] * 1e-6 for m in case["members"]])
    ftol = REL * fscale

    # 1. constitutive law: N = EA/L * (elongation - free elongation)
    for m in case["members"]:
        (xi, yi), (xj, yj) = xy[m["i"]], xy[m["j"]]
        L = math.hypot(xj - xi, yj - yi)
        c, s = (xj - xi) / L, (yj - yi) / L
        e = c * (u[m["j"]][0] - u[m["i"]][0]) + s * (u[m["j"]][1] - u[m["i"]][1])
        e0 = m.get("alpha", 0.0) * m.get("dT", 0.0) * L + m.get("lack_of_fit", 0.0)
        expected = m["E"] * m["A"] / L * (e - e0)
        assert abs(N[m["id"]] - expected) <= ftol, (
            f"member {m['id']}: force {N[m['id']]:.6g} N does not match its elongation ({expected:.6g} N)")

    # 2. equilibrium at every node
    resid = {k: [loads[k][0], loads[k][1]] for k in xy}
    for k, (rx, ry) in R.items():
        resid[k][0] += rx
        resid[k][1] += ry
    for m in case["members"]:
        (xi, yi), (xj, yj) = xy[m["i"]], xy[m["j"]]
        L = math.hypot(xj - xi, yj - yi)
        c, s = (xj - xi) / L, (yj - yi) / L
        n = N[m["id"]]  # tension pulls each end towards the other
        resid[m["i"]][0] += n * c
        resid[m["i"]][1] += n * s
        resid[m["j"]][0] -= n * c
        resid[m["j"]][1] -= n * s
    for k, (fx, fy) in resid.items():
        assert math.hypot(fx, fy) <= ftol, f"node {k} is not in equilibrium (residual {math.hypot(fx, fy):.4g} N)"

    # 3. support conditions
    dscale = max([math.hypot(*v) for v in u.values()] + [1e-9])
    for k, s in supports.items():
        ux, uy = u[k]
        rx, ry = R[k]
        if s["type"] == "pin":
            sx, sy = s.get("settlement", (0.0, 0.0))
            assert abs(ux - sx) <= 1e-9 + 1e-9 * dscale and abs(uy - sy) <= 1e-9 + 1e-9 * dscale, \
                f"pin {k}: displacement must equal its settlement"
        else:
            a = math.radians(s["angle_deg"])
            nx, ny = math.cos(a), math.sin(a)
            assert abs(ux * nx + uy * ny) <= 1e-9 + 1e-9 * dscale, f"roller {k}: moves through its support surface"
            assert abs(-rx * ny + ry * nx) <= ftol, f"roller {k}: reaction has a component along the surface"


def test_solver_exists():
    assert SOLVER.exists(), "expected the solver at /app/truss.py"


def test_solver_does_not_touch_grader_files():
    src = SOLVER.read_text(errors="replace")
    assert not re.search(r"/tests\b|/logs\b|reward\.txt", src), "the solver must not read or write grader files"


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_case(name, tmp_path):
    case = json.loads((CASES / f"{name}.json").read_text())
    res = run_solver(name, tmp_path)
    if not EXPECTED[name]["stable"]:
        assert res.get("stable") is False, "this truss is a mechanism; stable must be false"
        return
    assert res.get("stable") is True, "this truss is stable; stable must be true"
    check_solution(case, res)
