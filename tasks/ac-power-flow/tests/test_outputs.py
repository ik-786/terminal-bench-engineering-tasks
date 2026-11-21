"""Verifier for the AC power flow task.

The agent's solver is run on hidden networks. A result is accepted only if the
reported bus voltages satisfy the power flow equations, which this file checks
from scratch with complex matrix arithmetic. Nothing is compared with stored
answers, so the only way to pass is to actually solve the network.
"""

import cmath
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
SOLVER = APP / "powerflow.py"
EXPECTED = json.loads((CASES / "expected.json").read_text())

TOL_MW = 1e-3     # power mismatch, MW or MVAr (1e-5 per unit on a 100 MVA base)
TOL_V = 1e-6      # voltage set points, per unit
TOL_DEG = 1e-6    # slack angle, degrees
TIME_LIMIT = 120  # seconds per case


def run_solver(name, tmp_path):
    out = tmp_path / f"{name}_result.json"
    proc = subprocess.run([sys.executable, str(SOLVER), str(CASES / f"{name}.json"), str(out)],
                          capture_output=True, text=True, timeout=TIME_LIMIT, cwd=tmp_path)
    assert proc.returncode == 0, f"solver exited with {proc.returncode}: {proc.stderr[-800:]}"
    assert out.exists(), "solver did not write the result file"
    return json.loads(out.read_text())


def admittance(case):
    """Bus admittance matrix as a dict of dicts, keyed by bus id."""
    base = case["base_mva"]
    Y = {b["id"]: {} for b in case["buses"]}

    def add(i, k, y):
        Y[i][k] = Y[i].get(k, 0j) + y

    for b in case["buses"]:
        add(b["id"], b["id"], complex(b.get("gs", 0.0), b.get("bs", 0.0)) / base)
    for br in case["branches"]:
        if br.get("status", 1) == 0:
            continue
        ys = 1 / complex(br["r"], br["x"])
        a = (br.get("tap", 1.0) or 1.0) * cmath.exp(1j * math.radians(br.get("shift_deg", 0.0)))
        sh = 1j * br.get("b", 0.0) / 2
        add(br["from"], br["from"], (ys + sh) / (a * a.conjugate()))
        add(br["from"], br["to"], -ys / a.conjugate())
        add(br["to"], br["from"], -ys / a)
        add(br["to"], br["to"], ys + sh)
    return Y


def branch_flows(case, br, Vc):
    base = case["base_mva"]
    ys = 1 / complex(br["r"], br["x"])
    a = (br.get("tap", 1.0) or 1.0) * cmath.exp(1j * math.radians(br.get("shift_deg", 0.0)))
    sh = 1j * br.get("b", 0.0) / 2
    vf, vt = Vc[br["from"]], Vc[br["to"]]
    i_from = (ys + sh) / (a * a.conjugate()) * vf - ys / a.conjugate() * vt
    i_to = -ys / a * vf + (ys + sh) * vt
    return vf * i_from.conjugate() * base, vt * i_to.conjugate() * base


def check_solution(case, res):
    base = case["base_mva"]
    buses = {b["id"]: b for b in case["buses"]}
    got = {b["id"]: b for b in res["buses"]}
    assert set(got) == set(buses), "result must list every bus exactly once"
    for b in res["buses"]:
        assert 0.5 <= b["v"] <= 1.5, f"bus {b['id']}: voltage {b['v']} outside 0.5-1.5 pu"

    Vc = {i: cmath.rect(got[i]["v"], math.radians(got[i]["theta_deg"])) for i in buses}
    Y = admittance(case)
    S = {i: Vc[i] * sum((y * Vc[k] for k, y in Y[i].items()), 0j).conjugate() * base for i in buses}

    limited = set(res.get("q_limited", []))
    for i, b in buses.items():
        p_calc, q_calc = S[i].real, S[i].imag
        r = got[i]
        assert abs(r["p_inj_mw"] - p_calc) < TOL_MW, f"bus {i}: reported P injection disagrees with the voltages"
        assert abs(r["q_inj_mvar"] - q_calc) < TOL_MW, f"bus {i}: reported Q injection disagrees with the voltages"
        pd, qd = b.get("pd", 0.0), b.get("qd", 0.0)
        if b["type"] == "slack":
            assert abs(r["v"] - b["v"]) < TOL_V, f"slack bus {i}: voltage magnitude changed"
            assert abs(r["theta_deg"] - b.get("theta_deg", 0.0)) < TOL_DEG, f"slack bus {i}: angle changed"
        elif b["type"] == "pq":
            assert abs(p_calc - (b.get("pg", 0.0) - pd)) < TOL_MW, f"PQ bus {i}: active power mismatch"
            assert abs(q_calc - (b.get("qg", 0.0) - qd)) < TOL_MW, f"PQ bus {i}: reactive power mismatch"
        else:  # PV bus with reactive limits
            assert abs(p_calc - (b["pg"] - pd)) < TOL_MW, f"PV bus {i}: active power mismatch"
            qg = q_calc + qd
            at_set = abs(r["v"] - b["v"]) < TOL_V and b["qmin"] - TOL_MW <= qg <= b["qmax"] + TOL_MW
            at_max = abs(qg - b["qmax"]) < TOL_MW and r["v"] <= b["v"] + TOL_V
            at_min = abs(qg - b["qmin"]) < TOL_MW and r["v"] >= b["v"] - TOL_V
            assert at_set or at_max or at_min, (
                f"PV bus {i}: Qg={qg:.4f} MVAr, V={r['v']:.6f} pu violates the limit rules "
                f"(set point {b['v']}, limits [{b['qmin']}, {b['qmax']}])")
            held = (at_max or at_min) and not at_set
            assert (i in limited) == held, f"PV bus {i}: q_limited list does not match its state"

    flows = {(f["from"], f["to"]): f for f in res["branches"]}
    live = [br for br in case["branches"] if br.get("status", 1) != 0]
    assert len(flows) == len(live), "result must list every in-service branch once"
    for br in live:
        f = flows[(br["from"], br["to"])]
        s_from, s_to = branch_flows(case, br, Vc)
        assert abs(f["p_from_mw"] - s_from.real) < TOL_MW and abs(f["q_from_mvar"] - s_from.imag) < TOL_MW, \
            f"branch {br['from']}-{br['to']}: from-end flow is wrong"
        assert abs(f["p_to_mw"] - s_to.real) < TOL_MW and abs(f["q_to_mvar"] - s_to.imag) < TOL_MW, \
            f"branch {br['from']}-{br['to']}: to-end flow is wrong"
    assert abs(res["losses_mw"] - sum(s.real for s in S.values())) < TOL_MW, "losses_mw is wrong"


def test_solver_exists():
    assert SOLVER.exists(), "expected the solver at /app/powerflow.py"


def test_solver_does_not_touch_grader_files():
    src = SOLVER.read_text(errors="replace")
    assert not re.search(r"/tests\b|/logs\b|reward\.txt", src), "the solver must not read or write grader files"


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_case(name, tmp_path):
    case = json.loads((CASES / f"{name}.json").read_text())
    res = run_solver(name, tmp_path)
    if not EXPECTED[name]["converged"]:
        assert res.get("converged") is False, "this case has no solution; converged must be false"
        return
    assert res.get("converged") is True, "this case has a solution; converged must be true"
    check_solution(case, res)
