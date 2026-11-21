"""Builds the example problems (visible to the agent) and the hidden test problems.

Run from the task folder:  python dev/make_cases.py
"""

import importlib.util
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE / "solution"))
import lambert  # noqa: E402  (reference solver)

MU_EARTH = 398600.4418
MU_SUN = 1.32712440018e11
AU = 149597870.7


def verifier():
    spec = importlib.util.spec_from_file_location("verifier", HERE / "tests" / "test_outputs.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def circ(r, deg, inc_deg=0.0, raan_deg=0.0):
    """Position on a circle of radius r at argument deg, in an inclined plane."""
    u, i, O = math.radians(deg), math.radians(inc_deg), math.radians(raan_deg)
    x, y = r * math.cos(u), r * math.sin(u)
    return [x * math.cos(O) - y * math.cos(i) * math.sin(O),
            x * math.sin(O) + y * math.cos(i) * math.cos(O),
            y * math.sin(i)]


def period(a, mu):
    return 2 * math.pi * math.sqrt(a ** 3 / mu)


def hidden():
    R1, R2 = 7000.0, 42164.0
    t_h = period((R1 + R2) / 2, MU_EARTH) / 2
    P = [
        {"id": "leo_to_geo_short", "r1": circ(R1, 0), "r2": circ(R2, 120), "tof": 0.7 * t_h, "revs": 0, "direction": "prograde"},
        {"id": "long_way_prograde", "r1": circ(R1, 0), "r2": circ(R2, 250), "tof": 1.3 * t_h, "revs": 0, "direction": "prograde"},
        {"id": "retrograde", "r1": circ(9000, 10), "r2": circ(15000, 100), "tof": 9000.0, "revs": 0, "direction": "retrograde"},
        {"id": "hyperbolic_fast", "r1": circ(R1, 0), "r2": circ(30000, 60), "tof": 1800.0, "revs": 0, "direction": "prograde"},
        {"id": "near_parabolic", "r1": circ(8000, 0), "r2": circ(20000, 90), "tof": None, "revs": 0, "direction": "prograde"},
        {"id": "half_turn_with_normal", "r1": [R1, 0.0, 0.0], "r2": [-R2, 0.0, 0.0], "tof": t_h, "revs": 0,
         "normal": [0.0, 0.6, 0.8]},
        {"id": "tiny_angle_inclined", "r1": circ(7200, 30, 51.6, 40), "r2": circ(7300, 31, 51.6, 40), "tof": 160.0,
         "revs": 0, "direction": "prograde"},
        {"id": "one_rev_low", "r1": circ(R1, 0), "r2": circ(12000, 70), "tof": 3.2 * period(9500, MU_EARTH),
         "revs": 1, "branch": "low", "direction": "prograde"},
        {"id": "one_rev_high", "r1": circ(R1, 0), "r2": circ(12000, 70), "tof": 3.2 * period(9500, MU_EARTH),
         "revs": 1, "branch": "high", "direction": "prograde"},
        {"id": "two_rev_low_retro", "r1": circ(R1, 20), "r2": circ(10000, 205), "tof": 4.1 * period(8500, MU_EARTH),
         "revs": 2, "branch": "low", "direction": "retrograde"},
        {"id": "two_rev_high_retro", "r1": circ(R1, 20), "r2": circ(10000, 205), "tof": 4.1 * period(8500, MU_EARTH),
         "revs": 2, "branch": "high", "direction": "retrograde"},
        {"id": "three_rev_too_short", "r1": circ(R1, 0), "r2": circ(12000, 70), "tof": 2.5 * period(9500, MU_EARTH),
         "revs": 3, "branch": "low", "direction": "prograde"},
        {"id": "earth_mars_heliocentric", "r1": circ(1.0 * AU, 0), "r2": circ(1.524 * AU, 150, 1.85, 49.6),
         "tof": 210 * 86400.0, "revs": 0, "direction": "prograde", "mu": MU_SUN},
    ]
    return P


def main():
    out = HERE / "tests" / "cases"
    out.mkdir(parents=True, exist_ok=True)
    for f in ("problems.json", "expected.json"):  # the verifier reads these at import
        if not (out / f).exists():
            (out / f).write_text('{"mu": 1, "problems": []}' if f == "problems.json" else "{}")
    V = verifier()
    env = HERE / "environment" / "data"
    env.mkdir(parents=True, exist_ok=True)

    problems = hidden()

    # near-parabolic: time of flight of the parabola through the two points (Euler/Lambert equation)
    p = next(p for p in problems if p["id"] == "near_parabolic")
    r1, r2 = lambert.norm(p["r1"]), lambert.norm(p["r2"])
    c = lambert.norm(lambert.sub(p["r2"], p["r1"]))
    s = (r1 + r2 + c) / 2
    t_par = math.sqrt(2) / 3 * math.sqrt(s ** 3 / MU_EARTH) * (1 - ((s - c) / s) ** 1.5)
    p["tof"] = t_par * (1 + 5e-4)

    data = {"mu": MU_EARTH, "problems": problems}
    expected = {}
    for p in problems:
        mu = p.get("mu", MU_EARTH)
        sols = lambert.lambert(p["r1"], p["r2"], p["tof"], mu, p.get("revs", 0),
                               p.get("direction", "prograde") == "prograde", p.get("normal"))
        ref = lambert.solve_problem(p, mu)
        if ref is None:
            expected[p["id"]] = {"feasible": False}
            print(f"{p['id']:26s} no solution")
            continue
        a = 1 / (2 / lambert.norm(p["r1"]) - sum(x * x for x in ref["v1"]) / mu)
        expected[p["id"]] = {"feasible": True, "a_km": a}
        r_end, _ = V.propagate(p["r1"], ref["v1"], p["tof"], mu)
        err = lambert.norm(lambert.sub(r_end, p["r2"])) / lambert.norm(p["r2"])
        revs = V.full_revolutions(p["r1"], ref["v1"], p["tof"], mu)
        kind = "ellipse" if a > 0 else "hyperbola"
        other = [round(s[0], 3) for s in sols]
        print(f"{p['id']:26s} a={a:12.3f} km ({kind}) revs={revs} pos_err={err:.1e} branches_a={other}")
        assert err < 1e-9 and revs == p.get("revs", 0), p["id"]
    (out / "problems.json").write_text(json.dumps(data, indent=1))
    (out / "expected.json").write_text(json.dumps(expected, indent=1))

    # Visible examples: a different geometry, with reference answers so the agent can self-check.
    ex = {"mu": MU_EARTH, "problems": [
        {"id": "ex_zero_rev", "r1": circ(6800, 0), "r2": circ(26000, 140), "tof": 15000.0, "revs": 0, "direction": "prograde"},
        {"id": "ex_one_rev_low", "r1": circ(6800, 0), "r2": circ(9000, 45), "tof": 2.6 * period(8000, MU_EARTH),
         "revs": 1, "branch": "low", "direction": "prograde"},
        {"id": "ex_one_rev_high", "r1": circ(6800, 0), "r2": circ(9000, 45), "tof": 2.6 * period(8000, MU_EARTH),
         "revs": 1, "branch": "high", "direction": "prograde"},
    ]}
    (env / "example_problems.json").write_text(json.dumps(ex, indent=1))
    ex_sol = {p["id"]: lambert.solve_problem(p, MU_EARTH) for p in ex["problems"]}
    (env / "example_solutions.json").write_text(json.dumps({"solutions": ex_sol}, indent=1))
    print("examples:", {k: v is not None for k, v in ex_sol.items()})


if __name__ == "__main__":
    main()
