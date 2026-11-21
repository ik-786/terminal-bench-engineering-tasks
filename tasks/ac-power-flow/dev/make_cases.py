"""Generates the example case (visible to the agent) and the hidden test cases.

Run from the task folder:  python dev/make_cases.py

Networks are random but seeded, so the files are reproducible. Each hidden case
is checked against the reference solver for the property it is meant to test
(a reactive limit binds, a limit is released, the case is beyond collapse...).
"""

import copy
import json
import math
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE / "solution"))
import powerflow  # noqa: E402  (reference solver)


def network(seed, n, n_pv, extra_edges, n_xfmr=1, shifters=0, shunts=0, outages=0,
            load_scale=1.0, q_cap=None, ids=None, slack_angle=0.0):
    rng = random.Random(seed)
    ids = ids or list(range(1, n + 1))
    buses = []
    for k in range(n):
        buses.append({"id": ids[k], "type": "pq", "pd": round(rng.uniform(5, 40) * load_scale, 2),
                      "qd": round(rng.uniform(1, 15) * load_scale, 2)})
    buses[0].update({"type": "slack", "v": 1.05, "theta_deg": slack_angle, "pd": 0.0, "qd": 0.0})
    for k in rng.sample(range(1, n), n_pv):
        qmax = q_cap if q_cap is not None else round(rng.uniform(25, 60), 1)
        buses[k].update({"type": "pv", "v": round(rng.uniform(1.00, 1.045), 3), "pg": round(rng.uniform(30, 80), 1),
                         "qmin": round(-rng.uniform(10, 30), 1), "qmax": qmax})
    for k in rng.sample(range(1, n), shunts):
        buses[k]["bs"] = round(rng.choice([-1, 1]) * rng.uniform(5, 20), 1)
        buses[k]["gs"] = round(rng.uniform(0, 2), 2)

    branches, edges = [], set()

    def add(a, b):
        key = (min(a, b), max(a, b))
        if a == b or key in edges:
            return
        edges.add(key)
        x = round(rng.uniform(0.04, 0.25), 4)
        branches.append({"from": ids[a], "to": ids[b], "r": round(x * rng.uniform(0.15, 0.4), 4), "x": x,
                         "b": round(rng.uniform(0.0, 0.06), 4), "tap": 1.0, "shift_deg": 0.0, "status": 1})

    order = list(range(n))
    rng.shuffle(order[1:])
    for k in range(1, n):  # random spanning tree keeps the network connected
        add(order[k], order[rng.randrange(k)])
    tree = len(branches)
    while len(branches) < tree + extra_edges:
        add(rng.randrange(n), rng.randrange(n))
    for br in rng.sample(branches, n_xfmr + shifters):
        br.update({"tap": round(rng.uniform(0.93, 1.06), 4), "b": 0.0})
    for br in rng.sample([b for b in branches if b["tap"] != 1.0], shifters):
        br["shift_deg"] = round(rng.uniform(-8, 8), 2)
    # outages only on meshed (non-tree) branches so the network stays connected
    for br in rng.sample(branches[tree:], outages):
        br["status"] = 0
    return {"base_mva": 100.0, "buses": buses, "branches": branches}


def solve(case):
    return powerflow.run(copy.deepcopy(case))


def scaled(case, lam):
    c = copy.deepcopy(case)
    for b in c["buses"]:
        b["pd"] = b.get("pd", 0.0) * lam
        b["qd"] = b.get("qd", 0.0) * lam
        if b["type"] == "pv":
            b["pg"] = b["pg"] * lam
    return c


def collapse_factor(case):
    """Largest load scale with a solution, by bisection on reference convergence."""
    lo, hi = 1.0, 1.0
    while solve(scaled(case, hi))["converged"]:
        lo, hi = hi, hi * 1.5
    for _ in range(30):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if solve(scaled(case, mid))["converged"] else (lo, mid)
    return lo


def main():
    out = HERE / "tests" / "cases"
    out.mkdir(parents=True, exist_ok=True)
    env = HERE / "environment" / "data"
    env.mkdir(parents=True, exist_ok=True)

    example = network(1, 6, 2, 2, n_xfmr=1)
    assert solve(example)["converged"]
    (env / "example_case.json").write_text(json.dumps(example, indent=1))

    cases = {}
    # 1. small, transformer tap, slack angle not zero
    cases["tap_and_slack_angle"] = (network(11, 5, 1, 2, n_xfmr=1, slack_angle=7.5), {"converged": True})
    # 2. reactive limits bind at several PV buses
    for s in range(100, 300):
        c = network(s, 14, 4, 6, n_xfmr=2, q_cap=12.0, load_scale=1.3)
        r = solve(c)
        if r["converged"] and len(r["q_limited"]) >= 2:
            cases["reactive_limits"] = (c, {"converged": True, "min_limited": 2})
            break
    # 3. phase shifter, bus shunts, out-of-service branches, scrambled bus ids
    ids = random.Random(5).sample(range(2, 999), 24)
    cases["shifter_shunt_outage_ids"] = (network(21, 24, 5, 10, n_xfmr=2, shifters=1, shunts=3, outages=2, ids=ids),
                                         {"converged": True})
    # 4. a limit has to be released again (one-way switching gives a wrong answer)
    for s in range(300, 2000):
        c = network(s, 12, 4, 5, n_xfmr=2, shunts=2, load_scale=1.1)
        for b in c["buses"]:
            if b["type"] == "pv":
                b["qmax"] = round(random.Random(s).uniform(8, 30), 1)
        r = solve(c)
        if not r["converged"]:
            continue
        # detect release: rerun with one-way switching and compare
        one_way = copy.deepcopy(c)
        res_one = one_way_switching(one_way)
        if res_one is not None and res_one != sorted(r["q_limited"]):
            cases["limit_release"] = (c, {"converged": True})
            break
    # 5. larger network (time budget)
    for s in range(57, 400):
        c = network(s, 57, 9, 24, n_xfmr=5, shifters=1, shunts=4, outages=3, load_scale=0.8)
        r = solve(c)
        if r["converged"] and r["q_limited"]:
            cases["large_57_bus"] = (c, {"converged": True})
            break
    # 6. light load with heavy line charging: PV buses pushed to their lower limit
    for s in range(2000, 2400):
        c = network(s, 18, 4, 10, n_xfmr=2, load_scale=0.15)
        for br in c["branches"]:
            if br["tap"] == 1.0:
                br["b"] = round(br["b"] * 3 + 0.05, 4)
        r = solve(c)
        if (r["converged"] and any(lim_side(c, r, b) == "min" for b in r["q_limited"])
                and all(0.9 <= x["v"] <= 1.15 for x in r["buses"])):
            cases["lower_limit"] = (c, {"converged": True})
            break
    # 7. beyond voltage collapse: no solution exists
    base = network(77, 10, 2, 4, n_xfmr=1)
    lam = collapse_factor(base)
    cases["beyond_collapse"] = (scaled(base, round(lam * 1.6, 3)), {"converged": False})

    meta = {}
    for name, (case, expect) in cases.items():
        (out / f"{name}.json").write_text(json.dumps(case, indent=1))
        r = solve(case)
        assert r["converged"] == expect["converged"], name
        meta[name] = {**expect, "n_bus": len(case["buses"]),
                      "reference_q_limited": r.get("q_limited", [])}
        print(f"{name:28s} buses={len(case['buses']):3d} converged={r['converged']} "
              f"q_limited={r.get('q_limited', '-')}")
    (out / "expected.json").write_text(json.dumps(meta, indent=1))
    print(f"collapse factor of base case for beyond_collapse: {lam:.3f}")


def lim_side(case, res, bus_id):
    b = next(x for x in case["buses"] if x["id"] == bus_id)
    rb = next(x for x in res["buses"] if x["id"] == bus_id)
    qg = rb["q_inj_mvar"] + b.get("qd", 0.0)
    return "max" if abs(qg - b["qmax"]) < abs(qg - b["qmin"]) else "min"


def one_way_switching(case):
    """Q-limit handling without ever releasing a limit, used to find release cases."""
    base = case["base_mva"]
    buses = case["buses"]
    idx = {b["id"]: i for i, b in enumerate(buses)}
    n = len(buses)
    Y = powerflow.build_ybus(case, idx)
    kind = [b["type"] for b in buses]
    vset = [b.get("v", 1.0) for b in buses]
    V = [vset[i] if kind[i] in ("slack", "pv") else 1.0 for i in range(n)]
    th = [0.0] * n
    Psp = [(b.get("pg", 0.0) - b.get("pd", 0.0)) / base for b in buses]
    Qsp = [-b.get("qd", 0.0) / base if b["type"] == "pq" else 0.0 for b in buses]
    limited = []
    for _ in range(20):
        if not powerflow.newton(Y, kind, V, th, Psp, Qsp):
            return None
        P, Q = powerflow.injections(Y, V, th)
        changed = False
        for i, b in enumerate(buses):
            if kind[i] != "pv":
                continue
            qg = Q[i] + b.get("qd", 0.0) / base
            if qg > b["qmax"] / base:
                kind[i], Qsp[i] = "pq", (b["qmax"] - b.get("qd", 0.0)) / base
                limited.append(b["id"]); changed = True
            elif qg < b["qmin"] / base:
                kind[i], Qsp[i] = "pq", (b["qmin"] - b.get("qd", 0.0)) / base
                limited.append(b["id"]); changed = True
        if not changed:
            break
    return sorted(limited)


if __name__ == "__main__":
    main()
