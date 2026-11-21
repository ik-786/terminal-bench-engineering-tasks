"""Builds the example truss (visible to the agent) and the hidden test trusses.

Run from the task folder:  python dev/make_cases.py
"""

import json
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE / "solution"))
import truss  # noqa: E402  (reference solver)

STEEL = {"E": 200e9, "alpha": 1.2e-5}


def member(mid, i, j, A, mat=STEEL, **extra):
    return {"id": mid, "i": i, "j": j, "E": mat["E"], "A": A, "alpha": mat.get("alpha", 0.0), **extra}


def pratt(panels, w, h, x_brace=False, drop_diag=None):
    nodes = [{"id": f"L{k}", "x": k * w, "y": 0.0} for k in range(panels + 1)]
    nodes += [{"id": f"U{k}", "x": k * w, "y": h} for k in range(1, panels)]
    mem, n = [], 0

    def add(i, j, A):
        nonlocal n
        n += 1
        mem.append(member(f"m{n}", i, j, A))

    for k in range(panels):
        add(f"L{k}", f"L{k + 1}", 0.003)
    for k in range(1, panels - 1):
        add(f"U{k}", f"U{k + 1}", 0.004)
    add("L0", "U1", 0.004)
    add(f"U{panels - 1}", f"L{panels}", 0.004)
    for k in range(1, panels):
        add(f"L{k}", f"U{k}", 0.002)
    for k in range(1, panels - 1):
        if drop_diag == k:
            continue
        if k < panels / 2:
            add(f"U{k}", f"L{k + 1}", 0.002)
            if x_brace:
                add(f"L{k}", f"U{k + 1}", 0.002)
        else:
            add(f"L{k}", f"U{k + 1}", 0.002)
            if x_brace:
                add(f"U{k}", f"L{k + 1}", 0.002)
    return nodes, mem


def tower(levels, width, h, seed):
    rng = random.Random(seed)
    rubber = {"E": 2e5, "alpha": 0.0}
    nodes, mem = [], []
    for k in range(levels + 1):
        taper = 1 - 0.4 * k / levels
        nodes.append({"id": f"A{k}", "x": -width * taper / 2, "y": k * h})
        nodes.append({"id": f"B{k}", "x": width * taper / 2, "y": k * h})
    n = 0
    for k in range(levels):
        for i, j in ((f"A{k}", f"A{k + 1}"), (f"B{k}", f"B{k + 1}"), (f"A{k + 1}", f"B{k + 1}"),
                     (f"A{k}", f"B{k + 1}"), (f"B{k}", f"A{k + 1}")):
            n += 1
            soft = rng.random() < 0.15 and i[0] != j[0]  # some bracing members are very soft
            mem.append(member(f"t{n}", i, j, rng.choice([0.002, 0.004, 0.008]) if not soft else 0.01,
                              mat=rubber if soft else STEEL))
    loads = [{"node": f"B{k}", "fx": rng.uniform(2e3, 8e3), "fy": -rng.uniform(1e4, 3e4)} for k in range(2, levels + 1, 3)]
    supports = [{"node": "A0", "type": "pin"}, {"node": "B0", "type": "pin"}]
    return {"nodes": nodes, "members": mem, "supports": supports, "loads": loads}


def bottom_loads(panels, P):
    return [{"node": f"L{k}", "fx": 0.0, "fy": -P} for k in range(1, panels)]


def main():
    env = HERE / "environment" / "data"
    out = HERE / "tests" / "cases"
    env.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)

    nodes, mem = pratt(4, 3.0, 2.5)
    example = {"nodes": nodes, "members": mem, "loads": bottom_loads(4, 20e3),
               "supports": [{"node": "L0", "type": "pin"}, {"node": "L4", "type": "roller", "angle_deg": 90}]}
    (env / "example_truss.json").write_text(json.dumps(example, indent=1))

    cases = {}
    nodes, mem = pratt(6, 3.0, 2.5)
    cases["pratt_bridge"] = ({"nodes": nodes, "members": mem, "loads": bottom_loads(6, 40e3),
                              "supports": [{"node": "L0", "type": "pin"},
                                           {"node": "L6", "type": "roller", "angle_deg": 90}]}, True)

    nodes, mem = pratt(6, 3.0, 2.5, x_brace=True)
    for m in mem:
        if m["i"].startswith("L") and m["j"].startswith("L"):
            m["dT"] = 40.0  # bottom chord heated by 40 K
    cases["thermal_indeterminate"] = ({"nodes": nodes, "members": mem, "loads": [],
                                       "supports": [{"node": "L0", "type": "pin"}, {"node": "L6", "type": "pin"}]}, True)

    nodes, mem = pratt(5, 4.0, 3.0)
    cases["inclined_roller_settlement"] = ({"nodes": nodes, "members": mem,
                                            "loads": bottom_loads(5, 25e3) + [{"node": "U2", "fx": 15e3, "fy": 0.0}],
                                            "supports": [{"node": "L0", "type": "pin", "settlement": [0.0, -0.012]},
                                                         {"node": "L5", "type": "roller", "angle_deg": 120}]}, True)

    nodes, mem = pratt(4, 3.0, 3.0, x_brace=True)
    next(m for m in mem if (m["i"], m["j"]) == ("L1", "U2"))["lack_of_fit"] = 0.003
    cases["lack_of_fit"] = ({"nodes": nodes, "members": mem, "loads": [{"node": "U2", "fx": 0.0, "fy": -10e3}],
                             "supports": [{"node": "L0", "type": "pin"},
                                          {"node": "L4", "type": "roller", "angle_deg": 90}]}, True)

    nodes, mem = pratt(6, 3.0, 2.5, drop_diag=2)
    cases["missing_diagonal"] = ({"nodes": nodes, "members": mem, "loads": bottom_loads(6, 40e3),
                                  "supports": [{"node": "L0", "type": "pin"},
                                               {"node": "L6", "type": "roller", "angle_deg": 90}]}, False)

    nodes, mem = pratt(4, 3.0, 2.5)
    nodes.append({"id": "C", "x": 4.5, "y": 2.5})  # midpoint of the top chord U1-U2
    mem = [m for m in mem if (m["i"], m["j"]) != ("U1", "U2")]
    mem += [member("c1", "U1", "C", 0.004), member("c2", "C", "U2", 0.004)]
    cases["collinear_node"] = ({"nodes": nodes, "members": mem,
                                "loads": bottom_loads(4, 20e3) + [{"node": "C", "fx": 0.0, "fy": -5e3}],
                                "supports": [{"node": "L0", "type": "pin"},
                                             {"node": "L4", "type": "roller", "angle_deg": 90}]}, False)

    nodes, mem = pratt(4, 3.0, 2.5)
    cases["parallel_rollers"] = ({"nodes": nodes, "members": mem, "loads": bottom_loads(4, 20e3),
                                  "supports": [{"node": "L0", "type": "roller", "angle_deg": 90},
                                               {"node": "L4", "type": "roller", "angle_deg": 90}]}, False)

    cases["stiff_soft_tower"] = (tower(40, 6.0, 2.0, seed=3), True)

    meta = {}
    for name, (case, stable) in cases.items():
        (out / f"{name}.json").write_text(json.dumps(case, indent=1))
        r = truss.solve(json.loads(json.dumps(case)))
        assert r["stable"] == stable, name
        meta[name] = {"stable": stable, "n_nodes": len(case["nodes"]), "n_members": len(case["members"])}
        extra = ""
        if stable:
            f = r["member_forces"]
            extra = f"max |N| = {max(abs(v) for v in f.values()) / 1e3:.1f} kN"
        print(f"{name:28s} nodes={len(case['nodes']):3d} members={len(case['members']):3d} stable={stable} {extra}")
    (out / "expected.json").write_text(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
