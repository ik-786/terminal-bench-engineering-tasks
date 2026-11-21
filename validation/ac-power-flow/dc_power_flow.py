"""WRONG ON PURPOSE (dc_power_flow): lossless DC approximation, V = 1 everywhere.
Expected to fail the verifier."""
import json, math, sys

def solve(A, b):
    n = len(b); M = [r[:] + [b[i]] for i, r in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c])); M[c], M[p] = M[p], M[c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            for j in range(c, n + 1): M[r][j] -= f * M[c][j]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (M[r][n] - sum(M[r][j] * x[j] for j in range(r + 1, n))) / M[r][r]
    return x

case = json.load(open(sys.argv[1])); base = case["base_mva"]
buses = case["buses"]; idx = {b["id"]: i for i, b in enumerate(buses)}; n = len(buses)
B = [[0.0] * n for _ in range(n)]
for br in case["branches"]:
    if br.get("status", 1) == 0: continue
    f, t = idx[br["from"]], idx[br["to"]]; y = 1 / br["x"]
    B[f][f] += y; B[t][t] += y; B[f][t] -= y; B[t][f] -= y
s = next(i for i, b in enumerate(buses) if b["type"] == "slack")
keep = [i for i in range(n) if i != s]
P = [(buses[i].get("pg", 0.0) - buses[i].get("pd", 0.0)) / base for i in keep]
th = solve([[B[i][j] for j in keep] for i in keep], P)
theta = [0.0] * n
for k, i in enumerate(keep): theta[i] = th[k]
res = {"converged": True, "q_limited": [], "losses_mw": 0.0,
       "buses": [{"id": b["id"], "v": b.get("v", 1.0) if b["type"] != "pq" else 1.0,
                  "theta_deg": math.degrees(theta[i]) + buses[s].get("theta_deg", 0.0),
                  "p_inj_mw": 0.0, "q_inj_mvar": 0.0} for i, b in enumerate(buses)],
       "branches": [{"from": br["from"], "to": br["to"], "p_from_mw": 0.0, "q_from_mvar": 0.0,
                     "p_to_mw": 0.0, "q_to_mvar": 0.0} for br in case["branches"] if br.get("status", 1)]}
json.dump(res, open(sys.argv[2], "w"))
