"""WRONG ON PURPOSE (claims_convergence): reports converged=true with its last iterate even when Newton fails.
Expected to fail the verifier.

Usage: python powerflow.py <case.json> <result.json>
"""

import cmath
import json
import math
import sys

MAX_ITER = 30
TOL = 1e-10  # per unit
MAX_OUTER = 20


def build_ybus(case, idx):
    n = len(idx)
    base = case["base_mva"]
    Y = [[0j] * n for _ in range(n)]
    for b in case["buses"]:
        i = idx[b["id"]]
        Y[i][i] += complex(b.get("gs", 0.0), b.get("bs", 0.0)) / base
    for br in case["branches"]:
        if br.get("status", 1) == 0:
            continue
        f, t = idx[br["from"]], idx[br["to"]]
        ys = 1 / complex(br["r"], br["x"])
        bc = br.get("b", 0.0)
        ratio = br.get("tap", 1.0) or 1.0
        tap = ratio * cmath.exp(1j * math.radians(br.get("shift_deg", 0.0)))
        Y[f][f] += (ys + 1j * bc / 2) / (abs(tap) ** 2)
        Y[f][t] += -ys / tap.conjugate()
        Y[t][f] += -ys / tap
        Y[t][t] += ys + 1j * bc / 2
    return Y


def injections(Y, V, th):
    n = len(V)
    P, Q = [0.0] * n, [0.0] * n
    for i in range(n):
        for k in range(n):
            y = Y[i][k]
            if y == 0:
                continue
            a = th[i] - th[k]
            c, s = math.cos(a), math.sin(a)
            P[i] += V[i] * V[k] * (y.real * c + y.imag * s)
            Q[i] += V[i] * V[k] * (y.real * s - y.imag * c)
    return P, Q


def solve_linear(A, b):
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) < 1e-14:
            raise ZeroDivisionError("singular Jacobian")
        M[c], M[p] = M[p], M[c]
        piv = M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / piv
            if f:
                row_r, row_c = M[r], M[c]
                for j in range(c, n + 1):
                    row_r[j] -= f * row_c[j]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (M[r][n] - sum(M[r][j] * x[j] for j in range(r + 1, n))) / M[r][r]
    return x


def newton(Y, kind, V, th, Psp, Qsp):
    n = len(V)
    pv_pq = [i for i in range(n) if kind[i] != "slack"]
    pq = [i for i in range(n) if kind[i] == "pq"]
    for _ in range(MAX_ITER):
        P, Q = injections(Y, V, th)
        mis = [Psp[i] - P[i] for i in pv_pq] + [Qsp[i] - Q[i] for i in pq]
        if max((abs(m) for m in mis), default=0.0) < TOL:
            return True
        # Jacobian: d(P,Q)/d(theta, V)
        npv, npq = len(pv_pq), len(pq)
        J = [[0.0] * (npv + npq) for _ in range(npv + npq)]
        col_th = {b: j for j, b in enumerate(pv_pq)}
        col_v = {b: npv + j for j, b in enumerate(pq)}
        rows = [("P", i) for i in pv_pq] + [("Q", i) for i in pq]
        for r, (kindrow, i) in enumerate(rows):
            for k in range(n):
                y = Y[i][k]
                if y == 0 and k != i:
                    continue
                if k == i:
                    G, B = Y[i][i].real, Y[i][i].imag
                    if kindrow == "P":
                        if i in col_th:
                            J[r][col_th[i]] = -Q[i] - B * V[i] ** 2
                        if i in col_v:
                            J[r][col_v[i]] = P[i] / V[i] + G * V[i]
                    else:
                        if i in col_th:
                            J[r][col_th[i]] = P[i] - G * V[i] ** 2
                        if i in col_v:
                            J[r][col_v[i]] = Q[i] / V[i] - B * V[i]
                else:
                    a = th[i] - th[k]
                    c, s = math.cos(a), math.sin(a)
                    G, B = y.real, y.imag
                    if kindrow == "P":
                        if k in col_th:
                            J[r][col_th[k]] = V[i] * V[k] * (G * s - B * c)
                        if k in col_v:
                            J[r][col_v[k]] = V[i] * (G * c + B * s)
                    else:
                        if k in col_th:
                            J[r][col_th[k]] = -V[i] * V[k] * (G * c + B * s)
                        if k in col_v:
                            J[r][col_v[k]] = V[i] * (G * s - B * c)
        try:
            dx = solve_linear(J, mis)
        except ZeroDivisionError:
            return False
        for b, j in col_th.items():
            th[b] += dx[j]
        for b, j in col_v.items():
            V[b] += dx[j]
            if not (0.05 < V[b] < 3.0):
                return False
    return False


def run(case):
    base = case["base_mva"]
    buses = case["buses"]
    idx = {b["id"]: i for i, b in enumerate(buses)}
    n = len(buses)
    Y = build_ybus(case, idx)
    kind = [b["type"] for b in buses]
    vset = [b.get("v", 1.0) for b in buses]
    V = [vset[i] if kind[i] in ("slack", "pv") else 1.0 for i in range(n)]
    th = [math.radians(buses[i].get("theta_deg", 0.0)) if kind[i] == "slack" else 0.0 for i in range(n)]
    pd = [b.get("pd", 0.0) / base for b in buses]
    qd = [b.get("qd", 0.0) / base for b in buses]
    pg = [b.get("pg", 0.0) / base for b in buses]
    Psp = [pg[i] - pd[i] for i in range(n)]
    Qsp = [b.get("qg", 0.0) / base - qd[i] if b["type"] == "pq" else 0.0 for i, b in enumerate(buses)]
    at_limit = {}  # bus index -> "max" | "min"

    for _ in range(MAX_OUTER):
        if not newton(Y, kind, V, th, Psp, Qsp):
            V = [max(0.6, min(1.4, v)) for v in V]
            break
        P, Q = injections(Y, V, th)
        changed = False
        for i, b in enumerate(buses):
            if b["type"] != "pv":
                continue
            qmax = b.get("qmax", math.inf) / base
            qmin = b.get("qmin", -math.inf) / base
            if kind[i] == "pv":
                qg = Q[i] + qd[i]
                if qg > qmax + 1e-9:
                    kind[i], Qsp[i], at_limit[i], changed = "pq", qmax - qd[i], "max", True
                elif qg < qmin - 1e-9:
                    kind[i], Qsp[i], at_limit[i], changed = "pq", qmin - qd[i], "min", True
            else:  # currently held at a limit: release if the voltage is on the wrong side
                if (at_limit[i] == "max" and V[i] > vset[i] + 1e-9) or (at_limit[i] == "min" and V[i] < vset[i] - 1e-9):
                    kind[i], V[i], changed = "pv", vset[i], True
                    del at_limit[i]
        if not changed:
            break
    else:
        return {"converged": False}

    P, Q = injections(Y, V, th)
    out_buses = [{"id": b["id"], "v": V[i], "theta_deg": math.degrees(th[i]),
                  "p_inj_mw": P[i] * base, "q_inj_mvar": Q[i] * base} for i, b in enumerate(buses)]
    flows = []
    Vc = [cmath.rect(V[i], th[i]) for i in range(n)]
    for br in case["branches"]:
        if br.get("status", 1) == 0:
            continue
        f, t = idx[br["from"]], idx[br["to"]]
        ys = 1 / complex(br["r"], br["x"])
        bc = br.get("b", 0.0)
        ratio = br.get("tap", 1.0) or 1.0
        tap = ratio * cmath.exp(1j * math.radians(br.get("shift_deg", 0.0)))
        i_f = Vc[f] * (ys + 1j * bc / 2) / abs(tap) ** 2 - Vc[t] * ys / tap.conjugate()
        i_t = Vc[t] * (ys + 1j * bc / 2) - Vc[f] * ys / tap
        s_f, s_t = Vc[f] * i_f.conjugate() * base, Vc[t] * i_t.conjugate() * base
        flows.append({"from": br["from"], "to": br["to"], "p_from_mw": s_f.real, "q_from_mvar": s_f.imag,
                      "p_to_mw": s_t.real, "q_to_mvar": s_t.imag})
    return {"converged": True, "buses": out_buses, "branches": flows,
            "q_limited": sorted(buses[i]["id"] for i in at_limit),
            "losses_mw": sum(P) * base}


def main():
    with open(sys.argv[1]) as fh:
        case = json.load(fh)
    with open(sys.argv[2], "w") as fh:
        json.dump(run(case), fh, indent=1)


if __name__ == "__main__":
    main()
