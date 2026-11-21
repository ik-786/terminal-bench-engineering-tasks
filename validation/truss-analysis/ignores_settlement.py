"""WRONG ON PURPOSE (ignores_settlement): fixes pins at zero even when a settlement is given.
Expected to fail the verifier.

Usage: python truss.py <truss.json> <result.json>
"""

import json
import math
import sys

PIVOT_TOL = 1e-10  # relative to the largest diagonal stiffness


def member_geometry(nodes, m):
    xi, yi = nodes[m["i"]]
    xj, yj = nodes[m["j"]]
    L = math.hypot(xj - xi, yj - yi)
    return L, (xj - xi) / L, (yj - yi) / L


def free_elongation(m, L):
    return m.get("alpha", 0.0) * m.get("dT", 0.0) * L + m.get("lack_of_fit", 0.0)


def solve(data):
    nodes = {n["id"]: (n["x"], n["y"]) for n in data["nodes"]}
    order = [n["id"] for n in data["nodes"]]
    dof = {nid: 2 * k for k, nid in enumerate(order)}
    ndof = 2 * len(order)

    # Each node gets a local frame; rollers use (tangent, normal) so the normal DOF can be fixed.
    frame = {nid: ((1.0, 0.0), (0.0, 1.0)) for nid in order}
    fixed = {}  # global dof index (in local frame) -> prescribed value
    supports = {s["node"]: s for s in data["supports"]}
    for nid, s in supports.items():
        d = dof[nid]
        if s["type"] == "pin":
            sx, sy = 0.0, 0.0
            fixed[d], fixed[d + 1] = sx, sy
        elif s["type"] == "roller":
            a = math.radians(s["angle_deg"])
            n = (math.cos(a), math.sin(a))
            t = (-n[1], n[0])
            frame[nid] = (t, n)
            fixed[d + 1] = 0.0
        else:
            raise ValueError(f"unknown support type {s['type']}")

    # Global stiffness and fixed-end forces in the global frame.
    K = [[0.0] * ndof for _ in range(ndof)]
    F0 = [0.0] * ndof
    for m in data["members"]:
        L, c, s = member_geometry(nodes, m)
        k = m["E"] * m["A"] / L
        d = (c, s)
        di, dj = dof[m["i"]], dof[m["j"]]
        idx = [di, di + 1, dj, dj + 1]
        vec = [-c, -s, c, s]
        for a in range(4):
            for b in range(4):
                K[idx[a]][idx[b]] += k * vec[a] * vec[b]
        p0 = k * free_elongation(m, L)
        for a in range(4):
            F0[idx[a]] += p0 * vec[a]

    P = [0.0] * ndof
    for ld in data.get("loads", []):
        d = dof[ld["node"]]
        P[d] += ld.get("fx", 0.0)
        P[d + 1] += ld.get("fy", 0.0)

    # Rotate roller nodes into their (tangent, normal) frame: KL = T^T K T, FL = T^T F,
    # where T is block diagonal with one 2x2 rotation per node (identity elsewhere).
    KL = [row[:] for row in K]
    FL = [P[i] + F0[i] for i in range(ndof)]
    for nid in order:
        (tx, ty), (nx, ny) = frame[nid]
        if (tx, ty) == (1.0, 0.0):
            continue
        d = dof[nid]
        for row in KL:  # columns: K T
            a, b = row[d], row[d + 1]
            row[d], row[d + 1] = a * tx + b * ty, a * nx + b * ny
        ra, rb = KL[d], KL[d + 1]  # rows: T^T (K T)
        KL[d] = [a * tx + b * ty for a, b in zip(ra, rb)]
        KL[d + 1] = [a * nx + b * ny for a, b in zip(ra, rb)]
        fa, fb = FL[d], FL[d + 1]
        FL[d], FL[d + 1] = fa * tx + fb * ty, fa * nx + fb * ny

    free = [i for i in range(ndof) if i not in fixed]
    pres = sorted(fixed)
    up = [fixed[i] for i in pres]
    A = [[KL[i][j] for j in free] for i in free]
    b = [FL[i] - sum(KL[i][j] * u for j, u in zip(pres, up)) for i in free]

    # Gaussian elimination with partial pivoting; a tiny pivot means a mechanism.
    n = len(free)
    scale = max((abs(KL[i][i]) for i in range(ndof)), default=0.0)
    M = [A[i][:] + [b[i]] for i in range(n)]
    for col in range(n):
        p = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[p][col]) <= PIVOT_TOL * scale:
            return {"stable": False}
        M[col], M[p] = M[p], M[col]
        piv = M[col][col]
        for r in range(col + 1, n):
            f = M[r][col] / piv
            if f:
                rr, rc = M[r], M[col]
                for j in range(col, n + 1):
                    rr[j] -= f * rc[j]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (M[r][n] - sum(M[r][j] * x[j] for j in range(r + 1, n))) / M[r][r]

    uL = [0.0] * ndof
    for i, v in zip(free, x):
        uL[i] = v
    for i, v in zip(pres, up):
        uL[i] = v
    u = uL[:]
    for nid in order:
        (tx, ty), (nx, ny) = frame[nid]
        d = dof[nid]
        u[d], u[d + 1] = tx * uL[d] + nx * uL[d + 1], ty * uL[d] + ny * uL[d + 1]

    forces = {}
    for m in data["members"]:
        L, c, s = member_geometry(nodes, m)
        di, dj = dof[m["i"]], dof[m["j"]]
        e = c * (u[dj] - u[di]) + s * (u[dj + 1] - u[di + 1])
        forces[m["id"]] = m["E"] * m["A"] / L * (e - free_elongation(m, L))

    Ku = [sum(K[i][j] * u[j] for j in range(ndof) if K[i][j]) for i in range(ndof)]
    reactions = {}
    for nid in supports:
        d = dof[nid]
        reactions[nid] = [Ku[d] - P[d] - F0[d], Ku[d + 1] - P[d + 1] - F0[d + 1]]

    return {"stable": True,
            "displacements": {nid: [u[dof[nid]], u[dof[nid] + 1]] for nid in order},
            "member_forces": forces,
            "reactions": reactions}


def main():
    with open(sys.argv[1]) as fh:
        data = json.load(fh)
    with open(sys.argv[2], "w") as fh:
        json.dump(solve(data), fh, indent=1)


if __name__ == "__main__":
    main()
