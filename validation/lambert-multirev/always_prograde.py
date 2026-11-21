"""WRONG ON PURPOSE (always_prograde): ignores the direction field.
Expected to fail the verifier.

Lancaster-Blanchard time-of-flight equation in Izzo's form (x variable), solved
by bracketing rather than iteration from a guess, so both multi-revolution
branches and the infeasible cases are found reliably.

Usage: python lambert.py <problems.json> <solutions.json>
"""

import json
import math
import sys


def sub(a, b): return [a[i] - b[i] for i in range(3)]
def norm(a): return math.sqrt(sum(x * x for x in a))
def scale(a, k): return [x * k for x in a]
def unit(a): return scale(a, 1 / norm(a))
def cross(a, b): return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
def add(*vs): return [sum(v[i] for v in vs) for i in range(3)]


def hyp2f1b(z, tol=1e-15):
    """Battin's hypergeometric series 2F1(3, 1, 5/2, z)."""
    s, c, j = 1.0, 1.0, 0
    while True:
        c = c * (3.0 + j) * (1.0 + j) / (2.5 + j) * z / (j + 1)
        s += c
        j += 1
        if abs(c) < tol or j > 10000:
            return s


def tof_lagrange(x, N, lam):
    a = 1.0 / (1.0 - x * x)
    if a > 0:
        alfa = 2.0 * math.acos(x)
        beta = 2.0 * math.asin(math.sqrt(lam * lam / a))
        if lam < 0:
            beta = -beta
        return a * math.sqrt(a) * ((alfa - math.sin(alfa)) - (beta - math.sin(beta)) + 2.0 * math.pi * N) / 2.0
    alfa = 2.0 * math.acosh(x)
    beta = 2.0 * math.asinh(math.sqrt(-lam * lam / a))
    if lam < 0:
        beta = -beta
    return -a * math.sqrt(-a) * ((beta - math.sinh(beta)) - (alfa - math.sinh(alfa))) / 2.0


def tof(x, N, lam):
    """Non-dimensional time of flight T(x) for N full revolutions."""
    dist = abs(x - 1.0)
    if 0.01 < dist < 0.2:
        return tof_lagrange(x, N, lam)
    K, E = lam * lam, x * x - 1.0
    rho = abs(E)
    z = math.sqrt(1 + K * E)
    if dist <= 0.01:
        eta = z - lam * x
        S1 = 0.5 * (1.0 - lam - x * eta)
        Q = 4.0 / 3.0 * hyp2f1b(S1)
        return (eta ** 3 * Q + 4.0 * lam * eta) / 2.0 + N * math.pi / rho ** 1.5
    y = math.sqrt(rho)
    g = x * z - lam * E
    if E < 0:
        d = N * math.pi + math.acos(max(-1.0, min(1.0, g)))
    else:
        d = math.log(y * (z - lam * x) + g)
    return (x - lam * z - d / y) / E


def bisect(f, lo, hi, iters=200):
    flo = f(lo)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
        if hi - lo < 1e-15 * max(1.0, abs(mid)):
            break
    return 0.5 * (lo + hi)


def golden_min(f, lo, hi, iters=300):
    g = (math.sqrt(5) - 1) / 2
    a, b = lo, hi
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(iters):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = f(d)
        if b - a < 1e-14:
            break
    return 0.5 * (a + b)


def solve_x(T, N, lam):
    """All roots x of T(x) = T. N = 0: one root on (-1, inf). N > 0: zero or two roots on (-1, 1)."""
    eps = 1e-12
    if N == 0:
        hi = 1.5
        while tof(hi, 0, lam) > T:
            hi = 1 + 2 * (hi - 1) + 1.0
        return [bisect(lambda x: tof(x, 0, lam) - T, -1 + eps, hi)]
    xm = golden_min(lambda x: tof(x, N, lam), -1 + eps, 1 - eps)
    if tof(xm, N, lam) > T:
        return []
    left = bisect(lambda x: tof(x, N, lam) - T, -1 + eps, xm)
    right = bisect(lambda x: tof(x, N, lam) - T, xm, 1 - eps)
    return [left, right]


def lambert(r1, r2, t, mu, N=0, prograde=True, normal=None):
    c = norm(sub(r2, r1))
    R1, R2 = norm(r1), norm(r2)
    s = (R1 + R2 + c) / 2
    ir1, ir2 = unit(r1), unit(r2)
    h = cross(r1, r2)
    collinear = norm(h) < 1e-10 * R1 * R2
    lam = math.sqrt(max(0.0, 1 - c / s))
    if collinear:
        ih = unit(normal)
        it1, it2 = cross(ih, ir1), cross(ih, ir2)
    else:
        ih = unit(h)
        if ih[2] < 0:
            lam = -lam
            it1, it2 = cross(ir1, ih), cross(ir2, ih)
        else:
            it1, it2 = cross(ih, ir1), cross(ih, ir2)
        if not prograde:
            lam = -lam
            it1, it2 = scale(it1, -1), scale(it2, -1)
    T = math.sqrt(2 * mu / s ** 3) * t
    out = []
    for x in solve_x(T, N, lam):
        y = math.sqrt(1 - lam * lam * (1 - x * x))
        gamma = math.sqrt(mu * s / 2)
        rho = (R1 - R2) / c
        sigma = math.sqrt(max(0.0, 1 - rho * rho))
        vr1 = gamma * ((lam * y - x) - rho * (lam * y + x)) / R1
        vr2 = -gamma * ((lam * y - x) + rho * (lam * y + x)) / R2
        vt = gamma * sigma * (y + lam * x)
        v1 = add(scale(ir1, vr1), scale(it1, vt / R1))
        v2 = add(scale(ir2, vr2), scale(it2, vt / R2))
        a = s / 2 / (1 - x * x)
        out.append((a, v1, v2))
    return out


def solve_problem(p, mu):
    N = p.get("revs", 0)
    sols = lambert(p["r1"], p["r2"], p["tof"], mu, N, True, p.get("normal"))
    if not sols:
        return None
    if N > 0:
        sols.sort(key=lambda s: s[0])
        a, v1, v2 = sols[0] if p["branch"] == "low" else sols[-1]
    else:
        a, v1, v2 = sols[0]
    return {"v1": v1, "v2": v2}


def main():
    with open(sys.argv[1]) as fh:
        data = json.load(fh)
    res = {p["id"]: solve_problem(p, p.get("mu", data["mu"])) for p in data["problems"]}
    with open(sys.argv[2], "w") as fh:
        json.dump({"solutions": res}, fh, indent=1)


if __name__ == "__main__":
    main()
