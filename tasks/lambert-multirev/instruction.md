# Multi-revolution Lambert solver

Write a solver for Lambert's problem in the two-body model and save it as `/app/lambert.py`.
It is run as

```
python /app/lambert.py <problems.json> <solutions.json>
```

For each problem, find the velocity `v1` at position `r1` such that a spacecraft moving under
gravity alone reaches position `r2` exactly `tof` seconds later, and the velocity `v2` it arrives
with. Example problems and their answers are in `/app/data/example_problems.json` and
`/app/data/example_solutions.json`. Your solver will be tested on other problems that you do not
have.

## Input

```json
{
  "mu": 398600.4418,
  "problems": [
    {"id": "p1", "r1": [7000, 0, 0], "r2": [-21000, 36000, 0], "tof": 12000.0,
     "revs": 0, "direction": "prograde"},
    {"id": "p2", "r1": [7000, 0, 0], "r2": [4100, 11200, 0], "tof": 28000.0,
     "revs": 1, "branch": "low", "direction": "retrograde"},
    {"id": "p3", "r1": [7000, 0, 0], "r2": [-42164, 0, 0], "tof": 19000.0,
     "revs": 0, "normal": [0.0, 0.6, 0.8]},
    {"id": "p4", "r1": [1.5e8, 0, 0], "r2": [-1.9e8, 1.1e8, 7e6], "tof": 1.8e7,
     "revs": 0, "direction": "prograde", "mu": 1.32712440018e11}
  ]
}
```

- Units: km, s, km/s, km^3/s^2. `mu` is the gravitational parameter. A problem may give its own `mu`,
  which then replaces the top-level one for that problem.
- `direction`: `prograde` means the transfer's angular momentum `r1 x v1` has a positive z component,
  `retrograde` means negative. The transfer angle is measured in that direction of motion, so it can
  be more than 180 degrees.
- `revs`: the number of complete revolutions made before arriving, 0 or more. The spacecraft sweeps
  the transfer angle plus `revs` full turns.
- `branch`: when `revs` is 1 or more there are generally two transfers. `low` asks for the one with
  the smaller semi-major axis, `high` for the one with the larger.
- If `r1` and `r2` point in exactly opposite directions, the plane of motion is not defined by them.
  Such problems give a `normal` vector instead of `direction`, and the transfer's angular momentum must
  point along that normal.
- The time of flight can be short enough that the transfer is hyperbolic. It can also be close to the
  parabolic transfer time, or long enough for several revolutions.

## Output

```json
{
  "solutions": {
    "p1": {"v1": [vx, vy, vz], "v2": [vx, vy, vz]},
    "p2": null
  }
}
```

Every problem id must appear once. If no transfer exists with the requested number of revolutions
(the time of flight is shorter than the minimum for that many revolutions), the answer is `null`.

The verifier flies your `v1` from `r1` for `tof` seconds with its own Kepler propagator. The result
must land on `r2` within 1e-6 of `|r2|`, and the velocity there must match your `v2` to 1e-6 relative.
It also checks the direction of motion, the number of complete revolutions and, for multi-revolution
problems, that you picked the requested branch.

The whole file must be solved within 120 seconds. Python 3.11 with NumPy is available; no network
access.
