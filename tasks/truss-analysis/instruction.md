# Linear-elastic analysis of plane trusses

Write a solver for pin-jointed plane trusses and save it as `/app/truss.py`. It is run as

```
python /app/truss.py <truss.json> <result.json>
```

It reads a truss, analyses it with small-displacement linear elasticity, writes `<result.json>` and exits
with status 0. An example is in `/app/data/example_truss.json`. Your solver will be tested on other
trusses that you do not have, up to about 100 nodes and 250 members.

## Input

SI units throughout (m, N, Pa, K).

```json
{
  "nodes":   [{"id": "L0", "x": 0.0, "y": 0.0}, {"id": "U1", "x": 3.0, "y": 2.5}],
  "members": [{"id": "m1", "i": "L0", "j": "U1", "E": 2e11, "A": 0.003,
               "alpha": 1.2e-5, "dT": 40.0, "lack_of_fit": 0.002}],
  "supports": [{"node": "L0", "type": "pin", "settlement": [0.0, -0.01]},
               {"node": "L4", "type": "roller", "angle_deg": 120}],
  "loads":   [{"node": "U1", "fx": 15000.0, "fy": -20000.0}]
}
```

- Node and member ids are strings.
- `alpha` (thermal expansion coefficient), `dT` (temperature change) and `lack_of_fit` (the member was
  made this much too long, in metres) are optional and default to 0. A member's free (stress-free)
  elongation is `alpha * dT * L + lack_of_fit`.
- A **pin** fixes both displacement components. If `settlement` is given, the node is moved to that
  displacement `[ux, uy]` instead of being held at zero.
- A **roller** stops motion along one direction only. `angle_deg` is the direction of that direction,
  the normal to the rolling surface, measured counter-clockwise from +x. The node can move freely along
  the surface, and the reaction acts only along the normal. `angle_deg = 90` is an ordinary roller on
  horizontal ground.
- A node has at most one support. Loads act at nodes.

## Output

```json
{
  "stable": true,
  "displacements": {"L0": [0.0, -0.01], "U1": [0.00123, -0.0201]},
  "member_forces": {"m1": -15234.7},
  "reactions": {"L0": [-15000.0, 41200.0], "L4": [0.0, 0.0]}
}
```

- `displacements`: `[ux, uy]` in metres for every node.
- `member_forces`: axial force for every member, **tension positive**, in N.
- `reactions`: `[rx, ry]` in global axes for every supported node, and only those.

The verifier checks your numbers against the physics: each member force must match its elongation,
every node must be in equilibrium, and every support condition must hold, all to a relative tolerance
of about 1e-6.

If the truss is not stable, because it is a mechanism or has a rigid-body motion that the supports do
not stop (including infinitesimal mechanisms such as a node held only by two collinear members),
write `{"stable": false}` and exit with status 0. Do not reject a truss just because its members have
very different stiffnesses.

Each case must finish within 120 seconds. Python 3.11 with NumPy is available; no network access.
