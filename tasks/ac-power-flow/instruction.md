# AC power flow with transformer taps and generator reactive limits

Write a steady-state AC power flow solver for transmission networks and save it as
`/app/powerflow.py`. It is run as

```
python /app/powerflow.py <case.json> <result.json>
```

It reads a network from `<case.json>`, solves it, writes `<result.json>` and exits with
status 0. An example network is in `/app/data/example_case.json`. Your solver will be
tested on other networks, up to 60 buses, that you do not have.

## Input

```json
{
  "base_mva": 100.0,
  "buses": [
    {"id": 1, "type": "slack", "v": 1.05, "theta_deg": 0.0, "pd": 0.0, "qd": 0.0},
    {"id": 2, "type": "pv", "v": 1.02, "pg": 60.0, "qmin": -20.0, "qmax": 35.0, "pd": 12.0, "qd": 5.0},
    {"id": 3, "type": "pq", "pd": 30.0, "qd": 10.0, "gs": 0.5, "bs": 12.0}
  ],
  "branches": [
    {"from": 1, "to": 2, "r": 0.02, "x": 0.08, "b": 0.03, "tap": 1.0, "shift_deg": 0.0, "status": 1}
  ]
}
```

- Bus ids are arbitrary integers, not necessarily consecutive or sorted. There is exactly one slack bus.
- `pd`, `qd` (load), `pg` (PV generation), `qmin`, `qmax` (PV reactive limits) are in MW or MVAr.
  Missing `pd`, `qd`, `gs`, `bs` mean 0.
- `gs` + j`bs` is a shunt admittance at the bus, given as the MW and MVAr it draws at 1.0 pu voltage
  (`bs` > 0 is capacitive, it injects reactive power).
- `r`, `x` (series impedance) and `b` (total line charging susceptance) are per unit on `base_mva`.
- `tap` is the off-nominal turns ratio on the `from` side (0 or missing means 1.0) and `shift_deg` a phase
  shift. Branches with `status` 0 are out of service and must be ignored.

Use the standard pi model with the transformer on the `from` side. With `ys = 1/(r + jx)` and the complex
ratio `a = tap * exp(j * shift)`:

```
I_from = (ys + j b/2) / |a|^2 * V_from  -  ys / conj(a) * V_to
I_to   = -ys / a * V_from              +  (ys + j b/2) * V_to
```

## Bus types and reactive limits

- **Slack**: voltage magnitude `v` and angle `theta_deg` are fixed.
- **PQ**: net injection is fixed at `P = -pd`, `Q = -qd`.
- **PV**: `P = pg - pd` is fixed. Let `Qg = Q + qd` be the generator's reactive output. In the solution each
  PV bus must be in exactly one of these states:
  1. at its voltage set point `v`, with `qmin <= Qg <= qmax`;
  2. at `Qg = qmax`, with voltage at or below the set point;
  3. at `Qg = qmin`, with voltage at or above the set point.

  A bus that hits a limit during the iterations may have to come off it again later.

## Output

```json
{
  "converged": true,
  "buses": [{"id": 1, "v": 1.05, "theta_deg": 0.0, "p_inj_mw": 95.1, "q_inj_mvar": 12.3}],
  "branches": [{"from": 1, "to": 2, "p_from_mw": 40.2, "q_from_mvar": 3.1, "p_to_mw": -39.6, "q_to_mvar": -2.0}],
  "q_limited": [5, 9],
  "losses_mw": 3.4
}
```

- `buses`: every bus once, with voltage magnitude (pu), angle (degrees) and the net injection into the
  network (MW, MVAr), which is generation minus load. Bus shunts count as part of the network, like line
  charging, so they are not subtracted here.
- `branches`: every in-service branch once, with the complex power entering it at each end.
- `q_limited`: ids of PV buses held at a reactive limit (states 2 and 3 above).
- `losses_mw`: sum of `p_inj_mw` over all buses.

Write full-precision numbers. The verifier recomputes every injection and flow from your voltages,
so power mismatches must be below 1e-3 MW/MVAr and PV voltages must match their set points to 1e-6 pu.
All voltages must be between 0.5 and 1.5 pu.

If a network has no solution, for example because the load is beyond the point of voltage collapse,
write `{"converged": false}` and exit with status 0.

Each case must finish within 120 seconds. Python 3.11 with NumPy is available; no network access.
