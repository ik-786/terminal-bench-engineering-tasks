# Terminal-Bench engineering tasks

[![CI](https://github.com/ik-786/terminal-bench-engineering-tasks/actions/workflows/ci.yml/badge.svg)](https://github.com/ik-786/terminal-bench-engineering-tasks/actions/workflows/ci.yml)

Three agent benchmark tasks from different engineering disciplines, in the
[Terminal-Bench](https://www.tbench.ai/) / Harbor task format. Each one asks an agent to write a real
engineering solver inside a container, and grades it on hidden cases with a verifier that checks
physics rather than comparing against stored answers.

| task | discipline | the agent writes | what makes it hard |
|---|---|---|---|
| [`ac-power-flow`](tasks/ac-power-flow) | electrical | Newton-Raphson AC power flow | off-nominal taps and phase shifters, generator reactive limits that bind and must sometimes be released again, unsorted bus ids, out-of-service lines, a case beyond voltage collapse |
| [`truss-analysis`](tasks/truss-analysis) | civil | stiffness-method plane truss solver | inclined rollers, support settlement, temperature change and lack of fit in indeterminate trusses, mechanisms (including an infinitesimal one) to detect, member stiffnesses 10^6 apart |
| [`lambert-multirev`](tasks/lambert-multirev) | aerospace | multi-revolution Lambert solver | both branches of 1- and 2-revolution transfers, retrograde and 180-degree transfers, hyperbolic and near-parabolic times of flight, an infeasible case, a per-problem gravitational parameter |

A fourth task in the same format, a chemical reactor safety controller, lives in its own repository:
[semi-batch-reactor-safety-control](https://github.com/ik-786/semi-batch-reactor-safety-control).

## Task layout

```text
tasks/<task>/
  instruction.md        the task statement the agent sees
  task.toml             metadata and time limits
  environment/          Dockerfile and the files placed in /app (one example input)
  solution/solve.sh     reference solution (oracle), installs the solver in /app
  tests/test.sh         runs the verifier and writes /logs/verifier/reward.txt (1 = pass, 0 = fail)
  tests/test_outputs.py the verifier
  tests/cases/          hidden inputs, only copied in at grading time
  dev/make_cases.py     how the hidden cases were generated and checked
validation/<task>/      deliberately wrong solutions that must score 0
scripts/run_local.py    runs any task without Docker
```

## How the verifiers work

None of the verifiers store expected outputs to compare against. They recompute the governing
equations from what the solver returns:

- **Power flow**: rebuilds the bus admittance matrix and requires every power mismatch below 1e-3 MW,
  every branch flow and the losses to be consistent with the reported voltages, and every PV bus to be
  in one of the three valid states (at its set point within limits, or at a limit with the voltage on
  the correct side).
- **Truss**: checks the member constitutive law, equilibrium at every node and every support condition.
  For a stable truss those three determine the solution uniquely.
- **Lambert**: flies the reported departure velocity with its own universal-variable Kepler propagator
  and requires it to land on the target, then checks the direction of motion, the number of full
  revolutions and the branch (through the semi-major axis).

Cases that have no solution (beyond voltage collapse, mechanisms, too few seconds for three
revolutions) must be reported as such. A solver that claims a solution there fails, because the claimed
numbers cannot satisfy the equations.

## Validation

`python scripts/run_local.py --all` runs every task with its reference solution, with no solution, and
with each deliberately wrong solution in `validation/`. CI runs it on every push, together with the
oracle and null checks inside Docker.

| task | reference | no solution | wrong solutions rejected |
|---|---|---|---|
| ac-power-flow | 1 | 0 | 7 of 7: DC approximation, ignores reactive limits, never releases a limit, ignores transformer taps, claims convergence when Newton fails, stops at a loose tolerance, replays the example |
| truss-analysis | 1 | 0 | 6 of 6: ignores temperature and lack of fit, treats every roller as vertical, ignores settlement, a stability threshold that rejects stiff/soft mixes, no mechanism check, compression-positive sign |
| lambert-multirev | 1 | 0 | 6 of 6: zero revolutions only, swapped branches, always prograde, answers the infeasible case, ignores a per-problem mu, builds the plane from r1 x r2 in a 180-degree transfer |

Each wrong solution fails on the hidden cases aimed at its mistake. For example, `one_way_limits` fails
only the two power flow networks where a reactive limit has to be released, and `pivot_tol_too_loose`
fails only the tower with very soft bracing.

One variant was dropped: a Lambert solver without the series expansion near the parabolic case still
passed. Even very close to the parabolic time of flight, the closed-form expression stays accurate
enough for the 1e-6 position tolerance. So the near-parabolic case checks that a solver does not break
there; it does not reward the series.

## Running a task

With Docker, from a task folder:

```bash
docker build -t task-env environment
# reference solution, expect 1
docker run --rm -v "$PWD/tests:/tests" -v "$PWD/solution:/solution" task-env \
  bash -c "bash /solution/solve.sh && bash /tests/test.sh; cat /logs/verifier/reward.txt"
# no solution, expect 0
docker run --rm -v "$PWD/tests:/tests" task-env bash -c "bash /tests/test.sh; cat /logs/verifier/reward.txt"
```

Without Docker (Python 3.10+ and pytest), from the repository root:

```bash
python scripts/run_local.py truss-analysis --oracle
python scripts/run_local.py truss-analysis --solution my_truss.py
python scripts/run_local.py --all
```

To regenerate a task's hidden cases: `cd tasks/<task> && python dev/make_cases.py`.

## Notes

- The hidden cases are public in this repository, so these tasks suit demonstration and agent
  development, not contamination-free evaluation. For a private benchmark, keep `tests/cases/` out of
  any public copy.
- All data is synthetic and generated by the scripts in `dev/`. The networks and trusses are not copies
  of standard test systems.
- The reference solvers use only the Python standard library, so that any numerical library the agent
  chooses can be checked against an independent implementation.

## Related

- [reward-hacking-detector](https://github.com/ik-786/reward-hacking-detector): static, differential
  and trajectory checks for agents that pass tests without solving the task.

## License

MIT
