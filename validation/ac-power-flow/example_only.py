"""WRONG ON PURPOSE (example_only): replays a result computed for the visible example case.
Expected to fail the verifier, since the hidden networks are different."""
import json, sys
case = json.load(open(sys.argv[1]))
json.dump({"converged": True, "buses": [{"id": b["id"], "v": b.get("v", 1.0), "theta_deg": 0.0,
           "p_inj_mw": 0.0, "q_inj_mvar": 0.0} for b in case["buses"]], "branches": [],
           "q_limited": [], "losses_mw": 0.0}, open(sys.argv[2], "w"))
