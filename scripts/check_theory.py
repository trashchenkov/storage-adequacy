#!/usr/bin/env python3
"""Primary finite-case checks for heterogeneous-storage adequacy.

The output file must not exist. Predefined cases are in data/cases.json.
The physical MILP retains explicit dispatch, external-only charging and freely
chosen cyclic states. Finite computations complement the manuscript proofs.
"""
from __future__ import annotations

if not __debug__:
    raise SystemExit("Run without -O or PYTHONOPTIMIZE: assertions validate the results.")

import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np
import scipy
from scipy.optimize import Bounds, LinearConstraint, linprog, milp
from scipy.sparse import lil_matrix

CASE_HASH = "746a5a603b557717cf906b5d68f89ec5d3ceebdb3ad0b4d39ff7f6f5510a0d3f"
TOL = 5e-7
LOSS_TOL = 1e-7


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def numeric(case):
    return {k: np.array([float(F(v)) for v in case[k]])
            for k in ("x", "E", "Pc", "Pd", "eta")}


def physical(a, r, f, e, c, d, loss, spill):
    net = a["x"] + f
    charge, discharge = c.sum(axis=0), d.sum(axis=0)
    residuals = {
        "state_balance": np.max(np.abs(np.diff(e, axis=1) - a["eta"][:, None] * c + d)),
        "cyclic_state": np.max(np.abs(e[:, 0] - e[:, -1])),
        "node_balance": np.max(np.abs(net + discharge - charge + loss - spill)),
        "negative_variable": max(0.0, *[-np.min(v) for v in (f, e, c, d, loss, spill)]),
        "firm_cap": max(0.0, float(np.max(f-r))),
        "energy_cap": max(0.0, float(np.max(e-a["E"][:, None]))),
        "charge_cap": max(0.0, float(np.max(c-a["Pc"][:, None]))),
        "discharge_cap": max(0.0, float(np.max(d-a["Pd"][:, None]))),
        "external_charge": max(0.0, float(np.max(charge-np.maximum(net, 0)))),
        "useful_discharge": max(0.0, float(np.max(discharge-np.maximum(-net, 0)))),
        "loss_definition": np.max(np.abs(loss-(np.maximum(-net, 0)-discharge))),
        "spill_definition": np.max(np.abs(spill-(np.maximum(net, 0)-charge))),
        "simultaneous_fleet_charge_discharge": np.max(np.minimum(charge, discharge)),
    }
    residuals = {k: float(v) for k, v in residuals.items()}
    return {"max_violation": max(residuals.values()), "residuals": residuals,
            "ENS": float(sum(loss)), "LOLH_numeric": int(sum(loss > LOSS_TOL)),
            "loss_count_tolerance": LOSS_TOL,
            "note": "The theorem uses strict positivity; the tolerance is only for numerical reporting."}


def solve(a, Ebar=None, Hbar=None, fixed_r=None, objective="r", full_firm=False):
    """Exact finite-dimensional physical MILP, to the recorded solver tolerances."""
    m, T = len(a["E"]), len(a["x"])
    at = 0
    ix = {}
    for name, shape in (("r", (1,)), ("f", (T,)), ("e", (m, T+1)),
                        ("c", (m, T)), ("d", (m, T)), ("loss", (T,)),
                        ("spill", (T,)), ("mode", (T,)), ("event", (T,))):
        size = int(np.prod(shape))
        ix[name] = np.arange(at, at+size).reshape(shape)
        at += size
    lb, ub = np.zeros(at), np.full(at, np.inf)
    integrality = np.zeros(at, dtype=int)
    rmax = max(0.0, float(max(-a["x"])))
    ub[ix["r"]] = rmax
    if fixed_r is not None:
        lb[ix["r"]] = ub[ix["r"]] = fixed_r
    ub[ix["f"]] = rmax
    if full_firm:
        assert fixed_r is not None
        lb[ix["f"]] = ub[ix["f"]] = fixed_r
    for name, caps in (("e", a["E"]), ("c", a["Pc"]), ("d", a["Pd"])):
        ub[ix[name]] = np.broadcast_to(caps[:, None], ix[name].shape)
    loss_cap, spill_cap = np.maximum(-a["x"], 0), np.maximum(a["x"]+rmax, 0)
    ub[ix["loss"]], ub[ix["spill"]] = loss_cap, spill_cap
    ub[ix["mode"]] = ub[ix["event"]] = 1
    integrality[ix["mode"]] = 1
    if Hbar is not None:
        integrality[ix["event"]] = 1
    else:
        lb[ix["event"]] = ub[ix["event"]] = (loss_cap > 0).astype(float)
    for t in range(T):
        if a["x"][t] >= 0:
            lb[ix["mode"][t]] = ub[ix["mode"][t]] = 1
        if loss_cap[t] == 0:
            lb[ix["event"][t]] = ub[ix["event"][t]] = 0
    rows, lo, hi = [], [], []

    def row(entries, low=-np.inf, high=np.inf):
        rows.append(entries)
        lo.append(low)
        hi.append(high)

    for i in range(m):
        row({int(ix["e"][i, 0]): 1, int(ix["e"][i, T]): -1}, 0, 0)
        for t in range(T):
            row({int(ix["e"][i, t+1]): 1, int(ix["e"][i, t]): -1,
                 int(ix["c"][i, t]): -a["eta"][i], int(ix["d"][i, t]): 1}, 0, 0)
            row({int(ix["c"][i, t]): 1, int(ix["mode"][t]): -a["Pc"][i]}, high=0)
            row({int(ix["d"][i, t]): 1, int(ix["mode"][t]): a["Pd"][i]}, high=a["Pd"][i])
    for t in range(T):
        row({int(ix["f"][t]): 1, int(ix["r"][0]): -1}, high=0)
        balance = {int(ix["f"][t]): 1, int(ix["loss"][t]): 1, int(ix["spill"][t]): -1}
        for i in range(m):
            balance[int(ix["c"][i, t])] = -1
            balance[int(ix["d"][i, t])] = 1
        row(balance, -a["x"][t], -a["x"][t])
        row({int(ix["loss"][t]): 1, int(ix["mode"][t]): loss_cap[t]}, high=loss_cap[t])
        row({int(ix["spill"][t]): 1, int(ix["mode"][t]): -spill_cap[t]}, high=0)
        row({int(ix["loss"][t]): 1, int(ix["event"][t]): -loss_cap[t]}, high=0)
    if Ebar is not None:
        row({int(j): 1 for j in ix["loss"]}, high=Ebar)
    if Hbar is not None:
        row({int(j): 1 for j in ix["event"]}, high=Hbar)
    matrix = lil_matrix((len(rows), at))
    for k, entries in enumerate(rows):
        for j, val in entries.items():
            matrix[k, j] += val
    cost = np.zeros(at)
    cost[ix["r"] if objective == "r" else ix["loss"]] = 1
    result = milp(cost, integrality=integrality, bounds=Bounds(lb, ub),
                  constraints=LinearConstraint(matrix.tocsr(), np.array(lo), np.array(hi)),
                  options={"presolve": True, "mip_rel_gap": 1e-10, "time_limit": 30})
    report = {"status": int(result.status), "success": bool(result.success),
              "message": result.message, "objective_kind": objective}
    for key in ("fun", "mip_gap", "mip_dual_bound", "mip_node_count"):
        value = getattr(result, key, None)
        report[key] = float(value) if value is not None else None
    if not result.success:
        return report
    values = {name: result.x[index] for name, index in ix.items()}
    r = float(values["r"][0])
    check = physical(a, r, *(values[k] for k in ("f", "e", "c", "d", "loss", "spill")))
    report.update({"r": r, "values": {k: v.tolist() for k, v in values.items()}, "physical": check})
    if Ebar is not None:
        report["ENS_cap_violation"] = max(0.0, check["ENS"]-Ebar)
    if Hbar is not None:
        report["LOLH_cap_violation"] = max(0, check["LOLH_numeric"]-Hbar)
    report["passed"] = check["max_violation"] <= TOL and report.get("ENS_cap_violation", 0) <= TOL and report.get("LOLH_cap_violation", 0) == 0
    if report["mip_dual_bound"] is not None:
        report["objective_minus_dual_bound"] = report["fun"]-report["mip_dual_bound"]
        report["passed"] &= abs(report["objective_minus_dual_bound"]) <= TOL
    return report


def bound(case):
    eta, E, Pc, Pd = ([F(v) for v in case[k]] for k in ("eta", "E", "Pc", "Pd"))
    beta = max(eta)
    w = [1-v/beta for v in eta]
    D = [min(p, e) for p, e in zip(Pd, E)]
    C = [min(p, e/v) for p, e, v in zip(Pc, E, eta)]
    breaks = sorted({F(0), F(1), *(c/(c+d) for c, d in zip(C, D) if c+d > 0)})
    dual = [(t, sum(v*min(t*d, (1-t)*c) for v, d, c in zip(w, D, C))) for t in breaks]
    dual_t, exact = max(dual, key=lambda pair: pair[1])
    m = len(w)
    A = np.array([[*[float(x) for x in D], -1], [*[-float(x) for x in C], -1]])
    b = np.array([0, -float(sum(v*c for v, c in zip(w, C)))])
    result = linprog([0]*m+[1], A_ub=A, b_ub=b,
                     bounds=[(0, float(x)) for x in w]+[(0, None)], method="highs")
    assert result.success and abs(result.fun-float(exact)) < 1e-8
    delta = sum(v*p for v, p in zip(w, Pd))
    Dsum, Csum = sum(v*d for v, d in zip(w, D)), sum(v*c for v, c in zip(w, C))
    harmonic = Dsum*Csum/(Dsum+Csum) if Dsum+Csum else F(0)
    return {"beta": str(beta), "weights": list(map(str, w)), "D": list(map(str, D)), "C": list(map(str, C)),
            "Delta": str(delta), "Delta_float": float(delta), "B_exact": str(exact), "B_float": float(exact),
            "harmonic_bound": str(harmonic), "dual_maximizer": str(dual_t),
            "dual_breakpoints": [{"theta": str(t), "value": str(v)} for t, v in dual],
            "primal_u": result.x[:m].tolist(), "primal_objective": float(result.fun)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("Output already exists; choose a new path.")
    here = Path(__file__).resolve().parent
    case_path = here.parent / "data" / "cases.json"
    assert sha(case_path) == CASE_HASH
    frozen = json.loads(case_path.read_text())
    report = {"started_at_utc": datetime.now(timezone.utc).isoformat(),
              "attribution": "Primary implementation; see PROVENANCE.md for origin and adaptations.",
              "code_sha256": sha(__file__), "case_sha256": sha(case_path),
              "environment": {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__, "scipy": scipy.__version__, "solver": "scipy.optimize.milp (HiGHS)", "options": {"mip_rel_gap": 1e-10, "time_limit": 30, "presolve": True}}, "cases": []}
    all_passed = True
    for case in frozen["cases"]:
        a, bnd = numeric(case), bound(case)
        Ebar, Hbar = float(F(case["Ebar"])), case["Hbar"]
        out = {"name": case["name"], "parameters": case, "bound": bnd}
        out["energy_only"] = solve(a, Ebar=Ebar)
        out["hours_only"] = solve(a, Hbar=Hbar)
        out["joint"] = solve(a, Ebar=Ebar, Hbar=Hbar)
        successful = all(out[k].get("passed", False) for k in ("energy_only", "hours_only", "joint"))
        if successful:
            r0 = max(out["energy_only"]["r"], out["hours_only"]["r"])
            gamma = out["joint"]["r"]-r0
            out.update({"r0": r0, "Gamma": gamma, "bound_holds_numeric": -TOL <= gamma <= bnd["B_float"]+TOL})
            alpha = a["eta"]/max(a["eta"])
            hom = {k: v.copy() for k, v in a.items()}
            hom["eta"][:] = max(a["eta"])
            hom["Pc"] *= alpha
            auxiliary = solve(hom, Hbar=Hbar, fixed_r=r0, objective="ENS", full_firm=True)
            out["homogeneous_auxiliary"] = auxiliary
            if auxiliary.get("passed"):
                hv = {k: np.array(v) for k, v in auxiliary["values"].items()}
                u = np.array(bnd["primal_u"])
                z = 1-u
                e, c, d = z[:, None]*hv["e"], (z/alpha)[:, None]*hv["c"], z[:, None]*hv["d"]
                charge_extra = (c-hv["c"]).sum(axis=0)
                discharge_extra = (hv["d"]-d).sum(axis=0)
                f = r0+charge_extra+discharge_extra
                r_witness = r0+bnd["B_float"]
                check = physical(a, r_witness, f, e, c, d, hv["loss"], hv["spill"])
                mapped_ok = check["max_violation"] <= TOL and check["ENS"] <= Ebar+TOL and check["LOLH_numeric"] <= Hbar
                out["mapped_witness"] = {"r": r_witness, "f": f.tolist(), "e": e.tolist(), "c": c.tolist(), "d": d.tolist(), "loss": hv["loss"].tolist(), "spill": hv["spill"].tolist(), "charge_extra": charge_extra.tolist(), "discharge_extra": discharge_extra.tolist(), "physical": check, "passed": mapped_ok}
            else:
                mapped_ok = False
            positive = a["x"]+r0
            cap_sum = sum(float(F(v)) for v in bnd["C"])
            satisfies_B2 = all(v >= cap_sum-TOL for v in positive if v > LOSS_TOL)
            out["B2_condition_numeric"] = satisfies_B2
            out["B2_implication_passed"] = not satisfies_B2 or abs(gamma) <= TOL
            exact = case.get("expected_gamma")
            out["expected_gamma_check"] = exact is None or abs(gamma-float(F(exact))) <= TOL
            successful &= out["bound_holds_numeric"] and mapped_ok and out["B2_implication_passed"] and out["expected_gamma_check"]
        out["passed"] = bool(successful)
        all_passed &= bool(successful)
        report["cases"].append(out)
        print(json.dumps({"case": case["name"], "passed": bool(successful), "Gamma": out.get("Gamma"), "B": bnd["B_float"]}), flush=True)
    report["sharp_family_exact_control"] = []
    for N in (1, 4, 16, 64, 1024):
        gamma, delta = F(15*N, 95*N+156), F(3, 19)
        B = delta*F(N)/(F(N)+F(4, 5))
        report["sharp_family_exact_control"].append({"N": N, "Gamma": str(gamma), "Delta": str(delta), "B": str(B), "Gamma_over_B": str(gamma/B), "checked_exactly": gamma < B < delta})
    report["all_passed"] = all_passed and all(x["checked_exactly"] for x in report["sharp_family_exact_control"])
    report["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n")
    if not report["all_passed"]:
        raise SystemExit("At least one finite check failed; see the saved complete results.")


if __name__ == "__main__":
    main()
