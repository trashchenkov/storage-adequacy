#!/usr/bin/env python3
"""Independent targeted numerical checks for B1."""
from __future__ import annotations

if not __debug__:
    raise SystemExit("Run without -O or PYTHONOPTIMIZE: assertions validate the results.")

from fractions import Fraction as F
import argparse
import json
from pathlib import Path
import platform
import sys

import numpy as np
import scipy
from scipy.optimize import linprog


TOL = 1e-8
LOLH_TOL = 1e-8


def frac_list(values):
    return [F(str(value)) for value in values]


def bound_check(name, eta, E, Pc, Pd):
    eta_f, E_f, Pc_f, Pd_f = map(frac_list, (eta, E, Pc, Pd))
    beta = max(eta_f)
    w = [F(1) - e / beta for e in eta_f]
    D = [min(p, e) for p, e in zip(Pd_f, E_f)]
    C = [min(p, e / a) for p, e, a in zip(Pc_f, E_f, eta_f)]
    n = len(w)
    # variables u_i and B
    c_obj = [0.0] * n + [1.0]
    a_ub = [
        [*[float(d) for d in D], -1.0],
        [*[-float(c) for c in C], -1.0],
    ]
    b_ub = [0.0, -float(sum(wi * ci for wi, ci in zip(w, C)))]
    result = linprog(c_obj, A_ub=a_ub, b_ub=b_ub, bounds=[(0.0, float(wi)) for wi in w] + [(0.0, None)], method="highs")
    breakpoints = sorted({F(0), F(1), *(ci / (ci + di) for ci, di in zip(C, D) if ci + di > 0)})
    dual_values = []
    for lam in breakpoints:
        value = sum(wi * min(lam * di, (1 - lam) * ci) for wi, di, ci in zip(w, D, C))
        dual_values.append((lam, value))
    dual_lam, dual_value = max(dual_values, key=lambda item: item[1])
    per_device_harmonic = sum((wi * ci * di / (ci + di)) for wi, di, ci in zip(w, D, C) if ci + di > 0)
    delta = sum(wi * pd for wi, pd in zip(w, Pd_f))
    switch_groups = {}
    for idx, (ci, di) in enumerate(zip(C, D)):
        if ci + di > 0:
            switch_groups.setdefault(str(ci / (ci + di)), []).append(idx)
    degeneracies = {k: v for k, v in switch_groups.items() if len(v) > 1}
    return {
        "name": name,
        "success": bool(result.success),
        "status": int(result.status),
        "weights": [str(x) for x in w],
        "D": [str(x) for x in D],
        "C": [str(x) for x in C],
        "Delta_exact": str(delta),
        "Delta": float(delta),
        "B_primal": None if not result.success else float(result.fun),
        "B_dual_exact": str(dual_value),
        "B_dual": float(dual_value),
        "dual_lambda": str(dual_lam),
        "dual_breakpoints": [{"lambda": str(lam), "value": str(value)} for lam, value in dual_values],
        "primal_u": None if not result.success else result.x[:n].tolist(),
        "primal_dual_gap": None if not result.success else float(result.fun - float(dual_value)),
        "per_device_harmonic_sum": float(per_device_harmonic),
        "optimized_improves_per_device_harmonic_sum": result.success and float(result.fun) <= float(per_device_harmonic) + TOL,
        "switch_degeneracies": degeneracies,
    }


def build_fixed_sign_lp(x, E, Pc, Pd, eta, Ebar=None, require_small_zero=False, r_fixed=None, objective="r"):
    x = np.array(x, dtype=float)
    E = np.array(E, dtype=float)
    Pc = np.array(Pc, dtype=float)
    Pd = np.array(Pd, dtype=float)
    eta = np.array(eta, dtype=float)
    m, T = len(E), len(x)
    positive = x > 0
    negative = x < 0
    small_deficit = np.isclose(x, -1.5)
    r_cap = 0.25
    offset = 0
    idx = {}
    for name, shape in [("r", (1,)), ("f", (T,)), ("e", (m, T + 1)), ("c", (m, T)), ("d", (m, T)), ("loss", (T,)), ("spill", (T,))]:
        size = int(np.prod(shape))
        idx[name] = np.arange(offset, offset + size).reshape(shape)
        offset += size
    n = offset
    bounds = [(0.0, r_cap if r_fixed is None else r_fixed)]
    if r_fixed is not None:
        bounds[0] = (r_fixed, r_fixed)
    bounds += [(0.0, r_cap if r_fixed is None else r_fixed) for _ in range(T)]
    bounds += [(0.0, float(E[i])) for i in range(m) for _ in range(T + 1)]
    for i in range(m):
        for t in range(T):
            bounds.append((0.0, float(Pc[i]) if positive[t] else 0.0))
    for i in range(m):
        for t in range(T):
            bounds.append((0.0, float(Pd[i]) if negative[t] else 0.0))
    for t in range(T):
        hi = float(-x[t]) if negative[t] else 0.0
        if require_small_zero and small_deficit[t]:
            hi = 0.0
        bounds.append((0.0, hi))
    for t in range(T):
        bounds.append((0.0, float(x[t] + r_cap) if positive[t] else 0.0))

    eq_rows = []
    eq_rhs = []
    ub_rows = []
    ub_rhs = []
    for i in range(m):
        row = np.zeros(n)
        row[idx["e"][i, 0]] = 1.0
        row[idx["e"][i, T]] = -1.0
        eq_rows.append(row)
        eq_rhs.append(0.0)
        for t in range(T):
            row = np.zeros(n)
            row[idx["e"][i, t + 1]] = 1.0
            row[idx["e"][i, t]] = -1.0
            row[idx["c"][i, t]] = -eta[i]
            row[idx["d"][i, t]] = 1.0
            eq_rows.append(row)
            eq_rhs.append(0.0)
    for t in range(T):
        row = np.zeros(n)
        row[idx["f"][t]] = 1.0
        row[idx["loss"][t]] = 1.0
        row[idx["spill"][t]] = -1.0
        for i in range(m):
            row[idx["d"][i, t]] = 1.0
            row[idx["c"][i, t]] = -1.0
        eq_rows.append(row)
        eq_rhs.append(-x[t])
        row = np.zeros(n)
        row[idx["f"][t]] = 1.0
        row[idx["r"][0]] = -1.0
        ub_rows.append(row)
        ub_rhs.append(0.0)
    if Ebar is not None:
        row = np.zeros(n)
        for t in range(T):
            row[idx["loss"][t]] = 1.0
        ub_rows.append(row)
        ub_rhs.append(Ebar)
    c_obj = np.zeros(n)
    if objective == "r":
        c_obj[idx["r"][0]] = 1.0
    elif objective == "ens":
        for t in range(T):
            c_obj[idx["loss"][t]] = 1.0
    else:
        raise ValueError(objective)
    result = linprog(
        c_obj,
        A_ub=np.vstack(ub_rows) if ub_rows else None,
        b_ub=np.array(ub_rhs) if ub_rows else None,
        A_eq=np.vstack(eq_rows),
        b_eq=np.array(eq_rhs),
        bounds=bounds,
        method="highs",
    )
    return result, idx, bounds, np.vstack(eq_rows), np.array(eq_rhs), (np.vstack(ub_rows), np.array(ub_rhs)) if ub_rows else (None, None)


def unpack(result, idx):
    x = result.x
    return {name: x[index] for name, index in idx.items()}


def physical_check(params, values, r_value=None):
    x = np.array(params["x"], dtype=float)
    E = np.array(params["E"], dtype=float)
    Pc = np.array(params["Pc"], dtype=float)
    Pd = np.array(params["Pd"], dtype=float)
    eta = np.array(params["eta"], dtype=float)
    r = float(values["r"][0] if r_value is None else r_value)
    f, e, c, d, loss, spill = (np.array(values[k], dtype=float) for k in ("f", "e", "c", "d", "loss", "spill"))
    net = x + f
    charge = c.sum(axis=0)
    discharge = d.sum(axis=0)
    residuals = {
        "state": float(np.max(np.abs(np.diff(e, axis=1) - eta[:, None] * c + d), initial=0.0)),
        "cyclic": float(np.max(np.abs(e[:, 0] - e[:, -1]), initial=0.0)),
        "balance": float(np.max(np.abs(net + discharge - charge + loss - spill), initial=0.0)),
        "firm": float(max(0.0, np.max(f - r, initial=0.0), np.max(-f, initial=0.0))),
        "energy_bounds": float(max(0.0, np.max(e - E[:, None], initial=0.0), np.max(-e, initial=0.0))),
        "charge_bounds": float(max(0.0, np.max(c - Pc[:, None], initial=0.0), np.max(-c, initial=0.0))),
        "discharge_bounds": float(max(0.0, np.max(d - Pd[:, None], initial=0.0), np.max(-d, initial=0.0))),
        "external_charge": float(max(0.0, np.max(charge - np.maximum(net, 0.0), initial=0.0))),
        "useful_discharge": float(max(0.0, np.max(discharge - np.maximum(-net, 0.0), initial=0.0))),
        "loss_definition": float(np.max(np.abs(loss - (np.maximum(-net, 0.0) - discharge)), initial=0.0)),
        "spill_definition": float(np.max(np.abs(spill - (np.maximum(net, 0.0) - charge)), initial=0.0)),
        "simultaneous_charge_discharge": float(np.max(np.minimum(charge, discharge), initial=0.0)),
    }
    return {
        "max_residual": max(residuals.values()),
        "residuals": residuals,
        "ENS": float(loss.sum()),
        "LOLH": int(np.sum(loss > LOLH_TOL)),
        "net": net.tolist(),
        "charge": charge.tolist(),
        "discharge": discharge.tolist(),
    }


def solve_task(name, params, Ebar=None, Hbar=None, r_fixed=None, objective="r"):
    require_small_zero = Hbar == 1
    result, idx, bounds, aeq, beq, aub = build_fixed_sign_lp(
        params["x"], params["E"], params["Pc"], params["Pd"], params["eta"],
        Ebar=Ebar, require_small_zero=require_small_zero, r_fixed=r_fixed, objective=objective,
    )
    item = {
        "name": name,
        "status": int(result.status),
        "success": bool(result.success),
        "message": result.message,
        "objective": None if result.fun is None else float(result.fun),
    }
    if result.success:
        values = unpack(result, idx)
        check = physical_check(params, values)
        a_ub, b_ub = aub
        item.update({
            "r": float(values["r"][0]),
            "values": {k: np.asarray(v).tolist() for k, v in values.items()},
            "physical": check,
            "eq_residual": float(np.max(np.abs(aeq @ result.x - beq), initial=0.0)),
            "ub_violation": 0.0 if a_ub is None else float(np.max(np.maximum(a_ub @ result.x - b_ub, 0.0), initial=0.0)),
        })
    return item


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    output = parser.parse_args().output
    output.parent.mkdir(parents=True, exist_ok=True)
    main_params = {
        "x": [1.0, 1.0, -1.5, -3.0],
        "E": [2.0, 2.0],
        "Pc": [1.0, 1.0],
        "Pd": [1.0, 1.0],
        "eta": [0.8, 0.95],
    }
    Ebar = float(F(209, 80))
    bound_main = bound_check("main_two_store", main_params["eta"], main_params["E"], main_params["Pc"], main_params["Pd"])
    energy = solve_task("energy_only", main_params, Ebar=Ebar)
    hours = solve_task("hours_only", main_params, Hbar=1)
    joint = solve_task("joint", main_params, Ebar=Ebar, Hbar=1)

    alpha = np.array(main_params["eta"]) / max(main_params["eta"])
    auxiliary_params = dict(main_params)
    auxiliary_params["eta"] = [0.95, 0.95]
    auxiliary_params["Pc"] = alpha.tolist()
    aux = solve_task("homogeneous_auxiliary_r0", auxiliary_params, Ebar=Ebar, Hbar=1, r_fixed=0.0, objective="ens")

    transformed = None
    if aux["success"]:
        hv = {k: np.array(v, dtype=float) for k, v in aux["values"].items()}
        u = np.array(bound_main["primal_u"], dtype=float)
        z = 1.0 - u
        e = z[:, None] * hv["e"]
        d = z[:, None] * hv["d"]
        c = (z / alpha)[:, None] * hv["c"]
        charge_extra = (c - hv["c"]).sum(axis=0)
        discharge_extra = (hv["d"] - d).sum(axis=0)
        f = charge_extra + discharge_extra
        values = {"r": np.array([bound_main["B_primal"]]), "f": f, "e": e, "c": c, "d": d, "loss": hv["loss"], "spill": hv["spill"]}
        transformed = {
            "z": z.tolist(),
            "firm_charge_extra": charge_extra.tolist(),
            "firm_discharge_extra": discharge_extra.tolist(),
            "physical": physical_check(main_params, values, r_value=bound_main["B_primal"]),
            "values": {k: np.asarray(v).tolist() for k, v in values.items()},
        }

    bound_cases = [
        bound_check("zero_energy", [0.8, 0.95], [0, 2], [1, 1], [1, 1]),
        bound_check("zero_charge", [0.8, 0.95], [2, 2], [0, 1], [1, 1]),
        bound_check("zero_discharge", [0.8, 0.95], [2, 2], [1, 1], [0, 1]),
        bound_check("equal_eta", [0.95, 0.95], [2, 2], [1, 1], [1, 1]),
        bound_check("three_switching_devices", [0.6, 0.75, 0.95], [2, 3, 4], [0.4, 1.2, 3.0], [2.0, 0.7, 1.5]),
    ]

    expected = {
        "rE": 0.0,
        "rH": 0.0,
        "rjoint": float(F(13, 654)),
        "Delta": float(F(3, 19)),
        "B": float(F(3, 38)),
    }
    checks = {
        "bound_primal_dual_main": abs(bound_main["B_primal"] - bound_main["B_dual"]) <= TOL,
        "main_B_expected": abs(bound_main["B_primal"] - expected["B"]) <= TOL,
        "main_Delta_expected": abs(bound_main["Delta"] - expected["Delta"]) <= TOL,
        "rE_expected": energy["success"] and abs(energy["r"] - expected["rE"]) <= TOL,
        "rH_expected": hours["success"] and abs(hours["r"] - expected["rH"]) <= TOL,
        "rjoint_expected": joint["success"] and abs(joint["r"] - expected["rjoint"]) <= TOL,
        "joint_physical": joint["success"] and joint["physical"]["max_residual"] <= TOL and joint["physical"]["ENS"] <= Ebar + TOL and joint["physical"]["LOLH"] <= 1,
        "auxiliary_physical": aux["success"] and aux["physical"]["max_residual"] <= TOL and aux["physical"]["ENS"] <= Ebar + TOL and aux["physical"]["LOLH"] <= 1,
        "transformed_physical": transformed is not None and transformed["physical"]["max_residual"] <= TOL and transformed["physical"]["ENS"] <= Ebar + TOL and transformed["physical"]["LOLH"] <= 1,
        "all_bound_cases_primal_dual": all(case["success"] and abs(case["primal_dual_gap"]) <= TOL for case in bound_cases),
        "three_device_improvement": bound_cases[-1]["B_primal"] < bound_cases[-1]["per_device_harmonic_sum"] - TOL,
    }
    payload = {
        "environment": {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__, "scipy": scipy.__version__},
        "provenance": {
            "implementation": "independent targeted LP checker",
            "details": "See PROVENANCE.md.",
        },
        "main_bound": bound_main,
        "main_tasks": {"energy_only": energy, "hours_only": hours, "joint": joint},
        "homogeneous_auxiliary": aux,
        "mapped_transformation": transformed,
        "targeted_bound_cases": bound_cases,
        "expected": expected,
        "checks": checks,
        "all_passed": all(checks.values()),
    }
    output.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"all_passed": payload["all_passed"], "checks": checks, "output": str(output)}, indent=2))
    return 0 if payload["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
