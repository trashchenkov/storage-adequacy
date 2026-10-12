#!/usr/bin/env python3
"""Exact rational certificates for manuscript Appendices A and B.

Standard library only. No optimization solver is used.
"""
from __future__ import annotations

if not __debug__:
    raise SystemExit("Run without -O or PYTHONOPTIMIZE: assertions validate the results.")

from fractions import Fraction as F
import argparse
import json
from pathlib import Path




def s(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def f(value: F) -> float:
    return float(value)


def max_abs(values):
    return max(abs(v) for v in values) if values else F(0)


def check_schedule(x, eta, E, Pc, Pd, r, firm, charge, discharge, loss, spill=None):
    m, T = len(E), len(x)
    if spill is None:
        spill = [F(0) for _ in range(T)]
    states = [[F(0)] for _ in range(m)]
    residuals = []
    violations = []
    for t in range(T):
        net = x[t] + firm[t]
        total_c = sum(charge[i][t] for i in range(m))
        total_d = sum(discharge[i][t] for i in range(m))
        residuals.append(net + total_d - total_c + loss[t] - spill[t])
        violations.extend([r - firm[t], firm[t], r])
        violations.append(max(net, F(0)) - total_c)
        violations.append(max(-net, F(0)) - total_d)
        for i in range(m):
            next_e = states[i][-1] + eta[i] * charge[i][t] - discharge[i][t]
            states[i].append(next_e)
            violations.extend([
                charge[i][t],
                Pc[i] - charge[i][t],
                discharge[i][t],
                Pd[i] - discharge[i][t],
                next_e,
                E[i] - next_e,
            ])
        violations.extend([loss[t], spill[t]])
    for i in range(m):
        residuals.append(states[i][-1] - states[i][0])
    lolh = sum(1 for item in loss if item > 0)
    ens = sum(loss)
    return {
        "states": [[s(v) for v in row] for row in states],
        "ENS": s(ens),
        "LOLH": lolh,
        "max_abs_balance_or_cyclic_residual": s(max_abs(residuals)),
        "min_constraint_margin": s(min(violations)),
        "passed": max_abs(residuals) == 0 and min(violations) >= 0,
    }


def appendix_a():
    x = [F(1), F(1), -F(3, 2), -F(3)]
    eta = [F(4, 5), F(19, 20)]
    E = [F(2), F(2)]
    Pc = [F(1), F(1)]
    Pd = [F(1), F(1)]
    zero = [F(0)] * 4
    empty = [[F(0)] * 4, [F(0)] * 4]

    energy_charge = [[F(0)] * 4, [F(1), F(1), F(0), F(0)]]
    energy_discharge = [[F(0)] * 4, [F(0), F(0), F(9, 10), F(1)]]
    energy_loss = [F(0), F(0), F(3, 5), F(2)]

    hours_charge = [[F(5, 16), F(5, 16), F(0), F(0)], [F(11, 16), F(11, 16), F(0), F(0)]]
    hours_discharge = [[F(0), F(0), F(1, 2), F(0)], [F(0), F(0), F(1), F(49, 160)]]
    hours_loss = [F(0), F(0), F(0), F(431, 160)]

    r = F(13, 654)
    cA = F(5, 8) * (F(1, 2) - r)
    cB = F(1) + r - cA
    dA_mid = F(1, 2) - r
    dB_last = F(19, 10) * cB - F(1)
    joint_charge = [[cA, cA, F(0), F(0)], [cB, cB, F(0), F(0)]]
    joint_discharge = [[F(0), F(0), dA_mid, F(0)], [F(0), F(0), F(1), dB_last]]
    joint_firm = [r, r, r, r]
    joint_loss = [F(0), F(0), F(0), F(209, 80)]

    lower_intercept = F(431, 160)
    lower_slope = F(327, 80)
    lower_at_r = lower_intercept - lower_slope * r

    # Explicit B-witness from the proof text.
    beta = F(19, 20)
    alpha = [F(16, 19), F(1)]
    u = [F(3, 38), F(0)]
    z = [F(1) - u[0], F(1)]
    chat = [[F(4, 19), F(16, 19), F(0), F(0)], [F(15, 19), F(3, 19), F(0), F(0)]]
    dhat = [[F(0), F(0), F(1), F(0)], [F(0), F(0), F(1, 2), F(2, 5)]]
    mapped_charge = [[z[i] / alpha[i] * chat[i][t] for t in range(4)] for i in range(2)]
    mapped_discharge = [[z[i] * dhat[i][t] for t in range(4)] for i in range(2)]
    mapped_firm = [F(3, 152), F(3, 38), F(3, 38), F(0)]
    mapped_loss = [F(0), F(0), F(0), F(13, 5)]
    B = F(3, 38)

    checks = {
        "energy_r0": check_schedule(x, eta, E, Pc, Pd, F(0), zero, energy_charge, energy_discharge, energy_loss),
        "hours_r0": check_schedule(x, eta, E, Pc, Pd, F(0), zero, hours_charge, hours_discharge, hours_loss),
        "joint_optimal": check_schedule(x, eta, E, Pc, Pd, r, joint_firm, joint_charge, joint_discharge, joint_loss),
        "explicit_B_witness": check_schedule(x, eta, E, Pc, Pd, B, mapped_firm, mapped_charge, mapped_discharge, mapped_loss),
    }
    return {
        "lower_bound_intercept": s(lower_intercept),
        "lower_bound_slope": s(lower_slope),
        "r_star": s(r),
        "lower_bound_at_r_star": s(lower_at_r),
        "Ebar": "209/80",
        "r_star_less_than_quarter": r < F(1, 4),
        "explicit_B_firm": [s(v) for v in mapped_firm],
        "explicit_B_loss": [s(v) for v in mapped_loss],
        "checks": checks,
        "passed": lower_at_r == F(209, 80) and all(item["passed"] for item in checks.values()),
    }


def family_case(N: int):
    a = F(4, 5)
    b = F(19, 20)
    p = F(1)
    k = b / a
    K = N * k + 1 + b
    L = k * N * p
    Q = L * (N + 1) / b
    r = N * p * (k - 1) / K
    Delta = p * (1 - a / b)
    B = Delta * F(N, 1) / (F(N, 1) + a)
    x = [Q] + [-(L + p)] * N + [-(L + 2 * p)]
    eta = [a, b]
    E = [F(N), L * (N + 1)]
    Pc = [F(N) / a, Q]
    Pd = [p, L]
    T = N + 2

    # Separate energy schedule.
    c_energy = [[F(0)] * T, [Q] + [F(0)] * (T - 1)]
    d_energy = [[F(0)] * T, [F(0)] + [L] * (N + 1)]
    loss_energy = [F(0)] + [p] * N + [2 * p]

    # Separate shortage-hour schedule.
    c_hours = [[F(N) / a] + [F(0)] * (T - 1), [Q - F(N) / a] + [F(0)] * (T - 1)]
    d_hours = [[F(0)] + [p] * N + [F(0)], [F(0)] + [L] * N + [F(0)]]
    loss_hours = [F(0)] + [F(0)] * N + [L + 2 * p]

    # Joint attaining schedule.
    cA = N * (p - r) / a
    cB = Q + r - cA
    y = (K - 1) * r
    c_joint = [[cA] + [F(0)] * (T - 1), [cB] + [F(0)] * (T - 1)]
    d_joint = [[F(0)] + [p - r] * N + [F(0)], [F(0)] + [L] * N + [y]]
    f_joint = [r] * T
    loss_joint = [F(0)] + [F(0)] * N + [F(N) + 2]

    zeros = [F(0)] * T
    checks = {
        "energy_r0": check_schedule(x, eta, E, Pc, Pd, F(0), zeros, c_energy, d_energy, loss_energy),
        "hours_r0": check_schedule(x, eta, E, Pc, Pd, F(0), zeros, c_hours, d_hours, loss_hours),
        "joint_optimal": check_schedule(x, eta, E, Pc, Pd, r, f_joint, c_joint, d_joint, loss_joint),
    }
    gamma_over_delta = r / Delta
    gamma_over_B = r / B
    return {
        "N": N,
        "K": s(K),
        "L": s(L),
        "Q": s(Q),
        "Gamma": s(r),
        "Delta": s(Delta),
        "B": s(B),
        "Gamma_over_Delta": s(gamma_over_delta),
        "Gamma_over_B": s(gamma_over_B),
        "formula_ratio_Gamma_over_B": s((N * k + b) / (N * k + 1 + b)),
        "bound_order": {
            "Gamma_lt_B": r < B,
            "B_le_Delta": B <= Delta,
        },
        "checks": checks,
        "passed": all(item["passed"] for item in checks.values()) and gamma_over_B == (N * k + b) / (N * k + 1 + b),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    output = parser.parse_args().output
    appendix_a_result = appendix_a()
    appendix_b = [family_case(N) for N in (1, 4, 16, 64, 1024)]
    payload = {
        "source": "Exact rational certificates for the four-hour example and sharpness family.",
        "appendix_A": appendix_a_result,
        "appendix_B": appendix_b,
        "all_passed": appendix_a_result["passed"] and all(item["passed"] for item in appendix_b),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"all_passed": payload["all_passed"], "output": str(output)}, indent=2))
    return 0 if payload["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
