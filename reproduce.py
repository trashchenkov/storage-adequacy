#!/usr/bin/env python3
"""Recalculate, validate against reference results, and generate numerical plots."""
from fractions import Fraction as F
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

sys.dont_write_bytecode = True
from verify_repository import verify

ROOT = Path(__file__).resolve().parent
REFERENCE = ROOT / "results" / "reference"
OUTPUT = ROOT / "generated"
TOL = 5e-7
TASKS = (
    ("check_theory.py", "primary_results.json"),
    ("independent_check.py", "independent_results.json"),
    ("exact_certificates.py", "exact_certificates.json"),
    ("check_family.py", "family_checks.json"),
)


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def close(actual, expected, label):
    if not math.isfinite(actual) or abs(actual - expected) > TOL:
        raise RuntimeError(f"{label}: {actual} differs from {expected}")


def validate(directory):
    """Check optima and recompute physical residuals; allow nonunique schedules."""
    import numpy as np

    spec = importlib.util.spec_from_file_location("primary", ROOT / "scripts/check_theory.py")
    primary = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(primary)
    fresh = {name: read(directory / name) for _, name in TASKS}
    if not all(item["all_passed"] is True for item in fresh.values()):
        raise RuntimeError("A calculation reported failure.")
    result = fresh["primary_results.json"]
    reference = read(REFERENCE / "primary_results.json")
    cases = read(ROOT / "data/cases.json")["cases"]
    assert len(result["cases"]) == len(reference["cases"]) == len(cases) == 15
    physical_max = mapped_max = gap_max = 0.0
    for case, saved, parameters in zip(result["cases"], reference["cases"], cases, strict=True):
        assert case["name"] == saved["name"] == parameters["name"]
        assert case["parameters"] == saved["parameters"] == parameters
        assert case["passed"] and case["bound_holds_numeric"]
        for key in ("r0", "Gamma"):
            close(case[key], saved[key], f"{case['name']}/{key}")
        for key in ("B_exact", "Delta", "D", "C", "weights", "harmonic_bound"):
            assert case["bound"][key] == saved["bound"][key]
        for name in ("energy_only", "hours_only", "joint", "homogeneous_auxiliary"):
            task = case[name]
            assert task["status"] == 0 and task["success"] and task["passed"]
            for key in ("r", "fun"):
                close(task[key], saved[name][key], f"{case['name']}/{name}/{key}")
            # HiGHS versions can omit MIP fields when presolve leaves an LP.
            for key in ("mip_gap", "objective_minus_dual_bound"):
                if task.get(key) is not None:
                    close(task[key], 0.0, f"{case['name']}/{name}/{key}")
                    gap_max = max(gap_max, abs(task[key]))
            a = primary.numeric(parameters)
            if name == "homogeneous_auxiliary":
                a["Pc"] *= a["eta"] / max(a["eta"])
                a["eta"][:] = max(a["eta"])
            values = task["values"]
            check = primary.physical(a, task["r"], *(
                np.asarray(values[k]) for k in ("f", "e", "c", "d", "loss", "spill")
            ))
            assert check["max_violation"] <= TOL
            if name != "hours_only":
                assert check["ENS"] <= float(F(parameters["Ebar"])) + TOL
            if name != "energy_only":
                assert check["LOLH_numeric"] <= parameters["Hbar"]
            physical_max = max(physical_max, check["max_violation"])
        witness = case["mapped_witness"]
        check = primary.physical(primary.numeric(parameters), witness["r"], *(
            np.asarray(witness[k]) for k in ("f", "e", "c", "d", "loss", "spill")
        ))
        assert witness["passed"] and check["max_violation"] <= TOL
        assert check["ENS"] <= float(F(parameters["Ebar"])) + TOL
        assert check["LOLH_numeric"] <= parameters["Hbar"]
        mapped_max = max(mapped_max, check["max_violation"])
    independent = fresh["independent_results.json"]
    saved = read(REFERENCE / "independent_results.json")
    assert all(independent["checks"].values())
    for name in ("energy_only", "hours_only", "joint"):
        task = independent["main_tasks"][name]
        assert task["success"] and task["physical"]["max_residual"] <= TOL
        close(task["r"], saved["main_tasks"][name]["r"], f"independent/{name}")
    for filename in ("exact_certificates.json", "family_checks.json"):
        assert fresh[filename] == read(REFERENCE / filename), filename
    assert result["sharp_family_exact_control"] == reference["sharp_family_exact_control"]
    four_hour = result["cases"][0]
    gamma = F(fresh["exact_certificates.json"]["appendix_A"]["r_star"])
    delta = F(four_hour["bound"]["Delta"])
    bound = F(four_hour["bound"]["B_exact"])
    assert gamma == F(13, 654) and bound == F(3, 38) and delta == F(3, 19)
    assert gamma / delta == F(247, 1962) and bound / delta == F(1, 2)
    return {
        "all_passed": True,
        "environment": result["environment"],
        "physical_optimizations": 60,
        "equipment_bound_LPs": 15,
        "predefined_cases": len(cases),
        "general_family_rational_cases": fresh["family_checks.json"]["family_cases"],
        "tutorial_rational_cases": fresh["family_checks.json"]["tutorial_cases"],
        "max_physical_residual": physical_max,
        "max_mapped_residual": mapped_max,
        "max_reported_solver_gap": gap_max,
        "comparison_tolerance": TOL,
        "four_hour_MW": {"Gamma": str(gamma), "B": str(bound), "Delta": str(delta)},
        "restricted_class_example_Gamma_over_Delta": str(gamma / delta),
    }


def run(script, arguments):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    process = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / script), *map(str, arguments)],
        cwd=ROOT, env=env, text=True, capture_output=True,
    )
    if process.returncode:
        raise RuntimeError(f"{script} failed:\n{process.stdout}\n{process.stderr}")
    print(f"Passed: {script}", flush=True)


def main():
    if not __debug__:
        raise SystemExit("Run without -O or PYTHONOPTIMIZE: assertions validate the results.")
    count = verify()
    print(f"Integrity check passed: {count} files.", flush=True)
    OUTPUT.mkdir(exist_ok=True)
    with TemporaryDirectory(prefix="run-", dir=OUTPUT) as temporary:
        work = Path(temporary)
        for script, filename in TASKS:
            run(script, ["--output", work / filename])
        summary = validate(work)
        run("plot_results.py", ["--results", work, "--output", work / "figures"])
        (work / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
        # Only complete validated results replace the previous generated run.
        for path in sorted(work.rglob("*")):
            if path.is_file() and not any(part.startswith(".") for part in path.relative_to(work).parts):
                target = OUTPUT / path.relative_to(work)
                target.parent.mkdir(parents=True, exist_ok=True)
                path.replace(target)
    verify()
    print("All numerical checks passed. Results: generated/summary.json; plots: generated/figures/.")


if __name__ == "__main__":
    main()
