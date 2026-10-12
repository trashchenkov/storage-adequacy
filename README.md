# Joint ENS and LOLH limits in storage systems

Reproducibility code for **Additional firm capacity for joint energy-not-served and loss-of-load-hours limits in storage systems**, by Sergei Trashchenkov.

Source repository: [https://github.com/trashchenkov/storage-adequacy](https://github.com/trashchenkov/storage-adequacy).

The code reproduces the analytical examples and numerical checks of the manuscript: separate and joint minimum source capacities, equipment-dependent bounds, feasible schedules, and the sharpness family. [Русская инструкция](README.ru.md).

## Quick start

Use Python 3.12. From this directory:

```bash
python3 verify_repository.py
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python reproduce.py
```

On Windows, activate the environment with `.venv\Scripts\activate` instead. SciPy supplies the HiGHS solver; a commercial optimization solver is unnecessary. The pinned versions and test evidence are recorded in [VALIDATION.md](VALIDATION.md).

A successful run prints `All numerical checks passed` and creates:

- `generated/summary.json`: verification counts, maximum residuals, solver gaps and software versions;
- `generated/primary_results.json`: all primary optimizations and feasible schedules;
- `generated/independent_results.json`: independently formulated targeted LP checks;
- `generated/exact_certificates.json` and `generated/family_checks.json`: exact rational checks;
- `generated/figures/`: comparison and sharpness plots in PNG/SVG, with CSV data.

Repeat `python reproduce.py` to recalculate. Reference files remain unchanged. A failed calculation exits with a nonzero status; it does not replace the previous successful generated results. Run Python normally: `-O` and `PYTHONOPTIMIZE` disable assertions and are rejected.

## Contents and manuscript correspondence

| Material | Files | Result reproduced |
| --- | --- | --- |
| Predefined input cases | `data/cases.json` | Fifteen cases from Section 6.3 |
| Primary physical model | `scripts/check_theory.py` | 45 original-model and 15 auxiliary optimizations; 15 equipment-bound LPs; return-map feasibility |
| Independent formulation | `scripts/independent_check.py` | Four-hour optimum, auxiliary transformation, and targeted bound checks |
| Exact certificates | `scripts/exact_certificates.py` | Four-hour example and five selected sharpness-family members |
| General-family checks | `scripts/check_family.py` | 360 exact feasible-schedule checks and the flow-augmentation example |
| Numerical plots | `scripts/plot_results.py` | Four-hour comparison (Figure 4) and sharpness ratios (Figure 5), including the values in Tables 2 and 3 |
| Saved reference outputs | `results/reference/` | Archived results used for comparison |

The theoretical diagrams in Figures 1–3 and the manuscript typesetting are outside this code repository. Numerical checks support the examples; the general results are proved in the manuscript.

## Model and data conventions

The model uses one aggregate energy balance and a known finite time series. `x[t]` is generation minus demand before the additional source and storage. For each storage unit, `E`, `Pc`, `Pd` and `eta` specify usable energy capacity, charging power, discharging power and efficiency. The source has installed power `r` and dispatch `0 <= f[t] <= r`. Each time step lasts one hour; power is in MW and energy in MWh. Decimal and fractional input strings are parsed as rational numbers before conversion for optimization.

Storage energy satisfies `e[i,t+1] = e[i,t] + eta[i]*c[i,t] - d[i,t]`. Initial energy is freely chosen and must equal final energy for each unit. Charging uses external surplus; storage-to-storage transfers are excluded. The source is continuously available, with no fuel or ramp limit. Excess supply can be curtailed. The model describes energy adequacy without a spatial electrical-network model.

`Ebar` limits energy not served (ENS); `Hbar` limits the number of shortage hours (LOLH). `r_E` and `r_H` are the separate minimum source capacities, `r_joint` the joint minimum, and `Gamma = r_joint - max(r_E, r_H)`. The bounds `B` and `Delta` are sufficient additional capacities. The four-hour example gives:

```text
Gamma = 13/654 MW       = 19.877676 kW
B     = 3/38 MW         = 78.947368 kW
Delta = 3/19 MW         = 157.894737 kW
Gamma/Delta = 247/1962  = 0.125892
```

The stored primary formulation retains binary charging/discharging modes. At fixed source capacity, the normalization lemma permits a reduced formulation without those mode binaries; that computational reduction is explained in the manuscript. These scripts reproduce the original calculations.

The primary solver uses relative MIP gap `1e-10`, a 30-second limit per solve, and physical-check tolerance `5e-7`. LOLH is counted numerically using shortage power greater than `1e-7` MW (the independent checker uses `1e-8`). Mathematically, LOLH counts strictly positive shortage. The runner compares optimal values and rechecks physical feasibility, allowing different schedules at a nonunique optimum. A missing MIP diagnostic is accepted only when the solver reports optimal success; all available gaps are checked. Raw solver diagnostics are retained.

## Running individual checks

Each script accepts a destination file:

```bash
python scripts/check_theory.py --output generated/custom_primary.json
python scripts/independent_check.py --output generated/custom_independent.json
python scripts/exact_certificates.py --output generated/custom_exact.json
python scripts/check_family.py --output generated/custom_family.json
```

The primary script requires a new output filename. The other three replace their specified output file. The two exact-rational scripts use only the Python standard library. The full runner always verifies the delivered files first and writes only inside `generated/`.

## Provenance, integrity and citation

[PROVENANCE.md](PROVENANCE.md) identifies the implementations and packaging changes. `MANIFEST.json` records SHA-256 hashes of the delivered files, excluding itself and generated outputs. This checks accidental changes, not publisher authentication. Editing a delivered file makes the frozen integrity check fail; use individual scripts while developing a modified version.

Citation metadata are in [CITATION.cff](CITATION.cff). The Zenodo DOI and manuscript identifier will be added after registration. The code, included inputs and accompanying documentation are distributed under the [MIT license](LICENSE).
