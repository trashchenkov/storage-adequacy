#!/usr/bin/env python3
"""Generate publication numeric plots from saved reproducibility outputs."""

from __future__ import annotations

import argparse
import csv
from fractions import Fraction
import json
import os
from pathlib import Path


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_fraction(value: object) -> Fraction:
    if isinstance(value, int):
        return Fraction(value, 1)
    if isinstance(value, float):
        return Fraction(str(value))
    if isinstance(value, str):
        return Fraction(value)
    raise TypeError(f"cannot parse fraction from {value!r}")


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def import_matplotlib(output_dir: Path):
    os.environ["MPLCONFIGDIR"] = str(output_dir / ".matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.edgecolor": "0.25",
            "axes.labelcolor": "0.05",
            "xtick.color": "0.1",
            "ytick.color": "0.1",
        }
    )
    return plt


def curve_values(n: int) -> tuple[Fraction, Fraction, Fraction]:
    n_frac = Fraction(n, 1)
    gamma = Fraction(15 * n, 95 * n + 156)
    delta = Fraction(3, 19)
    bound_b = delta * n_frac / (n_frac + Fraction(4, 5))
    return gamma, delta, bound_b


def load_comparison(primary: dict, certificates: dict) -> list[dict[str, object]]:
    case = primary["cases"][0]
    bound = case["bound"]
    delta = parse_fraction(bound["Delta"])
    bound_b = parse_fraction(bound["B_exact"])
    gamma = parse_fraction(certificates["appendix_A"]["r_star"])
    return [
        {"label": "Simpler sufficient\nbound $\\Delta$", "metric": "Delta", "value": delta},
        {"label": "Refined sufficient\nbound $B$", "metric": "B", "value": bound_b},
        {"label": "Exact increment\n$\\Gamma$", "metric": "Gamma", "value": gamma},
    ]


def build_sharpness(certificates: dict) -> list[dict[str, object]]:
    rows = []
    for n in range(1, 257):
        gamma, delta, bound_b = curve_values(n)
        rows.append(
            {
                "N": n,
                "Gamma": gamma,
                "Delta": delta,
                "B": bound_b,
                "Gamma_over_Delta": gamma / delta,
                "Gamma_over_B": gamma / bound_b,
            }
        )

    certificate_rows = {int(row["N"]): row for row in certificates["appendix_B"]}
    for n in (1, 4, 16, 64):
        if n not in certificate_rows:
            raise ValueError(f"certificate row for N={n} is missing")
        gamma, delta, bound_b = curve_values(n)
        cert = certificate_rows[n]
        checks = {
            "Gamma": gamma,
            "Delta": delta,
            "B": bound_b,
            "Gamma_over_Delta": gamma / delta,
            "Gamma_over_B": gamma / bound_b,
        }
        for key, expected in checks.items():
            actual = parse_fraction(cert[key])
            if actual != expected:
                raise ValueError(
                    f"certificate mismatch for N={n} {key}: "
                    f"{fraction_text(actual)} != {fraction_text(expected)}"
                )
    return rows


def write_comparison_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "exact_MW", "value_MW", "value_kW"])
        writer.writeheader()
        for row in rows:
            value = row["value"]
            assert isinstance(value, Fraction)
            writer.writerow(
                {
                    "metric": row["metric"],
                    "exact_MW": fraction_text(value),
                    "value_MW": f"{float(value):.12g}",
                    "value_kW": f"{1000.0 * float(value):.12g}",
                }
            )


def write_sharpness_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = [
        "N",
        "Gamma",
        "Delta",
        "B",
        "Gamma_over_Delta",
        "Gamma_over_B",
        "Gamma_float",
        "B_float",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "N": row["N"],
                    "Gamma": fraction_text(row["Gamma"]),
                    "Delta": fraction_text(row["Delta"]),
                    "B": fraction_text(row["B"]),
                    "Gamma_over_Delta": fraction_text(row["Gamma_over_Delta"]),
                    "Gamma_over_B": fraction_text(row["Gamma_over_B"]),
                    "Gamma_float": f"{float(row['Gamma']):.12g}",
                    "B_float": f"{float(row['B']):.12g}",
                }
            )


def figure_comparison(plt, rows: list[dict[str, object]], output_dir: Path) -> None:
    labels = [str(row["label"]) for row in rows]
    values_kw = [1000.0 * float(row["value"]) for row in rows]
    colors = ["#64748b", "#0f766e", "#b91c1c"]

    fig, ax = plt.subplots(figsize=(6.6, 3.9), constrained_layout=True)
    bars = ax.bar(range(3), values_kw, color=colors, width=0.56)
    ax.set_xticks(range(3), labels)
    ax.set_ylabel("Power, kW")
    ax.set_ylim(0, 180)
    ax.grid(axis="y", color="#dde3e8", linewidth=0.8)
    ax.set_axisbelow(True)
    for bar, value in zip(bars, values_kw):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 4, f"{value:.1f}", ha="center", fontsize=9)
    ax.text(
        0.98,
        0.70,
        "Bounds are guaranteed margins;\n$\\Gamma$ is the actual optimum in this example.",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.3,
        bbox={"facecolor": "white", "edgecolor": "#d0d0d0", "boxstyle": "round,pad=.3"},
    )
    fig.savefig(output_dir / "comparison.png", dpi=300)
    fig.savefig(output_dir / "comparison.svg")
    plt.close(fig)


def figure_sharpness(plt, rows: list[dict[str, object]], output_dir: Path) -> None:
    xs = [row["N"] for row in rows]
    y_delta = [float(row["Gamma_over_Delta"]) for row in rows]
    y_b = [float(row["Gamma_over_B"]) for row in rows]
    marker_rows = [row for row in rows if row["N"] in {1, 4, 16, 64}]

    fig, ax = plt.subplots(figsize=(6.9, 4.3), constrained_layout=True)
    ax.plot(xs, y_b, color="#0f766e", lw=2.1, label="$\\Gamma_N/B_N$")
    ax.plot(xs, y_delta, color="#64748b", lw=1.8, ls="--", label="$\\Gamma_N/\\Delta$")
    ax.scatter(
        [row["N"] for row in marker_rows],
        [float(row["Gamma_over_B"]) for row in marker_rows],
        color="#b91c1c",
        edgecolor="white",
        zorder=3,
        s=44,
        label="certificate values",
    )
    offsets = {1: (5, 8), 4: (5, -13), 16: (5, 8), 64: (-21, -13)}
    for row in marker_rows:
        ax.annotate(
            str(row["N"]),
            (row["N"], float(row["Gamma_over_B"])),
            xytext=offsets[int(row["N"])],
            textcoords="offset points",
            fontsize=8,
        )
    ax.set_xscale("log", base=2)
    ax.set_xlim(0.85, 300)
    ax.set_ylim(0, 1.04)
    ax.set_xlabel("$N$ repeated deficit hours; $T=N+2$, equipment varies with $N$")
    ax.set_ylabel("ratio")
    ax.grid(True, which="both", color="#dde3e8", linewidth=0.75)
    ax.legend(frameon=False, loc="lower right")
    ax.text(
        0.02,
        0.29,
        "Changing-fleet family:\nratings, time series, and norm vary with $N$.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.3,
        bbox={"facecolor": "white", "edgecolor": "#d0d0d0", "boxstyle": "round,pad=.28"},
    )
    fig.savefig(output_dir / "sharpness.png", dpi=300)
    fig.savefig(output_dir / "sharpness.svg")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True, help="Directory containing saved JSON results.")
    parser.add_argument("--output", type=Path, required=True, help="Directory for plots, CSVs, and Matplotlib cache.")
    args = parser.parse_args()

    results_dir = args.results.resolve()
    output_dir = args.output.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    certificates = read_json(results_dir / "exact_certificates.json")
    primary = read_json(results_dir / "primary_results.json")

    comparison = load_comparison(primary, certificates)
    sharpness = build_sharpness(certificates)
    write_comparison_csv(output_dir / "comparison.csv", comparison)
    write_sharpness_csv(output_dir / "sharpness.csv", sharpness)

    plt = import_matplotlib(output_dir)
    figure_comparison(plt, comparison, output_dir)
    figure_sharpness(plt, sharpness, output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
