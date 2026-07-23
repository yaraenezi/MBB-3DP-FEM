"""Compare all accepted full-upright sparse-infill MBB FEM cases."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PATTERNS = ("cubic", "gyroid", "honeycomb")
DENSITIES = (15, 25)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--solid", type=Path, required=True)
    for pattern in PATTERNS:
        for density in DENSITIES:
            parser.add_argument(
                f"--{pattern}-{density}-fine",
                dest=f"{pattern}_{density}_fine",
                type=Path,
                required=True,
            )
            parser.add_argument(
                f"--{pattern}-{density}-coarse",
                dest=f"{pattern}_{density}_coarse",
                type=Path,
                required=True,
            )
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def stress_pa(result: dict) -> float:
    return float(
        result.get(
            "maximum_homogenized_von_mises_pa",
            result.get("maximum_von_mises_pa"),
        )
    )


def row(label: str, result: dict, mass_kg: float) -> dict:
    displacement_mm = -1000.0 * float(result["mean_load_displacement_m"])
    shell_core = result.get("shell_core", {})
    return {
        "case": label,
        "mesh_size_mm": float(result["mesh_size_mm"]),
        "tetrahedra": int(result["mesh"]["tetrahedra"]),
        "mean_load_displacement_mm": displacement_mm,
        "maximum_displacement_mm": 1000.0
        * float(result["maximum_displacement_m"]),
        "apparent_stiffness_n_per_mm": float(result["load_n"])
        / displacement_mm,
        "maximum_reported_von_mises_mpa": stress_pa(result) / 1.0e6,
        "estimated_mass_g": 1000.0 * mass_kg,
        "core_relative_stiffness": shell_core.get("core_relative_stiffness", 1.0),
        "equilibrium_relative_error": float(result["equilibrium_relative_error"]),
        "energy_identity_relative_error": float(
            result["energy_identity_relative_error"]
        ),
    }


def main() -> None:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    solid = load(args.solid)

    fine: dict[tuple[str, int], dict] = {}
    coarse: dict[tuple[str, int], dict] = {}
    for pattern in PATTERNS:
        for density in DENSITIES:
            fine[(pattern, density)] = load(
                getattr(args, f"{pattern}_{density}_fine")
            )
            coarse[(pattern, density)] = load(
                getattr(args, f"{pattern}_{density}_coarse")
            )

    accepted = []
    for density in DENSITIES:
        for pattern in PATTERNS:
            result = fine[(pattern, density)]
            accepted.append(
                row(
                    f"{pattern.title()} {density}%",
                    result,
                    result["estimated_mass_kg"],
                )
            )
    accepted.append(row("Dense 100%", solid, solid["mass_kg"]))

    convergence = {}
    for pattern in PATTERNS:
        for density in DENSITIES:
            fine_value = float(
                fine[(pattern, density)]["mean_load_displacement_m"]
            )
            coarse_value = float(
                coarse[(pattern, density)]["mean_load_displacement_m"]
            )
            convergence[f"{pattern}_{density}_mean_displacement_change_percent"] = (
                100.0 * abs(fine_value - coarse_value) / abs(fine_value)
            )
    convergence["coarse_mesh_size_mm"] = 3.0
    convergence["fine_mesh_size_mm"] = 2.4

    summary = {
        "geometry_policy": (
            "Every accepted case uses the same full upright original-mesh.stl."
        ),
        "load_n": float(solid["load_n"]),
        "screening_model": (
            "Volume-matched dense shell plus density-squared core stiffness, "
            "corrected by a G-code axial fourth-orientation moment normalized "
            "to cubic at each density."
        ),
        "accepted_fine_mesh_cases": accepted,
        "mesh_convergence": convergence,
        "interpretation": (
            "Pattern rankings are comparative homogenized screening results. "
            "Pattern-specific RVE or coupon calibration is required."
        ),
    }
    (args.output / "comparison.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    with (args.output / "comparison.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=accepted[0].keys())
        writer.writeheader()
        writer.writerows(accepted)

    labels = [item["case"].replace(" ", "\n", 1) for item in accepted]
    colors = [
        "#D55E00",
        "#009E73",
        "#0072B2",
        "#E69F00",
        "#56B4E9",
        "#CC79A7",
        "#4D4D4D",
    ]
    metrics = (
        ("mean_load_displacement_mm", "Mean load displacement", "mm"),
        ("apparent_stiffness_n_per_mm", "Apparent stiffness", "N/mm"),
        ("estimated_mass_g", "Estimated mass", "g"),
    )
    figure, axes = plt.subplots(1, 3, figsize=(15, 5.2), dpi=180)
    for axis, (key, title, unit) in zip(axes, metrics):
        values = [float(item[key]) for item in accepted]
        bars = axis.bar(labels, values, color=colors)
        axis.set_title(title)
        axis.set_ylabel(unit)
        axis.grid(axis="y", alpha=0.25)
        axis.tick_params(axis="x", labelsize=8)
        axis.bar_label(bars, fmt="%.2f", padding=3, fontsize=7)
    figure.suptitle(
        "Full upright MBB beam: six sparse-infill cases and dense baseline at 100 N"
    )
    figure.tight_layout()
    figure.savefig(args.output / "all_patterns_full_upright_comparison.png")
    plt.close(figure)


if __name__ == "__main__":
    main()
