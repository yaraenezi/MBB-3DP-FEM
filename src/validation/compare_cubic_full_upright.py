"""Compare accepted full-upright cubic MBB FEM cases."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--solid", type=Path, required=True)
    parser.add_argument("--cubic-15-fine", type=Path, required=True)
    parser.add_argument("--cubic-15-coarse", type=Path, required=True)
    parser.add_argument("--cubic-25-fine", type=Path, required=True)
    parser.add_argument("--cubic-25-coarse", type=Path, required=True)
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
        "equilibrium_relative_error": float(result["equilibrium_relative_error"]),
        "energy_identity_relative_error": float(
            result["energy_identity_relative_error"]
        ),
    }


def main() -> None:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    solid = load(args.solid)
    fifteen_fine = load(args.cubic_15_fine)
    fifteen_coarse = load(args.cubic_15_coarse)
    twenty_five_fine = load(args.cubic_25_fine)
    twenty_five_coarse = load(args.cubic_25_coarse)

    accepted = [
        row("Cubic 15%", fifteen_fine, fifteen_fine["estimated_mass_kg"]),
        row("Cubic 25%", twenty_five_fine, twenty_five_fine["estimated_mass_kg"]),
        row("Dense 100%", solid, solid["mass_kg"]),
    ]
    convergence = {
        "cubic_15_mean_displacement_change_percent": 100.0
        * abs(
            fifteen_fine["mean_load_displacement_m"]
            - fifteen_coarse["mean_load_displacement_m"]
        )
        / abs(fifteen_fine["mean_load_displacement_m"]),
        "cubic_25_mean_displacement_change_percent": 100.0
        * abs(
            twenty_five_fine["mean_load_displacement_m"]
            - twenty_five_coarse["mean_load_displacement_m"]
        )
        / abs(twenty_five_fine["mean_load_displacement_m"]),
        "coarse_mesh_size_mm": float(fifteen_coarse["mesh_size_mm"]),
        "fine_mesh_size_mm": float(fifteen_fine["mesh_size_mm"]),
    }
    summary = {
        "geometry_policy": (
            "Every accepted case uses the same full upright original-mesh.stl."
        ),
        "load_n": float(solid["load_n"]),
        "accepted_fine_mesh_cases": accepted,
        "mesh_convergence": convergence,
        "interpretation": (
            "Sparse cases are homogenized shell/core screening models. "
            "Coupon or RVE calibration is required before design allowables."
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

    labels = [item["case"] for item in accepted]
    displacement = [item["mean_load_displacement_mm"] for item in accepted]
    stiffness = [item["apparent_stiffness_n_per_mm"] for item in accepted]
    mass = [item["estimated_mass_g"] for item in accepted]
    colors = ["#d55e00", "#0072b2", "#4d4d4d"]
    figure, axes = plt.subplots(1, 3, figsize=(13, 4.5), dpi=180)
    for axis, values, title, unit in zip(
        axes,
        (displacement, stiffness, mass),
        ("Mean load displacement", "Apparent stiffness", "Estimated mass"),
        ("mm", "N/mm", "g"),
    ):
        bars = axis.bar(labels, values, color=colors)
        axis.set_title(title)
        axis.set_ylabel(unit)
        axis.grid(axis="y", alpha=0.25)
        axis.bar_label(bars, fmt="%.2f", padding=3, fontsize=9)
    figure.suptitle("Full upright MBB beam: cubic infill screening at 100 N")
    figure.tight_layout()
    figure.savefig(args.output / "cubic_full_upright_comparison.png")
    plt.close(figure)


if __name__ == "__main__":
    main()
