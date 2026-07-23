"""Inspect Bambu Studio .gcode.3mf archives and quantify cubic infill paths."""

from __future__ import annotations

import argparse
import json
import math
import re
import zipfile
from collections import defaultdict
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


SETTING_KEYS = (
    "sparse_infill_density",
    "sparse_infill_pattern",
    "wall_loops",
    "outer_wall_line_width",
    "inner_wall_line_width",
    "layer_height",
    "nozzle_diameter",
    "top_shell_layers",
    "bottom_shell_layers",
    "top_shell_thickness",
    "bottom_shell_thickness",
    "infill_direction",
    "infill_rotate_step",
    "infill_shift_step",
    "enable_support",
    "support_type",
    "filament_type",
)

MOVE_VALUE = re.compile(r"([XYZE])(-?(?:\d+(?:\.\d*)?|\.\d+))")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("archives", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def cluster_values(values: list[float], tolerance: float) -> list[float]:
    if not values:
        return []
    ordered = sorted(values)
    clusters: list[list[float]] = [[ordered[0]]]
    for value in ordered[1:]:
        if value - float(np.mean(clusters[-1])) <= tolerance:
            clusters[-1].append(value)
        else:
            clusters.append([value])
    return [float(np.mean(cluster)) for cluster in clusters]


def inspect_archive(path: Path) -> tuple[dict[str, object], dict[float, list[dict[str, float]]]]:
    settings: dict[str, str] = {}
    header: dict[str, str] = {}
    segments_by_z: dict[float, list[dict[str, float]]] = defaultdict(list)
    x = y = z = e_position = 0.0
    xyz_absolute = True
    extrusion_absolute = False
    feature = ""

    with zipfile.ZipFile(path) as archive:
        plate = json.loads(archive.read("Metadata/plate_1.json"))
        with archive.open("Metadata/plate_1.gcode") as raw:
            for raw_line in raw:
                line = raw_line.decode("utf-8", errors="replace").strip()
                if line.startswith("; "):
                    body = line[2:]
                    if " = " in body:
                        key, value = body.split(" = ", 1)
                        if key in SETTING_KEYS:
                            settings[key] = value
                    elif ": " in body:
                        key, value = body.split(": ", 1)
                        if key in {
                            "model printing time",
                            "total filament weight [g]",
                            "max_z_height",
                            "Z_HEIGHT",
                        }:
                            if key == "Z_HEIGHT":
                                z = float(value)
                            else:
                                header[key] = value
                if line.startswith("; FEATURE:"):
                    feature = line.split(":", 1)[1].strip()
                    continue
                if line == "G90":
                    xyz_absolute = True
                    continue
                if line == "G91":
                    xyz_absolute = False
                    continue
                if line == "M82":
                    extrusion_absolute = True
                    continue
                if line == "M83":
                    extrusion_absolute = False
                    continue
                if line.startswith("G92"):
                    values = {key: float(value) for key, value in MOVE_VALUE.findall(line)}
                    if "E" in values:
                        e_position = values["E"]
                    continue
                if not line.startswith(("G0 ", "G1 ")):
                    continue

                values = {key: float(value) for key, value in MOVE_VALUE.findall(line)}
                next_x = values.get("X", x)
                next_y = values.get("Y", y)
                next_z = values.get("Z", z)
                if not xyz_absolute:
                    next_x = x + values.get("X", 0.0)
                    next_y = y + values.get("Y", 0.0)
                    next_z = z + values.get("Z", 0.0)

                extrusion_delta = 0.0
                if "E" in values:
                    if extrusion_absolute:
                        extrusion_delta = values["E"] - e_position
                        e_position = values["E"]
                    else:
                        extrusion_delta = values["E"]
                        e_position += values["E"]

                dx = next_x - x
                dy = next_y - y
                length = math.hypot(dx, dy)
                if (
                    feature == "Sparse infill"
                    and extrusion_delta > 1.0e-7
                    and length >= 0.8
                ):
                    angle = math.degrees(math.atan2(dy, dx)) % 180.0
                    angle_bin = (round(angle / 5.0) * 5.0) % 180.0
                    segments_by_z[round(z, 6)].append(
                        {
                            "x0": x,
                            "y0": y,
                            "x1": next_x,
                            "y1": next_y,
                            "length": length,
                            "angle_deg": angle,
                            "angle_bin_deg": angle_bin,
                        }
                    )
                x, y, z = next_x, next_y, next_z

    orientation_weights: dict[float, float] = defaultdict(float)
    layer_dominant: list[dict[str, float]] = []
    spacing_candidates: list[float] = []
    for layer_z, segments in sorted(segments_by_z.items()):
        layer_weights: dict[float, float] = defaultdict(float)
        grouped_rho: dict[float, list[float]] = defaultdict(list)
        for segment in segments:
            angle_bin = segment["angle_bin_deg"]
            orientation_weights[angle_bin] += segment["length"]
            layer_weights[angle_bin] += segment["length"]
            theta = math.radians(angle_bin)
            midpoint_x = 0.5 * (segment["x0"] + segment["x1"])
            midpoint_y = 0.5 * (segment["y0"] + segment["y1"])
            rho = -math.sin(theta) * midpoint_x + math.cos(theta) * midpoint_y
            grouped_rho[angle_bin].append(rho)

        dominant_angle = max(layer_weights, key=layer_weights.get)
        layer_dominant.append(
            {
                "z_mm": layer_z,
                "dominant_angle_deg": dominant_angle,
                "sparse_path_length_mm": sum(layer_weights.values()),
            }
        )
        for rho_values in grouped_rho.values():
            centers = cluster_values(rho_values, tolerance=0.7)
            differences = np.diff(centers)
            spacing_candidates.extend(
                float(value) for value in differences if 1.0 <= value <= 40.0
            )

    if spacing_candidates:
        rounded = np.round(np.asarray(spacing_candidates) * 2.0) / 2.0
        unique, counts = np.unique(rounded, return_counts=True)
        pitch = float(unique[np.argmax(counts)])
        pitch_median = float(np.median(spacing_candidates))
    else:
        pitch = float("nan")
        pitch_median = float("nan")

    sorted_orientations = sorted(
        orientation_weights.items(), key=lambda item: item[1], reverse=True
    )
    result: dict[str, object] = {
        "archive": path.name,
        "settings": settings,
        "header": header,
        "plate_bbox_mm": plate.get("bbox_all"),
        "object_bbox": plate.get("bbox_objects"),
        "sparse_layer_count": len(segments_by_z),
        "sparse_segment_count": sum(len(value) for value in segments_by_z.values()),
        "orientation_length_fraction": {
            f"{angle:g}": weight / sum(orientation_weights.values())
            for angle, weight in sorted_orientations[:8]
        },
        "estimated_in_plane_line_pitch_mm": pitch,
        "median_spacing_candidate_mm": pitch_median,
        "layer_dominant_orientations": layer_dominant,
    }
    return result, segments_by_z


def plot_layers(
    archive_name: str,
    segments_by_z: dict[float, list[dict[str, float]]],
    output: Path,
) -> None:
    populated = sorted(segments_by_z)
    if not populated:
        return
    indices = np.linspace(0, len(populated) - 1, 4, dtype=int)
    selected = [populated[index] for index in indices]
    figure, axes = plt.subplots(2, 2, figsize=(13, 8), dpi=170)
    for axis, layer_z in zip(axes.flat, selected):
        segments = segments_by_z[layer_z]
        lines = [
            [(item["x0"], item["y0"]), (item["x1"], item["y1"])]
            for item in segments
        ]
        angles = np.asarray([item["angle_bin_deg"] for item in segments])
        collection = LineCollection(lines, array=angles, cmap="twilight", linewidths=1.2)
        axis.add_collection(collection)
        axis.autoscale()
        axis.set_aspect("equal")
        axis.set_title(f"Sparse infill at Z={layer_z:.2f} mm")
        axis.set_xlabel("Printer X (mm)")
        axis.set_ylabel("Printer Y (mm)")
    figure.suptitle(f"Cubic toolpath verification\n{archive_name}", fontsize=14)
    figure.tight_layout()
    figure.savefig(output)
    plt.close(figure)


def main() -> None:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    reports: list[dict[str, object]] = []
    for archive_path in args.archives:
        result, segments = inspect_archive(archive_path)
        reports.append(result)
        stem = archive_path.name.replace(" ", "_").replace(".gcode.3mf", "")
        plot_layers(
            archive_path.name,
            segments,
            args.output / f"{stem}_cubic_layers.png",
        )
    (args.output / "cubic_gcode_validation.json").write_text(
        json.dumps(reports, indent=2), encoding="utf-8"
    )
    print(json.dumps(reports, indent=2))


if __name__ == "__main__":
    main()
