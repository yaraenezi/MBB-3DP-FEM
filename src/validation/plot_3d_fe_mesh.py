from __future__ import annotations

import argparse
from pathlib import Path

import gmsh
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot the exterior facets of a Gmsh tetrahedral FE mesh."
    )
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--maximum-facets", type=int, default=32000)
    return parser.parse_args()


def load_exterior_facets(mesh_path: Path) -> tuple[np.ndarray, np.ndarray, int]:
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(mesh_path))
        node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
        _, tetrahedron_nodes = gmsh.model.mesh.getElementsByType(4)
    finally:
        gmsh.finalize()

    if len(tetrahedron_nodes) == 0:
        raise RuntimeError(f"No first-order tetrahedra found in {mesh_path}")

    coordinates = np.asarray(coordinates, dtype=float).reshape((-1, 3))
    node_tags = np.asarray(node_tags, dtype=np.int64)
    order = np.argsort(node_tags)
    sorted_tags = node_tags[order]
    sorted_coordinates = coordinates[order]

    tetrahedra = np.asarray(tetrahedron_nodes, dtype=np.int64).reshape((-1, 4))
    faces = np.concatenate(
        (
            tetrahedra[:, (0, 1, 2)],
            tetrahedra[:, (0, 1, 3)],
            tetrahedra[:, (0, 2, 3)],
            tetrahedra[:, (1, 2, 3)],
        )
    )
    faces.sort(axis=1)
    unique_faces, counts = np.unique(faces, axis=0, return_counts=True)
    exterior_tags = unique_faces[counts == 1]
    exterior_indices = np.searchsorted(sorted_tags, exterior_tags)

    if not np.array_equal(sorted_tags[exterior_indices], exterior_tags):
        raise RuntimeError("Mesh connectivity references unknown node tags.")
    return sorted_coordinates, exterior_indices, tetrahedra.shape[0]


def plot_mesh(
    coordinates: np.ndarray,
    exterior_facets: np.ndarray,
    tetrahedron_count: int,
    output: Path,
    title: str,
    maximum_facets: int,
) -> None:
    if exterior_facets.shape[0] > maximum_facets:
        sample = np.linspace(
            0, exterior_facets.shape[0] - 1, maximum_facets, dtype=int
        )
        plotted_facets = exterior_facets[sample]
    else:
        plotted_facets = exterior_facets

    points_mm = coordinates - coordinates.min(axis=0)
    triangles = points_mm[plotted_facets]

    figure = plt.figure(figsize=(14, 7), dpi=180)
    axes = figure.add_subplot(111, projection="3d")
    collection = Poly3DCollection(
        triangles,
        facecolor="#DDEAF2",
        edgecolor="#315A70",
        linewidth=0.08,
        alpha=0.88,
    )
    axes.add_collection3d(collection)
    axes.set_xlim(points_mm[:, 0].min(), points_mm[:, 0].max())
    axes.set_ylim(points_mm[:, 1].min(), points_mm[:, 1].max())
    axes.set_zlim(points_mm[:, 2].min(), points_mm[:, 2].max())
    axes.set_box_aspect(tuple(np.ptp(points_mm, axis=0)))
    axes.view_init(elev=22, azim=-62)
    axes.set_axis_off()
    figure.suptitle(
        f"{title}\n"
        f"{tetrahedron_count:,} first-order tetrahedra | "
        f"{exterior_facets.shape[0]:,} exterior facets",
        fontsize=17,
        y=0.93,
    )
    figure.subplots_adjust(left=0.01, right=0.99, bottom=0.02, top=0.86)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, facecolor="white")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    if not args.mesh.exists():
        raise FileNotFoundError(args.mesh)
    coordinates, exterior_facets, tetrahedron_count = load_exterior_facets(args.mesh)
    plot_mesh(
        coordinates,
        exterior_facets,
        tetrahedron_count,
        args.output,
        args.title,
        args.maximum_facets,
    )


if __name__ == "__main__":
    main()
