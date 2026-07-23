"""True 3D linear-elastic FEM solve of a dense topology-optimized MBB STL.

At 100 percent infill, a cubic slicer pattern has no internal void phase and
is represented by dense isotropic PLA. The source STL is tetrahedralized with
Gmsh and solved with DOLFINx using three displacement components.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime
from pathlib import Path

import gmsh
import matplotlib
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import ufl
from dolfinx import fem, io, mesh
from dolfinx.fem.petsc import LinearProblem, assemble_vector
from dolfinx.io import gmsh as dolfinx_gmsh


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stl", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, default=Path("outputs/fea"))
    parser.add_argument("--mesh-size-mm", type=float, default=3.0)
    parser.add_argument("--load-n", type=float, default=100.0)
    parser.add_argument("--youngs-modulus-pa", type=float, default=3.0e9)
    parser.add_argument("--poisson-ratio", type=float, default=0.35)
    parser.add_argument("--density-kg-m3", type=float, default=1240.0)
    parser.add_argument("--plot-scale", type=float, default=150.0)
    return parser.parse_args()


def make_run_directory(root: Path) -> Path:
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S") + "_mbb_cubic_100"
    run_dir = root / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    return run_dir


def tetrahedralize_stl(
    stl_path: Path, msh_path: Path, mesh_size_mm: float
) -> dict[str, int]:
    gmsh.initialize()
    try:
        gmsh.logger.start()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("topology_mbb_dense")
        gmsh.merge(str(stl_path))

        angle = math.radians(40.0)
        gmsh.model.mesh.classifySurfaces(angle, True, True, math.pi)
        gmsh.model.mesh.createGeometry()
        surfaces = gmsh.model.getEntities(2)
        if not surfaces:
            raise RuntimeError("Gmsh did not create any surfaces from the STL.")

        surface_loop = gmsh.model.geo.addSurfaceLoop([tag for _, tag in surfaces])
        volume = gmsh.model.geo.addVolume([surface_loop])
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(3, [volume], tag=1)
        gmsh.model.setPhysicalName(3, 1, "dense_mbb")

        gmsh.option.setNumber("Mesh.MeshSizeMin", 0.55 * mesh_size_mm)
        gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size_mm)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 12)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)
        gmsh.option.setNumber("Mesh.Optimize", 1)
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize("Netgen")

        element_types, element_tags, _ = gmsh.model.mesh.getElements(3)
        tetrahedra = sum(len(tags) for tags in element_tags)
        if not element_types or tetrahedra == 0:
            raise RuntimeError("Gmsh produced no volume elements.")
        gmsh.write(str(msh_path))
        messages = gmsh.logger.get()
        (msh_path.parent / "gmsh.log").write_text(
            "\n".join(messages) + "\n", encoding="utf-8"
        )
        return {
            "logged_warnings": sum(
                "Warning" in message or "invalid" in message.lower()
                for message in messages
            ),
            "generated_volume_elements_before_export": int(tetrahedra),
        }
    finally:
        gmsh.finalize()


def scalar_allreduce(domain: mesh.Mesh, value: float) -> float:
    return float(domain.comm.allreduce(value, op=MPI.SUM))


def locate_boundary_regions(
    domain: mesh.Mesh,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    coordinates = domain.geometry.x
    local_min = coordinates.min(axis=0)
    local_max = coordinates.max(axis=0)
    bounds_min = np.array(
        [domain.comm.allreduce(v, op=MPI.MIN) for v in local_min], dtype=float
    )
    bounds_max = np.array(
        [domain.comm.allreduce(v, op=MPI.MAX) for v in local_max], dtype=float
    )
    extents = bounds_max - bounds_min

    # STL convention: x=span, y=depth, z=vertical.
    span = extents[0]
    depth = extents[1]
    height = extents[2]
    support_width = max(0.004, 0.035 * span)
    support_height = max(0.004, 0.14 * height)
    load_half_width = max(0.004, 0.035 * span)
    load_half_depth = max(0.003, 0.20 * depth)
    load_height_band = max(0.0025, 0.08 * height)
    tolerance = 1.0e-9

    fdim = domain.topology.dim - 1
    left_facets = mesh.locate_entities_boundary(
        domain,
        fdim,
        lambda x: (
            (x[0] <= bounds_min[0] + support_width + tolerance)
            & (x[2] <= bounds_min[2] + support_height + tolerance)
        ),
    )
    right_facets = mesh.locate_entities_boundary(
        domain,
        fdim,
        lambda x: (
            (x[0] >= bounds_max[0] - support_width - tolerance)
            & (x[2] <= bounds_min[2] + support_height + tolerance)
        ),
    )
    load_facets = mesh.locate_entities_boundary(
        domain,
        fdim,
        lambda x: (
            (np.abs(x[0] - 0.5 * (bounds_min[0] + bounds_max[0])) <= load_half_width)
            & (np.abs(x[1] - 0.5 * (bounds_min[1] + bounds_max[1])) <= load_half_depth)
            & (x[2] >= bounds_max[2] - load_height_band - tolerance)
        ),
    )

    if min(len(left_facets), len(right_facets), len(load_facets)) == 0:
        raise RuntimeError(
            "A support or load region is empty: "
            f"left={len(left_facets)}, right={len(right_facets)}, "
            f"load={len(load_facets)}."
        )

    geometry = {
        "span_m": float(span),
        "depth_m": float(depth),
        "height_m": float(height),
        "support_width_m": support_width,
        "support_height_m": support_height,
        "load_half_width_m": load_half_width,
        "load_half_depth_m": load_half_depth,
        "load_height_band_m": load_height_band,
    }
    return left_facets, right_facets, load_facets, geometry


def create_facet_measure(
    domain: mesh.Mesh, facets: np.ndarray, tag: int
) -> ufl.Measure:
    sorted_facets = np.sort(np.unique(facets)).astype(np.int32)
    values = np.full(sorted_facets.size, tag, dtype=np.int32)
    tags = mesh.meshtags(
        domain, domain.topology.dim - 1, sorted_facets, values
    )
    return ufl.Measure("ds", domain=domain, subdomain_data=tags)


def component_dofs(
    vector_space: fem.FunctionSpace, facets: np.ndarray, component: int
) -> np.ndarray:
    return fem.locate_dofs_topological(
        vector_space.sub(component),
        vector_space.mesh.topology.dim - 1,
        facets,
    )


def mesh_quality(domain: mesh.Mesh) -> dict[str, float | int]:
    tdim = domain.topology.dim
    domain.topology.create_connectivity(tdim, 0)
    cells_to_vertices = domain.topology.connectivity(tdim, 0)
    index_map = domain.topology.index_map(tdim)
    cell_count = index_map.size_local
    coordinates = domain.geometry.x
    volumes = np.empty(cell_count, dtype=float)
    for cell in range(cell_count):
        vertices = cells_to_vertices.links(cell)
        points = coordinates[vertices]
        matrix = np.column_stack(
            (points[1] - points[0], points[2] - points[0], points[3] - points[0])
        )
        volumes[cell] = abs(np.linalg.det(matrix)) / 6.0
    degeneracy_tolerance = np.finfo(float).eps
    return {
        "tetrahedra": int(domain.comm.allreduce(cell_count, op=MPI.SUM)),
        "minimum_tetra_volume_m3": float(
            domain.comm.allreduce(float(volumes.min()), op=MPI.MIN)
        ),
        "maximum_tetra_volume_m3": float(
            domain.comm.allreduce(float(volumes.max()), op=MPI.MAX)
        ),
        "degenerate_tetrahedra": int(
            domain.comm.allreduce(
                int(np.count_nonzero(volumes <= degeneracy_tolerance)), op=MPI.SUM
            )
        ),
    }


def create_plot(
    domain: mesh.Mesh,
    displacement: fem.Function,
    output: Path,
    deformation_scale: float,
    title: str | None = None,
) -> None:
    if domain.comm.size != 1:
        return
    tdim = domain.topology.dim
    fdim = tdim - 1
    domain.topology.create_connectivity(fdim, tdim)
    domain.topology.create_connectivity(fdim, 0)
    exterior = mesh.exterior_facet_indices(domain.topology)
    facets_to_vertices = domain.topology.connectivity(fdim, 0)
    triangles = np.vstack([facets_to_vertices.links(facet) for facet in exterior])

    coordinates = domain.geometry.x.copy()
    nodal_displacement = displacement.x.array.reshape((-1, 3))
    if nodal_displacement.shape[0] != coordinates.shape[0]:
        raise RuntimeError("P1 displacement nodes do not match geometry nodes.")
    magnitude = np.linalg.norm(nodal_displacement, axis=1)
    deformed = coordinates + deformation_scale * nodal_displacement
    deformed_mm = (deformed - coordinates.min(axis=0)) * 1000.0

    maximum_triangles = 45000
    if triangles.shape[0] > maximum_triangles:
        indices = np.linspace(
            0, triangles.shape[0] - 1, maximum_triangles, dtype=int
        )
        triangles = triangles[indices]

    values = magnitude[triangles].mean(axis=1) * 1000.0
    normalization = plt.Normalize(vmin=float(values.min()), vmax=float(values.max()))
    colors = plt.cm.viridis(normalization(values))

    figure = plt.figure(figsize=(14, 8), dpi=180)
    axes = figure.add_subplot(111, projection="3d")
    collection = Poly3DCollection(
        deformed_mm[triangles],
        facecolors=colors,
        edgecolors="none",
        linewidths=0.0,
        alpha=1.0,
    )
    axes.add_collection3d(collection)
    axes.set_xlim(deformed_mm[:, 0].min(), deformed_mm[:, 0].max())
    axes.set_ylim(deformed_mm[:, 1].min(), deformed_mm[:, 1].max())
    axes.set_zlim(deformed_mm[:, 2].min(), deformed_mm[:, 2].max())
    axes.set_box_aspect(
        (
            np.ptp(deformed_mm[:, 0]),
            np.ptp(deformed_mm[:, 1]),
            np.ptp(deformed_mm[:, 2]),
        )
    )
    axes.view_init(elev=22, azim=-62)
    axes.set_axis_off()
    figure.suptitle(
        title
        or (
            "Full upright dense MBB true 3D FEM\n"
            f"100% cubic-equivalent solid | deformation x{deformation_scale:g}"
        ),
        fontsize=18,
        y=0.94,
    )
    scalar_map = plt.cm.ScalarMappable(norm=normalization, cmap="viridis")
    scalar_map.set_array([])
    color_axis = figure.add_axes((0.88, 0.20, 0.025, 0.60))
    colorbar = figure.colorbar(scalar_map, cax=color_axis)
    colorbar.set_label("Displacement magnitude (mm)")
    figure.subplots_adjust(left=0.01, right=0.84, bottom=0.03, top=0.88)
    figure.savefig(output)
    plt.close(figure)


def main() -> None:
    args = parse_args()
    if MPI.COMM_WORLD.size != 1:
        raise RuntimeError("This validation script currently requires one MPI rank.")
    if not args.stl.exists():
        raise FileNotFoundError(args.stl)
    if args.mesh_size_mm <= 0.0 or args.load_n <= 0.0:
        raise ValueError("Mesh size and load must be positive.")

    run_dir = make_run_directory(args.output_root)
    msh_path = run_dir / "dense_mbb_tetrahedra.msh"
    gmsh_report = tetrahedralize_stl(args.stl, msh_path, args.mesh_size_mm)

    mesh_data = dolfinx_gmsh.read_from_msh(
        msh_path, MPI.COMM_WORLD, rank=0, gdim=3
    )
    domain = mesh_data.mesh
    domain.name = "dense_mbb"
    domain.geometry.x[:] *= 0.001  # STL/Gmsh millimetres to SI metres.

    left_facets, right_facets, load_facets, geometry = locate_boundary_regions(domain)
    ds_load = create_facet_measure(domain, load_facets, tag=11)
    load_area = scalar_allreduce(
        domain, fem.assemble_scalar(fem.form(1.0 * ds_load(11)))
    )
    if load_area <= 0.0:
        raise RuntimeError("The integrated load area is zero.")

    vector_space = fem.functionspace(domain, ("Lagrange", 1, (3,)))
    zero_vector = np.zeros(3, dtype=PETSc.ScalarType)
    left_dofs = fem.locate_dofs_topological(
        vector_space, domain.topology.dim - 1, left_facets
    )
    left_bc = fem.dirichletbc(zero_vector, left_dofs, vector_space)
    right_depth_dofs = component_dofs(vector_space, right_facets, 1)
    right_vertical_dofs = component_dofs(vector_space, right_facets, 2)
    right_depth_bc = fem.dirichletbc(
        PETSc.ScalarType(0.0), right_depth_dofs, vector_space.sub(1)
    )
    right_vertical_bc = fem.dirichletbc(
        PETSc.ScalarType(0.0), right_vertical_dofs, vector_space.sub(2)
    )
    boundary_conditions = [left_bc, right_depth_bc, right_vertical_bc]

    traction = fem.Constant(
        domain,
        np.array(
            [0.0, 0.0, -args.load_n / load_area],
            dtype=PETSc.ScalarType,
        ),
    )
    trial = ufl.TrialFunction(vector_space)
    test = ufl.TestFunction(vector_space)
    lame_lambda = (
        args.youngs_modulus_pa
        * args.poisson_ratio
        / ((1.0 + args.poisson_ratio) * (1.0 - 2.0 * args.poisson_ratio))
    )
    lame_mu = args.youngs_modulus_pa / (2.0 * (1.0 + args.poisson_ratio))

    def strain(field):
        return ufl.sym(ufl.grad(field))

    def stress(field):
        epsilon = strain(field)
        return (
            lame_lambda * ufl.tr(epsilon) * ufl.Identity(3)
            + 2.0 * lame_mu * epsilon
        )

    bilinear = ufl.inner(stress(trial), strain(test)) * ufl.dx
    linear = ufl.dot(traction, test) * ds_load(11)
    problem = LinearProblem(
        bilinear,
        linear,
        bcs=boundary_conditions,
        petsc_options_prefix="dense_mbb_3d_",
        petsc_options={
            "ksp_type": "cg",
            "pc_type": "gamg",
            "ksp_rtol": 1.0e-10,
            "ksp_atol": 1.0e-12,
            "ksp_max_it": 10000,
            "ksp_error_if_not_converged": True,
        },
    )
    displacement = problem.solve()
    displacement.name = "displacement"
    displacement.x.scatter_forward()
    converged_reason = int(problem.solver.getConvergedReason())
    iterations = int(problem.solver.getIterationNumber())
    if converged_reason <= 0:
        raise RuntimeError(f"PETSc failed to converge: reason={converged_reason}.")

    compliance = scalar_allreduce(
        domain,
        fem.assemble_scalar(
            fem.form(ufl.dot(traction, displacement) * ds_load(11))
        ),
    )
    strain_energy = 0.5 * scalar_allreduce(
        domain,
        fem.assemble_scalar(
            fem.form(ufl.inner(stress(displacement), strain(displacement)) * ufl.dx)
        ),
    )
    energy_error = abs(strain_energy - 0.5 * compliance) / max(
        abs(strain_energy), np.finfo(float).eps
    )
    mean_load_displacement = scalar_allreduce(
        domain,
        fem.assemble_scalar(fem.form(displacement[2] * ds_load(11))),
    ) / load_area

    residual_form = fem.form(
        ufl.action(bilinear, displacement) - linear
    )
    residual = assemble_vector(residual_form)
    residual.ghostUpdate(
        addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE
    )
    left_vertical_dofs = component_dofs(vector_space, left_facets, 2)
    reaction_left_z = float(residual.array[left_vertical_dofs].sum())
    reaction_right_z = float(residual.array[right_vertical_dofs].sum())
    reaction_total_z = reaction_left_z + reaction_right_z
    equilibrium_error = abs(reaction_total_z - args.load_n) / args.load_n

    stress_tensor = stress(displacement)
    deviatoric = stress_tensor - (ufl.tr(stress_tensor) / 3.0) * ufl.Identity(3)
    von_mises_expression = ufl.sqrt(1.5 * ufl.inner(deviatoric, deviatoric))
    stress_space = fem.functionspace(domain, ("DG", 0))
    von_mises = fem.Function(stress_space)
    von_mises.name = "von_mises"
    interpolation = fem.Expression(
        von_mises_expression, stress_space.element.interpolation_points
    )
    von_mises.interpolate(interpolation)
    von_mises.x.scatter_forward()

    local_max_displacement = float(
        np.linalg.norm(displacement.x.array.reshape((-1, 3)), axis=1).max()
    )
    maximum_displacement = float(
        domain.comm.allreduce(local_max_displacement, op=MPI.MAX)
    )
    maximum_von_mises = float(
        domain.comm.allreduce(float(von_mises.x.array.max()), op=MPI.MAX)
    )
    volume = scalar_allreduce(
        domain, fem.assemble_scalar(fem.form(fem.Constant(domain, 1.0) * ufl.dx))
    )
    quality = mesh_quality(domain)

    with io.XDMFFile(domain.comm, run_dir / "dense_mbb_results.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(displacement)
        xdmf.write_function(von_mises)

    create_plot(
        domain,
        displacement,
        run_dir / "dense_mbb_3d_displacement.png",
        args.plot_scale,
    )

    results = {
        "case": "cubic_100_dense_equivalent",
        "geometry_policy": (
            "The complete original-mesh.stl is analyzed in the upright "
            "210 mm span by 40 mm height orientation."
        ),
        "interpretation": (
            "At 100% infill density, cubic infill is modeled as dense "
            "isotropic PLA without an internal void architecture."
        ),
        "source_stl": str(args.stl),
        "dolfinx_version": __import__("dolfinx").__version__,
        "mesh_size_mm": args.mesh_size_mm,
        "load_n": args.load_n,
        "youngs_modulus_pa": args.youngs_modulus_pa,
        "poisson_ratio": args.poisson_ratio,
        "density_kg_m3": args.density_kg_m3,
        "geometry": geometry,
        "mesh": quality,
        "gmsh": gmsh_report,
        "volume_m3": volume,
        "mass_kg": args.density_kg_m3 * volume,
        "load_area_m2": load_area,
        "support_facets": {
            "left": int(len(left_facets)),
            "right": int(len(right_facets)),
            "load": int(len(load_facets)),
        },
        "solver": {
            "converged_reason": converged_reason,
            "iterations": iterations,
        },
        "mean_load_displacement_m": mean_load_displacement,
        "maximum_displacement_m": maximum_displacement,
        "compliance_j": compliance,
        "strain_energy_j": strain_energy,
        "energy_identity_relative_error": energy_error,
        "reaction_left_z_n": reaction_left_z,
        "reaction_right_z_n": reaction_right_z,
        "reaction_total_z_n": reaction_total_z,
        "equilibrium_relative_error": equilibrium_error,
        "maximum_von_mises_pa": maximum_von_mises,
        "plot_deformation_scale": args.plot_scale,
    }
    (run_dir / "results.json").write_text(
        json.dumps(results, indent=2), encoding="utf-8"
    )
    (run_dir / "config.json").write_text(
        json.dumps(vars(args), indent=2, default=str), encoding="utf-8"
    )
    print(json.dumps(results, indent=2))
    print(f"run_directory={run_dir}")


if __name__ == "__main__":
    main()
