"""Full-upright MBB shell/core screening model for sparse FDM infill.

The exterior shell follows the topology-optimized STL. A dense-shell volume
is estimated from the measured surface area and two slicer perimeter widths.
The core is represented by a cellular-solid screening law corrected by an
axial road-orientation factor measured from archived Bambu Studio G-code.
This is a comparative model, not an explicit filament or failure simulation.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

import ufl
from dolfinx import fem, io, mesh
from dolfinx.fem.petsc import LinearProblem, assemble_vector
from dolfinx.io import gmsh as dolfinx_gmsh

from mbb_3d_dolfinx import (
    component_dofs,
    create_facet_measure,
    create_plot,
    locate_boundary_regions,
    mesh_quality,
    scalar_allreduce,
    tetrahedralize_stl,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stl", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, default=Path("outputs/fem"))
    parser.add_argument(
        "--pattern", choices=("cubic", "gyroid", "honeycomb"), required=True
    )
    parser.add_argument("--infill-density", type=float, required=True)
    parser.add_argument("--infill-pitch-mm", type=float, required=True)
    parser.add_argument(
        "--axial-orientation-factor",
        type=float,
        required=True,
        help="G-code fourth-moment ratio relative to planar isotropic roads.",
    )
    parser.add_argument("--mesh-size-mm", type=float, default=2.4)
    parser.add_argument(
        "--msh",
        type=Path,
        help="Reuse a tetrahedral mesh so pattern cases compare cell-for-cell.",
    )
    parser.add_argument("--load-n", type=float, default=100.0)
    parser.add_argument("--youngs-modulus-pa", type=float, default=3.0e9)
    parser.add_argument("--poisson-ratio", type=float, default=0.35)
    parser.add_argument("--density-kg-m3", type=float, default=1240.0)
    parser.add_argument("--wall-thickness-mm", type=float, default=0.87)
    parser.add_argument("--core-exponent", type=float, default=2.0)
    parser.add_argument("--plot-scale", type=float, default=15.0)
    return parser.parse_args()


def cell_volumes(domain: mesh.Mesh) -> np.ndarray:
    tdim = domain.topology.dim
    domain.topology.create_connectivity(tdim, 0)
    cells_to_vertices = domain.topology.connectivity(tdim, 0)
    count = domain.topology.index_map(tdim).size_local
    volumes = np.empty(count, dtype=float)
    for cell in range(count):
        points = domain.geometry.x[cells_to_vertices.links(cell)]
        matrix = np.column_stack(
            (points[1] - points[0], points[2] - points[0], points[3] - points[0])
        )
        volumes[cell] = abs(np.linalg.det(matrix)) / 6.0
    return volumes


def exterior_area_and_cells(domain: mesh.Mesh) -> tuple[float, np.ndarray]:
    tdim = domain.topology.dim
    fdim = tdim - 1
    domain.topology.create_connectivity(fdim, tdim)
    domain.topology.create_connectivity(fdim, 0)
    exterior = mesh.exterior_facet_indices(domain.topology)
    facets_to_vertices = domain.topology.connectivity(fdim, 0)
    facets_to_cells = domain.topology.connectivity(fdim, tdim)
    local_area = 0.0
    boundary_cells: set[int] = set()
    for facet in exterior:
        vertices = facets_to_vertices.links(facet)
        points = domain.geometry.x[vertices]
        local_area += 0.5 * np.linalg.norm(
            np.cross(points[1] - points[0], points[2] - points[0])
        )
        boundary_cells.update(int(cell) for cell in facets_to_cells.links(facet))
    area = float(domain.comm.allreduce(local_area, op=MPI.SUM))
    return area, np.array(sorted(boundary_cells), dtype=np.int32)


def create_stiffness_scale(
    domain: mesh.Mesh,
    volumes: np.ndarray,
    boundary_cells: np.ndarray,
    surface_area: float,
    wall_thickness_m: float,
    infill_density: float,
    exponent: float,
    axial_orientation_factor: float,
) -> tuple[fem.Function, dict[str, float]]:
    total_volume = scalar_allreduce(domain, float(volumes.sum()))
    target_shell_volume = min(surface_area * wall_thickness_m, total_volume)
    density_scale = infill_density**exponent
    core_scale = min(1.0, density_scale * axial_orientation_factor)

    tdim = domain.topology.dim
    fdim = tdim - 1
    domain.topology.create_connectivity(tdim, fdim)
    domain.topology.create_connectivity(fdim, tdim)
    cells_to_facets = domain.topology.connectivity(tdim, fdim)
    facets_to_cells = domain.topology.connectivity(fdim, tdim)

    shell_fraction = np.zeros_like(volumes)
    visited = np.zeros(volumes.size, dtype=bool)
    frontier = np.asarray(boundary_cells, dtype=np.int32)
    represented_shell_volume = 0.0
    shell_layers = 0
    while frontier.size and represented_shell_volume < target_shell_volume:
        frontier = np.unique(frontier[~visited[frontier]])
        if frontier.size == 0:
            break
        visited[frontier] = True
        layer_volume = float(volumes[frontier].sum())
        remaining = target_shell_volume - represented_shell_volume
        layer_fraction = min(1.0, remaining / layer_volume)
        shell_fraction[frontier] = layer_fraction
        represented_shell_volume += layer_fraction * layer_volume
        shell_layers += 1
        if layer_fraction < 1.0:
            break
        neighbours: list[int] = []
        for cell in frontier:
            for facet in cells_to_facets.links(int(cell)):
                neighbours.extend(int(item) for item in facets_to_cells.links(facet))
        frontier = np.asarray(neighbours, dtype=np.int32)

    scale_space = fem.functionspace(domain, ("DG", 0))
    scale = fem.Function(scale_space)
    scale.name = "relative_stiffness"
    for cell, local_shell_fraction in enumerate(shell_fraction):
        dof = scale_space.dofmap.cell_dofs(int(cell))[0]
        scale.x.array[dof] = (
            local_shell_fraction
            + (1.0 - local_shell_fraction) * core_scale
        )
    scale.x.scatter_forward()

    target_shell_fraction = target_shell_volume / total_volume
    material_fraction = (
        target_shell_fraction
        + (1.0 - target_shell_fraction) * infill_density
    )
    return scale, {
        "surface_area_m2": surface_area,
        "total_volume_m3": total_volume,
        "target_shell_volume_m3": target_shell_volume,
        "represented_shell_volume_m3": represented_shell_volume,
        "target_shell_volume_fraction": target_shell_fraction,
        "shell_cell_layers": shell_layers,
        "fractional_outermost_shell_layer": float(shell_fraction.max()),
        "fractional_innermost_shell_layer": float(
            shell_fraction[shell_fraction > 0.0].min()
        ),
        "core_relative_stiffness": core_scale,
        "density_only_core_relative_stiffness": density_scale,
        "axial_orientation_factor": axial_orientation_factor,
        "estimated_total_material_fraction": material_fraction,
    }


def main() -> None:
    args = parse_args()
    if MPI.COMM_WORLD.size != 1:
        raise RuntimeError("This validation script currently requires one MPI rank.")
    if not 0.0 < args.infill_density <= 1.0:
        raise ValueError("Infill density must be in (0, 1].")
    if args.axial_orientation_factor <= 0.0:
        raise ValueError("Axial orientation factor must be positive.")

    case_percent = int(round(100.0 * args.infill_density))
    run_dir = args.output_root / (
        datetime.now().strftime("%Y%m%d_%H%M%S")
        + f"_mbb_{args.pattern}_{case_percent:02d}_full_upright"
    )
    run_dir.mkdir(parents=True, exist_ok=False)
    msh_path = run_dir / "full_upright_tetrahedra.msh"
    if args.msh:
        if not args.msh.exists():
            raise FileNotFoundError(args.msh)
        msh_path.write_bytes(args.msh.read_bytes())
        gmsh_report = {
            "reused_mesh": True,
            "source_mesh": str(args.msh),
        }
    else:
        gmsh_report = tetrahedralize_stl(args.stl, msh_path, args.mesh_size_mm)
    mesh_data = dolfinx_gmsh.read_from_msh(
        msh_path, MPI.COMM_WORLD, rank=0, gdim=3
    )
    domain = mesh_data.mesh
    domain.name = f"full_upright_{args.pattern}_{case_percent:02d}"
    domain.geometry.x[:] *= 0.001

    left_facets, right_facets, load_facets, geometry = locate_boundary_regions(domain)
    ds_load = create_facet_measure(domain, load_facets, tag=11)
    load_area = scalar_allreduce(
        domain, fem.assemble_scalar(fem.form(1.0 * ds_load(11)))
    )

    volumes = cell_volumes(domain)
    surface_area, boundary_cells = exterior_area_and_cells(domain)
    stiffness_scale, shell_core = create_stiffness_scale(
        domain,
        volumes,
        boundary_cells,
        surface_area,
        args.wall_thickness_mm * 0.001,
        args.infill_density,
        args.core_exponent,
        args.axial_orientation_factor,
    )

    vector_space = fem.functionspace(domain, ("Lagrange", 1, (3,)))
    left_dofs = fem.locate_dofs_topological(
        vector_space, domain.topology.dim - 1, left_facets
    )
    left_bc = fem.dirichletbc(
        np.zeros(3, dtype=PETSc.ScalarType), left_dofs, vector_space
    )
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
            [0.0, 0.0, -args.load_n / load_area], dtype=PETSc.ScalarType
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

    def dense_stress(field):
        epsilon = strain(field)
        return (
            lame_lambda * ufl.tr(epsilon) * ufl.Identity(3)
            + 2.0 * lame_mu * epsilon
        )

    bilinear = (
        stiffness_scale
        * ufl.inner(dense_stress(trial), strain(test))
        * ufl.dx
    )
    linear = ufl.dot(traction, test) * ds_load(11)
    problem = LinearProblem(
        bilinear,
        linear,
        bcs=boundary_conditions,
        petsc_options_prefix=f"mbb_{args.pattern}_{case_percent:02d}_",
        petsc_options={
            "ksp_type": "cg",
            "pc_type": "gamg",
            "ksp_rtol": 1.0e-10,
            "ksp_atol": 1.0e-12,
            "ksp_max_it": 15000,
            "ksp_error_if_not_converged": True,
        },
    )
    displacement = problem.solve()
    displacement.name = "displacement"
    displacement.x.scatter_forward()

    compliance = scalar_allreduce(
        domain,
        fem.assemble_scalar(
            fem.form(ufl.dot(traction, displacement) * ds_load(11))
        ),
    )
    strain_energy = 0.5 * scalar_allreduce(
        domain,
        fem.assemble_scalar(
            fem.form(
                stiffness_scale
                * ufl.inner(dense_stress(displacement), strain(displacement))
                * ufl.dx
            )
        ),
    )
    mean_load_displacement = (
        scalar_allreduce(
            domain,
            fem.assemble_scalar(fem.form(displacement[2] * ds_load(11))),
        )
        / load_area
    )

    residual = assemble_vector(fem.form(ufl.action(bilinear, displacement) - linear))
    residual.ghostUpdate(
        addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE
    )
    left_vertical_dofs = component_dofs(vector_space, left_facets, 2)
    reaction_left = float(residual.array[left_vertical_dofs].sum())
    reaction_right = float(residual.array[right_vertical_dofs].sum())
    reaction_total = reaction_left + reaction_right

    macro_stress = stiffness_scale * dense_stress(displacement)
    deviatoric = macro_stress - (ufl.tr(macro_stress) / 3.0) * ufl.Identity(3)
    von_mises_expression = ufl.sqrt(1.5 * ufl.inner(deviatoric, deviatoric))
    stress_space = fem.functionspace(domain, ("DG", 0))
    von_mises = fem.Function(stress_space)
    von_mises.name = "homogenized_von_mises"
    von_mises.interpolate(
        fem.Expression(
            von_mises_expression, stress_space.element.interpolation_points
        )
    )

    with io.XDMFFile(domain.comm, run_dir / "results.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(displacement)
        xdmf.write_function(stiffness_scale)
        xdmf.write_function(von_mises)

    create_plot(
        domain,
        displacement,
        run_dir / "displacement_3d.png",
        args.plot_scale,
        (
            f"Full upright MBB {args.pattern} {case_percent}%\n"
            f"homogenized shell/core | deformation x{args.plot_scale:g}"
        ),
    )

    maximum_displacement = float(
        np.linalg.norm(displacement.x.array.reshape((-1, 3)), axis=1).max()
    )
    material_volume = (
        shell_core["estimated_total_material_fraction"]
        * shell_core["total_volume_m3"]
    )
    results = {
        "case": f"{args.pattern}_{case_percent:02d}_full_upright",
        "status": "homogenized_shell_core_screening_model",
        "geometry_policy": (
            "Full upright source STL used for every case. Half-beam archive "
            "geometry and orientation are excluded."
        ),
        "gcode_verified": {
            "pattern": args.pattern,
            "nominal_density": args.infill_density,
            "measured_characteristic_spacing_mm": args.infill_pitch_mm,
            "axial_orientation_factor_vs_planar_isotropic": (
                args.axial_orientation_factor
            ),
            "wall_loops": 2,
            "outer_wall_line_width_mm": 0.42,
            "inner_wall_line_width_mm": 0.45,
            "normalized_full_upright_top_layers": 5,
            "normalized_full_upright_bottom_layers": 3,
            "layer_height_mm": 0.2,
        },
        "model_assumptions": {
            "dense_pla_isotropic": True,
            "core_law": (
                "E_core/E_PLA = infill_density ** core_exponent * "
                "G-code axial orientation factor"
            ),
            "core_exponent": args.core_exponent,
            "wall_thickness_mm": args.wall_thickness_mm,
            "requires_rve_or_coupon_calibration": True,
        },
        "dolfinx_version": __import__("dolfinx").__version__,
        "load_n": args.load_n,
        "mesh_size_mm": args.mesh_size_mm,
        "mesh": mesh_quality(domain),
        "gmsh": gmsh_report,
        "geometry": geometry,
        "shell_core": shell_core,
        "estimated_material_volume_m3": material_volume,
        "estimated_mass_kg": args.density_kg_m3 * material_volume,
        "solver": {
            "converged_reason": int(problem.solver.getConvergedReason()),
            "iterations": int(problem.solver.getIterationNumber()),
        },
        "mean_load_displacement_m": mean_load_displacement,
        "maximum_displacement_m": maximum_displacement,
        "compliance_j": compliance,
        "strain_energy_j": strain_energy,
        "energy_identity_relative_error": abs(strain_energy - 0.5 * compliance)
        / abs(strain_energy),
        "reaction_left_z_n": reaction_left,
        "reaction_right_z_n": reaction_right,
        "reaction_total_z_n": reaction_total,
        "equilibrium_relative_error": abs(reaction_total - args.load_n)
        / args.load_n,
        "maximum_homogenized_von_mises_pa": float(von_mises.x.array.max()),
        "plot_deformation_scale": args.plot_scale,
    }
    (run_dir / "results.json").write_text(
        json.dumps(results, indent=2), encoding="utf-8"
    )
    print(json.dumps(results, indent=2))
    print(f"run_directory={run_dir}")


if __name__ == "__main__":
    main()
