"""DOLFINx 2D plane-stress benchmark for the MBB beam envelope."""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

import matplotlib
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

import ufl
from dolfinx import fem, io, mesh
from dolfinx.fem.petsc import LinearProblem, assemble_vector


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", type=Path, default=Path("outputs/fem"))
    parser.add_argument("--load-n", type=float, default=100.0)
    parser.add_argument("--nx", type=int, default=168)
    parser.add_argument("--ny", type=int, default=32)
    parser.add_argument("--youngs-modulus-pa", type=float, default=3.0e9)
    parser.add_argument("--poisson-ratio", type=float, default=0.35)
    parser.add_argument("--thickness-m", type=float, default=0.024)
    parser.add_argument("--plot-scale", type=float, default=50.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if MPI.COMM_WORLD.size != 1:
        raise RuntimeError("This validation script currently requires one MPI rank.")

    run_dir = args.output_root / (
        datetime.now().strftime("%Y%m%d_%H%M%S") + "_mbb_2d_plane_stress"
    )
    run_dir.mkdir(parents=True, exist_ok=False)

    length = 0.210
    height = 0.040
    domain = mesh.create_rectangle(
        MPI.COMM_WORLD,
        [np.array([0.0, 0.0]), np.array([length, height])],
        [args.nx, args.ny],
        cell_type=mesh.CellType.triangle,
    )
    vector_space = fem.functionspace(domain, ("Lagrange", 1, (2,)))
    fdim = domain.topology.dim - 1
    support_width = 0.0075
    load_half_width = 0.0075

    left_facets = mesh.locate_entities_boundary(
        domain,
        fdim,
        lambda x: np.isclose(x[1], 0.0) & (x[0] <= support_width),
    )
    right_facets = mesh.locate_entities_boundary(
        domain,
        fdim,
        lambda x: np.isclose(x[1], 0.0) & (x[0] >= length - support_width),
    )
    load_facets = mesh.locate_entities_boundary(
        domain,
        fdim,
        lambda x: np.isclose(x[1], height)
        & (np.abs(x[0] - 0.5 * length) <= load_half_width),
    )
    if min(len(left_facets), len(right_facets), len(load_facets)) == 0:
        raise RuntimeError("A support or load boundary is empty.")

    left_dofs = fem.locate_dofs_topological(vector_space, fdim, left_facets)
    left_bc = fem.dirichletbc(
        np.zeros(2, dtype=PETSc.ScalarType), left_dofs, vector_space
    )
    right_vertical_dofs = fem.locate_dofs_topological(
        vector_space.sub(1), fdim, right_facets
    )
    right_bc = fem.dirichletbc(
        PETSc.ScalarType(0.0), right_vertical_dofs, vector_space.sub(1)
    )
    boundary_conditions = [left_bc, right_bc]

    sorted_load_facets = np.sort(np.unique(load_facets)).astype(np.int32)
    load_tags = mesh.meshtags(
        domain,
        fdim,
        sorted_load_facets,
        np.full(sorted_load_facets.size, 1, dtype=np.int32),
    )
    ds = ufl.Measure("ds", domain=domain, subdomain_data=load_tags)
    local_length = fem.assemble_scalar(fem.form(1.0 * ds(1)))
    load_line_length = float(domain.comm.allreduce(local_length, op=MPI.SUM))
    traction = fem.Constant(
        domain,
        np.array(
            [0.0, -args.load_n / (args.thickness_m * load_line_length)],
            dtype=PETSc.ScalarType,
        ),
    )

    trial = ufl.TrialFunction(vector_space)
    test = ufl.TestFunction(vector_space)
    lame_mu = args.youngs_modulus_pa / (2.0 * (1.0 + args.poisson_ratio))
    plane_stress_lambda = (
        args.youngs_modulus_pa
        * args.poisson_ratio
        / (1.0 - args.poisson_ratio**2)
    )

    def strain(field):
        return ufl.sym(ufl.grad(field))

    def stress(field):
        epsilon = strain(field)
        return (
            plane_stress_lambda * ufl.tr(epsilon) * ufl.Identity(2)
            + 2.0 * lame_mu * epsilon
        )

    bilinear = (
        args.thickness_m * ufl.inner(stress(trial), strain(test)) * ufl.dx
    )
    linear = args.thickness_m * ufl.dot(traction, test) * ds(1)
    problem = LinearProblem(
        bilinear,
        linear,
        bcs=boundary_conditions,
        petsc_options_prefix="mbb_2d_plane_stress_",
        petsc_options={
            "ksp_type": "cg",
            "pc_type": "gamg",
            "ksp_rtol": 1.0e-11,
            "ksp_atol": 1.0e-13,
            "ksp_error_if_not_converged": True,
        },
    )
    displacement = problem.solve()
    displacement.name = "displacement"
    displacement.x.scatter_forward()

    compliance = float(
        domain.comm.allreduce(
            fem.assemble_scalar(
                fem.form(
                    args.thickness_m
                    * ufl.dot(traction, displacement)
                    * ds(1)
                )
            ),
            op=MPI.SUM,
        )
    )
    strain_energy = 0.5 * float(
        domain.comm.allreduce(
            fem.assemble_scalar(
                fem.form(
                    args.thickness_m
                    * ufl.inner(stress(displacement), strain(displacement))
                    * ufl.dx
                )
            ),
            op=MPI.SUM,
        )
    )
    mean_load_displacement = float(
        domain.comm.allreduce(
            fem.assemble_scalar(fem.form(displacement[1] * ds(1))),
            op=MPI.SUM,
        )
        / load_line_length
    )

    residual = assemble_vector(fem.form(ufl.action(bilinear, displacement) - linear))
    residual.ghostUpdate(
        addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE
    )
    left_vertical_dofs = fem.locate_dofs_topological(
        vector_space.sub(1), fdim, left_facets
    )
    reaction_left = float(residual.array[left_vertical_dofs].sum())
    reaction_right = float(residual.array[right_vertical_dofs].sum())

    sigma = stress(displacement)
    von_mises_expression = ufl.sqrt(
        sigma[0, 0] ** 2
        - sigma[0, 0] * sigma[1, 1]
        + sigma[1, 1] ** 2
        + 3.0 * sigma[0, 1] ** 2
    )
    stress_space = fem.functionspace(domain, ("DG", 0))
    von_mises = fem.Function(stress_space)
    von_mises.name = "von_mises"
    von_mises.interpolate(
        fem.Expression(
            von_mises_expression, stress_space.element.interpolation_points
        )
    )

    with io.XDMFFile(domain.comm, run_dir / "mbb_2d_results.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(displacement)
        xdmf.write_function(von_mises)

    coordinates = domain.geometry.x[:, :2]
    nodal_displacement = displacement.x.array.reshape((-1, 2))
    magnitude = np.linalg.norm(nodal_displacement, axis=1) * 1000.0
    deformed = (coordinates + args.plot_scale * nodal_displacement) * 1000.0
    domain.topology.create_connectivity(domain.topology.dim, 0)
    cell_vertices = domain.topology.connectivity(domain.topology.dim, 0)
    cells = np.vstack(
        [
            cell_vertices.links(cell)
            for cell in range(domain.topology.index_map(domain.topology.dim).size_local)
        ]
    )
    triangulation = mtri.Triangulation(deformed[:, 0], deformed[:, 1], cells)
    figure, axes = plt.subplots(figsize=(13, 4.5), dpi=180)
    result_plot = axes.tripcolor(
        triangulation,
        magnitude,
        shading="gouraud",
        cmap="viridis",
    )
    axes.set_aspect("equal")
    axes.set_axis_off()
    axes.set_title(
        "MBB 2D plane-stress FEM\n"
        f"100 N comparison load | deformation x{args.plot_scale:g}"
    )
    colorbar = figure.colorbar(result_plot, ax=axes, shrink=0.82)
    colorbar.set_label("Displacement magnitude (mm)")
    figure.tight_layout()
    figure.savefig(run_dir / "mbb_2d_displacement.png")
    plt.close(figure)

    maximum_displacement = float(magnitude.max() / 1000.0)
    maximum_von_mises = float(von_mises.x.array.max())
    reaction_total = reaction_left + reaction_right
    results = {
        "case": "mbb_2d_plane_stress_solid_rectangle",
        "load_n": args.load_n,
        "dimensions_m": [length, height, args.thickness_m],
        "mesh": {"nx": args.nx, "ny": args.ny, "triangles": int(cells.shape[0])},
        "material": {
            "youngs_modulus_pa": args.youngs_modulus_pa,
            "poisson_ratio": args.poisson_ratio,
        },
        "solver": {
            "converged_reason": int(problem.solver.getConvergedReason()),
            "iterations": int(problem.solver.getIterationNumber()),
        },
        "mean_load_displacement_m": mean_load_displacement,
        "maximum_displacement_m": maximum_displacement,
        "maximum_von_mises_pa": maximum_von_mises,
        "compliance_j": compliance,
        "strain_energy_j": strain_energy,
        "energy_identity_relative_error": abs(strain_energy - 0.5 * compliance)
        / abs(strain_energy),
        "reaction_left_y_n": reaction_left,
        "reaction_right_y_n": reaction_right,
        "reaction_total_y_n": reaction_total,
        "equilibrium_relative_error": abs(reaction_total - args.load_n) / args.load_n,
    }
    (run_dir / "results.json").write_text(
        json.dumps(results, indent=2), encoding="utf-8"
    )
    print(json.dumps(results, indent=2))
    print(f"run_directory={run_dir}")


if __name__ == "__main__":
    main()
