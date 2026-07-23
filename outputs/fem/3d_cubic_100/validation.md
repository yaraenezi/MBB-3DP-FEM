# Dense MBB True 3D FEM Validation

Case:

```text
100% cubic-equivalent infill density
```

At 100% density, cubic infill contains no intentional internal void phase and
is represented as dense isotropic PLA.

## Model

```text
Source geometry:       raw/meshes/original-mesh.stl
Formulation:           3D small-strain linear elasticity
Elements:              first-order tetrahedra
DOLFINx:               0.11.0.post0
Material modulus:      3.0 GPa
Poisson ratio:         0.35
Density:               1240 kg/m3
Applied load:          100 N downward at top midspan
Left support:          ux = uy = uz = 0 over bearing pad
Right support:         uy = uz = 0 over bearing pad; ux free
```

## Accepted fine result

```text
Nominal mesh size:              2.4 mm
Tetrahedra:                     134,566
Degenerate tetrahedra:          0
Gmsh logged warnings:           0
PETSc converged reason:         3
PETSc iterations:               158
Mean load-patch displacement:  -0.1243156 mm
Maximum displacement:           0.1274890 mm
Compliance:                     0.01243156 J
Strain energy:                  0.006215781 J
Energy identity error:          2.43e-11
Left vertical reaction:         50.21361 N
Right vertical reaction:        49.78639 N
Total vertical reaction:        99.9999984 N
Equilibrium relative error:     1.63e-8
Maximum von Mises stress:       2.45682 MPa
Calculated volume:              79,248.26 mm3
Calculated mass:                98.27 g
```

The FEM volume differs from the source STL volume by approximately 0.356%.

## Mesh convergence

| Mesh size | Tetrahedra | Mean load displacement | Maximum von Mises |
|---:|---:|---:|---:|
| 4.0 mm | 40,782 | -0.1208523 mm | 2.28575 MPa |
| 3.0 mm | 78,459 | -0.1229991 mm | 2.30872 MPa |
| 2.4 mm | 134,566 | -0.1243156 mm | 2.45682 MPa |

Displacement changed by 1.745% from the coarse to nominal mesh and 1.059%
from the nominal to fine mesh. The loading displacement therefore satisfies
the provisional 3% convergence criterion.

The absolute maximum stress changed by 6.028% between the two finest meshes.
It is mesh-sensitive and must not be treated as a converged failure stress.
Support/load-edge singularities and percentile stresses should be evaluated
before using stress to compare designs.

## Output files

```text
dense_mbb_3d_displacement.png
dense_mbb_results.xdmf
dense_mbb_results.h5
dense_mbb_tetrahedra.msh
results.json
config.json
gmsh.log
```

## Validation boundary

This is a numerical elastic stiffness result for the dense topology-optimized
beam under a 100 N comparison load. It does not establish physical failure
load, inter-layer strength, fatigue life, or a certified safety factor.
