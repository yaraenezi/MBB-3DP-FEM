# Six-Case Full-Upright Infill FEM Validation

## Geometry Policy

All accepted cases use the same complete upright
`raw/meshes/original-mesh.stl`, approximately
`210.27 x 24.15 x 40.22 mm`. The 25% slicer archives contain half-thickness
prints with a `12.00 mm` build height; they are used only to characterize the
pattern at that density. Their geometry, orientation, top-shell setting, and
mass are not used in the FEM comparison.

Every fine case uses the same `134,566`-tetrahedron mesh. Every coarse case
uses the same `78,459`-tetrahedron mesh.

## G-code Evidence

The generalized archive inspector selects only `FEATURE: Sparse infill`
extrusion moves and calculates the axial fourth-orientation moment
`<cos(theta)^4>`. Pattern factors are normalized against cubic at the same
density so the previously accepted cubic baseline is unchanged.

| Pattern | Density | Sparse layers | Sparse segments | Raw orientation factor | Factor vs. cubic |
|---|---:|---:|---:|---:|---:|
| Cubic | 15% | 193 | 9,207 | 1.05165 | 1.00000 |
| Gyroid | 15% | 193 | 16,613 | 1.13167 | 1.07609 |
| Honeycomb | 15% | 193 | 49,252 | 1.05693 | 1.00502 |
| Cubic | 25% | 57 | 14,451 | 1.12216 | 1.00000 |
| Gyroid | 25% | 57 | 26,322 | 1.29016 | 1.14971 |
| Honeycomb | 25% | 57 | 67,081 | 1.19053 | 1.06093 |

The orientation factor is an affine road-network screening descriptor. It is
not a periodic-cell homogenization result and does not include filament
contacts or inter-layer adhesion.

## Model

The model combines:

- a dense shell volume matched to the two slicer perimeter widths
  (`0.42 + 0.45 mm`);
- a homogenized core with density scaling
  `Ecore/EPLA = relative_density^2`;
- a G-code axial orientation correction normalized to cubic at the same
  density;
- dense PLA assumptions of `E = 3.0 GPa`, `nu = 0.35`, and
  `rho = 1240 kg/m3`;
- a common `100 N` top-midspan load and identical support patches.

## Accepted Fine-Mesh Results

| Case | Mean displacement | Maximum displacement | Apparent stiffness | Estimated mass |
|---|---:|---:|---:|---:|
| Cubic 15% | 0.39089 mm | 0.40766 mm | 255.83 N/mm | 40.93 g |
| Gyroid 15% | 0.38626 mm | 0.40255 mm | 258.89 N/mm | 40.93 g |
| Honeycomb 15% | 0.39057 mm | 0.40731 mm | 256.04 N/mm | 40.93 g |
| Cubic 25% | 0.32648 mm | 0.33746 mm | 306.30 N/mm | 47.68 g |
| Gyroid 25% | 0.31754 mm | 0.32790 mm | 314.92 N/mm | 47.68 g |
| Honeycomb 25% | 0.32270 mm | 0.33341 mm | 309.89 N/mm | 47.68 g |
| Dense 100% | 0.12432 mm | 0.12749 mm | 804.40 N/mm | 98.27 g |

Within this screening model, gyroid is the stiffest sparse pattern at both
densities. The differences between patterns are modest: approximately `1.2%`
across 15% cases and `2.8%` across 25% cases.

## Numerical Checks

Reaction-equilibrium errors are below `1.7e-8`, energy-identity errors are
below `2.5e-11`, and no degenerate tetrahedra were reported.

Changing nominal mesh size from `3.0` to `2.4 mm` changed mean displacement by:

| Pattern | 15% | 25% |
|---|---:|---:|
| Cubic | 3.53% | 2.29% |
| Gyroid | 3.43% | 2.15% |
| Honeycomb | 3.52% | 2.23% |

The 15% cases remain slightly above a provisional 3% displacement-convergence
target and should receive a finer mesh before being treated as mesh
independent.

## Validation Boundary

These are comparative linear-elastic screening results. They do not establish
printed-part failure load, strength, safety factor, fatigue life, buckling
capacity, or a universal ranking of infill patterns. Pattern-specific printed
coupons or representative-volume models and full-beam bending experiments are
required before physical design decisions.
