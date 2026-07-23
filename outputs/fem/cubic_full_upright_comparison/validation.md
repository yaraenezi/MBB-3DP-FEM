# Full-Upright Cubic FEM Validation

## Geometry decision

All accepted cases use `raw/meshes/original-mesh.stl` as the same complete
upright beam. Its FEM envelope is approximately `210.27 x 24.15 x 40.22 mm`.
The 25% archive contains a half-beam print with a `12.00 mm` build height; it
is used to verify the cubic toolpath only and is excluded from FEM geometry.

## G-code evidence

| Archive setting | Cubic 15% | Cubic 25% |
|---|---:|---:|
| Declared pattern | cubic | cubic |
| Declared density | 15% | 25% |
| Measured projected pitch | 8.142 mm | 4.885 mm |
| Dominant road families | 45, 105, 165 deg | 45, 105, 165 deg |

The pitch ratio is approximately `1.667`, matching the inverse density ratio
`25/15`. This supports a common cubic generator scaled by infill density.

## Accepted fine-mesh results

The load is `100 N`; dense PLA assumptions are `E = 3.0 GPa`, `nu = 0.35`,
and `rho = 1240 kg/m3`.

| Case | Mean displacement | Maximum displacement | Apparent stiffness | Estimated mass |
|---|---:|---:|---:|---:|
| Cubic 15% | 0.39089 mm | 0.40766 mm | 255.83 N/mm | 40.93 g |
| Cubic 25% | 0.32648 mm | 0.33746 mm | 306.30 N/mm | 47.68 g |
| Dense 100% | 0.12432 mm | 0.12749 mm | 804.40 N/mm | 98.27 g |

Both sparse cases use `134,566` tetrahedra at a nominal `2.4 mm` mesh size.
Reaction equilibrium errors are below `9e-9`, energy identity errors are
below `1.1e-11`, and no degenerate tetrahedra were found.

## Mesh check

Changing the nominal mesh size from `3.0` to `2.4 mm` changed mean load
displacement by:

- Cubic 15%: `3.53%`
- Cubic 25%: `2.29%`

The 15% case is close to, but slightly above, a provisional 3% convergence
target. A finer run is appropriate before treating its displacement as a
mesh-independent value.

## Model boundary

The sparse infill is a homogenized shell/core screening model, not an explicit
filament-contact mesh. The dense shell volume is matched to two perimeter
widths (`0.42 + 0.45 mm`), while the core uses
`E_core/E_PLA = relative_density^2`. Reported sparse-case stress is a
homogenized macro stress and must not be used as a printed PLA failure
criterion. Coupon or representative-volume calibration is required before
physical design decisions.
