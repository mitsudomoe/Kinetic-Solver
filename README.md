# Kinetic-Solver

A finite volume solver for the one-dimensional Saint-Venant (shallow water) system
with topography, implementing the **kinetic scheme with reflections** of
B. Perthame and C. Simeoni, *A kinetic scheme for the Saint-Venant system with a
source term*, Calcolo 38 (2001). The paper can be found
[here](https://hal.archives-ouvertes.fr/hal-00922664/document).

A Rusanov scheme with hydrostatic reconstruction (Audusse et al., 2004) is
included as a reference.

## Build and run

```sh
make                 # builds build/kinetic-solver (C99, needs only libm)
./build/kinetic-solver --list
./build/kinetic-solver --case dam-break-wet --scheme kinetic --cells 400 --snapshot-every 20
python3 scripts/plot.py output             # output/final.png
python3 scripts/plot.py output --animate   # output/animation.gif
```

`graph.m` animates the same output in MATLAB or GNU Octave: `graph('output')`.

Main options (see `--help`): `--case`, `--scheme`, `--cells`, `--t-end`, `--cfl`,
`--output-dir`, `--snapshot-every`, `--no-snapshots`.

Each run writes `summary.json` (iterations, mass balance, minimum water height and,
when an exact solution exists, L1 and max errors) and, unless `--no-snapshots` is
given, `snapshot_NNNNN.dat` files with the columns `x z h hu [h_exact hu_exact]`
plus `times.dat`.

## Test cases

| Case | Description | Exact solution |
|---|---|---|
| `lake-at-rest` | H = 2 over a parabolic bump | yes (steady) |
| `lake-at-rest-step` | H = 1 over a bottom step | yes (steady) |
| `lake-at-rest-island` | H = 0.1, the top of the bump is dry | yes (steady) |
| `dam-break-dry` | h = 1 \| 0 on a flat bottom | Ritter |
| `dam-break-wet` | h = 2 \| 1 on a flat bottom | Stoker |
| `dam-break-bump` | h = 2 \| 1 over the bump | no |
| `periodic-pulse` | moving pulse, periodic boundaries | no |
| `sloshing` | h = 2 \| 1 over the bump, closed tank | no |

## Tests

```sh
make test            # needs python3 with numpy and scipy
```

- `tests/test_fluxes.py` checks the kinetic fluxes against a reference that
  integrates the interface densities (3.7)-(3.8) of the paper directly in ξ,
  without any of the closed-form integrals of Section 4. It also checks the
  conservation property of Remark 3.1 and the hydrostatic balance of both schemes.
  It also checks the momentum interface law of M. Ersoy, *Perthame–Simeoni kinetic
  traces and nonconservative interface laws for DG approximations of the
  Saint–Venant system* (2026, [hal-05618693](https://hal.science/hal-05618693)),
  Proposition 4.1: the two one-sided momentum fluxes differ by exactly
  `g ΔZ h^jump`, with `h^jump ≥ 0` given by the reflected and transmitted particles.
- `tests/test_solver.py` runs the solver: lake at rest kept to 1e-12, errors
  against Ritter and Stoker, convergence, positivity and mass conservation.

The tests run on every push through GitHub Actions.

## Implementation notes

- **Grid.** `n` cells of width `dx = (xmax - xmin) / n`, cell `i` centred at
  `xmin + (i - 1/2) dx`, ghost cells `0` and `n + 1`. The bottom is stored as cell
  averages, as in Section 3.1 of the paper.
- **Sign convention.** The bottom jump is the one seen from the cell that owns the
  flux: `dZ_{i+1/2} = Z_{i+1} - Z_i` and `dZ_{i-1/2} = Z_{i-1} - Z_i`. This is what
  the lake-at-rest proof of Theorem 3.2 uses.
- **Interface fluxes.** Each interface is visited once and returns the flux seen by
  the left cell (`IF^-`) and by the right cell (`IF^+`). `IF^+` is obtained from the
  formula for `IF^-` by the mirror symmetry x → -x, u → -u.
- **The II3q integral** has no elementary primitive (Remark 4.1). The paper uses a
  midpoint rule; the integrand has square-root singularities at both ends, so that
  rule converges like N^(-3/2) and the lake at rest is only kept to that accuracy.
  This code uses a tanh-sinh rule (113 nodes), which handles the endpoint
  singularities and keeps the lake at rest to round-off.
- **A sign in the paper.** The printed closed form of IIh3 ends with
  `- Fr β √(1-β²) / 2`; integrating gives `+`, which is what the code uses and what
  matches the flux definition numerically.
- **Time step.** CFL condition (3.9), `dt max(|u| + √(2gh)) ≤ dx`, with a default
  CFL number of 0.95 for the kinetic scheme and 0.9 for Rusanov.

`Guidebook.pdf` documents the original 2013 version of the code (single-precision,
compile-time options); the structure has changed since but the scheme is the same.

## License

GPL v3, see `gpl License.pdf`.
