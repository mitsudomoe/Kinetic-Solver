# Getting started

## Build

The solver is plain C99 and only needs a C compiler and `libm`.

```sh
make            # builds build/kinetic-solver
make clean
```

The tests and plotting scripts need Python 3 with `numpy`, `scipy` and `matplotlib`.

## Run

```sh
./build/kinetic-solver --list          # cases, schemes and integrators
./build/kinetic-solver --case dam-break-dry --scheme kinetic --cells 400
```

| Option | Meaning | Default |
|---|---|---|
| `--case NAME` | test case (see below) | `lake-at-rest` |
| `--scheme NAME` | `kinetic` or `rusanov` | `kinetic` |
| `--integrator NAME` | `euler`, `ssp-rk2` or `ssp-rk3` | `euler` |
| `--cells N` | number of cells | per case (1000) |
| `--t-end T` | final time in seconds | per case |
| `--cfl C` | CFL number | 0.95 (kinetic), 0.9 (Rusanov) |
| `--dt DT` | constant time step instead of the CFL one | — |
| `--output-dir DIR` | where the files go | `output` |
| `--snapshot-every K` | also write a snapshot every K steps | initial and final only |
| `--no-snapshots` | only write `summary.json` | — |
| `--mirror` | start from the mirror image $x \to -x$ of the case | — |
| `--shift-cells K` | periodic cases: start from the case shifted by K cells | — |
| `--verbose` | print every time step | — |

`--mirror` and `--shift-cells` exist for the symmetry tests. With `--dt`, steps
where the time step exceeded the CFL limit are counted in `summary.json`; a run
whose solution stops being finite stops with an error.

## Test cases

| Case | Bottom | Initial state | Boundaries | Exact solution |
|---|---|---|---|---|
| `lake-at-rest` | parabolic bump | $h + Z = 2$, $u = 0$ | free | steady |
| `lake-at-rest-step` | step of 0.3 m | $h + Z = 1$, $u = 0$ | free | steady |
| `lake-at-rest-island` | bump | $h + Z = 0.1$: the top is dry | free | steady |
| `dam-break-dry` | flat | $h = 1 \mid 0$ | free | Ritter |
| `dam-break-wet` | flat | $h = 2 \mid 1$ | free | Stoker |
| `dam-break-bump` | bump | $h = 2 \mid 1$ | free | — |
| `periodic-pulse` | flat | moving Gaussian pulse | periodic | — |
| `sloshing` | bump | $h = 2 \mid 1$ | walls | — |
| `thacker` | parabolic bowl | tilted plane surface at rest, both shorelines move | walls | Thacker |
| `bump-subcritical` | bump | lake at rest, then inflow $q = 4.42$ | inflow / outflow | steady, Bernoulli |
| `bump-transcritical` | bump | lake at rest, then inflow $q = 1.53$ | inflow / outflow | steady, Bernoulli |
| `bump-hydraulic-jump` | bump | lake at rest, then inflow $q = 0.18$ | inflow / outflow | steady, with a jump |
| `beach-runup` | flat, then a 10% beach | wave travelling towards the dry beach | walls | — |

The domains are $[-10, 10]$ m, except $[-10, 15]$ m for the flows over the bump and the
beach. The bump is $Z(x) = 0.2 - 0.05\,x^2$ on $[-2, 2]$; the bowl is
$Z(x) = 0.5\,(x^2/64 - 1)$; the beach is $Z(x) = 0.1\,x$ for $x \ge 0$.

The flows over the bump take 100 to 200 s of simulated time to settle into their
steady state; `thacker` runs for one period of the oscillation (16.05 s). Runs longer
than the default only need `--t-end`.

## Output files

Each run writes to `--output-dir`:

- **`summary.json`**: case, scheme, integrator, number of cells and steps, CFL
  number, final time, wall time, initial and final mass, minimum water height over
  the run, CFL violations, and the L1 and maximum errors of $h$ and $hu$ when an
  exact solution exists (`null` otherwise).
- **`snapshot_NNNNN.dat`**: one line per cell with the columns
  `x z h hu`, plus `h_exact hu_exact` when an exact solution exists. Two header
  lines start with `#`.
- **`times.dat`**: `index time` for every snapshot.

## Plots

```sh
python3 scripts/plot.py output             # output/final.png
python3 scripts/plot.py output --all       # one PNG per snapshot
python3 scripts/plot.py output --animate   # output/animation.gif
```

In MATLAB or GNU Octave, `graph('output')` animates the same files.
