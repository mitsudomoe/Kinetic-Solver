# Contributing

## Running the tests

```sh
make test     # needs python3 with numpy and scipy
```

| File | What it checks |
|---|---|
| `tests/test_fluxes.py` | kinetic fluxes against the definitions of the paper; mass flux single-valued; momentum interface law of Ersoy (2026); hydrostatic balance; Rusanov on flat bottoms and dry states; Stoker's middle state |
| `tests/test_solver.py` | lake at rest on three topographies; errors against Ritter and Stoker; convergence; positivity; mass conservation; command line |
| `tests/test_integrators.py` | order in time of each integrator; SSP properties; Euler versus Runge-Kutta; mirror and translation symmetries; options |

`tests/flux_harness.c` exposes single interface fluxes to the Python tests: it
reads `scheme hL huL zL hR huR zR` lines and prints the two fluxes.

GitHub Actions builds the code with gcc and clang (warnings are errors) and runs the
tests on every push and pull request.

Before submitting a change to the numerics, it is worth also running the suite under
the sanitizers:

```sh
make clean
make CFLAGS="-std=c99 -O1 -g -fsanitize=address,undefined -Wall -Wextra"
make test CFLAGS="-std=c99 -O1 -g -fsanitize=address,undefined -Wall -Wextra"
```

## Adding a scheme

A scheme is a pair of functions and an entry in the table of `src/scheme.c`:

```c
/* Fluxes at the interface between `left` and `right`, seen by each cell. */
void my_interface_flux(CellState left, CellState right, Flux *seen_by_left, Flux *seen_by_right);

/* Bound on the particle or wave speeds, for the CFL condition. */
double my_wave_speed(double h, double hu);
```

```c
static const Scheme schemes[] = {
	...
	{"my-scheme", "one-line description", my_interface_flux, my_wave_speed, 0.9},
};
```

The two fluxes differ only when the scheme accounts for the bottom at the interface;
on a flat bottom they must be equal. Add the new scheme to the `SCHEMES` tuples of the
tests: the lake-at-rest, conservation and symmetry tests then cover it.

## Adding a test case

A case is an entry in the table of `src/cases.c`: domain, final time, default number
of cells, boundary condition, the **primitive** of the bottom (cell averages are
computed from it), the initial condition and, if one is known, the exact solution.
Keep the domain symmetric about $x = 0$ so that the mirror test applies.

## Adding a time integrator

An integrator is a list of Shu-Osher stages $(a_k, b_k)$ in `src/integrator.c`. Keep
$a_k, b_k \ge 0$ and $a_k + b_k = 1$ so that it stays strong-stability-preserving,
and add it to `INTEGRATORS` in `tests/test_integrators.py` with its expected order.

## Updating the documentation

The documentation lives in `docs/` and is built with
[MkDocs Material](https://squidfunk.github.io/mkdocs-material/):

```sh
pip install -r docs/requirements.txt
mkdocs serve            # http://127.0.0.1:8000, reloads on change
mkdocs build --strict   # what the CI runs
```

The pages also render directly on GitHub. Figures and tables of the
[Validation](validation.md) page come from

```sh
make && python3 scripts/make_doc_figures.py
```

which writes the PNG files to `docs/assets/` and prints the tables.

## Conventions

- C99, no dependency beyond `libm`; `-Wall -Wextra -Wpedantic -Wshadow` clean.
- Double precision everywhere; physical constants in `src/physics.h`.
- Bottom jumps are seen from the cell that owns the flux (see
  [Implementation](implementation.md#sign-convention-for-the-bottom-jump)).
- A numerical change comes with a test that fails without it.
