# Kinetic-Solver

A finite volume solver for the one-dimensional Saint-Venant (shallow water) system
with topography, implementing the **kinetic scheme with reflections** of
B. Perthame and C. Simeoni, *A kinetic scheme for the Saint-Venant system with a
source term*, Calcolo 38 (2001) ([hal-00922664](https://hal.science/hal-00922664)).

The scheme keeps the water height non-negative, preserves the lake at rest exactly
(dry areas included) and stays stable on discontinuous bottoms. A Rusanov scheme
with hydrostatic reconstruction (Audusse et al., 2004) is included for comparison.

![Dam break on a wet bed](docs/assets/dam-break-wet.png)

## Quick start

```sh
make                                     # C99, needs only libm
./build/kinetic-solver --list            # cases, schemes, integrators
./build/kinetic-solver --case dam-break-wet --cells 400 --snapshot-every 20
python3 scripts/plot.py output --animate # output/animation.gif
make test                                # needs python3 with numpy and scipy
```

## Documentation

The documentation is in [`docs/`](docs/index.md) and renders directly on GitHub:

- [Getting started](docs/getting-started.md): building, options, test cases, output files, plots
- [Theory](docs/theory.md): the equations, the kinetic formulation and the scheme
- [Implementation](docs/implementation.md): code organisation and numerical choices
- [Validation](docs/validation.md): exact solutions, convergence, well-balancing, symmetries
- [Contributing](docs/contributing.md): tests, adding a scheme, a case or an integrator
- [References](docs/references.md)

It can also be built as a website with MkDocs (`pip install -r docs/requirements.txt`,
then `mkdocs serve`).

`Guidebook.pdf` documents the original 2013 version of the code.

## License

GPL v3, see `gpl License.pdf`.
