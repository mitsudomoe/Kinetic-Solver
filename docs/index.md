# Kinetic-Solver

A finite volume solver for the one-dimensional **Saint-Venant (shallow water)
system with topography**, built around the *kinetic scheme with reflections* of
B. Perthame and C. Simeoni (Calcolo, 2001).

$$
\partial_t h + \partial_x (hu) = 0, \qquad
\partial_t (hu) + \partial_x\!\left(hu^2 + \tfrac{g}{2}h^2\right) = -g\,h\,\partial_x Z .
$$

The scheme treats the bottom $Z$ as a potential barrier for fictitious particles:
particles cross a step when they have enough kinetic energy and are reflected
otherwise. This gives a finite volume method that

- keeps the water height **non-negative**, dry areas included;
- preserves the **lake at rest** ($u = 0$, $h + Z$ constant) exactly;
- satisfies a discrete **entropy inequality**;
- stays stable on **discontinuous** bottoms.

A Rusanov scheme with hydrostatic reconstruction (Audusse et al., 2004) is included
for comparison.

![Dam break on a wet bed](assets/dam-break-wet.png)

## In thirty seconds

```sh
make
./build/kinetic-solver --case dam-break-wet --cells 400 --snapshot-every 20
python3 scripts/plot.py output --animate
make test
```

## Where to go next

| Page | What you will find |
|---|---|
| [Getting started](getting-started.md) | Building, running, options, output files, plots |
| [Theory](theory.md) | The equations, the kinetic formulation and the scheme of the paper |
| [Implementation](implementation.md) | How the code is organised and the numerical choices it makes |
| [Validation](validation.md) | Exact solutions, convergence, well-balancing, symmetries |
| [Contributing](contributing.md) | Tests, adding a scheme, a case or an integrator |
| [References](references.md) | Papers the code is based on |

## Status

The scheme of the paper is implemented and checked against its definitions to
round-off (see [Validation](validation.md)). The code was first written in 2013;
the 2026 revision fixed several bugs (including a truncated value of π that broke
the lake-at-rest property), restructured the code and added a test suite. The
original guide, `Guidebook.pdf` in the repository, describes the 2013 version.
