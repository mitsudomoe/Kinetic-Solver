# References

## The scheme

- B. Perthame and C. Simeoni, *A kinetic scheme for the Saint-Venant system with a
  source term*, Calcolo 38(4), 201–231 (2001).
  [hal-00922664](https://hal.science/hal-00922664).
  The scheme implemented here: kinetic formulation, interface densities with
  reflections, properties (Theorem 3.2) and closed-form fluxes (Section 4).

## Analysis

- M. Ersoy, *Perthame–Simeoni kinetic traces and nonconservative interface laws for
  DG approximations of the Saint–Venant system*, preprint (2026).
  [hal-05618693](https://hal.science/hal-05618693).
  Discontinuous Galerkin interpretation of the scheme; exact momentum interface
  law (Proposition 4.1), tested in this code; the finite volume scheme is the case of
  polynomial degree 0.
- M. Ersoy and C. Simeoni, *Consistency of the kinetic scheme with reflections*,
  research report, IMATH, Université de Toulon (2013).

## Reference scheme and numerical tools

- E. Audusse, F. Bouchut, M.-O. Bristeau, R. Klein and B. Perthame, *A fast and
  stable well-balanced scheme with hydrostatic reconstruction for shallow water
  flows*, SIAM J. Sci. Comput. 25(6), 2050–2065 (2004). The hydrostatic
  reconstruction used by the Rusanov scheme.
- C.-W. Shu and S. Osher, *Efficient implementation of essentially non-oscillatory
  shock-capturing schemes*, J. Comput. Phys. 77(2), 439–471 (1988). The SSP
  Runge-Kutta methods.
- H. Takahasi and M. Mori, *Double exponential formulas for numerical integration*,
  Publ. RIMS Kyoto Univ. 9(3), 721–741 (1974). The tanh-sinh quadrature used for
  $II^q_3$.

## Exact solutions

- A. Ritter, *Die Fortpflanzung der Wasserwellen*, Zeitschrift des Vereines
  Deutscher Ingenieure 36 (1892). Dam break on a dry bed.
- J. J. Stoker, *Water Waves: The Mathematical Theory with Applications*,
  Interscience (1957). Dam break on a wet bed.
