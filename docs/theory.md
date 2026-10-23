# Theory

This page summarises the mathematics the code implements. Equation numbers such as
(3.9) refer to Perthame & Simeoni (2001); see [References](references.md).

## The Saint-Venant system

With $h \ge 0$ the water height, $u$ the velocity, $q = hu$ the discharge, $Z(x)$
the bottom and $g$ the gravity,

$$
\begin{aligned}
\partial_t h + \partial_x (hu) &= 0, \\
\partial_t (hu) + \partial_x\!\left(hu^2 + \tfrac{g}{2} h^2\right) &= -g\,h\,\partial_x Z .
\end{aligned}
$$

The system is hyperbolic for $h > 0$, with wave speeds $u \pm \sqrt{gh}$; it loses
strict hyperbolicity at $h = 0$ (dry areas). Two properties matter for the numerics:

- **Steady states at rest.** Every state with $u = 0$ and $h + Z = H$ constant (a
  *lake at rest*) is steady: the pressure gradient exactly balances the bottom
  slope. A scheme that does not reproduce this balance creates spurious currents
  over any non-flat bottom.
- **Energy.** The energy $E = \tfrac12 h u^2 + \tfrac{g}{2} h^2 + g Z h$ is conserved
  by smooth solutions and decreases across shocks:
  $\partial_t E + \partial_x\!\big(u (E + \tfrac{g}{2} h^2)\big) \le 0$.

## Kinetic formulation

The kinetic approach describes the water as a density $f(t, x, \xi)$ of fictitious
particles moving at velocity $\xi$. The macroscopic quantities are its moments, and
at equilibrium the density is a *Maxwellian*

$$
M(h, \xi - u) = \sqrt{h}\;\chi\!\left(\frac{\xi - u}{\sqrt{h}}\right),
\qquad
\chi(\omega) = \frac{\sqrt{2}}{\pi \sqrt{g}} \left(1 - \frac{\omega^2}{2g}\right)_+^{1/2},
$$

which satisfies

$$
\int M \, d\xi = h, \qquad \int \xi M \, d\xi = hu, \qquad \int \xi^2 M \, d\xi = hu^2 + \tfrac{g}{2} h^2 .
$$

The Saint-Venant system is then equivalent to the single linear transport equation

$$
\partial_t f + \xi\, \partial_x f - g\, \partial_x Z\; \partial_\xi f = Q,
$$

where the collision term $Q$ has zero mass and momentum and relaxes $f$ towards $M$.
The bottom acts as a force on the particles: along a trajectory their mechanical
energy $\tfrac12 \xi^2 + g Z$ is conserved.

The half-ellipse $\chi$ is not arbitrary: it is the density that minimises the
kinetic energy $\int \big(\tfrac12 \xi^2 f + \tfrac{\pi^2 g^2}{6} f^3 + g Z f\big)\, d\xi$
under the constraints of given $h$ and $hu$ (Lemma 2.3). This choice is what gives
the scheme its entropy inequality. Its compact support, $|\xi - u| \le \sqrt{2gh}$,
also bounds the particle speeds, which is where the CFL condition comes from.

## The scheme with reflections

The domain is divided into cells $C_i$ of width $\Delta x$, with piecewise
constant bottom $Z_i$. At each time step the density is reset to the local
Maxwellian $M_i^n$ ("all collisions at once") and transported upwind:

$$
f_i^{n+1}(\xi) = M_i^n(\xi) - \frac{\Delta t}{\Delta x}\,\xi\left(M^-_{i+1/2}(\xi) - M^+_{i-1/2}(\xi)\right).
$$

The bottom jump at an interface is treated as a **potential barrier**. Write
$\Delta Z_{i+1/2} = Z_{i+1} - Z_i$ and $\Delta Z_{i-1/2} = Z_{i-1} - Z_i$: the jump
*seen from cell $i$*. Particles leaving cell $i$ through its right interface with
velocity $\xi > 0$ are counted directly. The particles arriving from the right
($\xi < 0$) are of two kinds:

- particles of cell $i$ that tried to climb the step without enough energy, and
  come back **reflected** with velocity $-\xi$ (when $\xi^2 \le 2g\,\Delta Z_{i+1/2}$);
- particles of cell $i+1$ that **crossed** the step, their speed changed by the
  energy balance $\tfrac12 \xi^2 = \tfrac12 \xi'^2 + g\,\Delta Z_{i+1/2}$.

This gives, for $\xi \le 0$,

$$
M_{i+1/2}(\xi) = M_i(-\xi)\,\mathbb{1}_{\xi^2 \le 2g\Delta Z_{i+1/2}}
+ M_{i+1}\!\left(-\sqrt{\xi^2 - 2g\Delta Z_{i+1/2}}\right)\mathbb{1}_{\xi^2 \ge 2g\Delta Z_{i+1/2}},
$$

and symmetrically at $x_{i-1/2}$ for $\xi \ge 0$, with $M_{i-1}$ and
$+\sqrt{\xi^2 - 2g\Delta Z_{i-1/2}}$. The interface densities are then

$$
M^-_{i+1/2} = M_i\,\mathbb{1}_{\xi \ge 0} + M_{i+1/2}\,\mathbb{1}_{\xi \le 0},
\qquad
M^+_{i-1/2} = M_{i-1/2}\,\mathbb{1}_{\xi \ge 0} + M_i\,\mathbb{1}_{\xi \le 0}.
$$

Taking moments gives the finite volume scheme

$$
U_i^{n+1} = U_i^n - \frac{\Delta t}{\Delta x}\left(F^-_{i+1/2} - F^+_{i-1/2}\right),
\qquad
F^\pm = \int \xi \begin{pmatrix} 1 \\ \xi \end{pmatrix} M^\pm(\xi)\, d\xi .
$$

The source term never appears as a volume term: it is entirely carried by the
difference between the flux seen from the left of an interface, $F^-$, and the flux
seen from its right, $F^+$ (interfacial upwind sources).

## Properties

Under the CFL condition (3.9)

$$
\Delta t \, \max_i \left(|u_i^n| + \sqrt{2 g h_i^n}\right) \le \Delta x,
$$

the scheme (Theorem 3.2)

1. keeps $h_i^n \ge 0$: $f_i^{n+1}$ is a convex combination of non-negative densities;
2. satisfies a discrete in-cell entropy (energy) inequality;
3. preserves every lake at rest: when $u_i = 0$ and $h_i + Z_i$ is constant,
   $M^-_{i+1/2} = M^+_{i-1/2} = M_i$, so nothing moves.

The **mass flux is single-valued**: $F^{h,-}_{i+1/2} = F^{h,+}_{i+1/2}$ (Remark 3.1),
so water is conserved exactly. The momentum fluxes differ across a step, which is
how the bottom force enters.

## Closed-form fluxes

With the dimensionless numbers

$$
Fr_i = \frac{u_i}{\sqrt{2 g h_i}}, \qquad
K^-_{i+1/2} = \frac{\Delta Z_{i+1/2}}{h_i}, \qquad
K^+_{i+1/2} = \frac{\Delta Z_{i+1/2}}{h_{i+1}},
$$

and $(a)_+ = \max(a, 0)$, $(a)_- = \max(-a, 0)$, the flux seen from cell $i$ is

$$
\begin{aligned}
F^{h,-}_{i+1/2} &= \frac{2\sqrt{2g}}{\pi}\left[h_i^{3/2}\left(II^h_1 + II^h_2\right) + h_{i+1}^{3/2}\, II^h_3\right], \\
F^{q,-}_{i+1/2} &= \frac{4g}{\pi}\left[h_i^{2}\left(II^q_1 + II^q_2\right) + h_{i+1}^{2}\, II^q_3\right],
\end{aligned}
$$

where, writing $s(\omega) = (1 - \omega^2)_+^{1/2}$,

| Integral | Particles | Definition |
|---|---|---|
| $II^h_1$, $II^q_1$ | of cell $i$, leaving to the right | $\int_{\omega \ge -Fr_i} (\omega + Fr_i)^{k}\, s(\omega)\, d\omega$ |
| $II^h_2$, $II^q_2$ | of cell $i$, reflected by an upward step | $\mp\int_{-Fr_i}^{\sqrt{(K^-)_+} - Fr_i} (\omega + Fr_i)^{k}\, s(\omega)\, d\omega$ |
| $II^h_3$ | of cell $i+1$, crossing to the left | $\int_{\omega \le -\sqrt{(K^+)_-} - Fr_{i+1}} (\omega + Fr_{i+1})\, s(\omega)\, d\omega$ |
| $II^q_3$ | of cell $i+1$, crossing to the left | $-\int_{\omega \le -\sqrt{(K^+)_-} - Fr_{i+1}} (\omega + Fr_{i+1})\sqrt{(\omega + Fr_{i+1})^2 + K^+}\; s(\omega)\, d\omega$ |

with $k = 1$ for the mass ($h$) and $k = 2$ for the momentum ($q$); the sign of
$II_2$ is $-$ for the mass and $+$ for the momentum. All of them have elementary
primitives except $II^q_3$ (Remark 4.1), which must be computed numerically (see
[Implementation](implementation.md#the-integral-ii3q)). The flux seen from the
right, $F^+$, follows from the same formulas by the symmetry $x \to -x$.

## The momentum interface law (Ersoy, 2026)

M. Ersoy revisited the scheme in a discontinuous Galerkin setting and identified
the macroscopic law it selects at a bottom jump. The defect between the two
one-sided momentum fluxes is exactly

$$
F^{q,-}_{i+1/2} - F^{q,+}_{i+1/2} = g\,\Delta Z_{i+1/2}\; h^{\text{jump}}_{i+1/2},
\qquad h^{\text{jump}}_{i+1/2} \ge 0,
$$

where $h^{\text{jump}}$ is an effective interface height built from three groups of
particles. With $\alpha = \sqrt{2g\lvert\Delta Z\rvert}$, $M^-$ the Maxwellian of the
left cell and $M^+$ that of the right cell, for an upward step ($\Delta Z > 0$):

$$
h^{\text{jump}} = \frac{2}{g\Delta Z}\int_0^\alpha \xi^2 M^-(\xi)\,d\xi
+ 2\int_\alpha^\infty \frac{\xi\, M^-(\xi)}{\xi + \sqrt{\xi^2 - \alpha^2}}\,d\xi
+ 2\int_0^\infty \frac{\xi\, M^+(-\xi)}{\sqrt{\xi^2 + \alpha^2} + \xi}\,d\xi :
$$

reflected particles, particles climbing the step, and particles coming down it. The
downward case swaps the roles of the two cells. At rest, $h^{\text{jump}}$ reduces to
the mean height $(h_i + h_{i+1})/2$, which is the hydrostatic pressure jump. For
polynomial degree 0 the discontinuous Galerkin scheme of that paper is exactly the
finite volume scheme above. The test suite checks this identity numerically.
