/* =============================================================================

	A finite volume code to solve the Saint-Venant system in an open channel.

   =============================================================================

	Copyright (C) 2013, 2026 Mathieu Besson

	This program is free software: you can redistribute
	it and/or modify it under the terms of the
	GNU General Public License as published by the Free Software Foundation, either
	version 3 of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty
	of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
	See the GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.*/

#include "kinetic.h"

#include <math.h>

#include "physics.h"
#include "quadrature.h"

/* (1 - w^2)^{1/2} and (1 - w^2)^{3/2} for w in [-1, 1]. */
static double root1(double w) { return sqrt(ks_max(0.0, 1.0 - w * w)); }
static double root3(double w) { double r = root1(w); return r * r * r; }

static double froude(double h, double hu)
{
	if (h <= KS_DRY_TOL)
	{
		return 0.0;
	}
	return hu / (h * sqrt(2.0 * KS_GRAVITY * h));
}

/* ---------------------------------------------------------------------------
 * Closed-form integrals. alpha = min(1, max(-1, -Fr)) everywhere; the bounds
 * beta are given with each integral.
 * ------------------------------------------------------------------------- */

/* Particles of the left cell leaving to the right: omega >= -Fr. */
double kinetic_IIh1(double fr)
{
	double a = ks_clamp(-fr, -1.0, 1.0);

	return root3(a) / 3.0 + 0.5 * fr * acos(a) - 0.5 * fr * a * root1(a);
}

/* Left-cell particles reflected by an upward step: -Fr <= omega <= sqrt((K-)+) - Fr. */
double kinetic_IIh2(double fr, double k_minus_pos)
{
	double a = ks_clamp(-fr, -1.0, 1.0);
	double b = ks_clamp(sqrt(k_minus_pos) - fr, -1.0, 1.0);

	return (root3(b) - root3(a)) / 3.0
		+ 0.5 * fr * (acos(b) - acos(a))
		+ 0.5 * fr * (a * root1(a) - b * root1(b));
}

/* Right-cell particles crossing to the left: omega <= -sqrt((K+)-) - Fr.
 * Note: the paper prints "- Fr beta sqrt(1 - beta^2) / 2" for the last term;
 * integrating (1 - w^2)^{1/2} over [-1, beta] gives "+", which is also what
 * matches the flux definition (3.5) numerically (tests/test_fluxes.py). */
double kinetic_IIh3(double fr, double k_plus_neg)
{
	double b = ks_clamp(-sqrt(k_plus_neg) - fr, -1.0, 1.0);

	return -root3(b) / 3.0 + 0.5 * fr * (KS_PI - acos(b)) + 0.5 * fr * b * root1(b);
}

double kinetic_IIq1(double fr)
{
	double a = ks_clamp(-fr, -1.0, 1.0);

	return root3(a) * (2.0 * fr + a) / 3.0
		+ 0.5 * acos(a) * (0.25 + fr * fr)
		+ 0.5 * a * root1(a) * (a * a / 6.0 - 5.0 / 12.0 - fr * fr);
}

double kinetic_IIq2(double fr, double k_minus_pos)
{
	double a = ks_clamp(-fr, -1.0, 1.0);
	double b = ks_clamp(sqrt(k_minus_pos) - fr, -1.0, 1.0);

	return root3(a) * (2.0 * fr + a) / 3.0 - root3(b) * (2.0 * fr + b) / 3.0
		+ 0.5 * (acos(a) - acos(b)) * (0.25 + fr * fr)
		+ 0.5 * b * root1(b) * (fr * fr - b * b / 6.0 + 5.0 / 12.0)
		- 0.5 * a * root1(a) * (fr * fr - a * a / 6.0 + 5.0 / 12.0);
}

/* ---------------------------------------------------------------------------
 * II3q has no elementary primitive (Remark 4.1 of the paper):
 *     II3q = - int_{-1}^{beta} (w + Fr) sqrt((w + Fr)^2 + K+) sqrt(1 - w^2) dw,
 *     beta = max(-1, min(-sqrt((K+)-) - Fr, 1)).
 * When the bottom steps down towards the left cell (K+ < 0) the radicand
 * (w + Fr)^2 + K+ vanishes exactly at w = beta, so the integrand has
 * square-root singularities at both ends. The midpoint rule of the paper then
 * converges like N^{-3/2}, which limits how well the lake at rest is kept;
 * tanh-sinh handles both singularities with ~100 nodes at full precision.
 * ------------------------------------------------------------------------- */
typedef struct
{
	double fr;
	double beta;
	double shift;    /* beta + Fr */
	double radicand; /* (beta + Fr)^2 + K+, exactly 0 when beta is the root */
} IIq3Context;

static double iiq3_integrand(double w, double from_a, double to_b, void *raw)
{
	const IIq3Context *ctx = raw;
	double one_plus_w = from_a;                 /* a = -1 */
	double one_minus_w = to_b + (1.0 - ctx->beta);
	double v = ctx->shift - to_b;               /* w + Fr */
	/* (w + Fr)^2 + K+ expanded around beta to avoid cancellation. */
	double radicand = ctx->radicand - 2.0 * ctx->shift * to_b + to_b * to_b;

	(void)w;
	return -v * sqrt(ks_max(0.0, radicand)) * sqrt(one_plus_w * one_minus_w);
}

double kinetic_IIq3(double fr, double k_plus_neg, double k_plus)
{
	double root = -sqrt(k_plus_neg) - fr;
	IIq3Context ctx;

	ctx.fr = fr;
	ctx.beta = ks_clamp(root, -1.0, 1.0);
	ctx.shift = ctx.beta + fr;
	if (k_plus_neg > 0.0 && ctx.beta == root)
	{
		/* beta is the zero of the radicand: (beta + Fr)^2 = (K+)- = -K+. */
		ctx.radicand = 0.0;
	}
	else
	{
		ctx.radicand = ks_max(0.0, ctx.shift * ctx.shift + k_plus);
	}
	return quad_tanh_sinh(iiq3_integrand, &ctx, -1.0, ctx.beta);
}

/* ---------------------------------------------------------------------------
 * Fluxes
 * ------------------------------------------------------------------------- */

/* IF^-_{i+1/2}: flux at the interface between `self` (cell i) and its
 * neighbour (cell i+1), as seen from `self`. */
static Flux one_sided_flux(CellState self, CellState other)
{
	const double c_h = 2.0 * sqrt(2.0 * KS_GRAVITY) / KS_PI;
	const double c_q = 4.0 * KS_GRAVITY / KS_PI;
	double dz = other.z - self.z;
	Flux flux = {0.0, 0.0};

	if (self.h > KS_DRY_TOL)
	{
		double fr = froude(self.h, self.hu);
		double k_minus_pos = ks_max(dz, 0.0) / self.h;

		flux.h += self.h * sqrt(self.h) * (kinetic_IIh1(fr) + kinetic_IIh2(fr, k_minus_pos));
		flux.hu += self.h * self.h * (kinetic_IIq1(fr) + kinetic_IIq2(fr, k_minus_pos));
	}
	if (other.h > KS_DRY_TOL)
	{
		double fr = froude(other.h, other.hu);
		double k_plus = dz / other.h;
		double k_plus_neg = ks_max(-dz, 0.0) / other.h;

		flux.h += other.h * sqrt(other.h) * kinetic_IIh3(fr, k_plus_neg);
		flux.hu += other.h * other.h * kinetic_IIq3(fr, k_plus_neg, k_plus);
	}
	flux.h *= c_h;
	flux.hu *= c_q;
	return flux;
}

void kinetic_interface_flux(CellState left, CellState right, Flux *seen_by_left, Flux *seen_by_right)
{
	CellState mirrored_right = {right.h, -right.hu, right.z};
	CellState mirrored_left = {left.h, -left.hu, left.z};
	Flux mirrored;

	*seen_by_left = one_sided_flux(left, right);

	/* IF^+_{i+1/2} seen from cell i+1 is IF^- of the mirrored problem
	 * (x -> -x, u -> -u, xi -> -xi): cell i+1 becomes the "self" cell with
	 * cell i on its right. The mass flux changes sign, the momentum flux
	 * (an integral of xi^2) does not. */
	mirrored = one_sided_flux(mirrored_right, mirrored_left);
	seen_by_right->h = -mirrored.h;
	seen_by_right->hu = mirrored.hu;
}

double kinetic_wave_speed(double h, double hu)
{
	if (h <= KS_DRY_TOL)
	{
		return 0.0;
	}
	return fabs(hu / h) + sqrt(2.0 * KS_GRAVITY * h);
}
