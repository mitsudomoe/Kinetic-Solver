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

/* Tanh-sinh (double exponential) quadrature on a finite interval.
 *
 * The rule converges exponentially even when the integrand has algebraic
 * endpoint singularities such as sqrt(x - a), which is exactly the
 * situation of the kinetic flux integral II3q. To avoid cancellation near
 * the ends, the integrand receives the distances to both endpoints along
 * with the abscissa. */
#ifndef KS_QUADRATURE_H
#define KS_QUADRATURE_H

/* Integrand: x is the abscissa, from_a = x - a and to_b = b - x, both
 * computed without cancellation. */
typedef double (*QuadIntegrand)(double x, double from_a, double to_b, void *ctx);

/* Integral of f over [a, b] (returns 0 when b <= a). */
double quad_tanh_sinh(QuadIntegrand f, void *ctx, double a, double b);

/* Number of nodes used by the rule (for documentation and tests). */
int quad_tanh_sinh_nodes(void);

#endif
