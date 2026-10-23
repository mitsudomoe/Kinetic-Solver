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

/* Physical constants and small numeric helpers shared by every module. */
#ifndef KS_PHYSICS_H
#define KS_PHYSICS_H

#include <math.h>

/* Full-precision pi. The closed-form flux integrals mix this constant with
 * acos(), which returns the exact value: any truncation here breaks the
 * lake-at-rest property of the kinetic scheme. */
#define KS_PI 3.14159265358979323846

/* Gravitational acceleration (m/s^2). */
#define KS_GRAVITY 9.81

/* Water heights at or below this value are treated as dry: zero velocity,
 * zero Froude number. */
#define KS_DRY_TOL 1e-10

static inline double ks_max(double a, double b) { return a > b ? a : b; }
static inline double ks_min(double a, double b) { return a < b ? a : b; }
static inline double ks_clamp(double x, double lo, double hi) { return ks_min(hi, ks_max(lo, x)); }

/* Velocity of a cell, zero when the cell is dry. */
static inline double ks_velocity(double h, double hu)
{
	return h > KS_DRY_TOL ? hu / h : 0.0;
}

#endif
