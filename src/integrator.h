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

/* Explicit time integrators for dU/dt = L(U), written in Shu-Osher form:
 *
 *     U^(0) = U^n
 *     U^(k) = a_k U^n + b_k (U^(k-1) + dt L(U^(k-1))),   k = 1..s
 *     U^{n+1} = U^(s)
 *
 * with a_k + b_k = 1 and a_k, b_k >= 0. Every stage is a convex combination
 * of forward Euler steps, so these strong-stability-preserving (SSP) methods
 * keep the properties the scheme has under forward Euler (positivity, lake
 * at rest, entropy inequality) under the same CFL condition (SSP coefficient
 * 1 for all three methods below). */
#ifndef KS_INTEGRATOR_H
#define KS_INTEGRATOR_H

#include "grid.h"

/* Evaluates D(U) = -dx L(U) for the current state of the grid, after filling
 * its ghost cells. dh and dhu receive D on the interior cells. */
typedef void (*BalanceFn)(Grid *grid, double *dh, double *dhu, void *ctx);

typedef struct
{
	double a; /* weight of U^n */
	double b; /* weight of the Euler step from the previous stage */
} ShuOsherStage;

typedef struct
{
	const char *name;
	const char *description;
	int order;                    /* formal order of accuracy in time */
	int stages;
	const ShuOsherStage *stage;
} Integrator;

typedef struct
{
	double *h0, *hu0; /* copy of U^n */
	double *dh, *dhu; /* D(U) of the current stage */
} IntegratorWorkspace;

/* Looks an integrator up by name; NULL if unknown. */
const Integrator *integrator_find(const char *name);
void integrator_list(void);

/* Advances the interior cells of the grid from t to t + dt. */
void integrator_step(const Integrator *integrator, Grid *grid, double dt,
	BalanceFn balance, void *ctx, IntegratorWorkspace *work);

#endif
