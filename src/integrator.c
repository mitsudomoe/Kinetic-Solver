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

#include "integrator.h"

#include <stdio.h>
#include <string.h>

#include "physics.h"

static const ShuOsherStage euler_stages[] = {{0.0, 1.0}};

/* Heun's method; SSP coefficient 1 (Shu & Osher, 1988). */
static const ShuOsherStage ssp_rk2_stages[] = {{0.0, 1.0}, {0.5, 0.5}};

/* Third-order SSP Runge-Kutta; SSP coefficient 1 (Shu & Osher, 1988). */
static const ShuOsherStage ssp_rk3_stages[] = {{0.0, 1.0}, {0.75, 0.25}, {1.0 / 3.0, 2.0 / 3.0}};

static const Integrator integrators[] = {
	{"euler", "forward Euler (the scheme of the paper)", 1, 1, euler_stages},
	{"ssp-rk2", "second-order SSP Runge-Kutta (Heun)", 2, 2, ssp_rk2_stages},
	{"ssp-rk3", "third-order SSP Runge-Kutta (Shu-Osher)", 3, 3, ssp_rk3_stages},
};

#define NB_INTEGRATORS ((int)(sizeof(integrators) / sizeof(integrators[0])))

const Integrator *integrator_find(const char *name)
{
	int i;

	for (i = 0; i < NB_INTEGRATORS; i++)
	{
		if (strcmp(integrators[i].name, name) == 0)
		{
			return &integrators[i];
		}
	}
	return NULL;
}

void integrator_list(void)
{
	int i;

	for (i = 0; i < NB_INTEGRATORS; i++)
	{
		printf("  %-10s %s\n", integrators[i].name, integrators[i].description);
	}
}

void integrator_step(const Integrator *integrator, Grid *grid, double dt,
	BalanceFn balance, void *ctx, IntegratorWorkspace *work)
{
	double ratio = dt / grid->dx;
	int n = grid->n;
	int k, i;

	if (integrator->stages > 1)
	{
		memcpy(work->h0, grid->h, (size_t)(n + 2) * sizeof(double));
		memcpy(work->hu0, grid->hu, (size_t)(n + 2) * sizeof(double));
	}

	for (k = 0; k < integrator->stages; k++)
	{
		const ShuOsherStage *stage = &integrator->stage[k];

		balance(grid, work->dh, work->dhu, ctx);
		for (i = 1; i <= n; i++)
		{
			double h = grid->h[i] - ratio * work->dh[i];
			double hu = grid->hu[i] - ratio * work->dhu[i];

			/* a = 0 (first stage, forward Euler) is kept exact on purpose. */
			if (stage->a != 0.0)
			{
				h = stage->a * work->h0[i] + stage->b * h;
				hu = stage->a * work->hu0[i] + stage->b * hu;
			}

			/* The scheme keeps h >= 0 under its CFL condition; this only
			 * removes round-off level negative values and stale momentum in
			 * cells that just dried out. */
			if (h <= KS_DRY_TOL)
			{
				h = ks_max(h, 0.0);
				hu = 0.0;
			}
			grid->h[i] = h;
			grid->hu[i] = hu;
		}
	}
}
