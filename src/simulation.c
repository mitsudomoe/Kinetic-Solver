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

#include "simulation.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "output.h"
#include "physics.h"

/* x -> -x, u -> -u. The domains are symmetric, so cell i maps exactly onto
 * cell n + 1 - i (ghost cells included). */
static void mirror_state(Grid *grid)
{
	int i, j;

	for (i = 0, j = grid->n + 1; i < j; i++, j--)
	{
		double h = grid->h[i], hu = grid->hu[i], z = grid->z[i];

		grid->h[i] = grid->h[j];
		grid->hu[i] = -grid->hu[j];
		grid->z[i] = grid->z[j];
		grid->h[j] = h;
		grid->hu[j] = -hu;
		grid->z[j] = z;
	}
	if (grid->n % 2 == 0)
	{
		return;
	}
	grid->hu[(grid->n + 1) / 2] = -grid->hu[(grid->n + 1) / 2];
}

/* Cyclic shift of the interior cells by k cells to the right. */
static int shift_state(Grid *grid, int k)
{
	int n = grid->n;
	double *copy = malloc((size_t)n * 3 * sizeof(double));
	int i;

	if (!copy)
	{
		return -1;
	}
	k = ((k % n) + n) % n;
	for (i = 0; i < n; i++)
	{
		copy[i] = grid->h[i + 1];
		copy[n + i] = grid->hu[i + 1];
		copy[2 * n + i] = grid->z[i + 1];
	}
	for (i = 0; i < n; i++)
	{
		int target = (i + k) % n + 1;

		grid->h[target] = copy[i];
		grid->hu[target] = copy[n + i];
		grid->z[target] = copy[2 * n + i];
	}
	free(copy);
	return 0;
}

static int state_is_finite(const Grid *grid)
{
	int i;

	for (i = 1; i <= grid->n; i++)
	{
		if (!isfinite(grid->h[i]) || !isfinite(grid->hu[i]))
		{
			return 0;
		}
	}
	return 1;
}

static void balance(Grid *grid, double *dh, double *dhu, void *ctx)
{
	Simulation *sim = ctx;

	case_apply_boundary(sim->config.test_case, grid);
	scheme_flux_balance(sim->config.scheme, grid, dh, dhu, sim->fluxes);
}

int simulation_create(Simulation *sim, const SimulationConfig *config)
{
	const TestCase *tc = config->test_case;
	size_t size;

	memset(sim, 0, sizeof(*sim));
	sim->config = *config;
	sim->t_end = config->t_end > 0.0 ? config->t_end : tc->t_end;
	sim->cfl = config->cfl > 0.0 ? config->cfl : config->scheme->default_cfl;

	if (config->mirror && tc->xmin != -tc->xmax)
	{
		fprintf(stderr, "--mirror needs a domain symmetric about x = 0\n");
		return -1;
	}
	if (config->shift_cells != 0 && tc->boundary != BC_PERIODIC)
	{
		fprintf(stderr, "--shift-cells needs periodic boundaries\n");
		return -1;
	}
	if (grid_init(&sim->grid, config->cells > 0 ? config->cells : tc->cells, tc->xmin, tc->xmax) != 0)
	{
		fprintf(stderr, "Out of memory\n");
		return -1;
	}

	size = (size_t)(sim->grid.n + 2);
	sim->fluxes = malloc(2 * (size_t)(sim->grid.n + 1) * sizeof(Flux));
	sim->work.h0 = calloc(size, sizeof(double));
	sim->work.hu0 = calloc(size, sizeof(double));
	sim->work.dh = calloc(size, sizeof(double));
	sim->work.dhu = calloc(size, sizeof(double));
	if (!sim->fluxes || !sim->work.h0 || !sim->work.hu0 || !sim->work.dh || !sim->work.dhu)
	{
		fprintf(stderr, "Out of memory\n");
		simulation_destroy(sim);
		return -1;
	}

	case_initialize(tc, &sim->grid);
	if (config->mirror)
	{
		mirror_state(&sim->grid);
	}
	if (config->shift_cells != 0)
	{
		if (shift_state(&sim->grid, config->shift_cells) != 0)
		{
			fprintf(stderr, "Out of memory\n");
			simulation_destroy(sim);
			return -1;
		}
		case_fill_ghost_bottom(tc, &sim->grid);
	}
	case_apply_boundary(tc, &sim->grid);

	sim->mass_initial = grid_mass(&sim->grid);
	sim->min_height = output_min_height(&sim->grid);
	return 0;
}

void simulation_destroy(Simulation *sim)
{
	grid_free(&sim->grid);
	free(sim->fluxes);
	free(sim->work.h0);
	free(sim->work.hu0);
	free(sim->work.dh);
	free(sim->work.dhu);
	sim->fluxes = NULL;
	sim->work.h0 = sim->work.hu0 = sim->work.dh = sim->work.dhu = NULL;
}

int simulation_step(Simulation *sim)
{
	double cfl_dt, dt;
	int done = 0;

	if (sim->t >= sim->t_end)
	{
		return 1;
	}
	case_apply_boundary(sim->config.test_case, &sim->grid);
	cfl_dt = scheme_time_step(sim->config.scheme, &sim->grid, sim->cfl);
	dt = cfl_dt;
	if (sim->config.fixed_dt > 0.0)
	{
		dt = sim->config.fixed_dt;
		if (dt > cfl_dt)
		{
			sim->cfl_violations++;
		}
	}
	if (sim->t + dt >= sim->t_end)
	{
		dt = sim->t_end - sim->t;
		done = 1;
	}

	integrator_step(sim->config.integrator, &sim->grid, dt, balance, sim, &sim->work);
	case_apply_boundary(sim->config.test_case, &sim->grid);

	sim->t = done ? sim->t_end : sim->t + dt;
	sim->iterations++;
	if (!state_is_finite(&sim->grid))
	{
		return -1;
	}
	sim->min_height = ks_min(sim->min_height, output_min_height(&sim->grid));
	return done;
}

int simulation_has_exact(const Simulation *sim)
{
	return sim->config.test_case->exact != NULL && !sim->config.mirror && sim->config.shift_cells == 0;
}
