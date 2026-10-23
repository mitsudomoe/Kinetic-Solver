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

#include "cases.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "physics.h"

/* ---------------------------------------------------------------------------
 * Topographies, given by their primitives
 * ------------------------------------------------------------------------- */

static double flat_bed(double x)
{
	(void)x;
	return 0.0;
}

/* z(x) = 0.2 - 0.05 x^2 on [-2, 2], 0 elsewhere. */
static double bump_bed(double x)
{
	double xc = ks_clamp(x, -2.0, 2.0);
	double inside = 0.2 * xc - xc * xc * xc / 60.0 + 0.8 / 3.0; /* 0 at x = -2 */

	return inside;
}

/* z(x) = 0.3 for x >= 0, 0 otherwise. */
static double step_bed(double x)
{
	return 0.3 * ks_max(x, 0.0);
}

/* ---------------------------------------------------------------------------
 * Initial conditions
 * ------------------------------------------------------------------------- */

static void lake_at_rest(const TestCase *tc, double x, double z, double *h, double *hu)
{
	(void)x;
	*h = ks_max(0.0, tc->level - z);
	*hu = 0.0;
}

static void dam_break(const TestCase *tc, double x, double z, double *h, double *hu)
{
	(void)z;
	*h = x < 0.0 ? tc->h_left : tc->h_right;
	*hu = 0.0;
}

static void gaussian_pulse(const TestCase *tc, double x, double z, double *h, double *hu)
{
	(void)z;
	*h = tc->level + 0.2 * exp(-x * x);
	*hu = 0.3 * *h;
}

/* ---------------------------------------------------------------------------
 * Exact solutions
 * ------------------------------------------------------------------------- */

static void lake_at_rest_exact(const TestCase *tc, double x, double t, double z, double *h, double *hu)
{
	(void)t;
	lake_at_rest(tc, x, z, h, hu);
}

/* Ritter's solution: dam break on a dry bed, dam at x = 0. */
static void ritter_exact(const TestCase *tc, double x, double t, double z, double *h, double *hu)
{
	double c_left = sqrt(KS_GRAVITY * tc->h_left);
	double xi = x / t;

	(void)z;
	if (xi <= -c_left)
	{
		*h = tc->h_left;
		*hu = 0.0;
	}
	else if (xi < 2.0 * c_left)
	{
		double c = (2.0 * c_left - xi) / 3.0;

		*h = c * c / KS_GRAVITY;
		*hu = *h * 2.0 * (xi + c_left) / 3.0;
	}
	else
	{
		*h = 0.0;
		*hu = 0.0;
	}
}

/* Middle height h_m of Stoker's solution: the rarefaction gives
 * u_m = 2 (c_l - c_m) and the shock (Rankine-Hugoniot) gives
 * u_m = (h_m - h_r) sqrt(g (h_m + h_r) / (2 h_m h_r)). Solved by bisection. */
double case_stoker_middle_height(double h_left, double h_right)
{
	double lo = h_right;
	double hi = h_left;
	int k;

	for (k = 0; k < 200; k++)
	{
		double hm = 0.5 * (lo + hi);
		double rarefaction = 2.0 * (sqrt(KS_GRAVITY * h_left) - sqrt(KS_GRAVITY * hm));
		double shock = (hm - h_right) * sqrt(KS_GRAVITY * (hm + h_right) / (2.0 * hm * h_right));

		if (rarefaction > shock)
		{
			lo = hm;
		}
		else
		{
			hi = hm;
		}
	}
	return 0.5 * (lo + hi);
}

/* Stoker's solution: dam break on a wet bed, dam at x = 0. */
static void stoker_exact(const TestCase *tc, double x, double t, double z, double *h, double *hu)
{
	double c_left = sqrt(KS_GRAVITY * tc->h_left);
	double h_mid = case_stoker_middle_height(tc->h_left, tc->h_right);
	double c_mid = sqrt(KS_GRAVITY * h_mid);
	double u_mid = 2.0 * (c_left - c_mid);
	double shock_speed = h_mid * u_mid / (h_mid - tc->h_right);
	double xi = x / t;

	(void)z;
	if (xi <= -c_left)
	{
		*h = tc->h_left;
		*hu = 0.0;
	}
	else if (xi < u_mid - c_mid)
	{
		double c = (2.0 * c_left - xi) / 3.0;

		*h = c * c / KS_GRAVITY;
		*hu = *h * 2.0 * (xi + c_left) / 3.0;
	}
	else if (xi < shock_speed)
	{
		*h = h_mid;
		*hu = h_mid * u_mid;
	}
	else
	{
		*h = tc->h_right;
		*hu = 0.0;
	}
}

/* ---------------------------------------------------------------------------
 * Registry
 * ------------------------------------------------------------------------- */

static const TestCase cases[] = {
	{"lake-at-rest", "lake at rest (H = 2) over a parabolic bump",
	 -10.0, 10.0, 1.5, 1000, BC_FREE, bump_bed, lake_at_rest, lake_at_rest_exact, 0.0, 0.0, 2.0},
	{"lake-at-rest-step", "lake at rest (H = 1) over a bottom step",
	 -10.0, 10.0, 1.5, 1000, BC_FREE, step_bed, lake_at_rest, lake_at_rest_exact, 0.0, 0.0, 1.0},
	{"lake-at-rest-island", "lake at rest (H = 0.1) with the top of the bump emerged (dry area)",
	 -10.0, 10.0, 1.5, 1000, BC_FREE, bump_bed, lake_at_rest, lake_at_rest_exact, 0.0, 0.0, 0.1},
	{"dam-break-dry", "dam break on a dry bed (h = 1 | 0), Ritter's solution",
	 -10.0, 10.0, 1.5, 1000, BC_FREE, flat_bed, dam_break, ritter_exact, 1.0, 0.0, 0.0},
	{"dam-break-wet", "dam break on a wet bed (h = 2 | 1), Stoker's solution",
	 -10.0, 10.0, 1.5, 1000, BC_FREE, flat_bed, dam_break, stoker_exact, 2.0, 1.0, 0.0},
	{"dam-break-bump", "dam break on a wet bed (h = 2 | 1) over the bump, no exact solution",
	 -10.0, 10.0, 1.5, 1000, BC_FREE, bump_bed, dam_break, NULL, 2.0, 1.0, 0.0},
	{"periodic-pulse", "moving Gaussian pulse with periodic boundaries (conservation check)",
	 -10.0, 10.0, 3.0, 1000, BC_PERIODIC, flat_bed, gaussian_pulse, NULL, 0.0, 0.0, 1.0},
	{"sloshing", "dam break in a closed tank (walls), no exact solution",
	 -10.0, 10.0, 3.0, 1000, BC_WALL, bump_bed, dam_break, NULL, 2.0, 1.0, 0.0},
};

#define NB_CASES ((int)(sizeof(cases) / sizeof(cases[0])))

const TestCase *case_find(const char *name)
{
	int i;

	for (i = 0; i < NB_CASES; i++)
	{
		if (strcmp(cases[i].name, name) == 0)
		{
			return &cases[i];
		}
	}
	return NULL;
}

void case_list(void)
{
	int i;

	for (i = 0; i < NB_CASES; i++)
	{
		printf("  %-20s %s\n", cases[i].name, cases[i].description);
	}
}

void case_initialize(const TestCase *tc, Grid *grid)
{
	int i;

	for (i = 1; i <= grid->n; i++)
	{
		double left = grid->x[i] - 0.5 * grid->dx;
		double right = grid->x[i] + 0.5 * grid->dx;

		grid->z[i] = (tc->topography(right) - tc->topography(left)) / grid->dx;
		tc->initial(tc, grid->x[i], grid->z[i], &grid->h[i], &grid->hu[i]);
	}

	/* The bottom of the ghost cells follows the boundary condition, so that
	 * a lake at rest stays at rest across the boundary. */
	if (tc->boundary == BC_PERIODIC)
	{
		grid->z[0] = grid->z[grid->n];
		grid->z[grid->n + 1] = grid->z[1];
	}
	else
	{
		grid->z[0] = grid->z[1];
		grid->z[grid->n + 1] = grid->z[grid->n];
	}
	case_apply_boundary(tc, grid);
}

void case_apply_boundary(const TestCase *tc, Grid *grid)
{
	int n = grid->n;

	switch (tc->boundary)
	{
	case BC_PERIODIC:
		grid->h[0] = grid->h[n];
		grid->hu[0] = grid->hu[n];
		grid->h[n + 1] = grid->h[1];
		grid->hu[n + 1] = grid->hu[1];
		break;
	case BC_WALL:
		grid->h[0] = grid->h[1];
		grid->hu[0] = -grid->hu[1];
		grid->h[n + 1] = grid->h[n];
		grid->hu[n + 1] = -grid->hu[n];
		break;
	case BC_FREE:
	default:
		grid->h[0] = grid->h[1];
		grid->hu[0] = grid->hu[1];
		grid->h[n + 1] = grid->h[n];
		grid->hu[n + 1] = grid->hu[n];
		break;
	}
}
