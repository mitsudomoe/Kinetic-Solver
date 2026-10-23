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

static double flat_bed(const TestCase *tc, double x)
{
	(void)tc;
	(void)x;
	return 0.0;
}

/* z(x) = 0.2 - 0.05 x^2 on [-2, 2], 0 elsewhere. */
static double bump_height(double x)
{
	return fabs(x) < 2.0 ? 0.2 - 0.05 * x * x : 0.0;
}

static double bump_bed(const TestCase *tc, double x)
{
	double xc = ks_clamp(x, -2.0, 2.0);
	double inside = 0.2 * xc - xc * xc * xc / 60.0 + 0.8 / 3.0; /* 0 at x = -2 */

	(void)tc;
	return inside;
}

/* z(x) = 0.3 for x >= 0, 0 otherwise. */
static double step_bed(const TestCase *tc, double x)
{
	(void)tc;
	return 0.3 * ks_max(x, 0.0);
}

/* Parabolic bowl of Thacker's solution: z(x) = h0 (x^2 / a^2 - 1). */
static double bowl_bed(const TestCase *tc, double x)
{
	return tc->depth * (x * x * x / (3.0 * tc->radius * tc->radius) - x);
}

/* Plane beach: z(x) = s x for x >= 0, 0 otherwise. */
#define BEACH_SLOPE 0.1

static double beach_bed(const TestCase *tc, double x)
{
	double xp = ks_max(x, 0.0);

	(void)tc;
	return 0.5 * BEACH_SLOPE * xp * xp;
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
 * Thacker's oscillation in a parabolic bowl
 *
 * W. C. Thacker, "Some exact solutions to the nonlinear shallow-water wave
 * equations", J. Fluid Mech. 107 (1981). In the bowl z = h0 (x^2/a^2 - 1)
 * the free surface stays a plane that rocks back and forth while the
 * velocity is uniform over the wet region:
 *     u(t) = B sin(w t),  eta(x, t) = -(B w / g) cos(w t) x + B^2 sin^2(w t) / (2 g),
 *     h = max(0, eta - z),  w = sqrt(2 g h0) / a.
 * Both shorelines move: a hard test for wet/dry fronts.
 * ------------------------------------------------------------------------- */

static void thacker_exact(const TestCase *tc, double x, double t, double z, double *h, double *hu)
{
	double omega = sqrt(2.0 * KS_GRAVITY * tc->depth) / tc->radius;
	double b = tc->amplitude;
	double s = sin(omega * t);
	double eta = -(b * omega / KS_GRAVITY) * cos(omega * t) * x + b * b * s * s / (2.0 * KS_GRAVITY);

	*h = ks_max(0.0, eta - z);
	*hu = *h * b * s;
}

static void thacker_initial(const TestCase *tc, double x, double z, double *h, double *hu)
{
	thacker_exact(tc, x, 0.0, z, h, hu);
}

/* ---------------------------------------------------------------------------
 * Steady flows over the bump
 *
 * Classic benchmarks (Goutal & Maurel, 1997; SWASHES, Delestre et al., 2013,
 * Section 3.1, with the bump centred at x = 0 instead of 10):
 * water enters on the left with discharge q and the flow settles to a steady
 * state where q is uniform and the Bernoulli head q^2/(2 g h^2) + h + z is
 * constant along each smooth branch. The critical height is
 * h_c = (q^2 / g)^(1/3).
 * ------------------------------------------------------------------------- */

static double bernoulli_residual(double q, double head, double z, double h)
{
	return q * q / (2.0 * KS_GRAVITY * h * h) + h + z - head;
}

double case_bernoulli_height(double q, double head, double z, int supercritical)
{
	double hc = cbrt(q * q / KS_GRAVITY);
	double lo, hi;
	int k;

	/* The residual is minimal at h_c: no root if the head is too low there
	 * (only at the crest of a transcritical flow, up to round-off). */
	if (bernoulli_residual(q, head, z, hc) >= 0.0)
	{
		return hc;
	}
	if (supercritical)
	{
		lo = 1e-12;
		hi = hc;
	}
	else
	{
		lo = hc;
		hi = head - z + 1.0;
	}
	for (k = 0; k < 200; k++)
	{
		double mid = 0.5 * (lo + hi);
		double r = bernoulli_residual(q, head, z, mid);

		/* The residual decreases on the supercritical branch, increases on
		 * the subcritical one. */
		if ((r > 0.0) == (supercritical != 0))
		{
			lo = mid;
		}
		else
		{
			hi = mid;
		}
	}
	return 0.5 * (lo + hi);
}

static double head_at_crest(double q)
{
	double hc = cbrt(q * q / KS_GRAVITY);

	return 1.5 * hc + bump_height(0.0);
}

static double head_at_outlet(const TestCase *tc)
{
	double q = tc->discharge;

	return q * q / (2.0 * KS_GRAVITY * tc->h_out * tc->h_out) + tc->h_out;
}

/* Fully subcritical flow: one Bernoulli head, fixed by the outlet. */
static void subcritical_exact(const TestCase *tc, double x, double t, double z, double *h, double *hu)
{
	(void)x;
	(void)t;
	*h = case_bernoulli_height(tc->discharge, head_at_outlet(tc), z, 0);
	*hu = tc->discharge;
}

/* Transcritical flow: subcritical upstream, critical at the crest,
 * supercritical downstream. */
static void transcritical_exact(const TestCase *tc, double x, double t, double z, double *h, double *hu)
{
	(void)t;
	*h = case_bernoulli_height(tc->discharge, head_at_crest(tc->discharge), z, x > 0.0);
	*hu = tc->discharge;
}

/* Conjugate (sequent) depth across a hydraulic jump. */
static double conjugate_height(double q, double h)
{
	double froude2 = q * q / (KS_GRAVITY * h * h * h);

	return 0.5 * h * (sqrt(1.0 + 8.0 * froude2) - 1.0);
}

/* Position of the hydraulic jump downstream of the crest: the supercritical
 * height jumps to the subcritical height fixed by the outlet where the
 * Rankine-Hugoniot condition (stationary shock) holds; x = 1.666 here. */
static double jump_position(const TestCase *tc)
{
	double q = tc->discharge;
	double crest = head_at_crest(q);
	double outlet = head_at_outlet(tc);
	double lo = 1e-9, hi = 2.0;
	int k;

	for (k = 0; k < 100; k++)
	{
		double mid = 0.5 * (lo + hi);
		double z = bump_height(mid);
		double upstream = case_bernoulli_height(q, crest, z, 1);
		double downstream = case_bernoulli_height(q, outlet, z, 0);

		/* Close to the crest the outlet head is too low for a subcritical
		 * flow to exist (no root): the jump must be further downstream.
		 * Beyond that, the conjugate depth of the supercritical flow
		 * decreases while the subcritical height increases; the jump sits
		 * where they meet. */
		if (bernoulli_residual(q, outlet, z, cbrt(q * q / KS_GRAVITY)) >= 0.0
			|| conjugate_height(q, upstream) > downstream)
		{
			lo = mid;
		}
		else
		{
			hi = mid;
		}
	}
	return 0.5 * (lo + hi);
}

static void hydraulic_jump_exact(const TestCase *tc, double x, double t, double z, double *h, double *hu)
{
	double q = tc->discharge;

	(void)t;
	if (x <= 0.0)
	{
		*h = case_bernoulli_height(q, head_at_crest(q), z, 0);
	}
	else if (x < jump_position(tc))
	{
		*h = case_bernoulli_height(q, head_at_crest(q), z, 1);
	}
	else
	{
		*h = case_bernoulli_height(q, head_at_outlet(tc), z, 0);
	}
	*hu = q;
}

/* ---------------------------------------------------------------------------
 * Wave running up a beach
 *
 * A Gaussian hump travelling to the right (simple wave: u = 2 (c - c0)) over
 * a flat bottom, which steepens into a bore, climbs the dry beach and
 * flows back. No exact solution.
 * ------------------------------------------------------------------------- */

static void runup_initial(const TestCase *tc, double x, double z, double *h, double *hu)
{
	double still = ks_max(0.0, tc->level - z);
	double hump = z > 0.0 ? 0.0 : tc->amplitude * exp(-(x + 5.0) * (x + 5.0));

	*h = still + hump;
	*hu = *h * 2.0 * (sqrt(KS_GRAVITY * *h) - sqrt(KS_GRAVITY * still));
}

/* ---------------------------------------------------------------------------
 * Registry
 * ------------------------------------------------------------------------- */

static const TestCase cases[] = {
	{.name = "lake-at-rest", .description = "lake at rest (H = 2) over a parabolic bump",
	 .xmin = -10.0, .xmax = 10.0, .t_end = 1.5, .cells = 1000, .boundary = BC_FREE,
	 .topography = bump_bed, .initial = lake_at_rest, .exact = lake_at_rest_exact, .level = 2.0},
	{.name = "lake-at-rest-step", .description = "lake at rest (H = 1) over a bottom step",
	 .xmin = -10.0, .xmax = 10.0, .t_end = 1.5, .cells = 1000, .boundary = BC_FREE,
	 .topography = step_bed, .initial = lake_at_rest, .exact = lake_at_rest_exact, .level = 1.0},
	{.name = "lake-at-rest-island",
	 .description = "lake at rest (H = 0.1) with the top of the bump emerged (dry area)",
	 .xmin = -10.0, .xmax = 10.0, .t_end = 1.5, .cells = 1000, .boundary = BC_FREE,
	 .topography = bump_bed, .initial = lake_at_rest, .exact = lake_at_rest_exact, .level = 0.1},
	{.name = "dam-break-dry", .description = "dam break on a dry bed (h = 1 | 0), Ritter's solution",
	 .xmin = -10.0, .xmax = 10.0, .t_end = 1.5, .cells = 1000, .boundary = BC_FREE,
	 .topography = flat_bed, .initial = dam_break, .exact = ritter_exact, .h_left = 1.0, .h_right = 0.0},
	{.name = "dam-break-wet", .description = "dam break on a wet bed (h = 2 | 1), Stoker's solution",
	 .xmin = -10.0, .xmax = 10.0, .t_end = 1.5, .cells = 1000, .boundary = BC_FREE,
	 .topography = flat_bed, .initial = dam_break, .exact = stoker_exact, .h_left = 2.0, .h_right = 1.0},
	{.name = "dam-break-bump", .description = "dam break on a wet bed (h = 2 | 1) over the bump, no exact solution",
	 .xmin = -10.0, .xmax = 10.0, .t_end = 1.5, .cells = 1000, .boundary = BC_FREE,
	 .topography = bump_bed, .initial = dam_break, .h_left = 2.0, .h_right = 1.0},
	{.name = "periodic-pulse", .description = "moving Gaussian pulse with periodic boundaries (conservation check)",
	 .xmin = -10.0, .xmax = 10.0, .t_end = 3.0, .cells = 1000, .boundary = BC_PERIODIC,
	 .topography = flat_bed, .initial = gaussian_pulse, .level = 1.0},
	{.name = "sloshing", .description = "dam break in a closed tank (walls), no exact solution",
	 .xmin = -10.0, .xmax = 10.0, .t_end = 3.0, .cells = 1000, .boundary = BC_WALL,
	 .topography = bump_bed, .initial = dam_break, .h_left = 2.0, .h_right = 1.0},
	{.name = "thacker", .description = "oscillation in a parabolic bowl, moving shorelines (Thacker's solution)",
	 .xmin = -10.0, .xmax = 10.0, .t_end = 16.05, .cells = 400, .boundary = BC_WALL,
	 .topography = bowl_bed, .initial = thacker_initial, .exact = thacker_exact,
	 .depth = 0.5, .radius = 8.0, .amplitude = 0.6},
	{.name = "bump-subcritical", .description = "steady subcritical flow over the bump (q = 4.42, h_out = 2)",
	 .xmin = -10.0, .xmax = 15.0, .t_end = 100.0, .cells = 400, .boundary = BC_INFLOW_OUTFLOW,
	 .topography = bump_bed, .initial = lake_at_rest, .exact = subcritical_exact,
	 .level = 2.0, .discharge = 4.42, .h_out = 2.0},
	{.name = "bump-transcritical",
	 .description = "steady transcritical flow over the bump, critical at the crest (q = 1.53)",
	 .xmin = -10.0, .xmax = 15.0, .t_end = 100.0, .cells = 400, .boundary = BC_INFLOW_OUTFLOW,
	 .topography = bump_bed, .initial = lake_at_rest, .exact = transcritical_exact,
	 .level = 0.66, .discharge = 1.53, .h_out = 0.66},
	{.name = "bump-hydraulic-jump",
	 .description = "steady transcritical flow over the bump with a hydraulic jump (q = 0.18, h_out = 0.33)",
	 .xmin = -10.0, .xmax = 15.0, .t_end = 200.0, .cells = 400, .boundary = BC_INFLOW_OUTFLOW,
	 .topography = bump_bed, .initial = lake_at_rest, .exact = hydraulic_jump_exact,
	 .level = 0.33, .discharge = 0.18, .h_out = 0.33},
	{.name = "beach-runup", .description = "wave steepening into a bore and running up a dry beach",
	 .xmin = -10.0, .xmax = 15.0, .t_end = 12.0, .cells = 1000, .boundary = BC_WALL,
	 .topography = beach_bed, .initial = runup_initial, .level = 0.5, .amplitude = 0.15},
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

		grid->z[i] = (tc->topography(tc, right) - tc->topography(tc, left)) / grid->dx;
		tc->initial(tc, grid->x[i], grid->z[i], &grid->h[i], &grid->hu[i]);
	}
	case_fill_ghost_bottom(tc, grid);
	case_apply_boundary(tc, grid);
}

void case_fill_ghost_bottom(const TestCase *tc, Grid *grid)
{
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
	case BC_INFLOW_OUTFLOW:
	{
		/* One characteristic enters at a subcritical outlet, none at a
		 * supercritical one: impose the height only in the first case. */
		double u_out = ks_velocity(grid->h[n], grid->hu[n]);
		int supercritical = u_out * u_out > KS_GRAVITY * grid->h[n];

		grid->h[0] = grid->h[1];
		grid->hu[0] = tc->discharge;
		grid->h[n + 1] = (tc->h_out > 0.0 && !supercritical) ? tc->h_out : grid->h[n];
		grid->hu[n + 1] = grid->hu[n];
		break;
	}
	case BC_FREE:
	default:
		grid->h[0] = grid->h[1];
		grid->hu[0] = grid->hu[1];
		grid->h[n + 1] = grid->h[n];
		grid->hu[n + 1] = grid->hu[n];
		break;
	}
}
