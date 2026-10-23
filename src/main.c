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

/* Command-line driver: runs one test case with one scheme.
 *
 *   kinetic-solver --case dam-break-dry --scheme kinetic --cells 400
 *
 * Run with --help for every option. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "cases.h"
#include "grid.h"
#include "output.h"
#include "physics.h"
#include "scheme.h"

typedef struct
{
	const char *case_name;
	const char *scheme_name;
	const char *output_dir;
	int cells;           /* 0: the case default */
	double t_end;        /* < 0: the case default */
	double cfl;          /* <= 0: the scheme default */
	int snapshot_every;  /* 0: initial and final states only */
	int snapshots;       /* write snapshot files at all */
	int verbose;
} Options;

static void usage(const char *program)
{
	printf("Usage: %s [options]\n\n", program);
	printf("  --case NAME          test case (default: lake-at-rest)\n");
	printf("  --scheme NAME        numerical scheme (default: kinetic)\n");
	printf("  --cells N            number of cells (default: per case)\n");
	printf("  --t-end T            final time in seconds (default: per case)\n");
	printf("  --cfl C              CFL number (default: per scheme)\n");
	printf("  --output-dir DIR     where snapshots and summary.json go (default: output)\n");
	printf("  --snapshot-every K   also write a snapshot every K time steps\n");
	printf("  --no-snapshots       only write summary.json\n");
	printf("  --verbose            print every time step\n");
	printf("  --list               list the available cases and schemes\n");
	printf("  --help               show this message\n");
}

static int parse_int(const char *text, int *value)
{
	char *end;
	long parsed = strtol(text, &end, 10);

	if (*text == '\0' || *end != '\0' || parsed < 0 || parsed > 100000000L)
	{
		return -1;
	}
	*value = (int)parsed;
	return 0;
}

static int parse_double(const char *text, double *value)
{
	char *end;

	*value = strtod(text, &end);
	return (*text == '\0' || *end != '\0') ? -1 : 0;
}

/* Returns 0 to run, 1 to exit successfully, -1 on error. */
static int parse_options(int argc, char **argv, Options *opt)
{
	int i;

	opt->case_name = "lake-at-rest";
	opt->scheme_name = "kinetic";
	opt->output_dir = "output";
	opt->cells = 0;
	opt->t_end = -1.0;
	opt->cfl = 0.0;
	opt->snapshot_every = 0;
	opt->snapshots = 1;
	opt->verbose = 0;

	for (i = 1; i < argc; i++)
	{
		const char *arg = argv[i];
		const char *value = i + 1 < argc ? argv[i + 1] : NULL;
		int bad = 0;

		if (strcmp(arg, "--help") == 0)
		{
			usage(argv[0]);
			return 1;
		}
		else if (strcmp(arg, "--list") == 0)
		{
			printf("Cases:\n");
			case_list();
			printf("Schemes:\n");
			scheme_list();
			return 1;
		}
		else if (strcmp(arg, "--no-snapshots") == 0)
		{
			opt->snapshots = 0;
			continue;
		}
		else if (strcmp(arg, "--verbose") == 0)
		{
			opt->verbose = 1;
			continue;
		}

		if (!value)
		{
			fprintf(stderr, "Unknown option or missing value: %s\n", arg);
			return -1;
		}
		if (strcmp(arg, "--case") == 0)
			opt->case_name = value;
		else if (strcmp(arg, "--scheme") == 0)
			opt->scheme_name = value;
		else if (strcmp(arg, "--output-dir") == 0)
			opt->output_dir = value;
		else if (strcmp(arg, "--cells") == 0)
			bad = parse_int(value, &opt->cells) != 0 || opt->cells < 2;
		else if (strcmp(arg, "--t-end") == 0)
			bad = parse_double(value, &opt->t_end) != 0 || opt->t_end <= 0.0;
		else if (strcmp(arg, "--cfl") == 0)
			bad = parse_double(value, &opt->cfl) != 0 || opt->cfl <= 0.0;
		else if (strcmp(arg, "--snapshot-every") == 0)
			bad = parse_int(value, &opt->snapshot_every) != 0;
		else
		{
			fprintf(stderr, "Unknown option: %s\n", arg);
			return -1;
		}
		if (bad)
		{
			fprintf(stderr, "Invalid value for %s: %s\n", arg, value);
			return -1;
		}
		i++;
	}
	return 0;
}

static void json_number(FILE *file, const char *key, double value, int has_value, int last)
{
	if (has_value)
		fprintf(file, "  \"%s\": %.17g%s\n", key, value, last ? "" : ",");
	else
		fprintf(file, "  \"%s\": null%s\n", key, last ? "" : ",");
}

int main(int argc, char **argv)
{
	Options opt;
	const TestCase *tc;
	const Scheme *scheme;
	Grid grid;
	Flux *workspace;
	ErrorNorms h_err, hu_err;
	double t = 0.0, t_end, cfl, mass_initial, mass_final, min_height, elapsed;
	int iterations = 0, snapshot_index = 0, has_exact, done = 0;
	clock_t start;
	char path[1024];
	FILE *summary;

	switch (parse_options(argc, argv, &opt))
	{
	case 1:
		return EXIT_SUCCESS;
	case -1:
		return EXIT_FAILURE;
	default:
		break;
	}

	tc = case_find(opt.case_name);
	scheme = scheme_find(opt.scheme_name);
	if (!tc || !scheme)
	{
		fprintf(stderr, "Unknown %s '%s' (see --list)\n", tc ? "scheme" : "case",
			tc ? opt.scheme_name : opt.case_name);
		return EXIT_FAILURE;
	}
	t_end = opt.t_end > 0.0 ? opt.t_end : tc->t_end;
	cfl = opt.cfl > 0.0 ? opt.cfl : scheme->default_cfl;

	if (grid_init(&grid, opt.cells > 0 ? opt.cells : tc->cells, tc->xmin, tc->xmax) != 0)
	{
		fprintf(stderr, "Out of memory\n");
		return EXIT_FAILURE;
	}
	workspace = malloc(2 * (size_t)(grid.n + 1) * sizeof(Flux));
	if (!workspace || output_prepare_dir(opt.output_dir) != 0)
	{
		grid_free(&grid);
		free(workspace);
		return EXIT_FAILURE;
	}

	case_initialize(tc, &grid);
	mass_initial = grid_mass(&grid);
	min_height = output_min_height(&grid);
	if (opt.snapshots && output_snapshot(opt.output_dir, snapshot_index++, t, &grid, tc) != 0)
	{
		return EXIT_FAILURE;
	}

	start = clock();
	while (!done)
	{
		double dt = scheme_time_step(scheme, &grid, cfl);

		if (t + dt >= t_end)
		{
			dt = t_end - t;
			done = 1;
		}
		scheme_step(scheme, &grid, dt, workspace);
		case_apply_boundary(tc, &grid);
		t = done ? t_end : t + dt;
		iterations++;
		min_height = ks_min(min_height, output_min_height(&grid));

		if (opt.verbose)
		{
			printf("step %6d  t = %.6f  dt = %.3e\n", iterations, t, dt);
		}
		if (opt.snapshots && (done || (opt.snapshot_every > 0 && iterations % opt.snapshot_every == 0)))
		{
			if (output_snapshot(opt.output_dir, snapshot_index++, t, &grid, tc) != 0)
			{
				return EXIT_FAILURE;
			}
		}
	}
	elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;
	mass_final = grid_mass(&grid);
	has_exact = output_errors(&grid, tc, t, &h_err, &hu_err);

	printf("case %s, scheme %s, %d cells, CFL %.2f\n", tc->name, scheme->name, grid.n, cfl);
	printf("t = %.6f after %d steps (%.2f s)\n", t, iterations, elapsed);
	printf("mass: %.15g -> %.15g (relative change %.2e)\n", mass_initial, mass_final,
		(mass_final - mass_initial) / mass_initial);
	printf("min h over the run: %.3e\n", min_height);
	if (has_exact)
	{
		printf("error h : L1 %.3e  Linf %.3e\n", h_err.l1, h_err.linf);
		printf("error hu: L1 %.3e  Linf %.3e\n", hu_err.l1, hu_err.linf);
	}

	snprintf(path, sizeof(path), "%s/summary.json", opt.output_dir);
	summary = fopen(path, "w");
	if (!summary)
	{
		perror(path);
		return EXIT_FAILURE;
	}
	fprintf(summary, "{\n  \"case\": \"%s\",\n  \"scheme\": \"%s\",\n  \"cells\": %d,\n  \"iterations\": %d,\n",
		tc->name, scheme->name, grid.n, iterations);
	json_number(summary, "cfl", cfl, 1, 0);
	json_number(summary, "t", t, 1, 0);
	json_number(summary, "wall_time", elapsed, 1, 0);
	json_number(summary, "mass_initial", mass_initial, 1, 0);
	json_number(summary, "mass_final", mass_final, 1, 0);
	json_number(summary, "min_h", min_height, 1, 0);
	json_number(summary, "h_l1", h_err.l1, has_exact, 0);
	json_number(summary, "h_linf", h_err.linf, has_exact, 0);
	json_number(summary, "hu_l1", hu_err.l1, has_exact, 0);
	json_number(summary, "hu_linf", hu_err.linf, has_exact, 1);
	fprintf(summary, "}\n");
	fclose(summary);

	free(workspace);
	grid_free(&grid);
	return EXIT_SUCCESS;
}
