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

/* Command-line driver: runs one test case with one scheme and one time
 * integrator.
 *
 *   kinetic-solver --case dam-break-dry --scheme kinetic --integrator ssp-rk2
 *
 * Run with --help for every option. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "cases.h"
#include "integrator.h"
#include "output.h"
#include "scheme.h"
#include "simulation.h"

typedef struct
{
	const char *case_name;
	const char *scheme_name;
	const char *integrator_name;
	const char *output_dir;
	SimulationConfig sim;
	int snapshot_every;  /* 0: initial and final states only */
	int snapshots;       /* write snapshot files at all */
	int verbose;
} Options;

static void usage(const char *program)
{
	printf("Usage: %s [options]\n\n", program);
	printf("  --case NAME          test case (default: lake-at-rest)\n");
	printf("  --scheme NAME        numerical scheme (default: kinetic)\n");
	printf("  --integrator NAME    time integrator (default: euler)\n");
	printf("  --cells N            number of cells (default: per case)\n");
	printf("  --t-end T            final time in seconds (default: per case)\n");
	printf("  --cfl C              CFL number (default: per scheme)\n");
	printf("  --dt DT              constant time step instead of the CFL one\n");
	printf("  --mirror             start from the mirror image x -> -x of the case\n");
	printf("  --shift-cells K      periodic cases: start from the case shifted by K cells\n");
	printf("  --output-dir DIR     where snapshots and summary.json go (default: output)\n");
	printf("  --snapshot-every K   also write a snapshot every K time steps\n");
	printf("  --no-snapshots       only write summary.json\n");
	printf("  --verbose            print every time step\n");
	printf("  --list               list the available cases, schemes and integrators\n");
	printf("  --help               show this message\n");
}

static int parse_int(const char *text, int *value, long min)
{
	char *end;
	long parsed = strtol(text, &end, 10);

	if (*text == '\0' || *end != '\0' || parsed < min || parsed > 100000000L)
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

	memset(opt, 0, sizeof(*opt));
	opt->case_name = "lake-at-rest";
	opt->scheme_name = "kinetic";
	opt->integrator_name = "euler";
	opt->output_dir = "output";
	opt->snapshots = 1;

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
		if (strcmp(arg, "--list") == 0)
		{
			printf("Cases:\n");
			case_list();
			printf("Schemes:\n");
			scheme_list();
			printf("Integrators:\n");
			integrator_list();
			return 1;
		}
		if (strcmp(arg, "--no-snapshots") == 0)
		{
			opt->snapshots = 0;
			continue;
		}
		if (strcmp(arg, "--verbose") == 0)
		{
			opt->verbose = 1;
			continue;
		}
		if (strcmp(arg, "--mirror") == 0)
		{
			opt->sim.mirror = 1;
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
		else if (strcmp(arg, "--integrator") == 0)
			opt->integrator_name = value;
		else if (strcmp(arg, "--output-dir") == 0)
			opt->output_dir = value;
		else if (strcmp(arg, "--cells") == 0)
			bad = parse_int(value, &opt->sim.cells, 2) != 0;
		else if (strcmp(arg, "--t-end") == 0)
			bad = parse_double(value, &opt->sim.t_end) != 0 || opt->sim.t_end <= 0.0;
		else if (strcmp(arg, "--cfl") == 0)
			bad = parse_double(value, &opt->sim.cfl) != 0 || opt->sim.cfl <= 0.0;
		else if (strcmp(arg, "--dt") == 0)
			bad = parse_double(value, &opt->sim.fixed_dt) != 0 || opt->sim.fixed_dt <= 0.0;
		else if (strcmp(arg, "--shift-cells") == 0)
			bad = parse_int(value, &opt->sim.shift_cells, -100000000L) != 0;
		else if (strcmp(arg, "--snapshot-every") == 0)
			bad = parse_int(value, &opt->snapshot_every, 0) != 0;
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

	opt->sim.test_case = case_find(opt->case_name);
	opt->sim.scheme = scheme_find(opt->scheme_name);
	opt->sim.integrator = integrator_find(opt->integrator_name);
	if (!opt->sim.test_case || !opt->sim.scheme || !opt->sim.integrator)
	{
		const char *kind = !opt->sim.test_case ? "case" : !opt->sim.scheme ? "scheme" : "integrator";
		const char *name = !opt->sim.test_case ? opt->case_name
			: !opt->sim.scheme ? opt->scheme_name : opt->integrator_name;

		fprintf(stderr, "Unknown %s '%s' (see --list)\n", kind, name);
		return -1;
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

static int write_summary(const Options *opt, const Simulation *sim, double elapsed)
{
	const TestCase *tc = sim->config.test_case;
	ErrorNorms h_err, hu_err;
	int has_exact = output_errors(&sim->grid, tc, sim->t, simulation_has_exact(sim), &h_err, &hu_err);
	double mass_final = grid_mass(&sim->grid);
	char path[1024];
	FILE *summary;

	printf("case %s, scheme %s, integrator %s, %d cells, CFL %.2f\n", tc->name,
		sim->config.scheme->name, sim->config.integrator->name, sim->grid.n, sim->cfl);
	printf("t = %.6f after %d steps (%.2f s)\n", sim->t, sim->iterations, elapsed);
	printf("mass: %.15g -> %.15g (relative change %.2e)\n", sim->mass_initial, mass_final,
		(mass_final - sim->mass_initial) / sim->mass_initial);
	printf("min h over the run: %.3e\n", sim->min_height);
	if (has_exact)
	{
		printf("error h : L1 %.3e  Linf %.3e\n", h_err.l1, h_err.linf);
		printf("error hu: L1 %.3e  Linf %.3e\n", hu_err.l1, hu_err.linf);
	}
	if (sim->cfl_violations > 0)
	{
		fprintf(stderr, "warning: --dt exceeded the CFL time step on %d steps\n", sim->cfl_violations);
	}

	snprintf(path, sizeof(path), "%s/summary.json", opt->output_dir);
	summary = fopen(path, "w");
	if (!summary)
	{
		perror(path);
		return -1;
	}
	fprintf(summary, "{\n  \"case\": \"%s\",\n  \"scheme\": \"%s\",\n  \"integrator\": \"%s\",\n"
		"  \"cells\": %d,\n  \"iterations\": %d,\n  \"cfl_violations\": %d,\n",
		tc->name, sim->config.scheme->name, sim->config.integrator->name, sim->grid.n,
		sim->iterations, sim->cfl_violations);
	json_number(summary, "cfl", sim->cfl, 1, 0);
	json_number(summary, "t", sim->t, 1, 0);
	json_number(summary, "wall_time", elapsed, 1, 0);
	json_number(summary, "mass_initial", sim->mass_initial, 1, 0);
	json_number(summary, "mass_final", mass_final, 1, 0);
	json_number(summary, "min_h", sim->min_height, 1, 0);
	json_number(summary, "h_l1", h_err.l1, has_exact, 0);
	json_number(summary, "h_linf", h_err.linf, has_exact, 0);
	json_number(summary, "hu_l1", hu_err.l1, has_exact, 0);
	json_number(summary, "hu_linf", hu_err.linf, has_exact, 1);
	fprintf(summary, "}\n");
	return fclose(summary) == 0 ? 0 : -1;
}

static int run(const Options *opt, Simulation *sim)
{
	int with_exact = simulation_has_exact(sim);
	int snapshot_index = 0;
	int done = 0;
	clock_t start;

	if (output_prepare_dir(opt->output_dir) != 0)
	{
		return -1;
	}
	if (opt->snapshots && output_snapshot(opt->output_dir, snapshot_index++, sim->t, &sim->grid,
			sim->config.test_case, with_exact) != 0)
	{
		return -1;
	}

	start = clock();
	while (!done)
	{
		double t_before = sim->t;

		done = simulation_step(sim);
		if (done < 0)
		{
			fprintf(stderr, "The solution became non-finite at t = %g after %d steps", sim->t, sim->iterations);
			if (sim->cfl_violations > 0)
			{
				fprintf(stderr, " (--dt exceeded the CFL time step on %d steps)", sim->cfl_violations);
			}
			fprintf(stderr, "\n");
			return -1;
		}
		if (opt->verbose)
		{
			printf("step %6d  t = %.6f  dt = %.3e\n", sim->iterations, sim->t, sim->t - t_before);
		}
		if (opt->snapshots && (done || (opt->snapshot_every > 0 && sim->iterations % opt->snapshot_every == 0)))
		{
			if (output_snapshot(opt->output_dir, snapshot_index++, sim->t, &sim->grid,
					sim->config.test_case, with_exact) != 0)
			{
				return -1;
			}
		}
	}
	return write_summary(opt, sim, (double)(clock() - start) / CLOCKS_PER_SEC);
}

int main(int argc, char **argv)
{
	Options opt;
	Simulation sim;
	int status;

	switch (parse_options(argc, argv, &opt))
	{
	case 1:
		return EXIT_SUCCESS;
	case -1:
		return EXIT_FAILURE;
	default:
		break;
	}
	if (simulation_create(&sim, &opt.sim) != 0)
	{
		return EXIT_FAILURE;
	}
	status = run(&opt, &sim);
	simulation_destroy(&sim);
	return status == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
