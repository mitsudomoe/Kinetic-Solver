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

/* A simulation: one test case advanced in time by one scheme and one time
 * integrator. The object owns every buffer it allocates; create it with
 * simulation_create and release it with simulation_destroy, which is safe to
 * call on a partially created simulation. */
#ifndef KS_SIMULATION_H
#define KS_SIMULATION_H

#include "cases.h"
#include "grid.h"
#include "integrator.h"
#include "scheme.h"

typedef struct
{
	const TestCase *test_case;
	const Scheme *scheme;
	const Integrator *integrator;
	int cells;          /* 0: the case default */
	double t_end;       /* <= 0: the case default */
	double cfl;         /* <= 0: the scheme default */
	double fixed_dt;    /* > 0: constant time step instead of the CFL one */
	int mirror;         /* start from the mirror image x -> -x of the case */
	int shift_cells;    /* periodic cases: start from the case shifted by k cells */
} SimulationConfig;

typedef struct
{
	SimulationConfig config;
	Grid grid;
	double t;
	double t_end;
	double cfl;
	int iterations;
	int cfl_violations;  /* steps where fixed_dt exceeded the CFL time step */
	double min_height;   /* over all the time steps so far */
	double mass_initial;
	Flux *fluxes;
	IntegratorWorkspace work;
} Simulation;

/* Returns 0 on success; on failure prints the reason, releases everything
 * and returns -1. */
int simulation_create(Simulation *sim, const SimulationConfig *config);
void simulation_destroy(Simulation *sim);

/* Advances by one time step. Returns 0 to continue, 1 once t_end has been
 * reached, and -1 if the solution stopped being finite (the simulation has
 * diverged, typically because --dt exceeded the CFL time step). */
int simulation_step(Simulation *sim);

/* True when the state can be compared with the case's exact solution (the
 * case has one and the initial state was not mirrored or shifted). */
int simulation_has_exact(const Simulation *sim);

#endif
