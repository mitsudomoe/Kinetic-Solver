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

/* Snapshot files and run summary. */
#ifndef KS_OUTPUT_H
#define KS_OUTPUT_H

#include <stdio.h>

#include "cases.h"
#include "grid.h"

typedef struct
{
	double l1;     /* dx * sum |e_i| */
	double linf;   /* max |e_i| */
} ErrorNorms;

/* Creates the directory if needed; returns 0 on success. */
int output_prepare_dir(const char *dir);

/* Writes <dir>/snapshot_<index>.dat with columns x z h hu, plus h_exact and
 * hu_exact when with_exact is set, and appends t to <dir>/times.dat.
 * Returns 0 on success. */
int output_snapshot(const char *dir, int index, double t, const Grid *grid, const TestCase *tc, int with_exact);

/* Errors of h and hu against the exact solution at time t, over the
 * interior cells. Returns 0 (and zero norms) when with_exact is not set or
 * the case has no exact solution. */
int output_errors(const Grid *grid, const TestCase *tc, double t, int with_exact, ErrorNorms *h_err, ErrorNorms *hu_err);

/* Smallest water height over the interior cells. */
double output_min_height(const Grid *grid);

#endif
