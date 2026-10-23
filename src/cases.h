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

/* Test cases: topography, initial condition, boundary condition and, when one
 * is known, the exact solution. */
#ifndef KS_CASES_H
#define KS_CASES_H

#include "grid.h"

typedef enum
{
	BC_FREE,      /* zero-gradient (transmissive) */
	BC_WALL,      /* reflective wall */
	BC_PERIODIC
} BoundaryKind;

typedef struct TestCase TestCase;

/* Primitive of the bottom height: the cell average of z over [a, b] is
 * (P(b) - P(a)) / (b - a), as in Section 3.1 of the paper. */
typedef double (*TopographyPrimitive)(double x);

/* Fills h and hu at position x, given the (cell-averaged) bottom z. */
typedef void (*InitialCondition)(const TestCase *tc, double x, double z, double *h, double *hu);

/* Exact solution at (x, t); NULL when none is available. */
typedef void (*ExactSolution)(const TestCase *tc, double x, double t, double z, double *h, double *hu);

struct TestCase
{
	const char *name;
	const char *description;
	double xmin, xmax;
	double t_end;
	int cells;
	BoundaryKind boundary;
	TopographyPrimitive topography;
	InitialCondition initial;
	ExactSolution exact;
	double h_left, h_right;   /* dam breaks */
	double level;             /* lake at rest: free surface h + z */
};

const TestCase *case_find(const char *name);
void case_list(void);

/* Sets z (including ghost cells), h and hu at t = 0. */
void case_initialize(const TestCase *tc, Grid *grid);

/* Fills the ghost cells of h and hu. */
void case_apply_boundary(const TestCase *tc, Grid *grid);

/* Stoker's middle state for a dam break on a wet bed (exposed for tests). */
double case_stoker_middle_height(double h_left, double h_right);

#endif
