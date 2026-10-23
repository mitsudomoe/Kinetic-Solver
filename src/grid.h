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

/* Uniform 1-D finite volume grid and the conserved state stored on it.
 *
 * Conventions (shared by every module, following Perthame & Simeoni):
 *   - the domain [xmin, xmax] is split into n cells of width dx = (xmax - xmin) / n;
 *   - interior cells are i = 1..n, cell C_i = [xmin + (i-1) dx, xmin + i dx),
 *     centred at x_i = xmin + (i - 1/2) dx;
 *   - cells 0 and n+1 are ghost cells filled by the boundary conditions;
 *   - z holds the bottom height as a cell average, h the water height and
 *     hu the discharge. */
#ifndef KS_GRID_H
#define KS_GRID_H

typedef struct
{
	int n;          /* number of interior cells */
	double xmin;    /* left end of the domain */
	double xmax;    /* right end of the domain */
	double dx;      /* cell width */
	double *x;      /* cell centres, size n + 2 */
	double *z;      /* bottom height, size n + 2 */
	double *h;      /* water height, size n + 2 */
	double *hu;     /* discharge, size n + 2 */
} Grid;

/* Allocates a grid; returns 0 on success, -1 on allocation failure. */
int grid_init(Grid *grid, int n, double xmin, double xmax);
void grid_free(Grid *grid);

/* Cell centre of cell i (valid for ghost cells too). */
double grid_cell_center(const Grid *grid, int i);

/* Total water volume over the interior cells. */
double grid_mass(const Grid *grid);

#endif
