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

#include "grid.h"

#include <stdlib.h>

int grid_init(Grid *grid, int n, double xmin, double xmax)
{
	int i;
	size_t size = (size_t)(n + 2);

	grid->n = n;
	grid->xmin = xmin;
	grid->xmax = xmax;
	grid->dx = (xmax - xmin) / n;
	grid->x = calloc(size, sizeof(double));
	grid->z = calloc(size, sizeof(double));
	grid->h = calloc(size, sizeof(double));
	grid->hu = calloc(size, sizeof(double));
	if (!grid->x || !grid->z || !grid->h || !grid->hu)
	{
		grid_free(grid);
		return -1;
	}
	for (i = 0; i <= n + 1; i++)
	{
		grid->x[i] = grid_cell_center(grid, i);
	}
	return 0;
}

void grid_free(Grid *grid)
{
	free(grid->x);
	free(grid->z);
	free(grid->h);
	free(grid->hu);
	grid->x = grid->z = grid->h = grid->hu = NULL;
}

double grid_cell_center(const Grid *grid, int i)
{
	return grid->xmin + (i - 0.5) * grid->dx;
}

double grid_mass(const Grid *grid)
{
	double mass = 0.0;
	int i;

	for (i = 1; i <= grid->n; i++)
	{
		mass += grid->h[i];
	}
	return mass * grid->dx;
}
