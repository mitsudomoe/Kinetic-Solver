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

#include "output.h"

#include <errno.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "physics.h"

int output_prepare_dir(const char *dir)
{
	if (mkdir(dir, 0755) != 0 && errno != EEXIST)
	{
		perror(dir);
		return -1;
	}
	return 0;
}

int output_snapshot(const char *dir, int index, double t, const Grid *grid, const TestCase *tc, int with_exact)
{
	char path[1024];
	FILE *file;
	int i;

	snprintf(path, sizeof(path), "%s/snapshot_%05d.dat", dir, index);
	file = fopen(path, "w");
	if (!file)
	{
		perror(path);
		return -1;
	}
	fprintf(file, "# case=%s t=%.17g\n", tc->name, t);
	fprintf(file, "# x z h hu%s\n", with_exact ? " h_exact hu_exact" : "");
	for (i = 1; i <= grid->n; i++)
	{
		fprintf(file, "%.17g %.17g %.17g %.17g", grid->x[i], grid->z[i], grid->h[i], grid->hu[i]);
		if (with_exact)
		{
			double h, hu;

			if (t > 0.0)
			{
				tc->exact(tc, grid->x[i], t, grid->z[i], &h, &hu);
			}
			else
			{
				tc->initial(tc, grid->x[i], grid->z[i], &h, &hu);
			}
			fprintf(file, " %.17g %.17g", h, hu);
		}
		fputc('\n', file);
	}
	if (fclose(file) != 0)
	{
		perror(path);
		return -1;
	}

	snprintf(path, sizeof(path), "%s/times.dat", dir);
	file = fopen(path, index == 0 ? "w" : "a");
	if (!file)
	{
		perror(path);
		return -1;
	}
	fprintf(file, "%d %.17g\n", index, t);
	return fclose(file) == 0 ? 0 : -1;
}

int output_errors(const Grid *grid, const TestCase *tc, double t, int with_exact, ErrorNorms *h_err, ErrorNorms *hu_err)
{
	int i;

	h_err->l1 = h_err->linf = 0.0;
	hu_err->l1 = hu_err->linf = 0.0;
	if (!with_exact || !tc->exact)
	{
		return 0;
	}
	for (i = 1; i <= grid->n; i++)
	{
		double h, hu, eh, ehu;

		if (t > 0.0)
		{
			tc->exact(tc, grid->x[i], t, grid->z[i], &h, &hu);
		}
		else
		{
			tc->initial(tc, grid->x[i], grid->z[i], &h, &hu);
		}
		eh = fabs(grid->h[i] - h);
		ehu = fabs(grid->hu[i] - hu);
		h_err->l1 += eh;
		hu_err->l1 += ehu;
		h_err->linf = ks_max(h_err->linf, eh);
		hu_err->linf = ks_max(hu_err->linf, ehu);
	}
	h_err->l1 *= grid->dx;
	hu_err->l1 *= grid->dx;
	return 1;
}

double output_min_height(const Grid *grid)
{
	double h_min = grid->h[1];
	int i;

	for (i = 2; i <= grid->n; i++)
	{
		h_min = ks_min(h_min, grid->h[i]);
	}
	return h_min;
}
