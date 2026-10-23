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

#include "scheme.h"

#include <stdio.h>
#include <string.h>

#include "kinetic.h"
#include "physics.h"
#include "rusanov.h"

static const Scheme schemes[] = {
	{"kinetic", "kinetic scheme with reflections (Perthame & Simeoni, 2001)",
	 kinetic_interface_flux, kinetic_wave_speed, 0.95},
	{"rusanov", "Rusanov flux with hydrostatic reconstruction (Audusse et al., 2004)",
	 rusanov_interface_flux, rusanov_wave_speed, 0.9},
};

#define NB_SCHEMES ((int)(sizeof(schemes) / sizeof(schemes[0])))

const Scheme *scheme_find(const char *name)
{
	int i;

	for (i = 0; i < NB_SCHEMES; i++)
	{
		if (strcmp(schemes[i].name, name) == 0)
		{
			return &schemes[i];
		}
	}
	return NULL;
}

void scheme_list(void)
{
	int i;

	for (i = 0; i < NB_SCHEMES; i++)
	{
		printf("  %-10s %s\n", schemes[i].name, schemes[i].description);
	}
}

double scheme_time_step(const Scheme *scheme, const Grid *grid, double cfl)
{
	double speed = 0.0;
	int i;

	/* Ghost cells take part in the interface fluxes, so they count too. */
	for (i = 0; i <= grid->n + 1; i++)
	{
		speed = ks_max(speed, scheme->wave_speed(grid->h[i], grid->hu[i]));
	}
	if (speed <= KS_DRY_TOL)
	{
		return cfl * grid->dx;
	}
	return cfl * grid->dx / speed;
}

void scheme_flux_balance(const Scheme *scheme, const Grid *grid, double *dh, double *dhu, Flux *workspace)
{
	Flux *seen_by_left = workspace;               /* F^-_{j+1/2}, j = 0..n */
	Flux *seen_by_right = workspace + grid->n + 1; /* F^+_{j+1/2}, j = 0..n */
	int i;

	/* One pass over the n + 1 interfaces: each flux is computed once. */
	for (i = 0; i <= grid->n; i++)
	{
		CellState left = {grid->h[i], grid->hu[i], grid->z[i]};
		CellState right = {grid->h[i + 1], grid->hu[i + 1], grid->z[i + 1]};

		scheme->interface_flux(left, right, &seen_by_left[i], &seen_by_right[i]);
	}

	for (i = 1; i <= grid->n; i++)
	{
		dh[i] = seen_by_left[i].h - seen_by_right[i - 1].h;
		dhu[i] = seen_by_left[i].hu - seen_by_right[i - 1].hu;
	}
}
