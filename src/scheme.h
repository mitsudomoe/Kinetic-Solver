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

/* Common interface of the numerical schemes.
 *
 * Both schemes in this code are "interfacial upwind source" schemes: at each
 * interface x_{j+1/2} they return two fluxes, one seen by the left cell j
 * (F^-_{j+1/2}) and one seen by the right cell j+1 (F^+_{j+1/2}). They differ
 * when the bottom jumps across the interface, which is how the topography
 * source term enters. The update is then
 *
 *     U_i^{n+1} = U_i^n - dt/dx (F^-_{i+1/2} - F^+_{i-1/2}). */
#ifndef KS_SCHEME_H
#define KS_SCHEME_H

#include "grid.h"

typedef struct
{
	double h;   /* mass component */
	double hu;  /* momentum component */
} Flux;

typedef struct
{
	double h;
	double hu;
	double z;
} CellState;

typedef void (*InterfaceFluxFn)(CellState left, CellState right, Flux *seen_by_left, Flux *seen_by_right);
typedef double (*WaveSpeedFn)(double h, double hu);

typedef struct
{
	const char *name;
	const char *description;
	InterfaceFluxFn interface_flux;
	WaveSpeedFn wave_speed;     /* bound on |xi| used by the CFL condition */
	double default_cfl;
} Scheme;

/* Looks a scheme up by name ("kinetic" or "rusanov"); NULL if unknown. */
const Scheme *scheme_find(const char *name);
void scheme_list(void);

/* Time step dx * cfl / max wave speed (dx * cfl if the domain is dry). */
double scheme_time_step(const Scheme *scheme, const Grid *grid, double cfl);

/* Advances the interior cells by dt. The ghost cells must be up to date;
 * workspace must hold 2 * (n + 1) Flux values. */
void scheme_step(const Scheme *scheme, Grid *grid, double dt, Flux *workspace);

#endif
