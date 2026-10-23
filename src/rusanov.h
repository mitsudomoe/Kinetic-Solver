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

/* Rusanov (local Lax-Friedrichs) flux combined with the hydrostatic
 * reconstruction of Audusse, Bouchut, Bristeau, Klein and Perthame,
 * "A fast and stable well-balanced scheme with hydrostatic reconstruction
 * for shallow water flows", SIAM J. Sci. Comput. 25 (2004).
 *
 * This is the reference scheme against which the kinetic scheme is compared:
 * it is also positive and preserves the lake at rest, including with dry
 * areas. */
#ifndef KS_RUSANOV_H
#define KS_RUSANOV_H

#include "scheme.h"

void rusanov_interface_flux(CellState left, CellState right, Flux *seen_by_left, Flux *seen_by_right);

/* |u| + sqrt(g h). */
double rusanov_wave_speed(double h, double hu);

#endif
