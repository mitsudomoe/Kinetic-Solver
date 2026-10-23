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

#include "rusanov.h"

#include <math.h>

#include "physics.h"

static Flux physical_flux(double h, double u)
{
	Flux flux;

	flux.h = h * u;
	flux.hu = h * u * u + 0.5 * KS_GRAVITY * h * h;
	return flux;
}

void rusanov_interface_flux(CellState left, CellState right, Flux *seen_by_left, Flux *seen_by_right)
{
	double u_left = ks_velocity(left.h, left.hu);
	double u_right = ks_velocity(right.h, right.hu);
	double z_star = ks_max(left.z, right.z);
	/* Hydrostatic reconstruction: water heights above the highest bottom. */
	double h_left = ks_max(0.0, left.h + left.z - z_star);
	double h_right = ks_max(0.0, right.h + right.z - z_star);
	double speed = ks_max(fabs(u_left) + sqrt(KS_GRAVITY * h_left),
		fabs(u_right) + sqrt(KS_GRAVITY * h_right));
	Flux f_left = physical_flux(h_left, u_left);
	Flux f_right = physical_flux(h_right, u_right);
	Flux flux;

	flux.h = 0.5 * (f_left.h + f_right.h) - 0.5 * speed * (h_right - h_left);
	flux.hu = 0.5 * (f_left.hu + f_right.hu) - 0.5 * speed * (h_right * u_right - h_left * u_left);

	/* Interface source terms: they balance the pressure difference between
	 * the cell value and the reconstructed value. */
	seen_by_left->h = flux.h;
	seen_by_left->hu = flux.hu + 0.5 * KS_GRAVITY * (left.h * left.h - h_left * h_left);
	seen_by_right->h = flux.h;
	seen_by_right->hu = flux.hu + 0.5 * KS_GRAVITY * (right.h * right.h - h_right * h_right);
}

double rusanov_wave_speed(double h, double hu)
{
	if (h <= KS_DRY_TOL)
	{
		return 0.0;
	}
	return fabs(hu / h) + sqrt(KS_GRAVITY * h);
}
