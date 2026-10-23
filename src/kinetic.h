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

/* Kinetic scheme with reflections of B. Perthame and C. Simeoni,
 * "A kinetic scheme for the Saint-Venant system with a source term",
 * Calcolo 38 (2001), https://hal.science/hal-00922664
 *
 * Notation follows Section 4 of the paper. For the interface x_{i+1/2}
 * seen from cell i (flux IF^-_{i+1/2}):
 *     dZ  = Z_{i+1} - Z_i              (bottom jump seen from cell i)
 *     Fr  = u / sqrt(2 g h)            (Froude-like number of each cell)
 *     K-  = dZ / h_i,  K+ = dZ / h_{i+1}
 *     IF^h = 2 sqrt(2g)/pi [h_i^{3/2} (IIh1 + IIh2) + h_{i+1}^{3/2} IIh3]
 *     IF^q = 4 g / pi      [h_i^2     (IIq1 + IIq2) + h_{i+1}^2     IIq3]
 * The flux IF^+_{i-1/2} seen from cell i at its left interface is obtained
 * from the same formula by the mirror symmetry x -> -x (see kinetic.c). */
#ifndef KS_KINETIC_H
#define KS_KINETIC_H

#include "scheme.h"

void kinetic_interface_flux(CellState left, CellState right, Flux *seen_by_left, Flux *seen_by_right);

/* |u| + sqrt(2 g h): the radius of the support of the kinetic density, used
 * by the CFL condition (3.9) of the paper. */
double kinetic_wave_speed(double h, double hu);

/* The integrals of Section 4.1, exposed for unit tests. */
double kinetic_IIh1(double fr);
double kinetic_IIh2(double fr, double k_minus_pos);
double kinetic_IIh3(double fr, double k_plus_neg);
double kinetic_IIq1(double fr);
double kinetic_IIq2(double fr, double k_minus_pos);
double kinetic_IIq3(double fr, double k_plus_neg, double k_plus);

#endif
