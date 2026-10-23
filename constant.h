/* =============================================================================

	A finite volume code to solve the Saint-Venant system in an open channel.

   =============================================================================

	Copyright (C) 2013 Mathieu Besson

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

/* =============================================================================

					All the constants needed for the program

   ============================================================================= */
#ifndef constant_h
	#define constant_h

	#define PI 3.14159265358979323846
	#define g 9.81 /* gravitational constant */
	#define bInf -10.0 /* the lower bound of the domain */
	#define bSup 10.0 /* the upper bound of the domain */
	#define nbCell 1000 /* Number of cell */
	#define freq 1 /* The frequence where you want to save data */
	#define T 1.5 /* time */
	#define epsilon 0.00001 /* the precision we want */
	#define N 1000 /* N is the number of iteration you want to do to approximate the integrals for each flux */
#endif
