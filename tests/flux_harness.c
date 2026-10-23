/* Test harness: evaluates one interface flux per input line.
 *
 * Input lines:  <scheme> hL huL zL hR huR zR
 * Output lines: F-.h F-.hu F+.h F+.hu   (fluxes seen by the left and right cell)
 *
 * Also accepts "stoker hL hR" and prints Stoker's middle height. */

#include <stdio.h>
#include <string.h>

#include "cases.h"
#include "scheme.h"

int main(void)
{
	char name[32];

	while (scanf("%31s", name) == 1)
	{
		if (strcmp(name, "stoker") == 0)
		{
			double hl, hr;

			if (scanf("%lf %lf", &hl, &hr) != 2)
				return 1;
			printf("%.17g\n", case_stoker_middle_height(hl, hr));
		}
		else
		{
			const Scheme *scheme = scheme_find(name);
			CellState left, right;
			Flux seen_by_left, seen_by_right;

			if (!scheme || scanf("%lf %lf %lf %lf %lf %lf", &left.h, &left.hu, &left.z,
					&right.h, &right.hu, &right.z) != 6)
				return 1;
			scheme->interface_flux(left, right, &seen_by_left, &seen_by_right);
			printf("%.17g %.17g %.17g %.17g\n", seen_by_left.h, seen_by_left.hu,
				seen_by_right.h, seen_by_right.hu);
		}
		fflush(stdout);
	}
	return 0;
}
