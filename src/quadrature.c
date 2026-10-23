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

#include "quadrature.h"

#include <math.h>

#include "physics.h"

/* Step size and truncation of the tanh-sinh rule. With h = 1/16 and
 * |t| <= 3.5 the rule uses 113 nodes and reaches about 1e-15 relative
 * accuracy on integrands with square-root endpoint singularities. */
#define TS_STEP (1.0 / 16.0)
#define TS_HALF_NODES 56

typedef struct
{
	double weight;   /* weight on [-1, 1] */
	double from_lo;  /* 1 + node, i.e. distance to -1 */
	double to_hi;    /* 1 - node, i.e. distance to +1 */
} Node;

static Node nodes[2 * TS_HALF_NODES + 1];
static int nodes_ready = 0;

static void build_nodes(void)
{
	int k;

	for (k = -TS_HALF_NODES; k <= TS_HALF_NODES; k++)
	{
		double t = k * TS_STEP;
		double s = 0.5 * KS_PI * sinh(t);
		double c = cosh(s);
		Node *node = &nodes[k + TS_HALF_NODES];

		node->weight = TS_STEP * 0.5 * KS_PI * cosh(t) / (c * c);
		/* 1 + tanh(s) = 2 / (1 + exp(-2s)) and 1 - tanh(s) = 2 / (1 + exp(2s)),
		 * both free of cancellation. */
		node->from_lo = 2.0 / (1.0 + exp(-2.0 * s));
		node->to_hi = 2.0 / (1.0 + exp(2.0 * s));
	}
	nodes_ready = 1;
}

double quad_tanh_sinh(QuadIntegrand f, void *ctx, double a, double b)
{
	double half = 0.5 * (b - a);
	double sum = 0.0;
	int k;

	if (!(b > a))
	{
		return 0.0;
	}
	if (!nodes_ready)
	{
		build_nodes();
	}
	for (k = 0; k < 2 * TS_HALF_NODES + 1; k++)
	{
		double from_a = half * nodes[k].from_lo;
		double to_b = half * nodes[k].to_hi;

		/* Nodes that collapse onto an endpoint carry negligible weight and
		 * could hit the singularity exactly: skip them. */
		if (from_a <= 0.0 || to_b <= 0.0)
		{
			continue;
		}
		sum += nodes[k].weight * f(a + from_a, from_a, to_b, ctx);
	}
	return half * sum;
}

int quad_tanh_sinh_nodes(void)
{
	return 2 * TS_HALF_NODES + 1;
}
