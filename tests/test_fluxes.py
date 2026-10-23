"""Checks the kinetic fluxes against the definitions of Perthame & Simeoni (2001).

The reference below does not use any of the closed-form integrals of Section 4:
it integrates in xi the interface densities (3.7)-(3.8), built from the
equilibrium M(h, xi - u) = sqrt(h) chi((xi - u) / sqrt(h)) with chi given by
(2.13), and takes the moments (3.5)-(3.6). Agreement therefore validates the
closed forms, the bounds, the sign conventions and the II3q quadrature.

Convention (from the lake-at-rest proof of Theorem 3.2): the bottom jump seen
from cell i is dZ_{i+1/2} = Z_{i+1} - Z_i and dZ_{i-1/2} = Z_{i-1} - Z_i.
"""

import math
import unittest
import warnings

import numpy as np
from scipy.integrate import IntegrationWarning, quad

from helpers import FluxHarness

G = 9.81


def equilibrium(h, u, xi):
    if h <= 0.0:
        return 0.0
    w = (xi - u) / math.sqrt(h)
    s = 1.0 - w * w / (2.0 * G)
    return math.sqrt(h) * math.sqrt(2.0) / (math.pi * math.sqrt(G)) * math.sqrt(s) if s > 0.0 else 0.0


def integrate(func, lo, hi, breakpoints):
    pts = sorted({p for p in breakpoints if lo < p < hi})
    with warnings.catch_warnings():
        # The integrands have square-root edges; quad still converges well
        # below the tolerance used in the tests but may warn about it.
        warnings.simplefilter("ignore", IntegrationWarning)
        return quad(func, lo, hi, points=pts or None, limit=500, epsabs=1e-14, epsrel=1e-13)[0]


def reference_fluxes(left, right):
    """Returns (F-.h, F-.hu, F+.h, F+.hu) at the interface between two cells."""
    (h_l, hu_l, z_l), (h_r, hu_r, z_r) = left, right
    u_l = hu_l / h_l if h_l > 0 else 0.0
    u_r = hu_r / h_r if h_r > 0 else 0.0
    dz_left = z_r - z_l   # jump seen from the left cell
    dz_right = z_l - z_r  # jump seen from the right cell

    def interface_left(xi):
        # M_{i+1/2}(xi) for xi <= 0, i = left cell (definition after (3.8))
        e = xi * xi - 2.0 * G * dz_left
        return equilibrium(h_l, u_l, -xi) if e <= 0 else equilibrium(h_r, u_r, -math.sqrt(e))

    def interface_right(xi):
        # M_{(i+1)-1/2}(xi) for xi >= 0, seen from the right cell
        e = xi * xi - 2.0 * G * dz_right
        return equilibrium(h_r, u_r, -xi) if e <= 0 else equilibrium(h_l, u_l, math.sqrt(e))

    breaks = []
    for h, u in ((h_l, u_l), (h_r, u_r)):
        c = math.sqrt(2.0 * G * max(h, 0.0))
        breaks += [u - c, u + c, -(u - c), -(u + c)]
    for dz in (dz_left, dz_right):
        if dz > 0:
            breaks += [math.sqrt(2 * G * dz), -math.sqrt(2 * G * dz)]
    # Images of the support edges through the energy shift xi^2 -> xi^2 + 2 g dZ.
    breaks += [s * math.sqrt(b * b + 2 * G * abs(dz))
               for b in list(breaks) for dz in (dz_left, dz_right) for s in (-1, 1)]
    limit = 60.0

    result = []
    for power in (1, 2):
        f_minus = (integrate(lambda x: x ** power * equilibrium(h_l, u_l, x), 0.0, limit, breaks)
                   + integrate(lambda x: x ** power * interface_left(x), -limit, 0.0, breaks))
        f_plus = (integrate(lambda x: x ** power * interface_right(x), 0.0, limit, breaks)
                  + integrate(lambda x: x ** power * equilibrium(h_r, u_r, x), -limit, 0.0, breaks))
        result.append((f_minus, f_plus))
    return result[0][0], result[1][0], result[0][1], result[1][1]


def random_states(seed, count):
    rng = np.random.default_rng(seed)
    states = [((1.0, 0.0, 0.0), (1.2, 0.0, -0.2)),   # lake at rest, step down
              ((0.8, 0.0, 0.2), (1.0, 0.0, 0.0))]    # lake at rest, step up
    for k in range(count):
        h = rng.uniform(0.05, 2.0, 2)
        u = rng.uniform(-3.0, 3.0, 2) * (3.0 if k % 5 == 0 else 1.0)  # some supercritical states
        z = rng.uniform(-0.5, 0.5, 2) * (k % 4 != 0)                 # some flat interfaces
        if k % 7 == 0:
            h[k % 2] = 0.0                                             # dry neighbour
        states.append(((h[0], h[0] * u[0], z[0]), (h[1], h[1] * u[1], z[1])))
    return states


class KineticFluxTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.harness = FluxHarness()

    @classmethod
    def tearDownClass(cls):
        cls.harness.close()

    def test_matches_kinetic_definition(self):
        for left, right in random_states(seed=2026, count=40):
            code = self.harness.flux("kinetic", left, right)
            ref = reference_fluxes(left, right)
            scale = 1.0 + max(abs(v) for v in ref)
            for name, c, r in zip(("F-.h", "F-.hu", "F+.h", "F+.hu"), code, ref):
                with self.subTest(left=left, right=right, component=name):
                    self.assertLess(abs(c - r) / scale, 1e-10)

    def test_mass_flux_is_conservative(self):
        # Remark 3.1: (IF_h)^- = (IF_h)^+ at every interface, even across a step.
        for left, right in random_states(seed=7, count=200):
            f_minus_h, _, f_plus_h, _ = self.harness.flux("kinetic", left, right)
            self.assertAlmostEqual(f_minus_h, f_plus_h, delta=1e-13 * (1 + abs(f_minus_h)))

    def test_flat_bottom_is_conservative(self):
        for left, right in random_states(seed=11, count=100):
            left, right = (left[0], left[1], 0.3), (right[0], right[1], 0.3)
            f = self.harness.flux("kinetic", left, right)
            self.assertAlmostEqual(f[0], f[2], delta=1e-13 * (1 + abs(f[0])))
            self.assertAlmostEqual(f[1], f[3], delta=1e-13 * (1 + abs(f[1])))

    def test_lake_at_rest_balance(self):
        # Hydrostatic balance: across a step, the momentum fluxes differ by
        # exactly the pressure jump g/2 (h_l^2 - h_r^2) when H is constant.
        for scheme in ("kinetic", "rusanov"):
            for h_left, z_left, z_right in ((1.0, 0.0, 0.3), (0.7, 0.3, 0.0), (2.0, 0.1, 0.25)):
                h_right = h_left + z_left - z_right
                f = self.harness.flux(scheme, (h_left, 0.0, z_left), (h_right, 0.0, z_right))
                with self.subTest(scheme=scheme, h_left=h_left):
                    self.assertAlmostEqual(f[0], 0.0, delta=1e-14)
                    self.assertAlmostEqual(f[2], 0.0, delta=1e-14)
                    self.assertAlmostEqual(f[1], 0.5 * G * h_left ** 2, delta=1e-13)
                    self.assertAlmostEqual(f[3], 0.5 * G * h_right ** 2, delta=1e-13)


class RusanovFluxTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.harness = FluxHarness()

    @classmethod
    def tearDownClass(cls):
        cls.harness.close()

    def test_dry_states_give_finite_fluxes(self):
        for left, right in (((1.0, 0.5, 0.0), (0.0, 0.0, 0.0)), ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0)),
                            ((0.0, 0.0, 0.5), (1.0, -0.3, 0.0))):
            f = self.harness.flux("rusanov", left, right)
            self.assertTrue(all(math.isfinite(v) for v in f))

    def test_flat_bottom_matches_rusanov(self):
        left, right = (1.5, 0.4, 0.0), (0.8, -0.2, 0.0)
        f = self.harness.flux("rusanov", left, right)
        u_l, u_r = left[1] / left[0], right[1] / right[0]
        c = max(abs(u_l) + math.sqrt(G * left[0]), abs(u_r) + math.sqrt(G * right[0]))
        flux_h = 0.5 * (left[1] + right[1]) - 0.5 * c * (right[0] - left[0])
        flux_q = (0.5 * (left[1] * u_l + 0.5 * G * left[0] ** 2 + right[1] * u_r + 0.5 * G * right[0] ** 2)
                  - 0.5 * c * (right[1] - left[1]))
        for got, expected in zip(f, (flux_h, flux_q, flux_h, flux_q)):
            self.assertAlmostEqual(got, expected, places=13)


def jump_height(left, right):
    """Effective interface height h^jump of Ersoy (2026), Proposition 4.1.

    M. Ersoy, "Perthame-Simeoni kinetic traces and nonconservative interface
    laws for DG approximations of the Saint-Venant system", hal-05618693.
    With dZ = z_right - z_left and alpha = sqrt(2 g |dZ|), the one-sided
    momentum fluxes satisfy F2^- - F2^+ = g dZ h^jump, where h^jump >= 0 is
    the sum of the reflected and transmitted particle contributions below.
    """
    (h_l, hu_l, z_l), (h_r, hu_r, z_r) = left, right
    u_l = hu_l / h_l if h_l > 0 else 0.0
    u_r = hu_r / h_r if h_r > 0 else 0.0
    dz = z_r - z_l
    alpha = math.sqrt(2.0 * G * abs(dz))

    def m_left(xi):
        return equilibrium(h_l, u_l, xi)

    def m_right_reversed(xi):
        return equilibrium(h_r, u_r, -xi)

    # Support edges of both Maxwellians, mirrored, plus alpha.
    breaks = [alpha]
    for h, u in ((h_l, u_l), (h_r, u_r)):
        c = math.sqrt(2.0 * G * h)
        breaks += [u - c, u + c, -(u - c), -(u + c)]

    def reflected(m):
        return 2.0 / (G * abs(dz)) * integrate(lambda x: x * x * m(x), 0.0, alpha, breaks)

    def transmitted_up(m):    # climbs the step: xi > alpha
        return 2.0 * integrate(lambda x: x / (x + math.sqrt(x * x - alpha * alpha)) * m(x),
                               alpha, 60.0, breaks)

    def transmitted_down(m):  # falls down the step
        return 2.0 * integrate(lambda x: x / (math.sqrt(x * x + alpha * alpha) + x) * m(x),
                               0.0, 60.0, breaks)

    if dz > 0:
        return reflected(m_left) + transmitted_up(m_left) + transmitted_down(m_right_reversed)
    return transmitted_down(m_left) + reflected(m_right_reversed) + transmitted_up(m_right_reversed)


class MomentumJumpTest(unittest.TestCase):
    """Ersoy (2026), Proposition 4.1: F2^- - F2^+ = g dZ h^jump with h^jump >= 0."""

    @classmethod
    def setUpClass(cls):
        cls.harness = FluxHarness()

    @classmethod
    def tearDownClass(cls):
        cls.harness.close()

    def test_momentum_defect_matches_jump_height(self):
        cases = [s for s in random_states(seed=2027, count=60) if s[0][2] != s[1][2]]
        self.assertGreater(len(cases), 30)
        for left, right in cases:
            _, f_minus_q, _, f_plus_q = self.harness.flux("kinetic", left, right)
            dz = right[2] - left[2]
            h_jump = jump_height(left, right)
            with self.subTest(left=left, right=right):
                self.assertGreaterEqual(h_jump, 0.0)
                self.assertAlmostEqual(f_minus_q - f_plus_q, G * dz * h_jump,
                                       delta=1e-10 * (1.0 + abs(f_minus_q) + abs(f_plus_q)))

    def test_lake_at_rest_selects_mean_height(self):
        # At rest (u = 0, h + z constant) the defect is the hydrostatic pressure
        # jump g/2 (h_l^2 - h_r^2) = g dZ (h_l + h_r)/2: h^jump is the mean height.
        for h_left, z_left, z_right in ((1.0, 0.0, 0.3), (0.7, 0.3, 0.0), (2.0, 0.1, 0.25)):
            h_right = h_left + z_left - z_right
            h_jump = jump_height((h_left, 0.0, z_left), (h_right, 0.0, z_right))
            self.assertAlmostEqual(h_jump, 0.5 * (h_left + h_right), places=10)


class StokerTest(unittest.TestCase):

    def test_middle_state_solves_jump_conditions(self):
        harness = FluxHarness()
        try:
            for h_left, h_right in ((2.0, 1.0), (1.0, 0.1), (5.0, 4.9)):
                h_mid = harness.stoker(h_left, h_right)
                u_rarefaction = 2.0 * (math.sqrt(G * h_left) - math.sqrt(G * h_mid))
                u_shock = (h_mid - h_right) * math.sqrt(G * (h_mid + h_right) / (2 * h_mid * h_right))
                self.assertGreater(h_mid, h_right)
                self.assertLess(h_mid, h_left)
                self.assertAlmostEqual(u_rarefaction, u_shock, places=12)
        finally:
            harness.close()


if __name__ == "__main__":
    unittest.main()
