"""Thacker's bowl, steady flows over the bump, and a wave running up a beach."""

import math
import unittest

import numpy as np

from helpers import run_final_state, run_snapshots, run_solver

G = 9.81
SCHEMES = ("kinetic", "rusanov")


def exact_columns(case, cells, t_end):
    """x, z, h_exact, hu_exact at t_end, without running the scheme for long."""
    _, data = run_final_state(case, "rusanov", "--cells", str(cells), "--t-end", repr(t_end))
    return data[:, 0], data[:, 1], data[:, 4], data[:, 5]


def bernoulli_head(q, h, z):
    return q * q / (2.0 * G * h * h) + h + z


class ExactSolutionTest(unittest.TestCase):
    """The exact solutions themselves, before comparing the schemes with them."""

    def test_subcritical_flow_has_one_bernoulli_head(self):
        _, z, h, hu = exact_columns("bump-subcritical", 1000, 1e-6)
        head = bernoulli_head(4.42, h, z)
        self.assertTrue(np.all(hu == 4.42))
        self.assertLess(np.ptp(head), 1e-12)
        self.assertAlmostEqual(head[0], bernoulli_head(4.42, 2.0, 0.0), places=12)
        self.assertTrue(np.all(h ** 3 > 4.42 ** 2 / G))  # subcritical everywhere

    def test_transcritical_flow_is_critical_at_the_crest(self):
        x, z, h, _ = exact_columns("bump-transcritical", 1000, 1e-6)
        q = 1.53
        froude = q / (h * np.sqrt(G * h))
        head = bernoulli_head(q, h, z)
        self.assertLess(np.ptp(head), 1e-6)
        self.assertTrue(np.all(froude[x < -0.1] < 1.0))
        self.assertTrue(np.all(froude[x > 0.1] > 1.0))
        self.assertAlmostEqual(froude[np.argmin(np.abs(x))], 1.0, delta=0.02)

    def test_hydraulic_jump_satisfies_rankine_hugoniot(self):
        x, z, h, _ = exact_columns("bump-hydraulic-jump", 5000, 1e-6)
        q = 0.18
        jump = np.argmax(np.diff(h))
        self.assertGreater(x[jump], 0.0)   # downstream of the crest
        self.assertLess(x[jump], 2.0)      # on the bump
        # Stationary shock: the heights on both sides are conjugate depths.
        h_left, h_right = h[jump], h[jump + 1]
        conjugate = 0.5 * h_left * (math.sqrt(1.0 + 8.0 * q * q / (G * h_left ** 3)) - 1.0)
        self.assertAlmostEqual(h_right, conjugate, delta=0.01)
        # Supercritical before the jump, subcritical with the outlet head after it.
        self.assertGreater(q / (h_left * math.sqrt(G * h_left)), 1.0)
        after = slice(jump + 1, None)
        self.assertLess(np.ptp(bernoulli_head(q, h[after], z[after])), 1e-6)
        self.assertAlmostEqual(bernoulli_head(q, h[-1], z[-1]), bernoulli_head(q, 0.33, 0.0), places=9)

    def test_numerical_jump_matches_the_exact_position(self):
        # Within two cells for the kinetic scheme, and still within two cells
        # for the (much faster) Rusanov scheme on a mesh eight times finer.
        for scheme, cells in (("kinetic", 400), ("rusanov", 3200)):
            _, data = run_final_state("bump-hydraulic-jump", scheme, "--cells", str(cells),
                                      "--snapshot-every", "100000000")
            x, h, h_exact = data[:, 0], data[:, 2], data[:, 4]
            numerical = x[np.argmax(np.diff(h))]
            exact = x[np.argmax(np.diff(h_exact))]
            with self.subTest(scheme=scheme, cells=cells):
                self.assertLessEqual(abs(numerical - exact), 2.0 * (x[1] - x[0]) + 1e-12)

    def test_thacker_solution_is_periodic_and_conserves_mass(self):
        x0, _, h0, hu0 = exact_columns("thacker", 2000, 1e-9)
        dx = x0[1] - x0[0]
        period = 2.0 * math.pi * 8.0 / math.sqrt(2.0 * G * 0.5)
        for fraction in (0.25, 0.5, 1.0):
            _, _, h, hu = exact_columns("thacker", 2000, fraction * period)
            with self.subTest(fraction=fraction):
                self.assertAlmostEqual(h.sum() * dx, h0.sum() * dx, delta=2e-4)
        self.assertLess(np.abs(h - h0).max(), 1e-12)
        self.assertLess(np.abs(hu).max(), 1e-12)


class ThackerTest(unittest.TestCase):

    BOUNDS = {"kinetic": 0.19, "rusanov": 0.026}  # L1 error of h, 200 cells, one period

    def test_error_positivity_and_mass(self):
        for scheme, bound in self.BOUNDS.items():
            with self.subTest(scheme=scheme):
                s = run_solver("thacker", scheme, "--cells", "200")
                self.assertLess(s["h_l1"], bound)
                self.assertGreaterEqual(s["min_h"], 0.0)
                self.assertLess(abs(s["mass_final"] / s["mass_initial"] - 1.0), 1e-12)

    def test_kinetic_converges(self):
        errors = [run_solver("thacker", "kinetic", "--cells", str(n))["h_l1"] for n in (100, 200, 400)]
        for k in range(2):
            self.assertGreater(math.log2(errors[k] / errors[k + 1]), 0.7)


class SteadyFlowTest(unittest.TestCase):

    # L1 error of h at the default final time, 200 cells (about 10% above the
    # values measured when the tests were written).
    BOUNDS = {
        ("bump-subcritical", "kinetic"): 0.064,
        ("bump-subcritical", "rusanov"): 0.125,
        ("bump-transcritical", "kinetic"): 0.105,
        ("bump-transcritical", "rusanov"): 0.66,
        ("bump-hydraulic-jump", "kinetic"): 0.0235,
        ("bump-hydraulic-jump", "rusanov"): 0.128,
    }

    def test_error_against_exact_steady_state(self):
        for (case, scheme), bound in self.BOUNDS.items():
            with self.subTest(case=case, scheme=scheme):
                s = run_solver(case, scheme, "--cells", "200")
                self.assertLess(s["h_l1"], bound)
                self.assertGreater(s["min_h"], 0.0)

    def test_kinetic_discharge_is_uniform_away_from_the_jump(self):
        for case, q in (("bump-subcritical", 4.42), ("bump-transcritical", 1.53)):
            with self.subTest(case=case):
                _, data = run_final_state(case, "kinetic", "--cells", "200")
                self.assertLess(np.abs(data[:, 3] - q).max() / q, 5e-3)


class BeachRunupTest(unittest.TestCase):

    @staticmethod
    def shorelines(scheme, cells):
        snapshots = run_snapshots("beach-runup", scheme, "--cells", str(cells), "--snapshot-every", "5")
        return [data[np.nonzero(data[:, 2] > 1e-3)[0][-1], 0] for data in snapshots]

    def test_conservation_and_positivity(self):
        for scheme in SCHEMES:
            with self.subTest(scheme=scheme):
                s = run_solver("beach-runup", scheme, "--cells", "500")
                self.assertGreaterEqual(s["min_h"], 0.0)
                self.assertLess(abs(s["mass_final"] / s["mass_initial"] - 1.0), 1e-12)

    def test_run_up_and_backwash(self):
        for scheme in SCHEMES:
            shore = self.shorelines(scheme, 500)
            peak = int(np.argmax(shore))
            run_up = 0.1 * shore[peak] - 0.5  # height above still water (beach slope 0.1)
            with self.subTest(scheme=scheme):
                self.assertAlmostEqual(shore[0], 5.0, delta=0.05)
                self.assertGreater(run_up, 0.30)
                self.assertLess(run_up, 0.40)
        # The kinetic scheme drains the beach during the backwash: the shoreline
        # goes back below its initial position.
        shore = self.shorelines("kinetic", 500)
        self.assertLess(min(shore[int(np.argmax(shore)):]), 4.5)


if __name__ == "__main__":
    unittest.main()
