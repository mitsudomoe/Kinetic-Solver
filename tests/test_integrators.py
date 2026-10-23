"""Time integrators and metamorphic properties of whole simulations."""

import math
import subprocess
import tempfile
import unittest

import numpy as np

from helpers import SOLVER, run_final_state, run_solver

INTEGRATORS = ("euler", "ssp-rk2", "ssp-rk3")
SCHEMES = ("kinetic", "rusanov")

# Symmetries hold exactly in IEEE arithmetic on x86-64. Compilers may fuse
# multiply-adds on other targets (e.g. arm64), so allow round-off.
SYMMETRY_TOL = 1e-12


class TemporalOrderTest(unittest.TestCase):
    """Order in time on a fixed grid, against a reference with a tiny time step.

    The scheme is first order in space, so the global error of every
    integrator is O(dx) on a converging grid; the time order can only be
    measured with dx fixed, which is what this test does.
    """

    CASE = ("periodic-pulse", "--cells", "100", "--t-end", "0.5")
    STEPS = (20, 40, 80)

    def final_height(self, integrator, steps):
        summary, data = run_final_state(self.CASE[0], "kinetic", *self.CASE[1:], "--integrator", integrator,
                                        "--dt", repr(0.5 / steps))
        self.assertEqual(summary["cfl_violations"], 0)
        return data[:, 2]

    def test_each_integrator_has_its_formal_order(self):
        reference = self.final_height("ssp-rk3", 2560)
        for integrator, order in zip(INTEGRATORS, (1, 2, 3)):
            errors = [np.abs(self.final_height(integrator, k) - reference).max() for k in self.STEPS]
            observed = [math.log2(errors[j] / errors[j + 1]) for j in range(len(errors) - 1)]
            with self.subTest(integrator=integrator, errors=errors):
                for value in observed:
                    self.assertGreater(value, order - 0.2)
                    self.assertLess(value, order + 0.3)


class SspPropertiesTest(unittest.TestCase):
    """SSP methods keep what the scheme has under forward Euler."""

    def test_lake_at_rest_is_preserved(self):
        for case in ("lake-at-rest", "lake-at-rest-step", "lake-at-rest-island"):
            for scheme in SCHEMES:
                for integrator in INTEGRATORS[1:]:
                    with self.subTest(case=case, scheme=scheme, integrator=integrator):
                        s = run_solver(case, scheme, "--integrator", integrator, "--cells", "400")
                        self.assertLess(s["h_linf"], 1e-12)
                        self.assertLess(s["hu_linf"], 1e-12)

    def test_positivity_and_mass(self):
        for case in ("dam-break-dry", "sloshing", "periodic-pulse"):
            for scheme in SCHEMES:
                for integrator in INTEGRATORS[1:]:
                    with self.subTest(case=case, scheme=scheme, integrator=integrator):
                        s = run_solver(case, scheme, "--integrator", integrator, "--cells", "400")
                        self.assertGreaterEqual(s["min_h"], 0.0)
                        self.assertLess(abs(s["mass_final"] / s["mass_initial"] - 1.0), 1e-12)

    def test_converges_to_the_semi_discrete_solution(self):
        # With first-order space, forward Euler at a large CFL number benefits
        # from a cancellation between time and space errors. As the CFL number
        # decreases, Euler must converge to the semi-discrete solution that the
        # Runge-Kutta methods already reach.
        def l1(integrator, cfl):
            return run_solver("dam-break-wet", "kinetic", "--cells", "200", "--cfl", cfl,
                              "--integrator", integrator)["h_l1"]

        semi_discrete = l1("ssp-rk3", "0.95")
        self.assertAlmostEqual(l1("ssp-rk2", "0.95"), semi_discrete, delta=5e-3 * semi_discrete)
        self.assertAlmostEqual(l1("euler", "0.02"), semi_discrete, delta=3e-2 * semi_discrete)


class MetamorphicTest(unittest.TestCase):
    """Transforming the input must transform the output in the same way."""

    def test_mirror_symmetry(self):
        # x -> -x, u -> -u: the mirrored run is the mirror of the original run.
        # An odd number of cells puts a cell centre on the axis.
        for case in ("dam-break-bump", "sloshing", "lake-at-rest-island", "dam-break-dry"):
            for scheme in SCHEMES:
                for integrator in ("euler", "ssp-rk3"):
                    args = ("--integrator", integrator, "--cells", "201")
                    _, direct = run_final_state(case, scheme, *args)
                    _, mirrored = run_final_state(case, scheme, *args, "--mirror")
                    mirrored = mirrored[::-1]
                    with self.subTest(case=case, scheme=scheme, integrator=integrator):
                        np.testing.assert_allclose(mirrored[:, 1], direct[:, 1], rtol=0, atol=SYMMETRY_TOL)
                        np.testing.assert_allclose(mirrored[:, 2], direct[:, 2], rtol=0, atol=SYMMETRY_TOL)
                        np.testing.assert_allclose(-mirrored[:, 3], direct[:, 3], rtol=0, atol=SYMMETRY_TOL)

    def test_translation_invariance_with_periodic_boundaries(self):
        for scheme in SCHEMES:
            for integrator in ("euler", "ssp-rk2"):
                for shift in (1, 37, -64):
                    args = ("--integrator", integrator, "--cells", "200", "--t-end", "1.0")
                    _, direct = run_final_state("periodic-pulse", scheme, *args)
                    _, shifted = run_final_state("periodic-pulse", scheme, *args, "--shift-cells", str(shift))
                    with self.subTest(scheme=scheme, integrator=integrator, shift=shift):
                        np.testing.assert_allclose(shifted[:, 2], np.roll(direct[:, 2], shift),
                                                   rtol=0, atol=SYMMETRY_TOL)
                        np.testing.assert_allclose(shifted[:, 3], np.roll(direct[:, 3], shift),
                                                   rtol=0, atol=SYMMETRY_TOL)

    def test_transformed_runs_do_not_report_exact_errors(self):
        summary, data = run_final_state("dam-break-dry", "kinetic", "--cells", "51", "--mirror")
        self.assertIsNone(summary["h_l1"])
        self.assertEqual(data.shape[1], 4)


class OptionsTest(unittest.TestCase):

    def test_rejects_invalid_combinations(self):
        with tempfile.TemporaryDirectory() as out:
            for args in (["--integrator", "rk4"], ["--case", "dam-break-dry", "--shift-cells", "3"],
                         ["--dt", "0"], ["--dt", "-1e-3"]):
                with self.subTest(args=args):
                    result = subprocess.run([SOLVER, "--output-dir", out, *args], capture_output=True, text=True)
                    self.assertNotEqual(result.returncode, 0)

    def test_reports_cfl_violations_of_fixed_time_step(self):
        # Slightly above the CFL time step: stays finite, but is reported.
        s = run_solver("dam-break-wet", "kinetic", "--cells", "100", "--dt", "0.033", "--t-end", "0.1")
        self.assertGreater(s["cfl_violations"], 0)

    def test_divergence_is_an_error(self):
        with tempfile.TemporaryDirectory() as out:
            result = subprocess.run([SOLVER, "--case", "dam-break-wet", "--cells", "100", "--dt", "0.05",
                                     "--no-snapshots", "--output-dir", out], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("non-finite", result.stderr)


if __name__ == "__main__":
    unittest.main()
