"""End-to-end regression tests: run the solver and check its summary."""

import math
import subprocess
import tempfile
import unittest

from helpers import SOLVER, run_solver

SCHEMES = ("kinetic", "rusanov")


class LakeAtRestTest(unittest.TestCase):
    """Theorem 3.2 (iii): steady states at rest are preserved."""

    def test_lake_at_rest_is_preserved(self):
        for case in ("lake-at-rest", "lake-at-rest-step", "lake-at-rest-island"):
            for scheme in SCHEMES:
                with self.subTest(case=case, scheme=scheme):
                    s = run_solver(case, scheme)
                    self.assertLess(s["h_linf"], 1e-12)
                    self.assertLess(s["hu_linf"], 1e-12)


class DamBreakTest(unittest.TestCase):

    # Upper bounds on the L1 error of h at the default resolution (1000 cells),
    # about 10% above the values measured when the tests were written.
    L1_BOUNDS = {
        ("dam-break-dry", "kinetic"): 0.041,
        ("dam-break-dry", "rusanov"): 0.080,
        ("dam-break-wet", "kinetic"): 0.065,
        ("dam-break-wet", "rusanov"): 0.073,
    }

    def test_error_against_exact_solution(self):
        for (case, scheme), bound in self.L1_BOUNDS.items():
            with self.subTest(case=case, scheme=scheme):
                s = run_solver(case, scheme)
                self.assertLess(s["h_l1"], bound)

    def test_positivity_and_mass(self):
        for case in ("dam-break-dry", "dam-break-wet", "dam-break-bump"):
            for scheme in SCHEMES:
                with self.subTest(case=case, scheme=scheme):
                    s = run_solver(case, scheme)
                    self.assertGreaterEqual(s["min_h"], 0.0)
                    # The waves do not reach the boundaries before t_end.
                    self.assertLess(abs(s["mass_final"] / s["mass_initial"] - 1.0), 1e-12)

    def test_kinetic_converges_on_dry_dam_break(self):
        errors = [run_solver("dam-break-dry", "kinetic", "--cells", str(n))["h_l1"] for n in (200, 400, 800)]
        orders = [math.log2(errors[k] / errors[k + 1]) for k in range(2)]
        for order in orders:
            self.assertGreater(order, 0.5)


class ConservationTest(unittest.TestCase):

    def test_closed_and_periodic_domains_conserve_mass(self):
        for case in ("periodic-pulse", "sloshing"):
            for scheme in SCHEMES:
                with self.subTest(case=case, scheme=scheme):
                    s = run_solver(case, scheme)
                    self.assertLess(abs(s["mass_final"] / s["mass_initial"] - 1.0), 1e-12)
                    self.assertGreaterEqual(s["min_h"], 0.0)


class CommandLineTest(unittest.TestCase):

    def test_rejects_unknown_case_and_bad_values(self):
        with tempfile.TemporaryDirectory() as out:
            for args in (["--case", "nope"], ["--scheme", "nope"], ["--cells", "-3"], ["--cfl", "abc"]):
                with self.subTest(args=args):
                    result = subprocess.run([SOLVER, "--output-dir", out, *args], capture_output=True, text=True)
                    self.assertNotEqual(result.returncode, 0)

    def test_snapshots_are_written(self):
        with tempfile.TemporaryDirectory() as out:
            subprocess.run([SOLVER, "--case", "dam-break-dry", "--cells", "50", "--snapshot-every", "10",
                            "--output-dir", out], check=True, capture_output=True)
            with open(f"{out}/times.dat") as handle:
                times = [float(line.split()[1]) for line in handle]
            self.assertEqual(times[0], 0.0)
            self.assertAlmostEqual(times[-1], 1.5)
            self.assertEqual(times, sorted(times))
            with open(f"{out}/snapshot_{len(times) - 1:05d}.dat") as handle:
                rows = [line.split() for line in handle if not line.startswith("#")]
            self.assertEqual(len(rows), 50)
            self.assertEqual(len(rows[0]), 6)


if __name__ == "__main__":
    unittest.main()
