"""Shared helpers for the test suite: paths, running the solver, the flux harness."""

import json
import os
import subprocess
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOLVER = os.path.join(ROOT, "build", "kinetic-solver")
HARNESS = os.path.join(ROOT, "build", "flux-harness")


def run_solver(case, scheme, *extra):
    """Runs one case and returns its summary.json as a dict."""
    with tempfile.TemporaryDirectory() as out:
        subprocess.run(
            [SOLVER, "--case", case, "--scheme", scheme, "--no-snapshots", "--output-dir", out, *extra],
            check=True, capture_output=True, text=True)
        with open(os.path.join(out, "summary.json")) as handle:
            return json.load(handle)


class FluxHarness:
    """Keeps one harness process open and evaluates interface fluxes through it."""

    def __init__(self):
        self.proc = subprocess.Popen([HARNESS], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)

    def flux(self, scheme, left, right):
        """left/right are (h, hu, z); returns (F-.h, F-.hu, F+.h, F+.hu)."""
        values = " ".join(repr(float(v)) for v in (*left, *right))
        self.proc.stdin.write(f"{scheme} {values}\n")
        self.proc.stdin.flush()
        return tuple(float(v) for v in self.proc.stdout.readline().split())

    def stoker(self, h_left, h_right):
        self.proc.stdin.write(f"stoker {h_left!r} {h_right!r}\n")
        self.proc.stdin.flush()
        return float(self.proc.stdout.readline())

    def close(self):
        self.proc.stdin.close()
        self.proc.wait()
        self.proc.stdout.close()
