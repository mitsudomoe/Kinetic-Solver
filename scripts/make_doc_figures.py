#!/usr/bin/env python3
"""Regenerates the figures and numbers of the documentation's validation page.

    make            # build the solver first
    python3 scripts/make_doc_figures.py

Writes PNG files to docs/assets/ and prints the tables used in docs/validation.md.
"""

import glob
import json
import math
import os
import subprocess
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOLVER = os.path.join(ROOT, "build", "kinetic-solver")
ASSETS = os.path.join(ROOT, "docs", "assets")

COLORS = {"kinetic": "#1f77b4", "rusanov": "#d62728", "exact": "black"}


def run(case, scheme="kinetic", *extra):
    """Runs the solver; returns (summary, final snapshot)."""
    with tempfile.TemporaryDirectory() as out:
        subprocess.run([SOLVER, "--case", case, "--scheme", scheme, "--output-dir", out, *extra],
                       check=True, capture_output=True)
        with open(os.path.join(out, "summary.json")) as handle:
            summary = json.load(handle)
        snapshots = sorted(glob.glob(os.path.join(out, "snapshot_*.dat")))
        data = np.loadtxt(snapshots[-1], ndmin=2) if snapshots else None
    return summary, data


def save(fig, name):
    path = os.path.join(ASSETS, name)
    fig.savefig(path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {os.path.relpath(path, ROOT)}")


def dam_break_figure(case, name, title):
    fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(8, 5.2))
    exact = None
    for scheme in ("kinetic", "rusanov"):
        _, data = run(case, scheme, "--cells", "200")
        x, z, h, hu = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
        top.plot(x, h + z, color=COLORS[scheme], label=scheme)
        bottom.plot(x, hu, color=COLORS[scheme], label=scheme)
        if data.shape[1] >= 6:
            exact = data
    if exact is not None:
        top.plot(exact[:, 0], exact[:, 4] + exact[:, 1], "--", color="black", linewidth=1, label="exact")
        bottom.plot(exact[:, 0], exact[:, 5], "--", color="black", linewidth=1, label="exact")
    top.set_title(f"{title} — 200 cells, t = 1.5 s")
    top.set_ylabel("h + z (m)")
    bottom.set_ylabel("hu (m²/s)")
    bottom.set_xlabel("x (m)")
    for ax in (top, bottom):
        ax.grid(alpha=0.3)
        ax.legend(fontsize="small")
    save(fig, name)


def topography_figure():
    fig, axes = plt.subplots(1, 3, figsize=(11, 3))
    for ax, case, title in zip(axes, ("lake-at-rest", "lake-at-rest-step", "lake-at-rest-island"),
                               ("bump, H = 2", "step, H = 1", "emerged bump, H = 0.1")):
        _, data = run(case, "kinetic", "--cells", "400")
        x, z, h = data[:, 0], data[:, 1], data[:, 2]
        ax.fill_between(x, z.min() - 0.05, z, color="0.75")
        ax.plot(x, h + z, color=COLORS["kinetic"])
        wet = h > 0
        ax.set_title(f"{title} ({np.count_nonzero(~wet)} dry cells)" if (~wet).any() else title, fontsize="medium")
        ax.set_xlabel("x (m)")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("elevation (m)")
    save(fig, "lake-at-rest.png")


def space_convergence():
    rows = []
    cells = (100, 200, 400, 800, 1600)
    fig, ax = plt.subplots(figsize=(5.5, 4))
    for case in ("dam-break-dry", "dam-break-wet"):
        for scheme in ("kinetic", "rusanov"):
            errors = [run(case, scheme, "--cells", str(n), "--no-snapshots")[0]["h_l1"] for n in cells]
            orders = [math.log2(errors[k] / errors[k + 1]) for k in range(len(errors) - 1)]
            rows.append((case, scheme, errors, orders))
            ax.loglog(cells, errors, "o-" if case.endswith("dry") else "s--", color=COLORS[scheme],
                      label=f"{case}, {scheme}")
    ref = np.array(cells, dtype=float)
    ax.loglog(ref, 0.5 * rows[0][2][0] * (ref / ref[0]) ** -1, ":", color="0.5", label="slope −1")
    ax.set_xlabel("cells")
    ax.set_ylabel("L1 error of h")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize="small")
    save(fig, "space-convergence.png")
    print("\nSpace convergence, L1 error of h at t = 1.5 s:")
    print("| Case | Scheme | " + " | ".join(str(n) for n in cells) + " | observed order |")
    for case, scheme, errors, orders in rows:
        print(f"| {case} | {scheme} | " + " | ".join(f"{e:.2e}" for e in errors)
              + f" | {orders[-2]:.2f}, {orders[-1]:.2f} |")


def time_convergence():
    base = ("--cells", "100", "--t-end", "0.5")
    _, ref = run("periodic-pulse", "kinetic", *base, "--integrator", "ssp-rk3", "--dt", repr(0.5 / 2560))
    steps = (20, 40, 80, 160)
    fig, ax = plt.subplots(figsize=(5.5, 4))
    print("\nTime convergence (periodic pulse, 100 cells, max error of h against dt -> 0):")
    print("| Integrator | " + " | ".join(f"dt = 0.5/{k}" for k in steps) + " | observed order |")
    for integrator, marker in (("euler", "o"), ("ssp-rk2", "s"), ("ssp-rk3", "^")):
        errors = []
        for k in steps:
            _, data = run("periodic-pulse", "kinetic", *base, "--integrator", integrator, "--dt", repr(0.5 / k))
            errors.append(np.abs(data[:, 2] - ref[:, 2]).max())
        orders = [math.log2(errors[j] / errors[j + 1]) for j in range(len(errors) - 1)]
        ax.loglog([0.5 / k for k in steps], errors, marker + "-", label=integrator)
        print(f"| {integrator} | " + " | ".join(f"{e:.2e}" for e in errors) + f" | {orders[-1]:.2f} |")
    ax.set_xlabel("dt (s)")
    ax.set_ylabel("max error of h")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize="small")
    save(fig, "time-convergence.png")


def euler_versus_rk():
    cfls = (0.95, 0.5, 0.25, 0.1, 0.05, 0.02)
    fig, ax = plt.subplots(figsize=(5.5, 4))
    print("\nWet dam break, 400 cells, L1 error of h as the CFL number decreases:")
    print("| Integrator | " + " | ".join(f"CFL {c}" for c in cfls) + " |")
    for integrator, marker in (("euler", "o"), ("ssp-rk2", "s"), ("ssp-rk3", "^")):
        errors = [run("dam-break-wet", "kinetic", "--cells", "400", "--cfl", str(c), "--integrator", integrator,
                      "--no-snapshots")[0]["h_l1"] for c in cfls]
        ax.semilogx(cfls, errors, marker + "-", label=integrator)
        print(f"| {integrator} | " + " | ".join(f"{e:.4f}" for e in errors) + " |")
    ax.set_xlabel("CFL number")
    ax.set_ylabel("L1 error of h")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize="small")
    save(fig, "euler-vs-rk.png")


def lake_at_rest_table():
    print("\nLake at rest, 1000 cells, t = 1.5 s (max errors):")
    print("| Case | Scheme | max abs(h - h0) | max abs(hu) |")
    for case in ("lake-at-rest", "lake-at-rest-step", "lake-at-rest-island"):
        for scheme in ("kinetic", "rusanov"):
            s, _ = run(case, scheme, "--no-snapshots")
            print(f"| {case} | {scheme} | {s['h_linf']:.1e} | {s['hu_linf']:.1e} |")


def wet_surface(h, z, threshold=1e-3):
    """Free surface h + z, with NaN on dry cells so that plots do not draw it there."""
    return np.where(h > threshold, h + z, np.nan)


def run_all_snapshots(case, scheme, *extra):
    """Runs the solver; returns (times, list of snapshots)."""
    with tempfile.TemporaryDirectory() as out:
        subprocess.run([SOLVER, "--case", case, "--scheme", scheme, "--output-dir", out, *extra],
                       check=True, capture_output=True)
        times = np.loadtxt(os.path.join(out, "times.dat"), ndmin=2)[:, 1]
        snapshots = [np.loadtxt(p, ndmin=2) for p in sorted(glob.glob(os.path.join(out, "snapshot_*.dat")))]
    return times, snapshots


def thacker_figure():
    period = 2.0 * math.pi * 8.0 / math.sqrt(2.0 * 9.81 * 0.5)
    fractions = (0.0, 0.25, 0.5, 1.0)
    fig, axes = plt.subplots(1, len(fractions), figsize=(13, 3.2), sharey=True)
    for ax, fraction in zip(axes, fractions):
        t_end = max(fraction * period, 1e-6)
        for scheme in ("kinetic", "rusanov"):
            _, data = run("thacker", scheme, "--cells", "200", "--t-end", repr(t_end))
            ax.plot(data[:, 0], wet_surface(data[:, 2], data[:, 1]), color=COLORS[scheme], label=scheme)
        wet = data[:, 4] > 0
        ax.plot(data[wet, 0], data[wet, 4] + data[wet, 1], "--", color="black", linewidth=1, label="exact")
        ax.fill_between(data[:, 0], data[:, 1].min() - 0.05, data[:, 1], color="0.8")
        ax.set_title(f"t = {fraction:g} T" if fraction else "t = 0")
        ax.set_xlabel("x (m)")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("elevation (m)")
    axes[-1].legend(fontsize="small")
    save(fig, "thacker.png")

    print("\nThacker, L1 error of h after one period:")
    print("| Scheme | 100 | 200 | 400 | 800 | 1600 |")
    for scheme in ("kinetic", "rusanov"):
        errors = [run("thacker", scheme, "--cells", str(n), "--no-snapshots")[0]["h_l1"] for n in (100, 200, 400, 800, 1600)]
        print(f"| {scheme} | " + " | ".join(f"{e:.2e}" for e in errors) + " |")

    print("\nThacker, 400 cells: tilt of the free surface (fit over |x| < 5) and mean velocity")
    g, h0, a, b = 9.81, 0.5, 8.0, 0.6
    omega = math.sqrt(2 * g * h0) / a
    for fraction in (0.5, 1.0):
        t = fraction * period
        line = f"t = {fraction:g} T: exact tilt {-(b * omega / g) * math.cos(omega * t):+.5f}, u {b * math.sin(omega * t):+.4f}"
        for scheme in ("kinetic", "rusanov"):
            _, data = run("thacker", scheme, "--cells", "400", "--t-end", repr(t))
            x, z, h, hu = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
            core = (h > 1e-3) & (np.abs(x) < 5)
            tilt = np.polyfit(x[core], (h + z)[core], 1)[0]
            line += f" | {scheme}: tilt {tilt:+.5f}, u {hu[core].sum() / h[core].sum():+.4f}"
        print(line)


def steady_flows_figure():
    cases = (("bump-subcritical", "subcritical"), ("bump-transcritical", "transcritical"),
             ("bump-hydraulic-jump", "transcritical with a hydraulic jump"))
    fig, axes = plt.subplots(2, 3, figsize=(13, 5.5), sharex=True)
    print("\nSteady flows over the bump, 400 cells, L1 error of h and max relative deviation of the discharge:")
    print("| Case | Scheme | L1 error of h | max abs(hu - q) / q |")
    for column, (case, title) in enumerate(cases):
        for scheme in ("kinetic", "rusanov"):
            summary, data = run(case, scheme, "--cells", "400")
            q = data[0, 5]
            axes[0, column].plot(data[:, 0], data[:, 2] + data[:, 1], color=COLORS[scheme], label=scheme)
            axes[1, column].plot(data[:, 0], data[:, 3], color=COLORS[scheme], label=scheme)
            print(f"| {case} | {scheme} | {summary['h_l1']:.2e} | {np.abs(data[:, 3] - q).max() / q:.1e} |")
        axes[0, column].plot(data[:, 0], data[:, 4] + data[:, 1], "--", color="black", linewidth=1, label="exact")
        axes[0, column].fill_between(data[:, 0], -0.02, data[:, 1], color="0.8")
        axes[1, column].plot(data[:, 0], data[:, 5], "--", color="black", linewidth=1, label="exact")
        axes[0, column].set_title(title)
        axes[1, column].set_xlabel("x (m)")
        for ax in axes[:, column]:
            ax.grid(alpha=0.3)
    axes[0, 0].set_ylabel("h + z (m)")
    axes[1, 0].set_ylabel("hu (m²/s)")
    axes[0, 0].legend(fontsize="small")
    save(fig, "steady-flows.png")


def runup_figure():
    times, snapshots = run_all_snapshots("beach-runup", "kinetic", "--snapshot-every", "10")
    moments = (0.0, 2.0, 3.5, 5.0, 6.5, 9.0)
    fig, axes = plt.subplots(2, 3, figsize=(13, 5.5), sharex=True, sharey=True)
    for ax, moment in zip(axes.flat, moments):
        data = snapshots[int(np.argmin(np.abs(times - moment)))]
        ax.fill_between(data[:, 0], -0.05, data[:, 1], color="0.8")
        ax.plot(data[:, 0], wet_surface(data[:, 2], data[:, 1]), color=COLORS["kinetic"])
        ax.axhline(0.5, color="0.5", linestyle=":", linewidth=1)
        ax.set_title(f"t = {moment:g} s")
        ax.set_ylim(-0.05, 1.2)
        ax.grid(alpha=0.3)
    for ax in axes[1]:
        ax.set_xlabel("x (m)")
    for ax in axes[:, 0]:
        ax.set_ylabel("elevation (m)")
    save(fig, "beach-runup.png")

    print("\nBeach run-up (shoreline = last cell with h > 1e-3): maximum run-up above still water, lowest shoreline after it")
    print("| Cells | Scheme | max run-up (m) | backwash shoreline x (m) |")
    for cells in (250, 500, 1000, 2000):
        for scheme in ("kinetic", "rusanov"):
            _, snaps = run_all_snapshots("beach-runup", scheme, "--cells", str(cells), "--snapshot-every", str(max(1, cells // 100)))
            shore = [d[np.nonzero(d[:, 2] > 1e-3)[0][-1], 0] for d in snaps]
            peak = int(np.argmax(shore))
            print(f"| {cells} | {scheme} | {0.1 * shore[peak] - 0.5:.3f} | {min(shore[peak:]):.2f} |")


def main():
    if not os.path.exists(SOLVER):
        raise SystemExit("build the solver first: make")
    os.makedirs(ASSETS, exist_ok=True)
    topography_figure()
    dam_break_figure("dam-break-dry", "dam-break-dry.png", "Dam break on a dry bed (Ritter)")
    dam_break_figure("dam-break-wet", "dam-break-wet.png", "Dam break on a wet bed (Stoker)")
    dam_break_figure("dam-break-bump", "dam-break-bump.png", "Dam break over the bump")
    lake_at_rest_table()
    space_convergence()
    time_convergence()
    euler_versus_rk()
    thacker_figure()
    steady_flows_figure()
    runup_figure()


if __name__ == "__main__":
    main()
