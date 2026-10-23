#!/usr/bin/env python3
"""Plots the snapshots written by kinetic-solver.

    python3 scripts/plot.py output              # last snapshot -> output/final.png
    python3 scripts/plot.py output --all        # one PNG per snapshot
    python3 scripts/plot.py output --animate    # output/animation.gif (needs Pillow)
"""

import argparse
import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def read_snapshot(path):
    with open(path) as handle:
        header = handle.readline()
    time = float(header.split("t=")[1])
    data = np.loadtxt(path, comments="#", ndmin=2)
    return time, data


def draw(axes, time, data, title):
    x, z, h, hu = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    top, bottom = axes
    top.clear()
    bottom.clear()
    top.fill_between(x, z.min() - 0.05, z, color="0.75", label="bottom")
    top.plot(x, h + z, label="free surface h + z")
    bottom.plot(x, hu, label="discharge hu")
    if data.shape[1] >= 6:
        top.plot(x, data[:, 4] + z, "k--", linewidth=1, label="exact")
        bottom.plot(x, data[:, 5], "k--", linewidth=1, label="exact")
    top.set_title(f"{title}   t = {time:.3f} s")
    top.set_ylabel("elevation (m)")
    bottom.set_ylabel("hu (m²/s)")
    bottom.set_xlabel("x (m)")
    for ax in axes:
        ax.legend(loc="upper right", fontsize="small")
        ax.grid(alpha=0.3)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("directory", help="output directory of kinetic-solver")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--all", action="store_true", help="write one PNG per snapshot")
    group.add_argument("--animate", action="store_true", help="write an animated GIF")
    args = parser.parse_args()

    paths = sorted(glob.glob(os.path.join(args.directory, "snapshot_*.dat")))
    if not paths:
        parser.error(f"no snapshot_*.dat file in {args.directory}")
    title = os.path.basename(os.path.normpath(args.directory))
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(9, 6))

    if args.animate:
        from matplotlib.animation import FuncAnimation, PillowWriter

        frames = [read_snapshot(p) for p in paths]
        animation = FuncAnimation(fig, lambda k: draw(axes, *frames[k], title), frames=len(frames))
        target = os.path.join(args.directory, "animation.gif")
        animation.save(target, writer=PillowWriter(fps=10))
    elif args.all:
        for path in paths:
            draw(axes, *read_snapshot(path), title)
            target = path.replace(".dat", ".png")
            fig.savefig(target, dpi=120, bbox_inches="tight")
    else:
        draw(axes, *read_snapshot(paths[-1]), title)
        target = os.path.join(args.directory, "final.png")
        fig.savefig(target, dpi=120, bbox_inches="tight")
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
