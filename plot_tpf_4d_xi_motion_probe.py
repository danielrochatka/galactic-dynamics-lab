#!/usr/bin/env python3
"""Render trajectory visualizations for tpf_4d_xi_motion_probe_benchmark outputs."""

from __future__ import annotations

import sys
from pathlib import Path

TRAJ_CSV = "tpf_4d_xi_motion_probe_trajectories.csv"

REQUIRED_PNGS = {
    "xy": "tpf_4d_xi_motion_probe_xy_trajectories.png",
    "xz": "tpf_4d_xi_motion_probe_xz_trajectories.png",
    "yz": "tpf_4d_xi_motion_probe_yz_trajectories.png",
    "radius": "tpf_4d_xi_motion_probe_radius_vs_time.png",
    "accel": "tpf_4d_xi_motion_probe_acceleration_norm_vs_time.png",
}


def _import_plot_stack():
    try:
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
    except Exception as exc:  # pragma: no cover - environment dependent
        print(
            "plot_tpf_4d_xi_motion_probe: warning: plotting dependencies unavailable "
            f"(numpy/pandas/matplotlib): {exc}",
            file=sys.stderr,
        )
        return None
    return plt, np, pd


def _title_prefix() -> str:
    return (
        "dynamic probe-motion benchmark using Xi-direct acceleration readout\n"
        "a=-K_xi*Xi_spatial | fixed-source Stage 7B benchmark"
    )


def _plot_plane(plt, df, out_path: Path, x_col: str, y_col: str, plane_name: str, mark_start_end: bool = False, equal_aspect: bool = False) -> bool:
    if x_col not in df.columns or y_col not in df.columns:
        print(
            f"plot_tpf_4d_xi_motion_probe: skipping {out_path.name}; missing columns {x_col}/{y_col}",
            file=sys.stderr,
        )
        return False

    fig, ax = plt.subplots(figsize=(7.5, 6), dpi=150)
    for probe_id, group in df.groupby("probe_id", sort=True):
        g = group.sort_values("time")
        ax.plot(g[x_col], g[y_col], linewidth=1.0, alpha=0.9, label=f"probe {probe_id}")
        if mark_start_end and not g.empty:
            ax.scatter([g.iloc[0][x_col]], [g.iloc[0][y_col]], marker="o", s=28, c="green", edgecolors="black", linewidths=0.4)
            ax.scatter([g.iloc[-1][x_col]], [g.iloc[-1][y_col]], marker="X", s=34, c="red", edgecolors="black", linewidths=0.4)

    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title(f"GravityXiMotionReadout_v1: {plane_name} probe trajectories\n{_title_prefix()}")
    if equal_aspect:
        ax.set_aspect("equal", adjustable="box")
    if df["probe_id"].nunique() <= 12:
        ax.legend(loc="best", fontsize=7)
    ax.grid(alpha=0.25)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return True


def _plot_time_series(plt, df, out_path: Path, value_col: str, y_label: str, title_label: str) -> bool:
    if "time" not in df.columns or value_col not in df.columns:
        print(
            f"plot_tpf_4d_xi_motion_probe: skipping {out_path.name}; missing columns time/{value_col}",
            file=sys.stderr,
        )
        return False

    fig, ax = plt.subplots(figsize=(8.0, 5.2), dpi=150)
    for probe_id, group in df.groupby("probe_id", sort=True):
        g = group.sort_values("time")
        ax.plot(g["time"], g[value_col], linewidth=1.0, alpha=0.9, label=f"probe {probe_id}")

    ax.set_xlabel("time")
    ax.set_ylabel(y_label)
    ax.set_title(f"{title_label}\n{_title_prefix()}")
    if df["probe_id"].nunique() <= 12:
        ax.legend(loc="best", fontsize=7)
    ax.grid(alpha=0.25)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return True


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: python3 plot_tpf_4d_xi_motion_probe.py <output_dir>", file=sys.stderr)
        return 2

    out_dir = Path(sys.argv[1]).resolve()
    stack = _import_plot_stack()
    if stack is None:
        return 0
    plt, _np, pd = stack

    traj_path = out_dir / TRAJ_CSV
    if not traj_path.is_file():
        print(f"plot_tpf_4d_xi_motion_probe: warning: missing {TRAJ_CSV} in {out_dir}", file=sys.stderr)
        return 0

    df = pd.read_csv(traj_path)
    if df.empty:
        print("plot_tpf_4d_xi_motion_probe: warning: trajectory CSV is empty; no PNGs generated", file=sys.stderr)
        return 0

    generated: list[str] = []

    targets = [
        (_plot_plane, (df, out_dir / REQUIRED_PNGS["xy"], "x", "y", "xy", True, True)),
        (_plot_plane, (df, out_dir / REQUIRED_PNGS["xz"], "x", "z", "xz", False, False)),
        (_plot_plane, (df, out_dir / REQUIRED_PNGS["yz"], "y", "z", "yz", False, False)),
    ]
    for fn, args in targets:
        out_path = args[1]
        if fn(plt, *args):
            generated.append(out_path.name)

    if _plot_time_series(
        plt,
        df,
        out_dir / REQUIRED_PNGS["radius"],
        "r_origin",
        "r_origin vs time",
        "GravityXiMotionReadout_v1: probe radius evolution",
    ):
        generated.append(REQUIRED_PNGS["radius"])

    if _plot_time_series(
        plt,
        df,
        out_dir / REQUIRED_PNGS["accel"],
        "a_norm",
        "a_norm",
        "GravityXiMotionReadout_v1: acceleration norm vs time",
    ):
        generated.append(REQUIRED_PNGS["accel"])

    if generated:
        print("plot_tpf_4d_xi_motion_probe: generated PNGs:")
        for name in generated:
            print(f"  {name}")
    else:
        print("plot_tpf_4d_xi_motion_probe: warning: no PNGs generated", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
