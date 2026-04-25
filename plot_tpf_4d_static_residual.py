#!/usr/bin/env python3
"""Post-process view-plane inspection artifacts for tpf_4d_static_residual_benchmark."""

from __future__ import annotations

import math
import sys
from pathlib import Path

CSV_BY_PLANE = {
    "xy": "tpf_4d_static_residual_slice_xy.csv",
    "xz": "tpf_4d_static_residual_slice_xz.csv",
    "yz": "tpf_4d_static_residual_slice_yz.csv",
}

SOURCES_CSV = "tpf_4d_static_residual_sources.csv"


def _import_plot_stack():
    try:
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
    except Exception as exc:  # pragma: no cover - environment-dependent
        print(
            "plot_tpf_4d_static_residual: warning: plotting dependencies unavailable "
            f"(numpy/pandas/matplotlib): {exc}",
            file=sys.stderr,
        )
        return None
    return plt, np, pd


def _safe_log10(arr, np):
    return np.log10(np.clip(np.abs(arr), 1e-300, None))


def _heatmap(plt, np, pd, df, x_col, y_col, value_col, title, out_path, sources_df=None, cmap="viridis"):
    if value_col not in df.columns:
        print(f"plot_tpf_4d_static_residual: skipping {out_path.name}; missing column {value_col}")
        return False

    pivot = df.pivot_table(index=y_col, columns=x_col, values=value_col, aggfunc="mean")
    if pivot.empty:
        return False

    z = pivot.to_numpy(dtype=float)
    if "residual" in value_col or "invariant" in value_col:
        z = _safe_log10(z, np)
        colorbar_label = f"log10(|{value_col}|)"
    else:
        colorbar_label = value_col

    fig, ax = plt.subplots(figsize=(7.5, 6), dpi=150)
    x_min, x_max = float(np.min(pivot.columns)), float(np.max(pivot.columns))
    y_min, y_max = float(np.min(pivot.index)), float(np.max(pivot.index))
    im = ax.imshow(z, extent=[x_min, x_max, y_min, y_max], origin="lower", aspect="equal", cmap=cmap)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label(colorbar_label)

    if "is_boundary" in df.columns:
        bd = df[df["is_boundary"] == 1]
        if not bd.empty:
            ax.scatter(bd[x_col], bd[y_col], s=5, c="white", alpha=0.25, marker="s", label="boundary")
    if "is_near_source" in df.columns:
        ns = df[df["is_near_source"] == 1]
        if not ns.empty:
            ax.scatter(ns[x_col], ns[y_col], s=7, c="red", alpha=0.4, marker="x", label="near-source")
    if sources_df is not None and x_col in sources_df.columns and y_col in sources_df.columns:
        ax.scatter(sources_df[x_col], sources_df[y_col], s=25, c="cyan", alpha=0.9, edgecolors="black", linewidths=0.5, marker="o",
                   label="sources")

    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title(title)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(loc="upper right", fontsize=8)
    fig.text(0.5, 0.02, "view-plane inspection artifact only; not full spatial-support physics domain", ha="center", fontsize=8)
    fig.subplots_adjust(bottom=0.12)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return True


def _quiver(plt, np, df, x_col, y_col, u_col, v_col, title, out_path):
    if u_col not in df.columns or v_col not in df.columns:
        return False
    fig, ax = plt.subplots(figsize=(7.5, 6), dpi=150)
    step = max(1, int(math.sqrt(max(len(df), 1)) // 24))
    q = df.iloc[::step]
    mag = np.sqrt(np.square(q[u_col].to_numpy(dtype=float)) + np.square(q[v_col].to_numpy(dtype=float)))
    ax.quiver(q[x_col], q[y_col], q[u_col], q[v_col], mag, angles="xy", scale_units="xy", scale=None, cmap="plasma")
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="box")
    fig.text(0.5, 0.02, "Xi projection quiver for view-plane inspection only", ha="center", fontsize=8)
    fig.subplots_adjust(bottom=0.12)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return True


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: plot_tpf_4d_static_residual.py <output_dir>", file=sys.stderr)
        return 2

    out_dir = Path(sys.argv[1]).resolve()
    stack = _import_plot_stack()
    if stack is None:
        return 0
    plt, np, pd = stack

    sources_path = out_dir / SOURCES_CSV
    sources_df = pd.read_csv(sources_path) if sources_path.is_file() else None

    generated = []
    for plane, csv_name in CSV_BY_PLANE.items():
        csv_path = out_dir / csv_name
        if not csv_path.is_file():
            print(f"plot_tpf_4d_static_residual: skipping plane {plane}; missing {csv_name}")
            continue
        df = pd.read_csv(csv_path)

        axis = {"xy": ("x", "y"), "xz": ("x", "z"), "yz": ("y", "z")}[plane]
        x_col, y_col = axis

        targets = [
            ("normalized_residual", f"tpf_4d_static_residual_{plane}_normalized_residual.png", "view-plane inspection: normalized residual"),
            ("residual_spatial_norm", f"tpf_4d_static_residual_{plane}_residual_spatial_norm.png", "view-plane inspection: residual spatial norm"),
            ("xi_spatial_norm", f"tpf_4d_static_residual_{plane}_xi_spatial_norm.png", "view-plane inspection: Xi spatial norm"),
            ("invariant_I_4d", f"tpf_4d_static_residual_{plane}_invariant_I.png", "view-plane inspection: invariant I"),
        ]
        for value_col, png_name, title in targets:
            out_path = out_dir / png_name
            if _heatmap(plt, np, pd, df, x_col, y_col, value_col, title, out_path, sources_df=sources_df):
                generated.append(out_path.name)

        q_map = {"xy": ("xi_x", "xi_y"), "xz": ("xi_x", "xi_z"), "yz": ("xi_y", "xi_z")}
        u_col, v_col = q_map[plane]
        q_out = out_dir / f"tpf_4d_static_residual_{plane}_xi_quiver.png"
        if _quiver(plt, np, df, x_col, y_col, u_col, v_col, f"view-plane inspection: Xi projected vector ({plane})", q_out):
            generated.append(q_out.name)

    if generated:
        print("plot_tpf_4d_static_residual: generated PNGs:")
        for p in generated:
            print(f"  {p}")
    else:
        print("plot_tpf_4d_static_residual: no PNGs generated (missing input slices or required columns)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
