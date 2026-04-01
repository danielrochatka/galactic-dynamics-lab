#!/usr/bin/env python3
"""
Visualization only: reads tpf_v11_earth_moon_line_correspondence_benchmark.csv (+ optional .plot_meta.json)
and writes PNGs. Does not change physics, calibration, or CSV data.

Requires: numpy, matplotlib (same stack as plot_cpp_run.py).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PNG_COMPARE = "tpf_v11_earth_moon_line_correspondence_benchmark_compare.png"
PNG_DIFF = "tpf_v11_earth_moon_line_correspondence_benchmark_difference.png"
PNG_EARTH_ZOOM = "tpf_v11_earth_moon_line_correspondence_benchmark_earth_zoom.png"
PNG_NORMALIZED = "tpf_v11_earth_moon_line_correspondence_benchmark_normalized_shape.png"

SUMMARY_NAME = "tpf_v11_earth_moon_line_correspondence_benchmark_summary.txt"
CSV_NAME = "tpf_v11_earth_moon_line_correspondence_benchmark.csv"
META_NAME = "tpf_v11_earth_moon_line_correspondence_benchmark.plot_meta.json"


def _subtitle() -> str:
    return (
        "v11 weak-field correspondence benchmark only — not full TPF dynamics\n"
        "correspondence / audit layer · ΔC omitted · no orbit integration"
    )


def _load_meta(out: Path) -> dict:
    p = out / META_NAME
    if not p.is_file():
        return {}
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}


def _vertical_guides(
    ax, re_m: float, xb_m: float, d_m: float, xlim_m: tuple[float, float]
) -> None:
    """Draw guides in Mm on x-axis; xlim_m is (min,max) in meters."""
    lo, hi = xlim_m
    if re_m > 0 and lo <= re_m <= hi:
        ax.axvline(re_m / 1e6, color="#2ca02c", lw=1.2, ls="-", zorder=0, label="Earth surface (calibration)")
    if xb_m > 0 and lo <= xb_m <= hi:
        ax.axvline(xb_m / 1e6, color="#ff7f0e", lw=1.2, ls="--", zorder=0, label="Barycenter")
    if d_m > 0 and lo <= d_m <= hi:
        ax.axvline(d_m / 1e6, color="#7f7f7f", lw=1.2, ls=":", zorder=0, label="Moon center (x = D)")


def _moon_beyond_note(ax, d_m: float, x_max_m: float) -> None:
    if d_m > x_max_m + 1e-3:
        ax.text(
            0.98,
            0.04,
            f"Moon center at D = {d_m / 1e6:.3f} Mm (past plotted x-range)",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
            color="#555555",
        )


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: plot_v11_earth_moon_line_benchmark.py <output_dir>", file=sys.stderr)
        return 2
    out = Path(sys.argv[1]).resolve()
    csv_path = out / CSV_NAME
    if not csv_path.is_file():
        print(f"plot_v11_earth_moon_line_benchmark: missing {csv_path}", file=sys.stderr)
        return 1

    df = pd.read_csv(csv_path)
    x_m = df["x_m"].to_numpy(dtype=float)
    a_tpf = df["a_tpf_correspondence_Eq45_m_s2"].to_numpy(dtype=float)
    a_new = df["a_newtonian_line_Eq46_m_s2"].to_numpy(dtype=float)
    diff = df["benchmark_difference_TPF_Eq45_minus_Newtonian_Eq46_m_s2"].to_numpy(dtype=float)

    d_row = float(df["D_m"].iloc[0])
    xb_row = float(df["x_b_m"].iloc[0])
    meta = _load_meta(out)
    re_m = float(meta.get("v11_em_calib_surface_radius_m", 0.0) or 0.0)
    d_m = float(meta.get("v11_em_mean_distance_m", d_row) or d_row)
    xb_m = float(meta.get("barycenter_m", xb_row) or xb_row)

    x_mm = x_m / 1e6

    # --- 1) Full comparison ---
    fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
    ax.plot(x_mm, a_tpf, color="#1f77b4", lw=2.0, label="TPF correspondence (Eq. 45)")
    ax.plot(x_mm, a_new, color="#d62728", lw=2.0, ls="--", label="Newtonian benchmark (Eq. 46)")
    ax.set_xlabel("Position from Earth center (Mm)")
    ax.set_ylabel(r"Acceleration $a_x$ (m/s²) — co-rotating frame")
    ax.set_title("Earth–Moon line of centers: correspondence vs Newtonian benchmark")
    fig.text(0.5, 0.02, _subtitle(), ha="center", fontsize=8, color="#444444")
    fig.subplots_adjust(bottom=0.18)
    xlim_m = (float(np.min(x_m)), float(np.max(x_m)))
    _vertical_guides(ax, re_m, xb_m, d_m, xlim_m)
    _moon_beyond_note(ax, d_m, xlim_m[1])
    ax.set_xlim(xlim_m[0] / 1e6, xlim_m[1] / 1e6)
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.savefig(out / PNG_COMPARE, bbox_inches="tight")
    plt.close(fig)

    # --- 2) Difference only ---
    fig, ax = plt.subplots(figsize=(10, 5), dpi=150)
    ax.plot(x_mm, diff, color="#9467bd", lw=2.0, label="Benchmark difference (algebra only; not new dynamics)")
    ax.set_xlabel("Position from Earth center (Mm)")
    ax.set_ylabel(r"$\Delta a = a_{\mathrm{TPF\,(45)}} - a_{\mathrm{Newtonian\,(46)}}$ (m/s²)")
    ax.set_title("Benchmark difference (TPF Eq. 45 minus Newtonian Eq. 46)")
    fig.text(0.5, 0.02, _subtitle(), ha="center", fontsize=8, color="#444444")
    fig.subplots_adjust(bottom=0.18)
    _vertical_guides(ax, re_m, xb_m, d_m, xlim_m)
    _moon_beyond_note(ax, d_m, xlim_m[1])
    ax.set_xlim(xlim_m[0] / 1e6, xlim_m[1] / 1e6)
    ax.legend(loc="best", fontsize=9)
    ax.axhline(0.0, color="#aaaaaa", lw=0.8)
    ax.grid(True, alpha=0.3)
    fig.savefig(out / PNG_DIFF, bbox_inches="tight")
    plt.close(fig)

    # --- 3) Earth-side zoom ---
    zoom_hi_m = min(5.0e7, float(np.max(x_m)), 0.15 * d_m)  # 50 Mm cap, data cap, or 15% of D
    zoom_hi_m = max(zoom_hi_m, float(np.min(x_m)) + 1e3)
    mask = x_m <= zoom_hi_m
    if np.count_nonzero(mask) < 2:
        mask = np.ones_like(x_m, dtype=bool)
    xz = x_m[mask]
    xzm = xz / 1e6
    fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
    ax.plot(xzm, a_tpf[mask], color="#1f77b4", lw=2.0, label="TPF correspondence (Eq. 45)")
    ax.plot(xzm, a_new[mask], color="#d62728", lw=2.0, ls="--", label="Newtonian benchmark (Eq. 46)")
    ax.set_xlabel("Position from Earth center (Mm)")
    ax.set_ylabel(r"$a_x$ (m/s²)")
    ax.set_title("Earth-side zoom (same curves; restricted x-range)")
    fig.text(0.5, 0.02, _subtitle(), ha="center", fontsize=8, color="#444444")
    fig.subplots_adjust(bottom=0.18)
    zlim = (0.0, float(np.max(xz)))
    _vertical_guides(ax, re_m, xb_m, d_m, zlim)
    _moon_beyond_note(ax, d_m, zlim[1])
    ax.set_xlim(0.0, zlim[1] / 1e6)
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.savefig(out / PNG_EARTH_ZOOM, bbox_inches="tight")
    plt.close(fig)

    # --- 4) Normalized shape (visualization only) ---
    scale = max(np.max(np.abs(a_tpf)), np.max(np.abs(a_new)), 1e-300)
    fig, ax = plt.subplots(figsize=(10, 5), dpi=150)
    ax.plot(x_mm, a_tpf / scale, color="#1f77b4", lw=2.0, label=r"TPF (Eq. 45) / max|$a$|")
    ax.plot(x_mm, a_new / scale, color="#d62728", lw=2.0, ls="--", label=r"Newtonian (Eq. 46) / max|$a$|")
    ax.set_xlabel("Position from Earth center (Mm)")
    ax.set_ylabel("Unitless (divide by max(|a_TPF|, |a_Newt|) on tabulated range)")
    ax.set_title("Shape comparison (normalized — visualization only)")
    fig.text(0.5, 0.02, _subtitle(), ha="center", fontsize=8, color="#444444")
    fig.subplots_adjust(bottom=0.18)
    _vertical_guides(ax, re_m, xb_m, d_m, xlim_m)
    _moon_beyond_note(ax, d_m, xlim_m[1])
    ax.set_xlim(xlim_m[0] / 1e6, xlim_m[1] / 1e6)
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.savefig(out / PNG_NORMALIZED, bbox_inches="tight")
    plt.close(fig)

    # Append summary
    summary_path = out / SUMMARY_NAME
    if summary_path.is_file():
        block = (
            "\n=== Auto-generated PNG plots (visualization only; same data as CSV) ===\n"
            f"  {out / PNG_COMPARE}\n"
            f"  {out / PNG_DIFF}\n"
            f"  {out / PNG_EARTH_ZOOM}\n"
            f"  {out / PNG_NORMALIZED}\n"
        )
        with open(summary_path, "a", encoding="utf-8") as sf:
            sf.write(block)

    print(f"plot_v11_earth_moon_line_benchmark: wrote {PNG_COMPARE}, {PNG_DIFF}, {PNG_EARTH_ZOOM}, {PNG_NORMALIZED}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
