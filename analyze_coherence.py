#!/usr/bin/env python3
"""
TPF stellar coherence / phase diagnostics from a snapshot CSV (SI: m, m/s).

Ring selection: the radial bin with maximum coherence C = 1 - σ_v/⟨|v|⟩ among bins
with at least 10 stars (not peak surface density).

Example:
  python analyze_coherence.py
  python analyze_coherence.py --snapshot cpp_sim/outputs/run1/snapshot_20000.csv
  # Default plot: <snapshot_dir>/coherence_plot.png unless --plot is given
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

G_SI = 6.6743e-11


def load_snapshot_csv(path: Path) -> pd.DataFrame:
    """C++ galaxy_sim: line0 # step,..., line1 header, then data."""
    lines = path.read_text().strip().splitlines()
    if len(lines) < 3:
        raise ValueError(f"Need at least 3 lines in {path}")

    header = [c.strip().lower() for c in lines[1].split(",")]
    rows = []
    for ln in lines[2:]:
        ln = ln.strip()
        if not ln:
            continue
        parts = [p.strip() for p in ln.split(",")]
        try:
            rows.append([float(p) for p in parts])
        except ValueError:
            continue
    if not rows:
        raise ValueError(f"No numeric rows in {path}")

    df = pd.DataFrame(rows, columns=header[: len(rows[0])])
    if "x" not in df.columns and len(df.columns) >= 5:
        df = pd.DataFrame(rows)
        df.columns = ["i", "x", "y", "vx", "vy", "mass"][: df.shape[1]]
    return df


def load_bh_mass_run_info(snapshot_path: Path) -> float | None:
    """Try snapshot_dir/run_info.txt for bh_mass (kg)."""
    run_info = snapshot_path.parent / "run_info.txt"
    if not run_info.exists():
        return None
    for line in run_info.read_text().splitlines():
        if line.startswith("bh_mass\t"):
            try:
                return float(line.split("\t", 1)[1].strip())
            except (IndexError, ValueError):
                return None
    return None


def annular_bin_edges(r_max: float, n_bins: int) -> np.ndarray:
    return np.linspace(0.0, r_max, n_bins + 1)


def bin_index(r: np.ndarray, edges: np.ndarray) -> np.ndarray:
    return np.clip(np.digitize(r, edges) - 1, 0, len(edges) - 2)


def main() -> int:
    parser = argparse.ArgumentParser(description="TPF coherence / ring analysis from snapshot CSV")
    parser.add_argument(
        "--snapshot",
        type=Path,
        default=Path("cpp_sim/outputs/20260329_025714/snapshot_20000.csv"),
        help="Path to snapshot_*.csv",
    )
    parser.add_argument(
        "--n-bins",
        type=int,
        default=48,
        help="Number of equal-width radial bins from 0 to max(r)",
    )
    parser.add_argument(
        "--bh-mass",
        type=float,
        default=None,
        help="Central mass for Newtonian comparison (kg); default from run_info or 1e41",
    )
    parser.add_argument(
        "--plot",
        type=Path,
        default=None,
        help="Output figure path (default: same directory as --snapshot, file coherence_plot.png)",
    )
    args = parser.parse_args()

    snap_path = args.snapshot.resolve()
    if not snap_path.is_file():
        raise SystemExit(f"Snapshot not found: {snap_path}")

    if args.plot is None:
        plot_out = snap_path.parent / "coherence_plot.png"
    else:
        plot_out = Path(args.plot).expanduser().resolve()

    df = load_snapshot_csv(snap_path)
    x = df["x"].to_numpy(dtype=float)
    y = df["y"].to_numpy(dtype=float)
    vx = df["vx"].to_numpy(dtype=float)
    vy = df["vy"].to_numpy(dtype=float)
    mass = df["mass"].to_numpy(dtype=float) if "mass" in df.columns else None

    r = np.hypot(x, y)
    v_mag = np.hypot(vx, vy)
    # Avoid r=0 for angles
    r_safe = np.maximum(r, 1e-300)
    rx, ry = x / r_safe, y / r_safe
    v_r = vx * rx + vy * ry
    v_theta = x * vy - y * vx
    v_theta /= r_safe

    r_max = float(np.max(r))
    edges = annular_bin_edges(r_max, args.n_bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)
    bin_id = bin_index(r, edges)

    # Most coherent ring: maximize C = 1 - σ_v/⟨|v|⟩ per bin (not peak density)
    coherence_per_bin = np.full(args.n_bins, np.nan)
    for b in range(args.n_bins):
        m = bin_id == b
        if np.sum(m) < 10:  # statistical threshold to avoid noise in sparse bins
            continue
        v_bin = v_mag[m]
        std_v = np.std(v_bin, ddof=1)
        mean_v_bin = np.mean(v_bin)
        if mean_v_bin > 1e-10:
            coherence_per_bin[b] = 1.0 - (std_v / mean_v_bin)

    if not np.any(np.isfinite(coherence_per_bin)):
        raise SystemExit("No bins meet the statistical threshold for coherence analysis.")

    peak_bin = int(np.nanargmax(coherence_per_bin))
    ring_r = float(centers[peak_bin])
    coherence_c = float(coherence_per_bin[peak_bin])

    in_ring = bin_id == peak_bin
    v_ring = v_mag[in_ring]
    sigma_v = float(np.std(v_ring, ddof=1))
    mean_v = float(np.mean(v_ring))

    ratio_rt = np.abs(v_r[in_ring]) / np.maximum(np.abs(v_theta[in_ring]), 1e-30)
    mean_ratio_rt = float(np.median(ratio_rt))  # robust to outliers
    mean_ratio_rt_mean = float(np.mean(ratio_rt))

    # Per-bin velocity dispersion for plot
    sigma_per_bin = np.full(args.n_bins, np.nan)
    for b in range(args.n_bins):
        m = bin_id == b
        if np.sum(m) < 2:
            continue
        sigma_per_bin[b] = np.std(v_mag[m], ddof=1)

    # Newtonian: Keplerian circular speed at ring radius (point mass)
    m_bh = args.bh_mass
    if m_bh is None:
        m_bh = load_bh_mass_run_info(snap_path)
    if m_bh is None:
        m_bh = 1.0e41
    v_newt_kepler = float(np.sqrt(G_SI * m_bh / max(ring_r, 1e-30)))
    if mass is not None:
        m_disk = float(np.sum(mass))
        m_enc = m_bh + m_disk  # crude enclosed mass for extended comparison
        v_newt_enc = float(np.sqrt(G_SI * m_enc / max(ring_r, 1e-30)))
    else:
        m_disk = float("nan")
        v_newt_enc = float("nan")

    # --- Text summary ---
    print("=== TPF coherence analysis (SI) ===")
    print(f"Snapshot: {snap_path}")
    print(f"Plot output: {plot_out}")
    print(f"Stars: {len(r)}")
    print(f"Ring radius (peak coherence C = 1 - σ_v/⟨|v|⟩): {ring_r:.6g} m")
    print(f"  bin index {peak_bin}/{args.n_bins - 1}, Δr ~ {widths[peak_bin]:.6g} m (≥10 stars/bin for C)")
    print(f"At ring: mean |v| = {mean_v:.6g} m/s, σ_v = {sigma_v:.6g} m/s")
    print(f"Coherence factor C = 1 - σ_v/⟨|v|⟩ = {coherence_c:.6f}")
    print(f"  (C → 1 solid-body-like; C → 0 disordered)")
    print(f"Phase |v_r|/|v_θ| in ring (median, mean): {mean_ratio_rt:.6g}, {mean_ratio_rt_mean:.6g}")
    print()
    print("--- Newtonian reference (same central mass, point-mass Kepler at R_ring) ---")
    print(f"M_bh used: {m_bh:.6g} kg")
    print(f"v_circ,Kepler(R_ring) = √(G M_bh / R) = {v_newt_kepler:.6g} m/s")
    print(f"Sim ⟨|v|⟩ at ring / v_Kepler = {mean_v / v_newt_kepler:.6g}")
    if mass is not None and np.isfinite(m_disk):
        print(f"Approx M_disk = Σ m_star = {m_disk:.6g} kg")
        print(f"v_circ(R_ring) if all mass enclosed ≈ √(G(M_bh+M_disk)/R) = {v_newt_enc:.6g} m/s")

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(9, 5))
    valid = np.isfinite(sigma_per_bin)
    ax.plot(centers[valid], sigma_per_bin[valid], "b.-", lw=1.2, markersize=4, label=r"$\sigma_v(|v|)$ per bin")
    ax.axvline(ring_r, color="crimson", ls="--", lw=1.0, label=f"Peak-coherence ring R = {ring_r:.3g} m")
    ax.set_xlabel("Radius r (m)")
    ax.set_ylabel(r"Velocity dispersion $\sigma_v$ (m/s)")
    ax.set_title(r"Coherence diagnostic: $\sigma_v(|v|)$ vs $r$; ring = argmax bin coherence ($\geq$10 stars/bin)")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    plot_out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_out, dpi=150)
    plt.close(fig)
    print()
    print(f"Saved coherence plot: {plot_out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
