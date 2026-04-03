#!/usr/bin/env python3
"""
Rotation curve from N-body snapshot CSV: compare simulated speeds to a Newtonian Keplerian baseline.

Primary glob: output_*.csv (highest step inferred from filename). Falls back to snapshot_*.csv
(C++ galaxy_sim format) if no output_*.csv is found.

`save_rotation_curve_png` is also used by plot_cpp_run.py.

The red curve is **pure** Newtonian v = sqrt(G*M_bh/r): it does **not** use tpf_kappa or any TPF ledger.
Only the **blue scatter** comes from the snapshot (so κ changes star speeds for TPFCore runs, not the red line).

By default, when run_info.txt sits next to the CSV, M_bh is taken from the run; the x-axis is capped
to 2× galaxy_radius (escapers omitted from scatter by default) so the disk stays readable. If
galaxy_radius is missing, x_max falls back to 2e20 m.

Usage:
  python plot_rotation_curve.py
  python plot_rotation_curve.py cpp_sim/outputs/20260327_120000
  python plot_rotation_curve.py --search-dir . --output rotation_curve.png
"""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from text_layout import add_fitted_footer, set_fitted_title

from display_units import apply_suppress_tick_offset, rotation_curve_display

G_SI = 6.6743e-11
# Fallback M_bh when run_info has no positive bh_mass (legacy astronomical benchmark only).
NEWTONIAN_REFERENCE_M_BH = 1.0e41
# X-axis fallback when run_info has no galaxy_radius (legacy / huge-r plots).
ROTATION_CURVE_X_FALLBACK_MAX = 2.0e20
# Kept for backward compatibility in auto branches inside save_rotation_curve_png.
DEFAULT_X_MAX = ROTATION_CURVE_X_FALLBACK_MAX


def load_run_info(run_dir: Path) -> dict[str, str | int | float]:
    """Parse run_info.txt (tab-separated key, value) into a dict."""
    path = Path(run_dir) / "run_info.txt"
    if not path.exists():
        return {}
    info: dict[str, str | int | float] = {}
    for line in path.read_text().strip().splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t", 1)
        if len(parts) != 2:
            continue
        key, val = parts[0].strip(), parts[1].strip()
        if key in (
            "git_commit_full",
            "git_commit_short",
            "git_branch",
            "git_tag",
            "code_version_label",
        ):
            info[key] = val
            continue
        try:
            if "." in val or "e" in val.lower():
                info[key] = float(val)
            else:
                info[key] = int(val)
        except ValueError:
            info[key] = val
    return info


def newtonian_m_bh_from_run_info(run_info: dict[str, str | int | float]) -> float:
    """Central mass for Keplerian overlay; matches simulation bh_mass when available."""
    if not run_info:
        return float(NEWTONIAN_REFERENCE_M_BH)
    m = run_info.get("bh_mass")
    if isinstance(m, (int, float)) and float(m) > 0:
        return float(m)
    return float(NEWTONIAN_REFERENCE_M_BH)


def resolve_newtonian_overlay_mass(
    run_info: dict[str, str | int | float],
    cli_m_bh: float | None,
) -> tuple[float, str, bool]:
    """
    Resolve Newtonian overlay mass source for CLI usage.
    Precedence:
      1) explicit --M-bh
      2) run_info bh_mass when present and positive
      3) fallback benchmark NEWTONIAN_REFERENCE_M_BH (explicitly reported)
    Returns: (mass_kg, source_label, is_fallback_benchmark).
    """
    if cli_m_bh is not None:
        return float(cli_m_bh), "cli(--M-bh)", False
    m = run_info.get("bh_mass") if run_info else None
    if isinstance(m, (int, float)) and float(m) > 0:
        return float(m), "run_info(bh_mass)", False
    return float(NEWTONIAN_REFERENCE_M_BH), "fallback_benchmark(NEWTONIAN_REFERENCE_M_BH)", True


def rotation_curve_x_max(
    run_info: dict[str, str | int | float],
    r_data: np.ndarray,
) -> float:
    """
    Horizontal plot limit (m): twice configured galaxy_radius so escaped outliers do not stretch the disk.
    Falls back to ROTATION_CURVE_X_FALLBACK_MAX (2e20 m) when galaxy_radius is missing.
    """
    _ = r_data  # kept for API compatibility with callers; x extent comes from run_info only
    gr = run_info.get("galaxy_radius") if run_info else None
    if isinstance(gr, (int, float)) and float(gr) > 0:
        return 2.0 * float(gr)
    return float(ROTATION_CURVE_X_FALLBACK_MAX)


def scatter_label_from_run_info(run_info: dict[str, str | int | float]) -> str:
    pkg = run_info.get("physics_package", "")
    if isinstance(pkg, str) and pkg == "TPFCore":
        return "Simulated stars (TPFCore)"
    if isinstance(pkg, str) and pkg:
        return f"Simulated stars ({pkg})"
    return "Simulated stars"


def _step_from_stem(stem: str) -> int:
    m = re.search(r"(?:snapshot_|output_)(\d+)", stem, re.IGNORECASE)
    if m:
        return int(m.group(1))
    nums = re.findall(r"\d+", stem)
    return int(max(nums, key=len)) if nums else -1


def find_latest_csv(search_dir: Path) -> Path:
    """
    Prefer output_*.csv with largest step number; if none exist, use snapshot_*.csv (e.g. C++ runs).
    """
    search_dir = search_dir.resolve()
    output_files = [Path(p) for p in glob.glob(str(search_dir / "output_*.csv"))]
    if output_files:
        return max(output_files, key=lambda p: _step_from_stem(p.stem))
    snapshot_files = [Path(p) for p in glob.glob(str(search_dir / "snapshot_*.csv"))]
    if snapshot_files:
        return max(snapshot_files, key=lambda p: _step_from_stem(p.stem))
    raise FileNotFoundError(
        f"No output_*.csv or snapshot_*.csv under {search_dir}. "
        "Pass --search-dir to your run folder (e.g. cpp_sim/outputs/...)."
    )


def load_particle_kinematics(path: Path) -> tuple[np.ndarray, np.ndarray, str]:
    """
    Load x,y,(z) and vx,vy,(vz). C++ snapshots: comment line, then i,x,y,vx,vy,mass.
    Returns r_3d, v_3d, title_suffix (filename + step if parsed).
    """
    lines = path.read_text().strip().splitlines()
    if len(lines) < 3:
        raise ValueError(f"Not enough lines in {path}")

    step_note = ""
    m = re.search(r"step\s*,\s*(\d+)", lines[0], re.IGNORECASE)
    if m:
        step_note = f" (step {m.group(1)})"

    header = [c.strip().lower() for c in lines[1].split(",")]
    rows: list[list[float]] = []
    for ln in lines[2:]:
        ln = ln.strip()
        if not ln:
            continue
        parts = [p.strip() for p in ln.split(",")]
        try:
            vals = [float(p) for p in parts]
        except ValueError:
            continue
        rows.append(vals)

    if not rows:
        raise ValueError(f"No numeric data rows in {path}")

    arr = np.array(rows)

    def idx(name: str) -> int | None:
        try:
            return header.index(name)
        except ValueError:
            return None

    ix, iy = idx("x"), idx("y")
    iz, ivx, ivy, ivz = idx("z"), idx("vx"), idx("vy"), idx("vz")

    if ix is None or iy is None or ivx is None or ivy is None:
        if arr.shape[1] >= 5:
            x = arr[:, 1]
            y = arr[:, 2]
            z = np.zeros_like(x)
            vx = arr[:, 3]
            vy = arr[:, 4]
            vz = np.zeros_like(vx)
        else:
            raise ValueError(f"Unrecognized columns in {path}: {header}")
    else:
        x = arr[:, ix]
        y = arr[:, iy]
        z = arr[:, iz] if iz is not None and iz < arr.shape[1] else np.zeros_like(x)
        vx = arr[:, ivx]
        vy = arr[:, ivy]
        vz = arr[:, ivz] if ivz is not None and ivz < arr.shape[1] else np.zeros_like(vx)

    r = np.sqrt(x * x + y * y + z * z)
    v = np.sqrt(vx * vx + vy * vy + vz * vz)
    title_suffix = f"{path.name}{step_note}"
    return r, v, title_suffix


def rv_from_plane_arrays(positions: np.ndarray, velocities: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    3D r and speed from 2D (n,2) or 3D (n,3) positions/velocities (C++ snapshot arrays).
    """
    x = positions[:, 0]
    y = positions[:, 1]
    z = positions[:, 2] if positions.shape[1] >= 3 else np.zeros_like(x)
    vx = velocities[:, 0]
    vy = velocities[:, 1]
    vz = velocities[:, 2] if velocities.shape[1] >= 3 else np.zeros_like(vx)
    r = np.sqrt(x * x + y * y + z * z)
    v = np.sqrt(vx * vx + vy * vy + vz * vz)
    return r, v


def save_rotation_curve_png(
    r: np.ndarray,
    v: np.ndarray,
    output_path: Path,
    title_suffix: str,
    *,
    M_bh: float = NEWTONIAN_REFERENCE_M_BH,
    x_max: float | None = None,
    scatter_label: str = "Simulated stars",
    newtonian_label: str | None = None,
    filter_scatter_by_xmax: bool = True,
    provenance_label: str | None = None,
    run_info: dict[str, str | int | float] | None = None,
) -> None:
    """
    Scatter r vs v and overlay Newtonian v = sqrt(G*M_bh/r). Arrays need not be pre-sorted.

    The overlay uses only G and M_bh (no κ, no TPF profile). x_max=None picks a capped scale from data.
    When filter_scatter_by_xmax is True, only points with r <= x_max are scattered (keeps v-scale disk-like).
    """
    order = np.argsort(r)
    r_s = r[order]
    v_s = v[order]

    if len(r_s) == 0:
        raise ValueError("No radius data for rotation curve")

    if x_max is None:
        r_hi = float(np.max(r_s))
        x_max = max(r_hi * 1.2, 1.0)
        if r_hi > 1e15:
            x_max = max(x_max, min(DEFAULT_X_MAX, r_hi * 1.2))

    x_max = float(x_max)
    if filter_scatter_by_xmax:
        in_disk = r_s <= x_max
        r_plot = r_s[in_disk]
        v_plot = v_s[in_disk]
    else:
        r_plot = r_s
        v_plot = v_s

    r_lo = float(np.min(r_plot)) if len(r_plot) else float(np.min(r_s))
    r_hi_plot = float(np.max(r_plot)) if len(r_plot) else float(np.max(r_s))
    r_line_lo = max(r_lo * 0.9, r_hi_plot * 1e-9, 1e-12)
    if r_line_lo >= x_max:
        r_line_lo = max(x_max * 1e-6, 1e-12)
    r_line = np.linspace(r_line_lo, x_max, 2000)
    v_newton = np.sqrt(np.maximum(0.0, (G_SI * M_bh) / r_line))

    if newtonian_label is None:
        newtonian_label = f"Newtonian √(GM_bh/r), M_bh={M_bh:.6g} kg (κ-independent)"

    rf, vf, xl, yl = rotation_curve_display(x_max, run_info)
    r_plot_d = r_plot * rf
    v_plot_d = v_plot * vf
    r_line_d = r_line * rf
    v_newton_d = v_newton * vf
    x_max_d = x_max * rf

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(
        r_plot_d,
        v_plot_d,
        s=5,
        alpha=0.6,
        c="blue",
        edgecolors="none",
        label=scatter_label,
    )
    ax.plot(
        r_line_d,
        v_newton_d,
        "r--",
        linewidth=2,
        label=newtonian_label,
    )
    ax.set_xlim(0.0, x_max_d)  # disk-focused cap from rotation_curve_x_max (e.g. 2× galaxy_radius)
    # Auto-fit y-axis to the disk-dominated bulk so a few extreme escapers do not flatten the plot.
    if len(v_plot_d):
        y_candidates = np.concatenate([v_plot_d, v_newton_d])
    else:
        y_candidates = np.concatenate([v_s * vf, v_newton_d])
    # Robust cap: keep almost all points while suppressing pathological outliers.
    y_hi = float(np.percentile(y_candidates, 99.5))
    if not np.isfinite(y_hi) or y_hi <= 0.0:
        y_hi = float(np.max(y_candidates)) if len(y_candidates) else 1.0
    ax.set_ylim(0.0, y_hi * 1.08)
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    apply_suppress_tick_offset(ax)
    set_fitted_title(
        ax,
        f"Rotation curve — {title_suffix}",
        fontsize=12,
        min_fontsize=6,
    )
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    disp_note = "Axes: display units; snapshot CSV SI."
    foot_txt = (provenance_label + " | " + disp_note) if provenance_label else disp_note
    add_fitted_footer(
        fig,
        foot_txt,
        x=0.01,
        y=0.02,
        ha="left",
        va="bottom",
        fontsize=7,
        min_fontsize=5,
        color="gray",
    ).set_alpha(0.9)
    output_path = output_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot rotation curve from latest snapshot/output CSV.")
    parser.add_argument(
        "search_dir",
        type=Path,
        nargs="?",
        default=Path("."),
        help="Directory to glob for output_*.csv / snapshot_*.csv (default: current directory)",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path("rotation_curve.png"),
        help="Output image path (default: rotation_curve.png)",
    )
    parser.add_argument(
        "--M-bh",
        type=float,
        default=None,
        help="Black hole mass for Newtonian baseline (kg). Default: bh_mass from run_info.txt if present, else 1e41",
    )
    parser.add_argument(
        "--no-filter-scatter",
        action="store_true",
        help="Plot all stars in scatter even if r > x_max (y-axis may stretch from escapers)",
    )
    args = parser.parse_args()

    csv_path = find_latest_csv(args.search_dir)
    run_dir = csv_path.parent
    run_info = load_run_info(run_dir)
    r, v, title_suffix = load_particle_kinematics(csv_path)
    M_bh, m_source, used_fallback_mass = resolve_newtonian_overlay_mass(run_info, args.M_bh)
    x_lim = rotation_curve_x_max(run_info, r)
    lbl = scatter_label_from_run_info(run_info)
    save_rotation_curve_png(
        r,
        v,
        args.output,
        title_suffix,
        M_bh=M_bh,
        x_max=x_lim,
        scatter_label=lbl,
        filter_scatter_by_xmax=not args.no_filter_scatter,
        run_info=run_info,
    )
    print(f"Loaded: {csv_path}")
    print(f"Newtonian overlay mass: {M_bh:g} kg (source={m_source})")
    if used_fallback_mass:
        print(
            "Warning: no run-derived bh_mass found and no --M-bh was provided; "
            "using fallback benchmark mass."
        )
    print(f"Saved: {args.output.resolve()}")


if __name__ == "__main__":
    main()
