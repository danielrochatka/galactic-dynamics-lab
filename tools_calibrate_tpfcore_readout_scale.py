#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np


@dataclass
class Metrics:
    scale: float
    output_dir: Path
    final_pair_sep_delta_m: float
    rms_pair_sep_error_m: float
    rms_rel_speed_error_m_s: float
    angular_momentum_drift: float
    energy_drift: float
    score: float


def run_cmd(cmd: list[str], cwd: Path, log_path: Path) -> None:
    cmd_str = " ".join(cmd)
    print(f"$ {cmd_str}")
    with log_path.open("a", encoding="utf-8") as f:
        f.write(f"$ {cmd_str}\n")
    res = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True)
    with log_path.open("a", encoding="utf-8") as f:
        f.write(res.stdout)
        if res.stderr:
            f.write("\n[stderr]\n")
            f.write(res.stderr)
        f.write("\n")
    if res.returncode != 0:
        print(res.stdout)
        print(res.stderr, file=sys.stderr)
        raise RuntimeError(f"Command failed: {cmd_str}")


def load_diag_csv(path: Path) -> dict[str, np.ndarray]:
    arr = np.genfromtxt(path, delimiter=",", names=True)
    cols = {name: np.asarray(arr[name], dtype=np.float64) for name in arr.dtype.names}
    return cols


def ensure_diag_csv(run_dir: Path, repo_root: Path, log_path: Path) -> Path:
    csv_path = run_dir / "diagnostic_two_body_timeseries.csv"
    if csv_path.exists():
        return csv_path
    run_cmd([sys.executable, "plot_cpp_run.py", str(run_dir), "--no-animation"], cwd=repo_root, log_path=log_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing diagnostic CSV after postprocess: {csv_path}")
    return csv_path


def score_against_reference(tpf: dict[str, np.ndarray], ref: dict[str, np.ndarray]) -> tuple[float, float, float, float, float, float]:
    if len(tpf["time"]) != len(ref["time"]):
        raise ValueError("TPF and reference have different row counts")
    if not np.allclose(tpf["time"], ref["time"], atol=1e-12, rtol=0.0):
        raise ValueError("TPF and reference time grids differ")

    d_sep = tpf["pair_separation_m"] - ref["pair_separation_m"]
    d_vel = tpf["pair_relative_speed_m_s"] - ref["pair_relative_speed_m_s"]
    d_lz = tpf["relative_angular_momentum_z_kg_m2_s"] - ref["relative_angular_momentum_z_kg_m2_s"]
    d_e = tpf["newtonian_specific_energy_J_per_kg"] - ref["newtonian_specific_energy_J_per_kg"]

    mse_sep = float(np.mean(d_sep * d_sep))
    mse_vel = float(np.mean(d_vel * d_vel))
    mse_lz = float(np.mean(d_lz * d_lz))
    mse_e = float(np.mean(d_e * d_e))

    scale_sep = float(np.var(ref["pair_separation_m"])) + 1e-30
    scale_vel = float(np.var(ref["pair_relative_speed_m_s"])) + 1e-30
    scale_lz = float(np.mean(ref["relative_angular_momentum_z_kg_m2_s"] ** 2)) + 1e-30
    scale_e = float(np.mean(ref["newtonian_specific_energy_J_per_kg"] ** 2)) + 1e-30

    # Primary emphasis on separation and speed; Lz and energy are light sanity terms.
    w_sep, w_vel, w_lz, w_e = 0.7, 0.25, 0.03, 0.02
    score = (
        w_sep * (mse_sep / scale_sep)
        + w_vel * (mse_vel / scale_vel)
        + w_lz * (mse_lz / scale_lz)
        + w_e * (mse_e / scale_e)
    )

    final_sep_delta = float(d_sep[-1])
    rms_sep = math.sqrt(mse_sep)
    rms_vel = math.sqrt(mse_vel)
    ang_drift = float(tpf["relative_angular_momentum_z_kg_m2_s"][-1] - tpf["relative_angular_momentum_z_kg_m2_s"][0])
    energy_drift = float(tpf["newtonian_specific_energy_J_per_kg"][-1] - tpf["newtonian_specific_energy_J_per_kg"][0])
    return final_sep_delta, rms_sep, rms_vel, ang_drift, energy_drift, float(score)


def make_output_dir(base: Path, name: str) -> Path:
    out = base / name
    out.mkdir(parents=True, exist_ok=False)
    return out


def run_sim(cpp_sim_dir: Path, output_dir: Path, physics_package: str, log_path: Path, scale: float | None = None) -> None:
    cmd = [
        "./galaxy_sim",
        "earth_moon_benchmark",
        f"--output_dir={output_dir}",
        f"--physics_package={physics_package}",
        "--dt=0.01",
        "--validation_n_steps=500000",
        "--validation_snapshot_every=10",
        "--save_snapshots=true",
        "--save_run_info=true",
        "--tpf_vdsg_coupling=0",
        "--tpf_cooling_fraction=0",
    ]
    if physics_package == "TPFCore":
        if scale is None:
            raise ValueError("scale is required for TPFCore runs")
        cmd.extend(
            [
                "--tpf_dynamics_mode=legacy_readout",
                "--tpfcore_enable_provisional_readout=true",
                "--tpfcore_readout_mode=experimental_radial_r_scaling",
                f"--tpfcore_readout_scale={scale:.17e}",
            ]
        )
    run_cmd(cmd, cwd=cpp_sim_dir, log_path=log_path)


def format_table(rows: Iterable[Metrics]) -> str:
    hdr = (
        "scale,final_pair_sep_delta_m,rms_pair_sep_error_m,rms_relative_speed_error_m_s,"
        "angular_momentum_drift,energy_drift,total_score,output_dir"
    )
    lines = [hdr]
    for r in rows:
        lines.append(
            f"{r.scale:.17e},{r.final_pair_sep_delta_m:.9e},{r.rms_pair_sep_error_m:.9e},"
            f"{r.rms_rel_speed_error_m_s:.9e},{r.angular_momentum_drift:.9e},"
            f"{r.energy_drift:.9e},{r.score:.12e},{r.output_dir}"
        )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="Calibrate tpfcore_readout_scale against Newtonian earth_moon_benchmark.")
    parser.add_argument("--lo", type=float, default=3.3615e-11)
    parser.add_argument("--hi", type=float, default=3.3630e-11)
    parser.add_argument("--coarse-points", type=int, default=7)
    parser.add_argument("--refine-points", type=int, default=5)
    parser.add_argument("--max-refinements", type=int, default=3)
    parser.add_argument("--tol-width", type=float, default=2.0e-16)
    parser.add_argument("--improvement-eps", type=float, default=1.0e-7)
    parser.add_argument("--outputs-root", default="cpp_sim/outputs")
    parser.add_argument("--session-name", default=None)
    parser.add_argument("--reuse-existing-reference", action="store_true", help="Reuse existing Newtonian reference if present.")
    parser.add_argument("--stop-when-interior", action="store_true", help="Continue refinements until best point is interior to current refinement bracket.")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent
    cpp_sim_dir = repo_root / "cpp_sim"
    outputs_root = (repo_root / args.outputs_root).resolve()
    stamp = args.session_name or datetime.now(timezone.utc).strftime("calibration_%Y%m%dT%H%M%SZ")
    session_root = outputs_root / stamp
    session_root.mkdir(parents=True, exist_ok=True)
    log_path = session_root / "commands.log"

    ref_out = session_root / "newtonian_reference"
    ref_out.mkdir(parents=True, exist_ok=True)
    ref_csv = ref_out / "diagnostic_two_body_timeseries.csv"
    if args.reuse_existing_reference and ref_csv.exists():
        pass
    else:
        run_cmd(["make"], cwd=cpp_sim_dir, log_path=log_path)
        run_sim(cpp_sim_dir, ref_out, "Newtonian", log_path=log_path)
        ref_csv = ensure_diag_csv(ref_out, repo_root, log_path=log_path)
    ref = load_diag_csv(ref_csv)

    evaluated: dict[float, Metrics] = {}

    summary_csv = session_root / "calibration_summary.csv"
    if summary_csv.exists():
        with summary_csv.open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                key = float(row["scale"])
                evaluated[key] = Metrics(
                    scale=key,
                    output_dir=Path(row["output_dir"]),
                    final_pair_sep_delta_m=float(row["final_pair_sep_delta_m"]),
                    rms_pair_sep_error_m=float(row["rms_pair_sep_error_m"]),
                    rms_rel_speed_error_m_s=float(row["rms_relative_speed_error_m_s"]),
                    angular_momentum_drift=float(row["angular_momentum_drift"]),
                    energy_drift=float(row["energy_drift"]),
                    score=float(row["total_score"]),
                )

    def eval_scale(scale: float) -> Metrics:
        key = float(f"{scale:.17e}")
        if key in evaluated:
            return evaluated[key]
        tag = f"tpf_scale_{key:.17e}".replace("+", "")
        out = session_root / tag
        out.mkdir(parents=True, exist_ok=True)
        csv_path = out / "diagnostic_two_body_timeseries.csv"
        if not csv_path.exists():
            run_sim(cpp_sim_dir, out, "TPFCore", log_path=log_path, scale=key)
            csv_path = ensure_diag_csv(out, repo_root, log_path=log_path)
        tpf = load_diag_csv(csv_path)
        final_sep_delta, rms_sep, rms_vel, lz_drift, e_drift, score = score_against_reference(tpf, ref)
        m = Metrics(
            scale=key,
            output_dir=out,
            final_pair_sep_delta_m=final_sep_delta,
            rms_pair_sep_error_m=rms_sep,
            rms_rel_speed_error_m_s=rms_vel,
            angular_momentum_drift=lz_drift,
            energy_drift=e_drift,
            score=score,
        )
        evaluated[key] = m
        return m

    lo, hi = args.lo, args.hi
    points = np.linspace(lo, hi, args.coarse_points)
    for s in points:
        eval_scale(float(s))

    best_prev = min(evaluated.values(), key=lambda x: x.score).score

    for _ in range(args.max_refinements):
        ordered = sorted(evaluated.values(), key=lambda x: x.scale)
        best = min(ordered, key=lambda x: x.score)
        i_best = next(i for i, r in enumerate(ordered) if r.scale == best.scale)

        left_bound = ordered[max(0, i_best - 1)].scale
        right_bound = ordered[min(len(ordered) - 1, i_best + 1)].scale
        if left_bound == right_bound:
            break

        lo, hi = left_bound, right_bound
        if (hi - lo) <= args.tol_width:
            break

        bracket_points = np.linspace(lo, hi, args.refine_points)
        for s in bracket_points:
            eval_scale(float(s))

        best_now = min(evaluated.values(), key=lambda x: x.score).score
        ordered = sorted(evaluated.values(), key=lambda x: x.scale)
        best = min(ordered, key=lambda x: x.score)
        i_best = next(i for i, r in enumerate(ordered) if r.scale == best.scale)
        on_boundary = (best.scale == lo) or (best.scale == hi) or (i_best == 0) or (i_best == len(ordered) - 1)

        if args.stop_when_interior and not on_boundary:
            break
        if (best_prev - best_now < args.improvement_eps * max(best_prev, 1e-30)) and (not args.stop_when_interior or not on_boundary):
            break
        best_prev = best_now

    all_rows = sorted(evaluated.values(), key=lambda x: x.score)

    with summary_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "scale",
                "final_pair_sep_delta_m",
                "rms_pair_sep_error_m",
                "rms_relative_speed_error_m_s",
                "angular_momentum_drift",
                "energy_drift",
                "total_score",
                "output_dir",
            ]
        )
        for r in sorted(evaluated.values(), key=lambda x: x.scale):
            writer.writerow(
                [
                    f"{r.scale:.17e}",
                    f"{r.final_pair_sep_delta_m:.9e}",
                    f"{r.rms_pair_sep_error_m:.9e}",
                    f"{r.rms_rel_speed_error_m_s:.9e}",
                    f"{r.angular_momentum_drift:.9e}",
                    f"{r.energy_drift:.9e}",
                    f"{r.score:.12e}",
                    str(r.output_dir),
                ]
            )

    report_txt = session_root / "calibration_report.txt"
    with report_txt.open("w", encoding="utf-8") as f:
        f.write(f"session_root: {session_root}\n")
        f.write(f"reference_csv: {ref_csv}\n")
        f.write(f"search_interval_initial: [{args.lo:.17e}, {args.hi:.17e}]\n")
        f.write(f"tested_scales: {len(evaluated)}\n\n")
        f.write(format_table(sorted(evaluated.values(), key=lambda x: x.scale)) + "\n\n")
        f.write("best_by_score:\n")
        for i, r in enumerate(all_rows[:2], start=1):
            f.write(f"  {i}. scale={r.scale:.17e}, score={r.score:.12e}, output={r.output_dir}\n")

    print("\n=== Calibration complete ===")
    print(f"Session root: {session_root}")
    print(f"Reference CSV: {ref_csv}")
    print(format_table(sorted(evaluated.values(), key=lambda x: x.scale)))
    print("Top 2 by score:")
    for i, r in enumerate(all_rows[:2], start=1):
        print(f"  {i}. scale={r.scale:.17e}, score={r.score:.12e}, output={r.output_dir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
