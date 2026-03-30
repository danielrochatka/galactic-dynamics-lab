"""
Audit overlays for galaxy scatter renders (plot_cpp_run / render.py).

Branch labels mirror cpp_sim/render_audit.cpp logic for consistency when run_info
omits explicit keys (older runs).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Optional

# matplotlib is imported inside draw_galaxy_render_overlay so tests can import branch helpers
# without requiring matplotlib at import time.

_DERIVED_READOUT_MODES = frozenset({"tr_coherence_readout", "derived_tpf_radial_readout"})


def load_render_manifest(run_dir: Path) -> Optional[dict[str, Any]]:
    p = run_dir / "render_manifest.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return None


def _f(ri: dict[str, Any], key: str, default: float = 0.0) -> float:
    v = ri.get(key, default)
    if v is None:
        return default
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def _i(ri: dict[str, Any], key: str, default: int = 0) -> int:
    v = ri.get(key, default)
    if v is None:
        return default
    try:
        return int(v)
    except (TypeError, ValueError):
        return default


def _pick_f(mf: dict[str, Any], ri: dict[str, Any], key: str, default: float = 0.0) -> float:
    if key in mf and mf[key] is not None:
        try:
            return float(mf[key])
        except (TypeError, ValueError):
            pass
    return _f(ri, key, default)


def _pick_i(mf: dict[str, Any], ri: dict[str, Any], key: str, default: int = 0) -> int:
    if key in mf and mf[key] is not None:
        try:
            return int(mf[key])
        except (TypeError, ValueError):
            pass
    return _i(ri, key, default)


def _pick_s(mf: dict[str, Any], ri: dict[str, Any], key: str, default: str = "") -> str:
    if key in mf and mf[key] is not None:
        return str(mf[key])
    v = ri.get(key, default)
    return default if v is None else str(v)


def _pick_boolish(mf: dict[str, Any], ri: dict[str, Any], key: str, default: bool = False) -> bool:
    """Merge manifest + run_info for 0/1, true/false, or string flags."""
    for src in (mf, ri):
        if key not in src:
            continue
        v = src[key]
        if isinstance(v, bool):
            return v
        if isinstance(v, int):
            return v != 0
        s = str(v).strip().lower()
        if s in ("1", "true", "yes", "on"):
            return True
        if s in ("0", "false", "no", "off"):
            return False
    return default


def infer_branches_from_run_info(ri: dict[str, Any]) -> tuple[str, str, str]:
    """Return (dynamics, metrics, acceleration_code_path) when explicit keys are missing."""
    pp = str(ri.get("physics_package", "Newtonian"))
    vdsg = _f(ri, "tpf_vdsg_coupling", 0.0)
    prov = _i(ri, "tpfcore_enable_provisional_readout", 0)
    rmode = str(ri.get("tpfcore_readout_mode", ""))
    if pp == "Newtonian":
        return (
            "Newtonian_pairwise_G_SI",
            "none",
            "NewtonianPackage::compute_accelerations",
        )
    if pp != "TPFCore":
        return (f"{pp} (non-TPFCore)", "unknown", "unknown_package")
    if prov == 0:
        return (
            "TPFCore_dynamics_DISABLED (provisional_readout off)",
            "TPFCore_metrics_n/a (provisional_readout off)",
            "TPFCorePackage::compute_accelerations (throws without provisional readout)",
        )
    dyn = f"TPF_readout_acceleration:{rmode}"
    met = f"tpfcore_readout:{rmode}"
    if rmode in _DERIVED_READOUT_MODES:
        base = "TPFCorePackage::compute_provisional_readout_acceleration + derived_tpf_radial_profile"
    else:
        base = f"TPFCorePackage::compute_provisional_readout_acceleration ({rmode})"
    shunt_on = bool(_pick_boolish({}, ri, "tpf_global_accel_shunt_enable", False))
    if shunt_on:
        acc = (
            base
            + " + accumulate_vdsg_velocity_modifier + apply_global_accel_magnitude_shunt (when tpf_global_accel_shunt_enable)"
        )
    else:
        acc = (
            base
            + " + accumulate_vdsg_velocity_modifier (global |a| shunt OFF — clean readout+VDSG path without velocity cap)"
        )
    return dyn, met, acc


def build_overlay_spec(
    run_dir: Path,
    run_info: dict[str, Any],
) -> dict[str, Any]:
    """Merge run_info + optional render_manifest.json into one dict for overlays."""
    mf = load_render_manifest(run_dir) or {}
    ri = dict(run_info)
    dyn = str(mf.get("active_dynamics_branch") or ri.get("active_dynamics_branch") or "")
    met = str(mf.get("active_metrics_branch") or ri.get("active_metrics_branch") or "")
    acc = str(mf.get("acceleration_code_path") or ri.get("acceleration_code_path") or "")
    if not dyn or not met or not acc:
        idyn, imet, iacc = infer_branches_from_run_info(ri)
        dyn = dyn or idyn
        met = met or imet
        acc = acc or iacc

    run_id = str(mf.get("run_id") or ri.get("run_id") or run_dir.name)
    n_steps = _pick_i(mf, ri, "n_steps", _i(ri, "n_steps", 0))
    spec: dict[str, Any] = {
        "run_id": run_id,
        "run_dir": str(run_dir),
        "active_dynamics_branch": dyn,
        "active_metrics_branch": met,
        "acceleration_code_path": acc,
        "physics_package": _pick_s(mf, ri, "physics_package", ""),
        "simulation_mode": _pick_s(mf, ri, "simulation_mode", ""),
        "tpf_vdsg_coupling": _pick_f(mf, ri, "tpf_vdsg_coupling", 0.0),
        "tpf_kappa": _pick_f(mf, ri, "tpf_kappa", 0.0),
        "tpf_cooling_fraction": _pick_f(mf, ri, "tpf_cooling_fraction", 0.0),
        "tpfcore_readout_mode": _pick_s(mf, ri, "tpfcore_readout_mode", ""),
        "galaxy_init_template": _pick_s(mf, ri, "galaxy_init_template", ""),
        "galaxy_init_seed": mf.get("galaxy_init_seed", ri.get("galaxy_init_seed", "")),
        "n_stars": _pick_i(mf, ri, "n_stars", _i(ri, "n_stars", 0)),
        "bh_mass": _pick_f(mf, ri, "bh_mass", 0.0),
        "softening": _pick_f(mf, ri, "softening", 0.0),
        "galaxy_radius": _pick_f(mf, ri, "galaxy_radius", 0.0),
        "enable_star_star_gravity": bool(_pick_i(mf, ri, "enable_star_star_gravity", 1)),
        "n_steps_total": n_steps,
        "galaxy_init_master_chaos": _pick_f(mf, ri, "galaxy_init_master_chaos", 1.0),
        "galaxy_init_position_noise": _pick_f(mf, ri, "galaxy_init_position_noise", 0.0),
        "galaxy_init_velocity_angle_noise": _pick_f(
            mf, ri, "galaxy_init_velocity_angle_noise", 0.0
        ),
        "galaxy_init_velocity_magnitude_noise": _pick_f(
            mf, ri, "galaxy_init_velocity_magnitude_noise", 0.0
        ),
        "galaxy_init_clumpiness": _pick_f(mf, ri, "galaxy_init_clumpiness", 0.0),
        "galaxy_init_num_clumps": _pick_i(mf, ri, "galaxy_init_num_clumps", 0),
        "galaxy_init_clump_radius_fraction": _pick_f(
            mf, ri, "galaxy_init_clump_radius_fraction", 0.0
        ),
        "galaxy_init_m2_amplitude": _pick_f(mf, ri, "galaxy_init_m2_amplitude", 0.0),
        "galaxy_init_m3_amplitude": _pick_f(mf, ri, "galaxy_init_m3_amplitude", 0.0),
        "galaxy_init_bar_amplitude": _pick_f(mf, ri, "galaxy_init_bar_amplitude", 0.0),
        "galaxy_init_bar_axis_ratio": _pick_f(mf, ri, "galaxy_init_bar_axis_ratio", 1.0),
        "galaxy_init_spiral_amplitude": _pick_f(mf, ri, "galaxy_init_spiral_amplitude", 0.0),
        "galaxy_init_spiral_winding": _pick_f(mf, ri, "galaxy_init_spiral_winding", 0.0),
        "galaxy_init_spiral_phase": _pick_f(mf, ri, "galaxy_init_spiral_phase", 0.0),
        "git_commit_full": _pick_s(mf, ri, "git_commit_full", "unknown"),
        "git_commit_short": _pick_s(mf, ri, "git_commit_short", "unknown"),
        "git_branch": _pick_s(mf, ri, "git_branch", ""),
        "git_tag": _pick_s(mf, ri, "git_tag", ""),
        "git_dirty": _pick_boolish(mf, ri, "git_dirty", False),
        "code_version_label": _pick_s(mf, ri, "code_version_label", "unknown"),
    }
    return spec


def resolve_overlay_mode(run_info: dict[str, Any], cli_override: Optional[str]) -> str:
    if cli_override is not None:
        return cli_override.strip().lower()
    raw = run_info.get("render_overlay_mode")
    if raw is None:
        return "none"
    return str(raw).strip().lower()


def provenance_overlay_watermark_text(mode: str, spec: dict[str, Any]) -> str:
    """
    Lower-left watermark for minimal / audit_full (empty for none or missing data).
    Does not duplicate the upper-left audit box; this is rev / branch / tag only.
    """
    if mode == "none":
        return ""
    lbl = str(spec.get("code_version_label") or "").strip()
    short = str(spec.get("git_commit_short") or "").strip()
    dirty_b = bool(spec.get("git_dirty"))
    br = str(spec.get("git_branch") or "").strip()
    tg = str(spec.get("git_tag") or "").strip()

    if mode == "minimal":
        if lbl and lbl != "unknown":
            return lbl
        if short and short != "unknown":
            s = f"rev: {short}"
            if dirty_b:
                s += " (dirty)"
            return s
        return "rev: unknown"

    if mode == "audit_full":
        lines: list[str] = []
        lines.append(f"rev: {short if short else 'unknown'}")
        if br:
            lines.append(f"branch: {br}")
        if tg:
            lines.append(f"tag: {tg}")
        lines.append(f"dirty: {'yes' if dirty_b else 'no'}")
        return "\n".join(lines)

    return ""


def _cooling_label(ri: dict[str, Any], spec: dict[str, Any]) -> str:
    frac = float(spec.get("tpf_cooling_fraction") or 0.0)
    pp = str(spec.get("physics_package", ""))
    sm = ri.get("simulation_mode", "")
    sm_s = str(sm).lower()
    is_galaxy = sm == 0 or sm_s == "0" or sm_s == "galaxy"
    if pp != "TPFCore" or not is_galaxy or frac <= 0:
        return "off"
    return f"{frac:g} (active first fraction of steps)"


def draw_galaxy_render_overlay(
    ax: Any,
    mode: str,
    spec: dict[str, Any],
    *,
    run_info: dict[str, Any],
    step: int,
    time_s: float,
) -> None:
    if mode == "none":
        return

    vdsg = float(spec.get("tpf_vdsg_coupling") or 0.0)
    coup = f"{vdsg:.6g}" if spec.get("physics_package") == "TPFCore" else "N/A"
    ntot = int(spec.get("n_steps_total") or 0)
    step_str = f"{step} / {ntot}" if ntot > 0 else str(step)

    lines_min = [
        f"Run ID: {spec.get('run_id', '?')}",
        f"Dynamics: {spec.get('active_dynamics_branch', '?')}",
        f"Metrics: {spec.get('active_metrics_branch', '?')}",
        f"Coupling (VDSG λ): {coup}",
        f"Cooling: {_cooling_label(run_info, spec)}",
        f"Template: {spec.get('galaxy_init_template', '?')}",
        f"Seed: {spec.get('galaxy_init_seed', '?')}",
        f"Stars: {spec.get('n_stars', '?')}",
        f"BH mass: {float(spec.get('bh_mass') or 0):.6g}",
        f"Step: {step_str}",
        f"t = {time_s:.6g}",
    ]

    if mode == "minimal":
        text = "\n".join(lines_min[:10])
        fontsize = 7
    elif mode == "audit_full":
        extra = [
            "",
            f"Physics: {spec.get('physics_package', '')}",
            f"Sim mode: {spec.get('simulation_mode', '')}",
            f"Readout mode: {spec.get('tpfcore_readout_mode', '')}",
            f"κ: {float(spec.get('tpf_kappa') or 0):.6g}",
            f"Softening: {float(spec.get('softening') or 0):.6g}",
            f"R_gal: {float(spec.get('galaxy_radius') or 0):.6g}",
            f"Star–star G: {spec.get('enable_star_star_gravity', '')}",
            f"Accel path: {spec.get('acceleration_code_path', '')}",
            f"master_chaos: {float(spec.get('galaxy_init_master_chaos') or 0):.6g}",
            f"pos_noise: {float(spec.get('galaxy_init_position_noise') or 0):.6g}",
            f"v_angle_noise: {float(spec.get('galaxy_init_velocity_angle_noise') or 0):.6g}",
            f"v_mag_noise: {float(spec.get('galaxy_init_velocity_magnitude_noise') or 0):.6g}",
            f"clumpiness: {float(spec.get('galaxy_init_clumpiness') or 0):.6g}  n_clumps: {spec.get('galaxy_init_num_clumps', '')}",
        ]
        amps = []
        for k, lab in (
            ("galaxy_init_m2_amplitude", "m2"),
            ("galaxy_init_m3_amplitude", "m3"),
            ("galaxy_init_bar_amplitude", "bar"),
            ("galaxy_init_spiral_amplitude", "spiral"),
        ):
            v = float(spec.get(k) or 0.0)
            if v != 0.0:
                amps.append(f"{lab}={v:.6g}")
        if amps:
            extra.append("modes: " + " ".join(amps))
        text = "\n".join(lines_min + extra)
        fontsize = 6
    else:
        return

    ax.text(
        0.02,
        0.98,
        text,
        transform=ax.transAxes,
        fontsize=fontsize,
        family="monospace",
        verticalalignment="top",
        horizontalalignment="left",
        color="white",
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "black",
            "edgecolor": "gray",
            "alpha": 0.55,
        },
        zorder=20,
    )

    wm = provenance_overlay_watermark_text(mode, spec)
    if wm:
        wm_fs = 5 if mode == "audit_full" else 6
        ax.text(
            0.01,
            0.01,
            wm,
            transform=ax.transAxes,
            fontsize=wm_fs,
            family="monospace",
            verticalalignment="bottom",
            horizontalalignment="left",
            color="white",
            bbox={
                "boxstyle": "round,pad=0.28",
                "facecolor": "black",
                "edgecolor": "gray",
                "alpha": 0.42,
            },
            zorder=21,
        )
