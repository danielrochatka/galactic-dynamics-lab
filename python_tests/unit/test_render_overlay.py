"""Unit tests for render_overlay helpers (manifest load, branch inference, overlay spec)."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

# Repo root (parent of python_tests/)
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from render_overlay import (  # noqa: E402
    build_overlay_spec,
    draw_galaxy_render_overlay,
    infer_branches_from_run_info,
    load_render_manifest,
    provenance_overlay_watermark_text,
    resolve_overlay_mode,
)


class TestRenderOverlay(unittest.TestCase):
    def test_resolve_overlay_mode_default_none(self) -> None:
        self.assertEqual(resolve_overlay_mode({}, None), "none")

    def test_resolve_overlay_mode_cli_wins(self) -> None:
        self.assertEqual(
            resolve_overlay_mode({"render_overlay_mode": "minimal"}, "audit_full"),
            "audit_full",
        )

    def test_infer_branches_newtonian(self) -> None:
        d, m, a = infer_branches_from_run_info({"physics_package": "Newtonian"})
        self.assertEqual(d, "Newtonian_pairwise_G_SI")
        self.assertEqual(m, "none")
        self.assertIn("NewtonianPackage", a)

    def test_infer_branches_tpfcore_vdsg(self) -> None:
        ri = {
            "physics_package": "TPFCore",
            "tpfcore_enable_provisional_readout": 1,
            "tpf_vdsg_coupling": 1e-20,
            "tpfcore_readout_mode": "derived_tpf_radial_readout",
        }
        d, m, acc = infer_branches_from_run_info(ri)
        self.assertEqual(d, "TPF_legacy_readout_plus_VDSG:derived_tpf_radial_readout")
        self.assertIn("derived_tpf_radial", m)
        self.assertIn("accumulate_vdsg_velocity_modifier", acc)
        self.assertIn("shunt OFF", acc)

    def test_infer_branches_fills_missing_explicit_keys(self) -> None:
        """Older run_info without active_* keys still yields consistent labels."""
        ri = {
            "physics_package": "TPFCore",
            "tpfcore_enable_provisional_readout": 1,
            "tpf_vdsg_coupling": 0.0,
            "tpfcore_readout_mode": "tensor_radial_projection",
        }
        d, m, acc = infer_branches_from_run_info(ri)
        self.assertIn("TPF_legacy_readout:", d)
        self.assertIn("tensor_radial_projection", m)

    def test_load_render_manifest_missing(self) -> None:
        self.assertIsNone(load_render_manifest(Path("/nonexistent/does/not/exist")))

    def test_build_overlay_spec_merges_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            mf = {
                "active_dynamics_branch": "TPF_legacy_readout_plus_VDSG:derived_tpf_radial_readout",
                "active_metrics_branch": "tpfcore_readout:derived_tpf_radial_readout",
                "acceleration_code_path": "TPFCorePackage::compute_provisional_readout_acceleration + "
                "derived_tpf_radial_profile + accumulate_vdsg_velocity_modifier + "
                "apply_global_accel_magnitude_shunt",
                "run_id": "test_run",
                "physics_package": "TPFCore",
                "tpf_vdsg_coupling": 1e-18,
            }
            (run_dir / "render_manifest.json").write_text(
                json.dumps(mf), encoding="utf-8"
            )
            ri = {"physics_package": "TPFCore", "simulation_mode": "galaxy"}
            spec = build_overlay_spec(run_dir, ri)
            self.assertEqual(
                spec["active_dynamics_branch"],
                "TPF_legacy_readout_plus_VDSG:derived_tpf_radial_readout",
            )
            self.assertEqual(spec["run_id"], "test_run")

    def test_build_overlay_spec_git_defaults(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            ri: dict = {"physics_package": "Newtonian"}
            spec = build_overlay_spec(run_dir, ri)
            self.assertEqual(spec["git_commit_full"], "unknown")
            self.assertEqual(spec["git_commit_short"], "unknown")
            self.assertEqual(spec["code_version_label"], "unknown")
            self.assertFalse(spec["git_dirty"])

    def test_build_overlay_spec_git_from_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            mf = {
                "git_commit_full": "abcd" * 10,
                "git_commit_short": "deadbeef",
                "git_branch": "main",
                "git_tag": "",
                "git_dirty": True,
                "code_version_label": "main@deadbeef-dirty",
            }
            (run_dir / "render_manifest.json").write_text(
                json.dumps(mf), encoding="utf-8"
            )
            spec = build_overlay_spec(run_dir, {"physics_package": "Newtonian"})
            self.assertEqual(spec["git_commit_short"], "deadbeef")
            self.assertTrue(spec["git_dirty"])
            self.assertEqual(spec["code_version_label"], "main@deadbeef-dirty")

    def test_provenance_watermark_text(self) -> None:
        spec = {
            "code_version_label": "main@abc1234",
            "git_commit_short": "abc1234",
            "git_dirty": False,
            "git_branch": "main",
            "git_tag": "",
        }
        self.assertEqual(
            provenance_overlay_watermark_text("minimal", spec), "main@abc1234"
        )
        spec["code_version_label"] = "unknown"
        self.assertEqual(
            provenance_overlay_watermark_text("minimal", spec), "rev: abc1234"
        )
        spec["git_dirty"] = True
        self.assertIn("dirty", provenance_overlay_watermark_text("minimal", spec))
        full = provenance_overlay_watermark_text("audit_full", spec)
        self.assertIn("branch: main", full)
        self.assertIn("dirty: yes", full)

    def test_minimal_overlay_includes_display_units_when_enabled(self) -> None:
        class _FakeAx:
            def __init__(self) -> None:
                self.transAxes = object()
                self.text_calls: list[tuple[float, float, str]] = []

            def text(self, x, y, text, **kwargs):  # type: ignore[no-untyped-def]
                _ = kwargs
                self.text_calls.append((x, y, text))

        ax = _FakeAx()
        spec = {
            "run_id": "r1",
            "active_dynamics_branch": "dyn",
            "active_metrics_branch": "met",
            "physics_package": "Newtonian",
            "tpf_vdsg_coupling": 0.0,
            "galaxy_init_template": "disk",
            "galaxy_init_seed": 7,
            "n_stars": 128,
            "bh_mass": 1.0,
            "n_steps_total": 100,
            "display_units_in_overlay": True,
            "active_display_distance_unit": "kpc",
            "active_display_time_unit": "Myr",
            "active_display_velocity_unit": "km/s",
            "code_version_label": "unknown",
            "git_commit_short": "abc1234",
            "git_dirty": False,
            "git_branch": "",
            "git_tag": "",
        }
        draw_galaxy_render_overlay(
            ax,
            "minimal",
            spec,
            run_info={"simulation_mode": "galaxy"},
            step=10,
            time_s=1.0,
        )
        overlay_text = ax.text_calls[0][2]
        self.assertIn("Display distance: kpc", overlay_text)
        self.assertIn("Display time: Myr", overlay_text)
        self.assertIn("Display velocity: km/s", overlay_text)


if __name__ == "__main__":
    unittest.main()
