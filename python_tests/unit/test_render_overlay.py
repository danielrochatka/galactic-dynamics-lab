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
        self.assertEqual(d, "TPF_readout_acceleration:derived_tpf_radial_readout")
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
        self.assertIn("TPF_readout_acceleration", d)
        self.assertIn("tensor_radial_projection", m)

    def test_load_render_manifest_missing(self) -> None:
        self.assertIsNone(load_render_manifest(Path("/nonexistent/does/not/exist")))

    def test_build_overlay_spec_merges_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            mf = {
                "active_dynamics_branch": "TPF_readout_acceleration:derived_tpf_radial_readout",
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
            self.assertEqual(spec["active_dynamics_branch"], "TPF_readout_acceleration:derived_tpf_radial_readout")
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


if __name__ == "__main__":
    unittest.main()
