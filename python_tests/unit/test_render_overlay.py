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
        self.assertEqual(d, "VDSG_centripetal_SI")
        self.assertIn("derived_tpf_radial", m)
        self.assertIn("centripetal", acc)

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
                "active_dynamics_branch": "VDSG_centripetal_SI",
                "active_metrics_branch": "tpfcore_readout:derived_tpf_radial_readout",
                "acceleration_code_path": "TPFCorePackage::accumulate_velocity_deformed_centripetal_gravity",
                "run_id": "test_run",
                "physics_package": "TPFCore",
                "tpf_vdsg_coupling": 1e-18,
            }
            (run_dir / "render_manifest.json").write_text(
                json.dumps(mf), encoding="utf-8"
            )
            ri = {"physics_package": "TPFCore", "simulation_mode": "galaxy"}
            spec = build_overlay_spec(run_dir, ri)
            self.assertEqual(spec["active_dynamics_branch"], "VDSG_centripetal_SI")
            self.assertEqual(spec["run_id"], "test_run")


if __name__ == "__main__":
    unittest.main()
