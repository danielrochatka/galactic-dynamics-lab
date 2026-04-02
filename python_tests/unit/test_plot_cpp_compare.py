from __future__ import annotations

import json
import sys
import tempfile
import unittest
import importlib.util
from pathlib import Path

# Repo root (parent of python_tests/)
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from plot_cpp_compare import (
    _resolve_side_run_dir,
    calculate_compare_smart_viewport,
    matched_steps_strict,
    render_compare,
)


def _write_snapshot(path: Path, step: int, t: float, x: float) -> None:
    path.write_text(
        "\n".join(
            [
                f"# step,{step},time,{t:.6e}",
                "i,x,y,vx,vy,mass",
                f"0,{x},0.0,0.0,1.0,1.0",
            ]
        ),
        encoding="utf-8",
    )


def _write_run_info(run_dir: Path, pkg: str) -> None:
    (run_dir / "run_info.txt").write_text(
        "\n".join(
            [
                f"physics_package\t{pkg}",
                "render_overlay_mode\tnone",
                "code_version_label\ttest@abc1234",
                "n_steps\t10",
                "dt\t0.1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


class TestPlotCppCompare(unittest.TestCase):
    def test_resolve_side_run_dir_cpp_sim_relative_manifest(self) -> None:
        """Manifest stores paths relative to cpp_sim; compare_parent is .../outputs/RUN."""
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            run = root / "cpp_sim" / "outputs" / "RUN1"
            (run / "left_X").mkdir(parents=True)
            m = "outputs/RUN1/left_X"
            got = _resolve_side_run_dir(m, run.resolve())
            self.assertEqual(got, (run / "left_X").resolve())

    def test_matched_steps_strict_pass(self) -> None:
        self.assertEqual(matched_steps_strict({0, 5}, {0, 5}), [0, 5])

    def test_matched_steps_strict_fail(self) -> None:
        with self.assertRaises(ValueError):
            matched_steps_strict({0, 5}, {0, 10})

    @unittest.skipUnless(
        importlib.util.find_spec("matplotlib") is not None and importlib.util.find_spec("numpy") is not None,
        "matplotlib/numpy not installed",
    )
    def test_render_compare_writes_static_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            parent = Path(td)
            left = parent / "left_TPFCore"
            right = parent / "right_Newtonian"
            left.mkdir()
            right.mkdir()
            _write_run_info(left, "TPFCore")
            _write_run_info(right, "Newtonian")
            _write_snapshot(left / "snapshot_00000.csv", 0, 0.0, 1.0)
            _write_snapshot(left / "snapshot_00010.csv", 10, 1.0, 2.0)
            _write_snapshot(right / "snapshot_00000.csv", 0, 0.0, 1.2)
            _write_snapshot(right / "snapshot_00010.csv", 10, 1.0, 2.2)
            (parent / "compare_manifest.json").write_text(
                json.dumps(
                    {
                        "compare_run_id": "r1",
                        "left_dir": str(left),
                        "right_dir": str(right),
                    }
                ),
                encoding="utf-8",
            )

            render_compare(parent, no_animation=True, overlay_mode="none")
            self.assertTrue((parent / "galaxy_initial_compare.png").exists())
            self.assertTrue((parent / "galaxy_final_compare.png").exists())

    def test_compare_animation_does_not_use_ema_radius_smoothing(self) -> None:
        """EMA on shared radius lagged behind per-frame max(r_l,r_r) and clipped stars in video."""
        import plot_cpp_compare

        text = Path(plot_cpp_compare.__file__).read_text(encoding="utf-8")
        self.assertNotIn("_smooth_alpha", text)
        self.assertNotIn("_smooth_state", text)

    def test_compare_smart_viewport_contains_both_extents(self) -> None:
        """Shared viewport is computed from the union of both sides' point clouds (+ origin)."""
        import numpy as np

        from plot_cpp_run import Snapshot

        pos_l = np.random.default_rng(0).normal(scale=1e19, size=(100, 2))
        vel_l = np.zeros_like(pos_l)
        th = np.linspace(0, 2 * np.pi, 100, endpoint=False)
        pos_r = np.column_stack([1e27 * np.cos(th), 1e27 * np.sin(th)])
        vel_r = np.column_stack([np.cos(th) * 1e5, np.sin(th) * 1e5])
        snap_l = Snapshot(0, 0.0, pos_l, vel_l)
        snap_r = Snapshot(0, 0.0, pos_r, vel_r)
        vp = calculate_compare_smart_viewport([snap_l], [snap_r], 150.0)
        self.assertGreater(vp.half_axis, 1e26)


if __name__ == "__main__":
    unittest.main()

