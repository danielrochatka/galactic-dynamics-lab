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

from plot_cpp_compare import matched_steps_strict, render_compare


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


if __name__ == "__main__":
    unittest.main()

