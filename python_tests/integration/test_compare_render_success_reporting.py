"""Integration test for compare render success/failure reporting invariants."""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def _have_plot_stack() -> bool:
    return importlib.util.find_spec("numpy") is not None and importlib.util.find_spec("matplotlib") is not None


def _write_snapshot(path: Path, step: int, time_s: float, radius: float) -> None:
    path.write_text(
        "\n".join(
            [
                f"# step,{step},time,{time_s}",
                "i,x,y,vx,vy,mass",
                f"0,{radius},0.0,0.0,0.0,1.0",
                "1,0.0,1.0,0.0,0.0,1.0",
            ]
        ),
        encoding="utf-8",
    )


@unittest.skipUnless(_have_plot_stack(), "requires numpy and matplotlib")
class TestCompareRenderSuccessReporting(unittest.TestCase):
    def test_render_exits_nonzero_when_required_pngs_missing(self) -> None:
        import plot_cpp_compare as pcc

        with tempfile.TemporaryDirectory() as td:
            parent = Path(td)
            left = parent / "left"
            right = parent / "right"
            left.mkdir()
            right.mkdir()
            (left / "run_info.json").write_text(json.dumps({"simulation_mode": "galaxy", "physics_package": "TPFCore"}), encoding="utf-8")
            (right / "run_info.json").write_text(json.dumps({"simulation_mode": "galaxy", "physics_package": "Newtonian"}), encoding="utf-8")
            _write_snapshot(left / "snapshot_00000.csv", 0, 0.0, 1.0)
            _write_snapshot(left / "snapshot_00001.csv", 1, 1.0, 2.0)
            _write_snapshot(right / "snapshot_00000.csv", 0, 0.0, 1.1)
            _write_snapshot(right / "snapshot_00001.csv", 1, 1.0, 2.1)
            (parent / "compare_manifest.json").write_text(
                json.dumps({"compare_run_id": "r_missing", "left_dir": str(left), "right_dir": str(right)}),
                encoding="utf-8",
            )

            original_copy2 = pcc.shutil.copy2
            pcc.shutil.copy2 = lambda *_args, **_kwargs: None
            try:
                with self.assertRaises(SystemExit):
                    pcc.render_compare(parent, no_animation=True, overlay_mode="none")
            finally:
                pcc.shutil.copy2 = original_copy2


if __name__ == "__main__":
    unittest.main()
