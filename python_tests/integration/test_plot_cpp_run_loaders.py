"""Integration-style tests for plot_cpp_run helpers (no matplotlib display)."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def _have_plot_stack() -> bool:
    return importlib.util.find_spec("numpy") is not None and importlib.util.find_spec(
        "matplotlib"
    ) is not None


@unittest.skipUnless(_have_plot_stack(), "requires numpy and matplotlib")
class TestPlotCppRunLoaders(unittest.TestCase):
    def test_load_snapshot_csv_roundtrip(self) -> None:
        import plot_cpp_run as pcr
        content = """# step,0,time,0.0e+00
i,x,y,vx,vy,mass
0,1.0,2.0,0.1,0.2,0.05
1,-3.0,4.0,-0.1,0.3,0.05
"""
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "snapshot_00000.csv"
            p.write_text(content, encoding="utf-8")
            snap = pcr.load_snapshot_csv(p)
            self.assertIsNotNone(snap)
            assert snap is not None
            self.assertEqual(snap.step, 0)
            self.assertEqual(snap.positions.shape, (2, 2))


if __name__ == "__main__":
    unittest.main()
