"""Integration tests for plot_cpp_run burn-in filtering and smart bounds."""

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


def _write_snapshot(path: Path, step: int, rows: list[tuple[float, float]]) -> None:
    lines = [f"# step,{step},time,{float(step):.1f}", "i,x,y,vx,vy,mass"]
    for i, (x, y) in enumerate(rows):
        lines.append(f"{i},{x},{y},0.0,0.0,1.0")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@unittest.skipUnless(_have_plot_stack(), "requires numpy and matplotlib")
class TestPlotCppRunBurninFilter(unittest.TestCase):
    def test_no_skip_keeps_all(self) -> None:
        import plot_cpp_run as pcr

        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            _write_snapshot(d / "snapshot_00000.csv", 0, [(1.0, 0.0)])
            _write_snapshot(d / "snapshot_00001.csv", 1, [(2.0, 0.0)])
            rec = pcr.load_all_snapshot_records(d)
            out = pcr.filter_snapshots_for_plotting(rec, 0, 0)
            self.assertEqual(len(out), 2)
            self.assertEqual(out[0][1].step, 0)
            self.assertEqual(out[-1][1].step, 1)

    def test_skip_initial_steps(self) -> None:
        import plot_cpp_run as pcr

        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            _write_snapshot(d / "snapshot_00000.csv", 0, [(1.0, 0.0)])
            _write_snapshot(d / "snapshot_00001.csv", 10, [(2.0, 0.0)])
            _write_snapshot(d / "snapshot_00002.csv", 20, [(3.0, 0.0)])
            rec = pcr.load_all_snapshot_records(d)
            out = pcr.filter_snapshots_for_plotting(rec, skip_initial_steps=10)
            self.assertEqual([s.step for _, s in out], [10, 20])

    def test_skip_initial_snapshots(self) -> None:
        import plot_cpp_run as pcr

        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            _write_snapshot(d / "snapshot_00000.csv", 0, [(1.0, 0.0)])
            _write_snapshot(d / "snapshot_00001.csv", 1, [(2.0, 0.0)])
            _write_snapshot(d / "snapshot_00002.csv", 2, [(3.0, 0.0)])
            rec = pcr.load_all_snapshot_records(d)
            out = pcr.filter_snapshots_for_plotting(rec, skip_initial_snapshots=2)
            self.assertEqual([s.step for _, s in out], [2])

    def test_filter_combines_step_and_count(self) -> None:
        import plot_cpp_run as pcr

        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            _write_snapshot(d / "snapshot_00000.csv", 0, [(1.0, 0.0)])
            _write_snapshot(d / "snapshot_00001.csv", 5, [(2.0, 0.0)])
            _write_snapshot(d / "snapshot_00002.csv", 10, [(3.0, 0.0)])
            _write_snapshot(d / "snapshot_00003.csv", 15, [(4.0, 0.0)])
            rec = pcr.load_all_snapshot_records(d)
            out = pcr.filter_snapshots_for_plotting(
                rec, skip_initial_steps=5, skip_initial_snapshots=1
            )
            self.assertEqual([s.step for _, s in out], [10, 15])

    def test_smart_bounds_uses_filtered_paths(self) -> None:
        import plot_cpp_run as pcr

        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            p0 = d / "snapshot_00000.csv"
            p1 = d / "snapshot_00001.csv"
            _write_snapshot(p0, 0, [(1.0, 0.0), (1.0, 0.0)])  # median r = 1
            _write_snapshot(p1, 1, [(100.0, 0.0), (100.0, 0.0)])  # median r = 100
            b_all = pcr.calculate_smart_bounds(d, fallback=1.0, snapshot_paths=[p0, p1])
            b_filtered = pcr.calculate_smart_bounds(d, fallback=1.0, snapshot_paths=[p1])
            self.assertAlmostEqual(b_all, 60.6, places=6)  # median([1,1,100,100]) * 1.2
            self.assertAlmostEqual(b_filtered, 120.0, places=6)


if __name__ == "__main__":
    unittest.main()
