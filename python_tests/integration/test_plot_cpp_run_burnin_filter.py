"""Integration tests for plot_cpp_run burn-in filtering and geometric framing."""

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

    def test_geometric_framing_respects_filtered_snapshots(self) -> None:
        """Burn-in should shrink the global viewport when outlier early snapshots are dropped."""
        import numpy as np

        import plot_cpp_run as pcr
        from framing import global_viewport_from_snapshots

        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            p0 = d / "snapshot_00000.csv"
            p1 = d / "snapshot_00001.csv"
            # Early frame pulls x extent to ~1000 m; late frame is compact around 10 m.
            _write_snapshot(p0, 0, [(1000.0, 0.0), (1000.0, 0.0)])
            _write_snapshot(p1, 1, [(10.0, 0.0), (10.0, 0.0)])
            rec = pcr.load_all_snapshot_records(d)
            snaps_all = [s for _, s in rec]
            snaps_filt = [s for _, s in pcr.filter_snapshots_for_plotting(rec, 0, 1)]
            extra = np.array([[0.0, 0.0]], dtype=np.float64)
            vp_all = global_viewport_from_snapshots(
                snaps_all, extra_xy=extra, fallback_half_axis=1.0
            )
            vp_filt = global_viewport_from_snapshots(
                snaps_filt, extra_xy=extra, fallback_half_axis=1.0
            )
            self.assertGreater(vp_all.half_axis, vp_filt.half_axis)
            # Filtered: x in {0, 10}; span ~10 after trim → half ≈ 0.5 * span * 1.15
            self.assertGreater(vp_filt.half_axis, 5.0)
            self.assertLess(vp_filt.half_axis, 6.0)

    def test_resolve_skip_settings_cli_overrides_run_info(self) -> None:
        import plot_cpp_run as pcr

        run_info = {
            "plot_skip_initial_steps": 10,
            "plot_skip_initial_snapshots": 3,
        }
        steps, snaps = pcr.resolve_burnin_skip_settings(
            run_info,
            cli_skip_initial_steps=5,
            cli_skip_initial_snapshots=None,
        )
        self.assertEqual(steps, 5)
        self.assertEqual(snaps, 3)

    def test_resolve_skip_settings_missing_defaults_to_zero(self) -> None:
        import plot_cpp_run as pcr

        steps, snaps = pcr.resolve_burnin_skip_settings({}, None, None)
        self.assertEqual(steps, 0)
        self.assertEqual(snaps, 0)


if __name__ == "__main__":
    unittest.main()
