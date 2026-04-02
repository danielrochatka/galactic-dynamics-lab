"""Unit tests for geometric square viewport framing (framing module)."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from framing import (
    DynamicViewportSmoother,
    SquareViewport,
    compute_square_viewport,
    global_viewport_from_snapshots,
)


class TestGeometricFraming(unittest.TestCase):
    def test_two_body_line_with_origin(self) -> None:
        """Earth at origin, Moon at +R: viewport spans ~R/2 to each side from center."""
        xy = np.array([[0.0, 0.0], [3.84e8, 0.0]], dtype=np.float64)
        vp = compute_square_viewport(xy, trim_fraction=0.0, margin=1.15, fallback_half_axis=1e6)
        self.assertAlmostEqual(vp.center_x, 1.92e8, delta=1e3)
        self.assertAlmostEqual(vp.half_axis, 0.5 * 3.84e8 * 1.15, delta=1e3)

    def test_quantile_trim_drops_outlier(self) -> None:
        x = np.concatenate([np.arange(1, 101, dtype=np.float64), [1.0e6]])
        pts = np.column_stack([x, np.zeros_like(x)])
        vp_strict = compute_square_viewport(
            pts, trim_fraction=0.0, margin=1.0, fallback_half_axis=1.0
        )
        self.assertGreater(vp_strict.half_axis, 1e5)
        vp_trim = compute_square_viewport(
            pts, trim_fraction=0.01, margin=1.0, fallback_half_axis=1.0
        )
        self.assertLess(vp_trim.half_axis, 100.0)

    def test_global_viewport_stacks_snapshots(self) -> None:
        class S:
            def __init__(self, pos: np.ndarray) -> None:
                self.positions = pos

        s0 = S(np.array([[0.0, 0.0], [2.0, 0.0]]))
        s1 = S(np.array([[0.0, 0.0], [2.0, 0.0]]))
        vp = global_viewport_from_snapshots(
            [s0, s1],
            extra_xy=None,
            trim_fraction=0.0,
            margin=1.0,
            fallback_half_axis=1.0,
        )
        self.assertAlmostEqual(vp.half_axis, 1.0)

    def test_dynamic_smoother_ema_log_half_axis(self) -> None:
        sm = DynamicViewportSmoother(alpha=0.5, min_half_axis=1.0)
        v0 = sm.step(SquareViewport(0.0, 0.0, 10.0))
        self.assertEqual(v0.half_axis, 10.0)
        v1 = sm.step(SquareViewport(10.0, 0.0, 40.0))
        self.assertAlmostEqual(v1.center_x, 5.0)
        self.assertAlmostEqual(v1.half_axis, 20.0)


if __name__ == "__main__":
    unittest.main()
