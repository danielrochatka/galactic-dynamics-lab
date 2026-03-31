"""Regression: static viewport must not blank plots when data has large radial spread / escapers."""

from __future__ import annotations

import unittest

import numpy as np

from plot_cpp_run import (
    bulk_chebyshev_extent_half_axis,
    fraction_stars_inside_square_viewport,
    galaxy_velocity_gated_target_limit,
    initial_snapshot_plot_title,
    static_viewport_radius_validated,
)


class TestViewportStatic(unittest.TestCase):
    def test_escapers_velocity_gate_collapses_static_viewport(self) -> None:
        """When stable-star branch fails, velocity gate uses 1.2*galaxy_radius; escapers can be far outside."""
        n = 1000
        th = np.linspace(0, 2 * np.pi, n, endpoint=False)
        r = 1e27
        pos = np.column_stack([r * np.cos(th), r * np.sin(th)])
        # Radial motion -> unstable mask -> fallback 1.2*galaxy_radius, far below star radius
        vel = np.column_stack([np.cos(th) * 1e5, np.sin(th) * 1e5])
        gr = 1e20
        fb = 150.0
        cand = galaxy_velocity_gated_target_limit(pos, vel, gr, fb)
        self.assertAlmostEqual(cand, 1.2e20, delta=1e10)
        self.assertLess(fraction_stars_inside_square_viewport(pos, cand), 0.5)

        out = static_viewport_radius_validated(pos, vel, gr, fb)
        self.assertGreaterEqual(fraction_stars_inside_square_viewport(pos, out), 0.5)
        self.assertGreater(out, cand)

    def test_weird_velocities_still_show_bulk(self) -> None:
        """Degenerate / odd velocities should not yield an empty frame when positions are spread."""
        n = 400
        rng = np.random.default_rng(1)
        pos = rng.normal(scale=1e19, size=(n, 2))
        vel = np.full((n, 2), np.nan)
        gr = 1e20
        fb = 150.0
        out = static_viewport_radius_validated(pos, vel, gr, fb)
        frac = fraction_stars_inside_square_viewport(pos, out)
        self.assertGreaterEqual(frac, 0.5)
        self.assertTrue(np.isfinite(out) and out > 0)

    def test_bulk_chebyshev_reasonable(self) -> None:
        pos = np.array([[1.0, 0.0], [0.0, 2.0], [100.0, 0.0]], dtype=np.float64)
        h = bulk_chebyshev_extent_half_axis(pos, 66.0, 1.0)
        self.assertGreaterEqual(h, 2.0)

    def test_initial_title_burn_in(self) -> None:
        self.assertIn("First plotted", initial_snapshot_plot_title(False, 2500, burn_in_plotting=True))
        self.assertEqual(
            initial_snapshot_plot_title(False, 0, burn_in_plotting=False),
            "Galaxy – Initial (C++ run)",
        )
        self.assertIn("cooling", initial_snapshot_plot_title(True, 10, burn_in_plotting=True).lower())


if __name__ == "__main__":
    unittest.main()
