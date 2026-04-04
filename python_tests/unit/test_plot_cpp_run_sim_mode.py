"""Regression: validation-mode runs use snapshot-based scale, not config galaxy_radius."""

from __future__ import annotations

import unittest

import numpy as np

from plot_cpp_run import (
    Snapshot,
    resolve_galaxy_radius_meters,
    should_draw_central_bh_marker,
    simulation_mode_name_from_run_info,
)


class TestSimulationModeFromRunInfo(unittest.TestCase):
    def test_int_from_cpp_enum_second_line_wins(self) -> None:
        # load_run_info overwrites duplicate keys; final simulation_mode is often int
        ri = {"simulation_mode": 1}
        self.assertEqual(simulation_mode_name_from_run_info(ri), "earth_moon_benchmark")

    def test_string_mode(self) -> None:
        ri = {"simulation_mode": "symmetric_pair"}
        self.assertEqual(simulation_mode_name_from_run_info(ri), "symmetric_pair")

    def test_empty_defaults_galaxy(self) -> None:
        self.assertEqual(simulation_mode_name_from_run_info({}), "galaxy")


class TestResolveGalaxyRadiusNonGalaxy(unittest.TestCase):
    def test_two_body_ignores_config_galaxy_radius(self) -> None:
        """Stale galaxy_radius (e.g. 50) must not shrink the velocity-gate reference for Earth–Moon SI."""
        pos = np.array([[0.0, 0.0], [3.844e8, 0.0]], dtype=float)
        vel = np.zeros((2, 2), dtype=float)
        snap = Snapshot(step=0, time=0.0, positions=pos, velocities=vel)
        ri = {"galaxy_radius": 50.0, "simulation_mode": 13}
        r = resolve_galaxy_radius_meters(ri, [snap], render_radius_arg=150.0)
        self.assertAlmostEqual(r, 3.844e8, delta=1.0)

    def test_bh_orbit_int_maps(self) -> None:
        self.assertEqual(simulation_mode_name_from_run_info({"simulation_mode": 14}), "bh_orbit_validation")

    def test_galaxy_mode_uses_config_radius(self) -> None:
        ri = {"galaxy_radius": 50.0, "simulation_mode": 0}
        pos = np.array([[10.0, 0.0]], dtype=float)
        vel = np.zeros((1, 2), dtype=float)
        snap = Snapshot(step=0, time=0.0, positions=pos, velocities=vel)
        r = resolve_galaxy_radius_meters(ri, [snap], render_radius_arg=150.0)
        self.assertEqual(r, 50.0)


class TestCentralBhMarkerDecision(unittest.TestCase):
    def test_earth_moon_benchmark_hides_marker(self) -> None:
        ri = {"simulation_mode": 13, "bh_mass": 1.0e20}
        self.assertFalse(should_draw_central_bh_marker(ri, "earth_moon_benchmark"))

    def test_bh_orbit_validation_shows_marker_when_mass_positive(self) -> None:
        ri = {"simulation_mode": 14, "bh_mass": 1.0e30}
        self.assertTrue(should_draw_central_bh_marker(ri, "bh_orbit_validation"))

    def test_symmetric_pair_respects_include_bh_flag(self) -> None:
        ri = {"simulation_mode": 2, "bh_mass": 1.0e30, "validation_symmetric_include_bh": 0}
        self.assertFalse(should_draw_central_bh_marker(ri, "symmetric_pair"))

    def test_galaxy_requires_positive_bh_mass(self) -> None:
        ri = {"simulation_mode": 0, "bh_mass": 0.0}
        self.assertFalse(should_draw_central_bh_marker(ri, "galaxy"))


if __name__ == "__main__":
    unittest.main()
