"""Pair/COM diagnostics for earth_moon_benchmark and bh_orbit_validation postprocess."""

from __future__ import annotations

import unittest

import numpy as np

from diagnostics import compute_two_body_pair_diagnostics, G_SI


class _Snap:
    __slots__ = ("time", "positions", "velocities")

    def __init__(self, time: float, positions: np.ndarray, velocities: np.ndarray) -> None:
        self.time = time
        self.positions = positions
        self.velocities = velocities


class TestTwoBodyPairDiagnostics(unittest.TestCase):
    def test_equal_mass_binary_com_origin(self) -> None:
        pos = np.array([[-1.0, 0.0], [1.0, 0.0]], dtype=float)
        vel = np.array([[0.0, -1.0], [0.0, 1.0]], dtype=float)
        snaps = [_Snap(0.0, pos, vel)]
        m = np.array([1.0, 1.0])
        d = compute_two_body_pair_diagnostics(snaps, m, "earth_moon_benchmark", 0.0)
        self.assertAlmostEqual(float(d["pair_separation"][0]), 2.0)
        self.assertAlmostEqual(float(d["center_of_mass_x"][0]), 0.0)
        self.assertAlmostEqual(float(d["center_of_mass_y"][0]), 0.0)
        self.assertAlmostEqual(float(d["center_of_mass_radius"][0]), 0.0)
        self.assertEqual(d["variant"], "two_mass")

    def test_bh_orbit_separation_matches_radius(self) -> None:
        pos = np.array([[3.0, 4.0]], dtype=float)
        vel = np.array([[0.0, 1.0]], dtype=float)
        snaps = [_Snap(0.0, pos, vel)]
        m = np.array([1.0])
        d = compute_two_body_pair_diagnostics(snaps, m, "bh_orbit_validation", 100.0)
        self.assertAlmostEqual(float(d["pair_separation"][0]), 5.0)
        self.assertEqual(d["variant"], "star_bh")
        # Newtonian specific: 0.5*v^2 - G*M/r
        r = 5.0
        M = 100.0
        expected_E = 0.5 * 1.0**2 - G_SI * M / r
        self.assertAlmostEqual(float(d["newtonian_specific_energy"][0]), expected_E, delta=1e-20)


if __name__ == "__main__":
    unittest.main()
