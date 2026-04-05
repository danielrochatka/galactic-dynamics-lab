from __future__ import annotations

import unittest

import numpy as np

from plot_rotation_curve import (
    NEWTONIAN_REFERENCE_M_BH,
    ROTATION_CURVE_X_FALLBACK_MAX,
    _resolve_display_unit_preferences,
    newtonian_m_bh_from_run_info,
    resolve_newtonian_overlay_mass,
    rotation_curve_x_max,
    scatter_label_from_run_info,
)


class TestRotationCurveRunInfoPrecedence(unittest.TestCase):
    def test_newtonian_mass_prefers_effective_then_configured_then_legacy_then_fallback(self) -> None:
        self.assertEqual(
            newtonian_m_bh_from_run_info(
                {
                    "effective_bh_mass": 3.0,
                    "configured_bh_mass": 2.0,
                    "bh_mass": 1.0,
                }
            ),
            3.0,
        )
        self.assertEqual(
            newtonian_m_bh_from_run_info({"configured_bh_mass": 2.0, "bh_mass": 1.0}),
            2.0,
        )
        self.assertEqual(newtonian_m_bh_from_run_info({"bh_mass": 1.0}), 1.0)
        self.assertEqual(newtonian_m_bh_from_run_info({}), NEWTONIAN_REFERENCE_M_BH)

    def test_resolve_overlay_mass_uses_same_precedence_and_cli_override(self) -> None:
        mass, source, used_fallback = resolve_newtonian_overlay_mass({"effective_bh_mass": 9.0}, None)
        self.assertEqual(mass, 9.0)
        self.assertEqual(source, "run_info(effective/configured/legacy bh_mass)")
        self.assertFalse(used_fallback)

        mass, source, used_fallback = resolve_newtonian_overlay_mass({"configured_bh_mass": 7.0}, None)
        self.assertEqual(mass, 7.0)
        self.assertEqual(source, "run_info(effective/configured/legacy bh_mass)")
        self.assertFalse(used_fallback)

        mass, source, used_fallback = resolve_newtonian_overlay_mass({"bh_mass": 5.0}, None)
        self.assertEqual(mass, 5.0)
        self.assertEqual(source, "run_info(effective/configured/legacy bh_mass)")
        self.assertFalse(used_fallback)

        mass, source, used_fallback = resolve_newtonian_overlay_mass({}, None)
        self.assertEqual(mass, NEWTONIAN_REFERENCE_M_BH)
        self.assertEqual(source, "fallback_benchmark(NEWTONIAN_REFERENCE_M_BH)")
        self.assertTrue(used_fallback)

        mass, source, used_fallback = resolve_newtonian_overlay_mass({"effective_bh_mass": 9.0}, 11.0)
        self.assertEqual(mass, 11.0)
        self.assertEqual(source, "cli(--M-bh)")
        self.assertFalse(used_fallback)

    def test_rotation_curve_x_max_prefers_effective_then_configured_then_legacy_then_fallback(self) -> None:
        r_data = np.array([1.0, 2.0])
        self.assertEqual(
            rotation_curve_x_max(
                {
                    "effective_galaxy_radius": 30.0,
                    "configured_galaxy_radius": 20.0,
                    "galaxy_radius": 10.0,
                },
                r_data,
            ),
            60.0,
        )
        self.assertEqual(
            rotation_curve_x_max({"configured_galaxy_radius": 20.0, "galaxy_radius": 10.0}, r_data),
            40.0,
        )
        self.assertEqual(rotation_curve_x_max({"galaxy_radius": 10.0}, r_data), 20.0)
        self.assertEqual(rotation_curve_x_max({}, r_data), ROTATION_CURVE_X_FALLBACK_MAX)

    def test_scatter_label_prefers_effective_then_configured_then_legacy(self) -> None:
        self.assertEqual(
            scatter_label_from_run_info(
                {
                    "effective_physics_package": "TPFCore",
                    "configured_physics_package": "Newtonian",
                    "physics_package": "Legacy",
                }
            ),
            "Simulated stars (TPFCore)",
        )
        self.assertEqual(
            scatter_label_from_run_info({"configured_physics_package": "Newtonian"}),
            "Simulated stars (Newtonian)",
        )
        self.assertEqual(
            scatter_label_from_run_info({"physics_package": "Legacy"}),
            "Simulated stars (Legacy)",
        )
        self.assertEqual(scatter_label_from_run_info({}), "Simulated stars")

    def test_display_unit_preference_effective_then_configured_then_legacy(self) -> None:
        self.assertEqual(
            _resolve_display_unit_preferences(
                {
                    "effective_display_distance_unit": "kpc",
                    "configured_display_distance_unit": "pc",
                    "display_distance_unit": "m",
                    "effective_display_velocity_unit": "km/s",
                    "configured_display_velocity_unit": "m/s",
                    "display_velocity_unit": "cm/s",
                }
            ),
            ("kpc", "km/s"),
        )
        self.assertEqual(
            _resolve_display_unit_preferences(
                {
                    "configured_display_distance_unit": "pc",
                    "display_distance_unit": "m",
                    "configured_display_velocity_unit": "m/s",
                    "display_velocity_unit": "cm/s",
                }
            ),
            ("pc", "m/s"),
        )
        self.assertEqual(
            _resolve_display_unit_preferences(
                {
                    "display_distance_unit": "m",
                    "display_velocity_unit": "cm/s",
                }
            ),
            ("m", "cm/s"),
        )
        self.assertEqual(_resolve_display_unit_preferences({}), ("auto", "auto"))


if __name__ == "__main__":
    unittest.main()
