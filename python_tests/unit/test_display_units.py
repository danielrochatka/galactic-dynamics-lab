from __future__ import annotations

import unittest

from display_units import format_animation_time_caption


class TestAnimationTimeCaptionUnits(unittest.TestCase):
    def test_auto_mode_stable_when_active_time_unit_is_provided(self) -> None:
        c0 = format_animation_time_caption(30.0, "galaxy", preferred_time_unit="auto", active_time_unit="day")
        c1 = format_animation_time_caption(30.0 * 86400.0, "galaxy", preferred_time_unit="auto", active_time_unit="day")
        c2 = format_animation_time_caption(30.0 * 365.25 * 86400.0, "galaxy", preferred_time_unit="auto", active_time_unit="day")
        self.assertTrue(c0.endswith(" day"))
        self.assertTrue(c1.endswith(" day"))
        self.assertTrue(c2.endswith(" day"))

    def test_explicit_forced_unit_remains_respected(self) -> None:
        c = format_animation_time_caption(7200.0, "galaxy", preferred_time_unit="hr")
        self.assertEqual(c, "t = 2.00 hr")

    def test_animation_time_caption_uses_fixed_two_decimal_precision(self) -> None:
        self.assertEqual(
            format_animation_time_caption(30.0, "galaxy", active_time_unit="day"),
            "t = 0.00 day",
        )
        self.assertEqual(
            format_animation_time_caption(86400.0, "galaxy", active_time_unit="day"),
            "t = 1.00 day",
        )
        self.assertEqual(
            format_animation_time_caption(1.2345 * 86400.0, "galaxy", active_time_unit="day"),
            "t = 1.23 day",
        )

    def test_active_time_unit_overrides_auto_resolution(self) -> None:
        c = format_animation_time_caption(
            0.5 * 365.25 * 86400.0,
            "galaxy",
            preferred_time_unit="auto",
            active_time_unit="min",
        )
        self.assertTrue(c.endswith(" min"))


if __name__ == "__main__":
    unittest.main()
