"""
Regression: default 6-digit scientific text collapses ~10–30 m motion at Earth–Moon distance scales.

The C++ snapshot writer used to rely on ostringstream defaults; that matched this failure mode.
Full double round-trip formatting is required in snapshot_*.csv (see cpp_sim/output.cpp).
"""

from __future__ import annotations

import unittest


class TestSnapshotCsvNumericPrecision(unittest.TestCase):
    def test_six_digit_scientific_hides_tens_of_meters_at_1e8_scale(self) -> None:
        x0 = 3.843663e8
        x1 = x0 + 30.0
        self.assertEqual(f"{x0:.6e}", f"{x1:.6e}")

    def test_full_precision_strings_differ_for_sub_meter_delta(self) -> None:
        x0 = 3.843663e8
        x1 = x0 + 0.25
        # 17 significant digits — same spirit as std::numeric_limits<double>::max_digits10
        self.assertNotEqual(f"{x0:.17e}", f"{x1:.17e}")


if __name__ == "__main__":
    unittest.main()
