"""Snapshot filename ordering for plot_cpp_run (numeric step, not lexicographic stem)."""

from __future__ import annotations

import unittest
from pathlib import Path

from plot_cpp_run import snapshot_csv_numeric_sort_key


class TestSnapshotCsvNumericSort(unittest.TestCase):
    def test_mixed_digit_widths_not_lex_order(self) -> None:
        # Lex stem sort would place snapshot_100000 before snapshot_99999 ('1' < '9').
        names = [
            Path("snapshot_100000.csv"),
            Path("snapshot_99999.csv"),
            Path("snapshot_10000.csv"),
        ]
        ordered = sorted(names, key=snapshot_csv_numeric_sort_key)
        self.assertEqual(
            [p.name for p in ordered],
            ["snapshot_10000.csv", "snapshot_99999.csv", "snapshot_100000.csv"],
        )

    def test_padded_five_digit_still_sorted(self) -> None:
        names = [
            Path("snapshot_00100.csv"),
            Path("snapshot_00010.csv"),
            Path("snapshot_01000.csv"),
        ]
        ordered = sorted(names, key=snapshot_csv_numeric_sort_key)
        self.assertEqual(
            [p.name for p in ordered],
            ["snapshot_00010.csv", "snapshot_00100.csv", "snapshot_01000.csv"],
        )


if __name__ == "__main__":
    unittest.main()
