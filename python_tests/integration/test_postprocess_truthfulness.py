"""Regression tests for post-processing cutoff and mass-source truthfulness."""

from __future__ import annotations

import contextlib
import importlib.util
import io
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


def _write_snapshot(path: Path, step: int = 0) -> None:
    content = "\n".join(
        [
            f"# step,{step},time,0.0e+00",
            "i,x,y,vx,vy,mass",
            "0,1.0,0.0,0.0,1.0,1.0",
            "1,0.0,1.0,-1.0,0.0,1.0",
        ]
    )
    path.write_text(content + "\n", encoding="utf-8")


@unittest.skipUnless(_have_plot_stack(), "requires numpy and matplotlib")
class TestPostprocessTruthfulness(unittest.TestCase):
    def test_diagnostic_cutoff_precedence_cli(self) -> None:
        import plot_cpp_run as pcr

        cutoff, src = pcr.resolve_diagnostic_cutoff_radius(
            {"diagnostic_cutoff_radius": 100.0, "galaxy_radius": 50.0},
            cli_cutoff_radius=12.5,
        )
        self.assertEqual(cutoff, 12.5)
        self.assertEqual(src, "cli(--diagnostic-cutoff-radius)")

    def test_diagnostic_cutoff_precedence_runinfo_override(self) -> None:
        import plot_cpp_run as pcr

        cutoff, src = pcr.resolve_diagnostic_cutoff_radius(
            {"diagnostic_cutoff_radius": 100.0, "galaxy_radius": 50.0},
            cli_cutoff_radius=None,
        )
        self.assertEqual(cutoff, 100.0)
        self.assertEqual(src, "run_info(diagnostic_cutoff_radius)")

    def test_diagnostic_cutoff_precedence_galaxy_radius_fallback(self) -> None:
        import plot_cpp_run as pcr

        cutoff, src = pcr.resolve_diagnostic_cutoff_radius(
            {"galaxy_radius": 77.0},
            cli_cutoff_radius=None,
        )
        self.assertEqual(cutoff, 77.0)
        self.assertEqual(src, "run_info(galaxy_radius)")

    def test_diagnostic_cutoff_error_when_undefined(self) -> None:
        import plot_cpp_run as pcr

        with self.assertRaises(SystemExit):
            pcr.resolve_diagnostic_cutoff_radius({}, cli_cutoff_radius=None)

    def test_rotation_mass_source_precedence(self) -> None:
        import plot_rotation_curve as prc

        m, src, fb = prc.resolve_newtonian_overlay_mass({"bh_mass": 42.0}, cli_m_bh=11.0)
        self.assertEqual(m, 11.0)
        self.assertEqual(src, "cli(--M-bh)")
        self.assertFalse(fb)

        m, src, fb = prc.resolve_newtonian_overlay_mass({"bh_mass": 42.0}, cli_m_bh=None)
        self.assertEqual(m, 42.0)
        self.assertEqual(src, "run_info(bh_mass)")
        self.assertFalse(fb)

    def test_rotation_fallback_mass_prints_warning(self) -> None:
        import plot_rotation_curve as prc

        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            _write_snapshot(d / "snapshot_00000.csv", step=0)
            out = d / "rotation_curve.png"
            buf = io.StringIO()
            argv_old = sys.argv[:]
            try:
                sys.argv = ["plot_rotation_curve.py", str(d), "--output", str(out)]
                with contextlib.redirect_stdout(buf):
                    prc.main()
            finally:
                sys.argv = argv_old
            txt = buf.getvalue()
            self.assertIn("source=fallback_benchmark", txt)
            self.assertIn("using fallback benchmark mass", txt)


if __name__ == "__main__":
    unittest.main()
