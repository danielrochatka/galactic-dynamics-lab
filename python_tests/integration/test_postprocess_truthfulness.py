"""Regression tests for post-processing cutoff and mass-source truthfulness."""

from __future__ import annotations

import contextlib
import importlib.util
import io
import os
import subprocess
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

    def test_cooled_run_title_label_helper(self) -> None:
        import plot_cpp_run as pcr

        self.assertEqual(
            pcr.initial_snapshot_plot_title(True, 100),
            "Galaxy – First saved snapshot after cooling (C++ run)",
        )
        self.assertEqual(
            pcr.initial_snapshot_plot_title(True, 0),
            "Galaxy – Initial (C++ run)",
        )

    def test_plot_cpp_run_warns_on_cooling_active(self) -> None:
        import plot_cpp_run as pcr

        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            _write_snapshot(d / "snapshot_00000.csv", step=100)
            (d / "run_info.txt").write_text(
                "\n".join(
                    [
                        "cooling_active\t1",
                        "cooling_steps\t50",
                        "first_saved_snapshot_step\t100",
                        "first_saved_snapshot_time\t1.0",
                        "bh_mass\t1.0e6",
                        "galaxy_radius\t50.0",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            argv_old = sys.argv[:]
            buf = io.StringIO()
            try:
                sys.argv = ["plot_cpp_run.py", str(d), "--no-animation", "--no-diagnostics"]
                with contextlib.redirect_stdout(buf):
                    pcr.main()
            finally:
                sys.argv = argv_old
            out = buf.getvalue()
            self.assertIn("WARNING: cooling phase was active", out)
            self.assertIn("first_saved_snapshot_step=100", out)

    def test_run_info_contains_explicit_cooling_fields(self) -> None:
        cpp_dir = REPO_ROOT / "cpp_sim"
        bin_path = cpp_dir / "galaxy_sim"
        if not bin_path.exists():
            self.skipTest("requires built cpp_sim/galaxy_sim")

        run_id = "test_cooling_meta_py"
        out_dir = f"outputs/{run_id}"
        cmd = [
            str(bin_path),
            "galaxy",
            "--physics_package=TPFCore",
            "--n_stars=8",
            "--n_steps=20",
            "--snapshot_every=10",
            "--tpf_cooling_fraction=0.5",
            "--save_snapshots=false",
            "--save_run_info=true",
            f"--output_dir={out_dir}",
        ]
        subprocess.run(
            cmd,
            cwd=str(cpp_dir),
            check=True,
            env={**os.environ, "GALAXY_RUN_CONFIG": ""},
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        run_info_path = cpp_dir / out_dir / "run_info.txt"
        txt = run_info_path.read_text(encoding="utf-8")
        self.assertIn("cooling_active\t1", txt)
        self.assertIn("cooling_steps\t10", txt)
        self.assertIn("cooling_end_step\t9", txt)
        self.assertIn("first_saved_snapshot_step\t10", txt)
        self.assertIn("first_saved_snapshot_time\t", txt)


if __name__ == "__main__":
    unittest.main()
