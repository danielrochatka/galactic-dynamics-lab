from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "engine/physics/TPFCore/plots/plot_tpf_4d_xi_motion_probe.py"


def _write_fake_csv(out_dir: Path) -> None:
    (out_dir / "tpf_4d_xi_motion_probe_trajectories.csv").write_text(
        "\n".join(
            [
                "step,time,probe_id,x,y,z,vx,vy,vz,ax,ay,az,a_norm,xi_x,xi_y,xi_z,xi_spatial_norm,theta_trace_4d,invariant_I_4d,r_origin,radial_alignment_to_origin_inward,transverse_fraction_origin,is_near_source,escaped,valid",
                "0,0.0,0,10.0,0.0,0.0,0.0,0.316,0.0,-0.01,0.0,0.0,0.01,1.0,0.0,0.0,1.0,0.1,0.2,10.0,1.0,0.0,0,0,1",
                "1,0.1,0,9.999,0.0316,0.0,-0.001,0.316,0.0,-0.00999,-0.00003,0.0,0.00999,0.999,0.003,0.0,0.999,0.1,0.2,9.99905,1.0,0.0,0,0,1",
                "0,0.0,1,-10.0,0.0,0.0,0.0,-0.316,0.0,0.01,0.0,0.0,0.01,-1.0,0.0,0.0,1.0,0.1,0.2,10.0,1.0,0.0,0,0,1",
                "1,0.1,1,-9.999,-0.0316,0.0,0.001,-0.316,0.0,0.00999,0.00003,0.0,0.00999,-0.999,-0.003,0.0,0.999,0.1,0.2,9.99905,1.0,0.0,0,0,1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


class TestPlotTPF4DXiMotionProbe(unittest.TestCase):
    def test_script_syntax_py_compile(self) -> None:
        subprocess.run([sys.executable, "-m", "py_compile", str(SCRIPT)], check=True)

    def test_script_smoke_generates_pngs_or_warns_if_missing_dependencies(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            out_dir = Path(td)
            _write_fake_csv(out_dir)

            proc = subprocess.run(
                [sys.executable, str(SCRIPT), str(out_dir)],
                check=False,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(proc.returncode, 0)

            try:
                import matplotlib.pyplot  # noqa: F401
                import numpy  # noqa: F401
                import pandas  # noqa: F401
                deps_available = True
            except Exception:
                deps_available = False
            if deps_available:
                self.assertTrue((out_dir / "tpf_4d_xi_motion_probe_xy_trajectories.png").exists())
                self.assertTrue((out_dir / "tpf_4d_xi_motion_probe_xz_trajectories.png").exists())
                self.assertTrue((out_dir / "tpf_4d_xi_motion_probe_yz_trajectories.png").exists())
                self.assertTrue((out_dir / "tpf_4d_xi_motion_probe_radius_vs_time.png").exists())
                self.assertTrue((out_dir / "tpf_4d_xi_motion_probe_acceleration_norm_vs_time.png").exists())
            else:
                self.assertIn("plotting dependencies unavailable", proc.stderr)


if __name__ == "__main__":
    unittest.main()
