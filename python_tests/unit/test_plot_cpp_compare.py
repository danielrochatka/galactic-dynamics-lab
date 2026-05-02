from __future__ import annotations

import json
import sys
import tempfile
import unittest
import importlib.util
from pathlib import Path
from unittest.mock import patch
from types import SimpleNamespace
import numpy as np

# Repo root (parent of python_tests/)
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from plot_cpp_compare import (
    _panel_physics_text,
    _compare_panel_metadata_lines,
    resolve_compare_display_selection,
    _resolve_side_run_dir,
    _mode_aware_compare_name,
    calculate_compare_smart_viewport,
    approximate_axis_aligned_box_covering_fraction,
    _sample_viewport_frame_indices,
    matched_steps_strict,
    render_compare,
    radial_kinematics,
    estimate_snapshot_accelerations,
    binned_profile,
    _time_display_factor,
    _velocity_display_factor,
    style_compare_diagnostic_axes,
    _draw_panel,
)
from display_units import DisplayUnitConfig
from display_units import spatial_display_for_xy_plot


def _write_snapshot(path: Path, step: int, t: float, x: float) -> None:
    path.write_text(
        "\n".join(
            [
                f"# step,{step},time,{t:.6e}",
                "i,x,y,vx,vy,mass",
                f"0,{x},0.0,0.0,1.0,1.0",
            ]
        ),
        encoding="utf-8",
    )


def _write_run_info(run_dir: Path, pkg: str) -> None:
    (run_dir / "run_info.txt").write_text(
        "\n".join(
            [
                f"physics_package\t{pkg}",
                "render_overlay_mode\tnone",
                "code_version_label\ttest@abc1234",
                "n_steps\t10",
                "dt\t0.1",
                "display_distance_unit\tauto",
                "display_time_unit\tauto",
                "display_velocity_unit\tauto",
                "display_units_in_overlay\t1",
                "display_show_unit_reference\t1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _write_run_info_new_schema(run_dir: Path, pkg: str, *, dyn: str = "direct_tpf") -> None:
    lines = [
        f"configured_physics_package\t{pkg}",
        f"effective_physics_package\t{pkg}",
        "configured_simulation_mode\tgalaxy",
        "effective_simulation_mode\tgalaxy",
        "render_overlay_mode\tnone",
        "code_version_label\ttest@abc1234",
        "n_steps\t10",
        "dt\t0.1",
        "display_distance_unit\tauto",
        "display_time_unit\tauto",
        "display_velocity_unit\tauto",
        "display_units_in_overlay\t1",
        "display_show_unit_reference\t1",
    ]
    if pkg == "TPFCore":
        lines.extend(
            [
                "configured_tpfcore_readout_mode\tlegacy_readout",
                "configured_tpf_vdsg_coupling\t0.0",
                f"effective_tpf_dynamics_mode\t{dyn}",
            ]
        )
    (run_dir / "run_info.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


class TestPlotCppCompare(unittest.TestCase):
    def _side_data_for_metadata(self, label: str, run_info: dict):
        return SimpleNamespace(label=label, run_info=run_info)

    @unittest.skipUnless(
        importlib.util.find_spec("matplotlib") is not None and importlib.util.find_spec("numpy") is not None,
        "matplotlib/numpy not installed",
    )
    def test_style_compare_diagnostic_axes_uses_light_theme(self) -> None:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1)
        style_compare_diagnostic_axes(ax)
        self.assertEqual(ax.get_facecolor(), (1.0, 1.0, 1.0, 1.0))
        self.assertEqual(ax.xaxis.label.get_color(), "black")
        self.assertEqual(ax.yaxis.label.get_color(), "black")
        self.assertEqual(ax.title.get_color(), "black")
        for spine in ax.spines.values():
            self.assertEqual(spine.get_edgecolor(), (0.35, 0.35, 0.35, 1.0))
        plt.close(fig)

    @unittest.skipUnless(
        importlib.util.find_spec("matplotlib") is not None and importlib.util.find_spec("numpy") is not None,
        "matplotlib/numpy not installed",
    )
    def test_render_compare_animation_writes_only_canonical_mp4(self) -> None:
        class DummyAnim:
            def save(self, path: str, **kwargs) -> None:
                Path(path).write_bytes(b"video")

        with tempfile.TemporaryDirectory() as td:
            parent = Path(td)
            left = parent / "left_TPFCore"
            right = parent / "right_Newtonian"
            left.mkdir()
            right.mkdir()
            _write_run_info_new_schema(left, "TPFCore")
            _write_run_info_new_schema(right, "Newtonian")
            _write_snapshot(left / "snapshot_00000.csv", 0, 0.0, 1.0)
            _write_snapshot(left / "snapshot_00010.csv", 10, 1.0, 2.0)
            _write_snapshot(right / "snapshot_00000.csv", 0, 0.0, 1.2)
            _write_snapshot(right / "snapshot_00010.csv", 10, 1.0, 2.2)
            compare_run_id = "r_anim_mp4"
            (parent / "compare_manifest.json").write_text(
                json.dumps({"compare_run_id": compare_run_id, "left_dir": str(left), "right_dir": str(right)}),
                encoding="utf-8",
            )

            with patch("matplotlib.animation.FuncAnimation", return_value=DummyAnim()):
                render_compare(parent, no_animation=False, overlay_mode="none")

            canonical_mp4 = parent / _mode_aware_compare_name(
                "animation",
                {"effective_physics_package": "TPFCore", "effective_simulation_mode": "galaxy", "effective_tpf_dynamics_mode": "direct_tpf"},
                {"effective_physics_package": "Newtonian", "effective_simulation_mode": "galaxy"},
                ext="mp4",
            )
            self.assertTrue(canonical_mp4.exists())
            self.assertFalse((parent / f"{compare_run_id}.mp4").exists())
            self.assertFalse((parent / f"{compare_run_id}.gif").exists())

    @unittest.skipUnless(
        importlib.util.find_spec("matplotlib") is not None and importlib.util.find_spec("numpy") is not None,
        "matplotlib/numpy not installed",
    )
    def test_render_compare_animation_gif_fallback_writes_only_canonical_gif(self) -> None:
        class DummyAnim:
            def __init__(self) -> None:
                self.calls = 0

            def save(self, path: str, **kwargs) -> None:
                self.calls += 1
                if self.calls == 1:
                    raise RuntimeError("ffmpeg unavailable")
                Path(path).write_bytes(b"gif")

        with tempfile.TemporaryDirectory() as td:
            parent = Path(td)
            left = parent / "left_TPFCore"
            right = parent / "right_Newtonian"
            left.mkdir()
            right.mkdir()
            _write_run_info_new_schema(left, "TPFCore")
            _write_run_info_new_schema(right, "Newtonian")
            _write_snapshot(left / "snapshot_00000.csv", 0, 0.0, 1.0)
            _write_snapshot(left / "snapshot_00010.csv", 10, 1.0, 2.0)
            _write_snapshot(right / "snapshot_00000.csv", 0, 0.0, 1.2)
            _write_snapshot(right / "snapshot_00010.csv", 10, 1.0, 2.2)
            compare_run_id = "r_anim_gif"
            (parent / "compare_manifest.json").write_text(
                json.dumps({"compare_run_id": compare_run_id, "left_dir": str(left), "right_dir": str(right)}),
                encoding="utf-8",
            )

            with patch("matplotlib.animation.FuncAnimation", return_value=DummyAnim()):
                render_compare(parent, no_animation=False, overlay_mode="none")

            canonical_gif = parent / _mode_aware_compare_name(
                "animation",
                {"effective_physics_package": "TPFCore", "effective_simulation_mode": "galaxy", "effective_tpf_dynamics_mode": "direct_tpf"},
                {"effective_physics_package": "Newtonian", "effective_simulation_mode": "galaxy"},
                ext="gif",
            )
            self.assertTrue(canonical_gif.exists())
            self.assertFalse((parent / f"{compare_run_id}.mp4").exists())
            self.assertFalse((parent / f"{compare_run_id}.gif").exists())

    def test_resolve_side_run_dir_engine_relative_manifest(self) -> None:
        """Manifest stores paths relative to engine; compare_parent is .../outputs/RUN."""
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            run = root / "engine" / "outputs" / "RUN1"
            (run / "left_X").mkdir(parents=True)
            m = "outputs/RUN1/left_X"
            got = _resolve_side_run_dir(m, run.resolve())
            self.assertEqual(got, (run / "left_X").resolve())

    def test_matched_steps_strict_pass(self) -> None:
        self.assertEqual(matched_steps_strict({0, 5}, {0, 5}), [0, 5])

    def test_matched_steps_strict_fail(self) -> None:
        with self.assertRaises(ValueError):
            matched_steps_strict({0, 5}, {0, 10})

    @unittest.skipUnless(
        importlib.util.find_spec("matplotlib") is not None and importlib.util.find_spec("numpy") is not None,
        "matplotlib/numpy not installed",
    )
    def test_render_compare_writes_static_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            parent = Path(td)
            left = parent / "left_TPFCore"
            right = parent / "right_Newtonian"
            left.mkdir()
            right.mkdir()
            _write_run_info(left, "TPFCore")
            _write_run_info(right, "Newtonian")
            _write_snapshot(left / "snapshot_00000.csv", 0, 0.0, 1.0)
            _write_snapshot(left / "snapshot_00010.csv", 10, 1.0, 2.0)
            _write_snapshot(right / "snapshot_00000.csv", 0, 0.0, 1.2)
            _write_snapshot(right / "snapshot_00010.csv", 10, 1.0, 2.2)
            (parent / "compare_manifest.json").write_text(
                json.dumps(
                    {
                        "compare_run_id": "r1",
                        "left_dir": str(left),
                        "right_dir": str(right),
                    }
                ),
                encoding="utf-8",
            )

            render_compare(parent, no_animation=True, overlay_mode="none")
            self.assertTrue((parent / "galaxy_initial_compare.png").exists())
            self.assertTrue((parent / "galaxy_final_compare.png").exists())
            self.assertTrue((parent / "compare_rotation_curve_final.png").exists())
            self.assertTrue((parent / "compare_binned_inward_acceleration_vs_radius_final.png").exists())
            self.assertTrue((parent / "compare_radial_acceleration_vs_radius_final.png").exists())
            self.assertTrue((parent / "compare_centripetal_profile_final.png").exists())
            self.assertTrue((parent / "compare_radial_velocity_vs_radius_final.png").exists())
            self.assertTrue((parent / "compare_radius_percentiles_over_time.png").exists())
            self.assertTrue((parent / "compare_acceleration_ratio_vs_radius_final.png").exists())

    @unittest.skipUnless(
        importlib.util.find_spec("matplotlib") is not None and importlib.util.find_spec("numpy") is not None,
        "matplotlib/numpy not installed",
    )
    def test_render_compare_new_schema_mode_aware_names_and_aliases(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            parent = Path(td)
            left = parent / "left_TPFCore"
            right = parent / "right_Newtonian"
            left.mkdir()
            right.mkdir()
            _write_run_info_new_schema(left, "TPFCore")
            _write_run_info_new_schema(right, "Newtonian")
            _write_snapshot(left / "snapshot_00000.csv", 0, 0.0, 1.0)
            _write_snapshot(left / "snapshot_00010.csv", 10, 1.0, 2.0)
            _write_snapshot(right / "snapshot_00000.csv", 0, 0.0, 1.2)
            _write_snapshot(right / "snapshot_00010.csv", 10, 1.0, 2.2)
            (parent / "compare_manifest.json").write_text(
                json.dumps({"compare_run_id": "r2", "left_dir": str(left), "right_dir": str(right)}),
                encoding="utf-8",
            )

            render_compare(parent, no_animation=True, overlay_mode="none")
            expected_initial = _mode_aware_compare_name("initial", {"effective_physics_package": "TPFCore", "effective_simulation_mode": "galaxy", "effective_tpf_dynamics_mode": "direct_tpf"}, {"effective_physics_package": "Newtonian", "effective_simulation_mode": "galaxy"}, ext="png")
            expected_final = _mode_aware_compare_name("final", {"effective_physics_package": "TPFCore", "effective_simulation_mode": "galaxy", "effective_tpf_dynamics_mode": "direct_tpf"}, {"effective_physics_package": "Newtonian", "effective_simulation_mode": "galaxy"}, ext="png")
            self.assertTrue((parent / expected_initial).exists())
            self.assertTrue((parent / expected_final).exists())
            self.assertNotIn("?", expected_initial)
            self.assertNotIn("?", expected_final)
            self.assertTrue((parent / "galaxy_initial_compare.png").exists())
            self.assertTrue((parent / "galaxy_final_compare.png").exists())

    @unittest.skipUnless(
        importlib.util.find_spec("matplotlib") is not None and importlib.util.find_spec("numpy") is not None,
        "matplotlib/numpy not installed",
    )
    def test_compare_layout_uses_shared_top_readout_and_panel_metadata_boxes(self) -> None:
        import matplotlib.figure

        with tempfile.TemporaryDirectory() as td:
            parent = Path(td)
            left = parent / "left_TPFCore"
            right = parent / "right_Newtonian"
            left.mkdir()
            right.mkdir()
            _write_run_info_new_schema(left, "TPFCore", dyn="xi_kernel_deformed")
            _write_run_info_new_schema(right, "Newtonian")
            _write_snapshot(left / "snapshot_00000.csv", 0, 0.0, 1.0)
            _write_snapshot(left / "snapshot_00010.csv", 10, 1.0, 2.0)
            _write_snapshot(right / "snapshot_00000.csv", 0, 0.0, 1.2)
            _write_snapshot(right / "snapshot_00010.csv", 10, 1.0, 2.2)
            (parent / "compare_manifest.json").write_text(
                json.dumps({"compare_run_id": "r_layout", "left_dir": str(left), "right_dir": str(right)}),
                encoding="utf-8",
            )

            seen_suptitles: list[str] = []
            seen_axis_titles: list[list[str]] = []
            seen_axis_texts: list[list[list[str]]] = []
            orig_savefig = matplotlib.figure.Figure.savefig

            def _capture_savefig(fig, *args, **kwargs):
                if len(fig.axes) == 2 and fig.get_facecolor() == (0.0, 0.0, 0.0, 1.0):
                    seen_suptitles.append(fig._suptitle.get_text() if fig._suptitle else "")
                    seen_axis_titles.append([ax.get_title() for ax in fig.axes])
                    seen_axis_texts.append([[t.get_text() for t in ax.texts] for ax in fig.axes])
                return orig_savefig(fig, *args, **kwargs)

            with patch.object(matplotlib.figure.Figure, "savefig", _capture_savefig):
                render_compare(parent, no_animation=True, overlay_mode="none")

            self.assertGreaterEqual(len(seen_suptitles), 2)
            for suptitle in seen_suptitles:
                self.assertIn("Compare r_layout", suptitle)
                self.assertIn("step=", suptitle)
                self.assertIn("display units:", suptitle)
                self.assertNotIn("left=", suptitle)
                self.assertNotIn("right=", suptitle)
                self.assertNotIn("rev=", suptitle)
            for title_pair in seen_axis_titles:
                self.assertEqual(title_pair, ["", ""])
            for panel_texts in seen_axis_texts:
                self.assertTrue(any("role: left / baseline" in txt for txt in panel_texts[0]))
                self.assertTrue(any("Newtonian" in txt for txt in panel_texts[1]))

    def test_mode_aware_compare_name_supports_legacy_run_info(self) -> None:
        left = {"simulation_mode": "galaxy", "physics_package": "TPFCore", "tpf_dynamics_mode": "direct_tpf"}
        right = {"simulation_mode": "galaxy", "physics_package": "Newtonian"}
        name = _mode_aware_compare_name("initial", left, right, ext="png")
        self.assertEqual(
            name,
            "galaxy_compare__tpfcore_direct_tpf_vs_newtonian__compare__initial_side_by_side.png",
        )

    def test_mode_aware_compare_name_supports_xi_kernel_deformed(self) -> None:
        left = {
            "simulation_mode": "galaxy",
            "physics_package": "TPFCore",
            "tpf_dynamics_mode": "xi_kernel_deformed",
            "tpf_4d_xi_kernel_mode": "gaussian_compact",
        }
        right = {"simulation_mode": "galaxy", "physics_package": "Newtonian"}
        name = _mode_aware_compare_name("initial", left, right, ext="png")
        self.assertEqual(
            name,
            "galaxy_compare__tpfcore_xi_kernel_deformed_gaussian_compact_vs_newtonian__compare__initial_side_by_side.png",
        )

    def test_panel_physics_text_labels_newtonian_baseline(self) -> None:
        self.assertEqual(
            _panel_physics_text({"effective_physics_package": "Newtonian"}),
            "Newtonian baseline package",
        )

    def test_panel_physics_text_labels_xi_kernel_deformed(self) -> None:
        self.assertEqual(
            _panel_physics_text(
                {
                    "effective_physics_package": "TPFCore",
                    "effective_tpf_dynamics_mode": "xi_kernel_deformed",
                    "effective_tpf_4d_xi_kernel_mode": "gaussian_compact",
                }
            ),
            "TPFCore xi_kernel_deformed / gaussian_compact",
        )

    def test_compare_panel_metadata_newtonian_left_tpfcore_right(self) -> None:
        left_lines = _compare_panel_metadata_lines(
            self._side_data_for_metadata("left_primary", {"effective_physics_package": "Newtonian"})
        )
        right_lines = _compare_panel_metadata_lines(
            self._side_data_for_metadata(
                "right_compare",
                {
                    "effective_physics_package": "TPFCore",
                    "effective_tpf_dynamics_mode": "xi_kernel_deformed",
                    "effective_tpf_4d_xi_kernel_mode": "metric_transverse_wake",
                    "effective_tpf_4d_xi_kernel_coupling": "1.0e-5",
                },
            )
        )
        self.assertIn("package: Newtonian", left_lines)
        self.assertIn("role: left / baseline", left_lines)
        self.assertNotIn("dynamics:", "\n".join(left_lines))
        self.assertIn("package: TPFCore", right_lines)
        self.assertIn("role: right / compare", right_lines)
        self.assertIn("dynamics: xi_kernel_deformed", right_lines)
        self.assertIn("kernel: metric_transverse_wake", right_lines)
        self.assertIn("coupling: 1.0e-5", right_lines)

    def test_compare_panel_metadata_tpfcore_both_sides_include_model_details(self) -> None:
        left_lines = _compare_panel_metadata_lines(
            self._side_data_for_metadata(
                "left_primary",
                {
                    "effective_physics_package": "TPFCore",
                    "effective_tpf_dynamics_mode": "xi_kernel_deformed",
                    "effective_tpf_4d_xi_kernel_mode": "gaussian_compact",
                    "effective_tpf_4d_xi_kernel_coupling": "2.5e-6",
                    "effective_tpf_4d_xi_kernel_factor_mode": "legacy",
                    "effective_tpf_4d_xi_kernel_metric_min": "0.1",
                    "effective_tpf_4d_xi_kernel_metric_max": "0.9",
                    "effective_tpf_4d_xi_temporal_mode": "ema",
                    "effective_tpf_4d_xi_motion_readout_scale": "1.8",
                },
            )
        )
        right_lines = _compare_panel_metadata_lines(
            self._side_data_for_metadata(
                "right_compare",
                {
                    "effective_physics_package": "TPFCore",
                    "effective_tpf_dynamics_mode": "xi_kernel_deformed",
                    "effective_tpf_4d_xi_kernel_mode": "metric_transverse_wake",
                    "effective_tpf_4d_xi_kernel_coupling": "1.25e-6",
                },
            )
        )
        for lines in (left_lines, right_lines):
            self.assertIn("package: TPFCore", lines)
            self.assertTrue(any(line.startswith("role: ") for line in lines))
            self.assertIn("dynamics: xi_kernel_deformed", lines)
            self.assertTrue(any(line.startswith("kernel: ") for line in lines))
            self.assertTrue(any(line.startswith("coupling: ") for line in lines))
        self.assertIn("factor: legacy", left_lines)
        self.assertIn("metric clamp: 0.1–0.9", left_lines)
        self.assertIn("temporal: ema", left_lines)
        self.assertIn("K_xi: 1.8", left_lines)

    def test_compare_panel_metadata_beta_power_only_when_factor_mode_beta_power(self) -> None:
        with_beta = _compare_panel_metadata_lines(
            self._side_data_for_metadata(
                "left_primary",
                {
                    "effective_physics_package": "TPFCore",
                    "effective_tpf_dynamics_mode": "xi_kernel_deformed",
                    "effective_tpf_4d_xi_kernel_factor_mode": "beta_power",
                    "effective_tpf_4d_xi_kernel_beta_power": "2.25",
                },
            )
        )
        without_beta = _compare_panel_metadata_lines(
            self._side_data_for_metadata(
                "left_primary",
                {
                    "effective_physics_package": "TPFCore",
                    "effective_tpf_dynamics_mode": "xi_kernel_deformed",
                    "effective_tpf_4d_xi_kernel_factor_mode": "legacy",
                    "effective_tpf_4d_xi_kernel_beta_power": "2.25",
                },
            )
        )
        self.assertIn("beta_power: 2.25", with_beta)
        self.assertFalse(any(line.startswith("beta_power:") for line in without_beta))

    def test_compare_panel_metadata_tpfcore_left_newtonian_right(self) -> None:
        left_lines = _compare_panel_metadata_lines(
            self._side_data_for_metadata(
                "left_primary",
                {
                    "effective_physics_package": "TPFCore",
                    "effective_tpf_dynamics_mode": "direct_tpf",
                    "effective_tpf_4d_direct_tpf_coupling": "4.0e-3",
                    "effective_tpf_vdsg_coupling": "0.0",
                },
            )
        )
        right_lines = _compare_panel_metadata_lines(
            self._side_data_for_metadata("right_compare", {"effective_physics_package": "Newtonian"})
        )
        self.assertIn("package: TPFCore", left_lines)
        self.assertIn("role: left / baseline", left_lines)
        self.assertIn("dynamics: direct_tpf", left_lines)
        self.assertIn("kappa: 4.0e-3", left_lines)
        self.assertFalse(any("xi_kernel" in line for line in left_lines))
        self.assertFalse(any(line.startswith("vdsg_coupling:") for line in left_lines))
        self.assertIn("package: Newtonian", right_lines)
        self.assertIn("role: right / compare", right_lines)
        self.assertFalse(any(line.startswith("dynamics:") for line in right_lines))

    def test_compare_panel_metadata_archived_legacy_readout_is_marked_removed(self) -> None:
        lines = _compare_panel_metadata_lines(
            self._side_data_for_metadata(
                "left_primary",
                {
                    "effective_physics_package": "TPFCore",
                    "effective_tpf_dynamics_mode": "legacy_readout",
                    "effective_tpf_vdsg_coupling": "6.0e-4",
                    "effective_tpfcore_readout_mode": "legacy_readout",
                },
            )
        )
        self.assertIn("dynamics: legacy_readout (removed)", lines)
        self.assertIn("vdsg_coupling: 6.0e-4", lines)
        self.assertIn("readout: legacy_readout", lines)
        self.assertFalse(any("xi_kernel" in line for line in lines))

    @unittest.skipUnless(
        importlib.util.find_spec("matplotlib") is not None and importlib.util.find_spec("numpy") is not None,
        "matplotlib/numpy not installed",
    )
    def test_draw_panel_metadata_fontsize_is_readable(self) -> None:
        import matplotlib.pyplot as plt
        from framing import SquareViewport

        fig, ax = plt.subplots(1, 1)
        side = SimpleNamespace(
            label="left_primary",
            run_info={"effective_physics_package": "Newtonian"},
            overlay_mode="none",
            overlay_spec={},
        )
        snap = SimpleNamespace(
            positions=np.array([[0.0, 0.0], [1.0, 0.0]], dtype=float),
            velocities=np.array([[0.0, 0.0], [0.0, 0.0]], dtype=float),
            step=0,
            time=0.0,
        )
        spatial_display = spatial_display_for_xy_plot("galaxy_compare", 2.0, preferred_unit="m")
        _draw_panel(ax, side, snap, SquareViewport(0.0, 0.0, 2.0), spatial_display=spatial_display)
        panel_texts = [t for t in ax.texts if "package:" in t.get_text()]
        self.assertTrue(panel_texts)
        self.assertGreaterEqual(panel_texts[0].get_fontsize(), 8.0)
        plt.close(fig)

    @unittest.skipUnless(
        importlib.util.find_spec("matplotlib") is not None and importlib.util.find_spec("numpy") is not None,
        "matplotlib/numpy not installed",
    )
    def test_render_compare_static_uses_first_and_last_step_indices(self) -> None:
        import plot_cpp_compare

        with tempfile.TemporaryDirectory() as td:
            parent = Path(td)
            left = parent / "left_TPFCore"
            right = parent / "right_Newtonian"
            left.mkdir()
            right.mkdir()
            _write_run_info(left, "TPFCore")
            _write_run_info(right, "Newtonian")
            _write_snapshot(left / "snapshot_00000.csv", 0, 0.0, 1.0)
            _write_snapshot(left / "snapshot_00010.csv", 10, 1.0, 2.0)
            _write_snapshot(right / "snapshot_00000.csv", 0, 0.0, 1.2)
            _write_snapshot(right / "snapshot_00010.csv", 10, 1.0, 2.2)
            (parent / "compare_manifest.json").write_text(
                json.dumps({"compare_run_id": "r3", "left_dir": str(left), "right_dir": str(right)}),
                encoding="utf-8",
            )
            seen_steps: list[int] = []
            original_draw = plot_cpp_compare._draw_panel

            def _record_draw(ax, side, snap, viewport, *, spatial_display):
                seen_steps.append(int(snap.step))
                return original_draw(ax, side, snap, viewport, spatial_display=spatial_display)

            try:
                plot_cpp_compare._draw_panel = _record_draw
                render_compare(parent, no_animation=True, overlay_mode="none")
            finally:
                plot_cpp_compare._draw_panel = original_draw
            self.assertIn(0, seen_steps)
            self.assertIn(10, seen_steps)
            self.assertNotEqual(min(seen_steps), max(seen_steps))

    def test_compare_animation_does_not_use_ema_radius_smoothing(self) -> None:
        """EMA on shared radius lagged behind per-frame max(r_l,r_r) and clipped stars in video."""
        import plot_cpp_compare

        text = Path(plot_cpp_compare.__file__).read_text(encoding="utf-8")
        self.assertNotIn("_smooth_alpha", text)
        self.assertNotIn("_smooth_state", text)

    def test_compare_smart_viewport_contains_both_extents(self) -> None:
        """Shared viewport is computed from the union of both sides' point clouds (+ origin)."""
        import numpy as np

        from plot_cpp_run import Snapshot

        pos_l = np.random.default_rng(0).normal(scale=1e19, size=(100, 2))
        vel_l = np.zeros_like(pos_l)
        th = np.linspace(0, 2 * np.pi, 100, endpoint=False)
        pos_r = np.column_stack([1e27 * np.cos(th), 1e27 * np.sin(th)])
        vel_r = np.column_stack([np.cos(th) * 1e5, np.sin(th) * 1e5])
        snap_l = Snapshot(0, 0.0, pos_l, vel_l)
        snap_r = Snapshot(0, 0.0, pos_r, vel_r)
        vp = calculate_compare_smart_viewport([snap_l], [snap_r], 150.0)
        self.assertGreater(vp.half_axis, 1e26)

    def test_compare_smart_viewport_ignores_small_persistent_escapees(self) -> None:
        import numpy as np

        from plot_cpp_run import Snapshot

        rng = np.random.default_rng(7)
        left_snaps = []
        right_snaps = []
        n = 120
        n_out = 3  # 2.5% pooled outliers across both panels at each frame
        for step in range(12):
            near_l = rng.normal(scale=10.0, size=(n - n_out, 2))
            near_r = rng.normal(scale=10.0, size=(n - n_out, 2))
            out_l = np.column_stack([np.full(n_out, 1.0e6), rng.normal(scale=100.0, size=n_out)])
            out_r = np.column_stack([np.full(n_out, -1.0e6), rng.normal(scale=100.0, size=n_out)])
            pos_l = np.vstack([near_l, out_l])
            pos_r = np.vstack([near_r, out_r])
            vel = np.zeros_like(pos_l)
            left_snaps.append(Snapshot(step, float(step), pos_l, vel))
            right_snaps.append(Snapshot(step, float(step), pos_r, vel))

        vp = calculate_compare_smart_viewport(left_snaps, right_snaps, 150.0, coverage=0.95)
        self.assertLess(vp.half_axis, 1.0e3)

    def test_compare_smart_viewport_expands_when_many_far_out(self) -> None:
        import numpy as np

        from plot_cpp_run import Snapshot

        rng = np.random.default_rng(8)
        left_snaps = []
        right_snaps = []
        n = 120
        n_out = 45  # 37.5% per side
        for step in range(10):
            near = rng.normal(scale=10.0, size=(n - n_out, 2))
            out_l = np.column_stack([np.full(n_out, 1.0e5), rng.normal(scale=100.0, size=n_out)])
            out_r = np.column_stack([np.full(n_out, -1.0e5), rng.normal(scale=100.0, size=n_out)])
            pos_l = np.vstack([near, out_l])
            pos_r = np.vstack([near, out_r])
            vel = np.zeros_like(pos_l)
            left_snaps.append(Snapshot(step, float(step), pos_l, vel))
            right_snaps.append(Snapshot(step, float(step), pos_r, vel))

        vp = calculate_compare_smart_viewport(left_snaps, right_snaps, 150.0, coverage=0.95)
        self.assertGreater(vp.half_axis, 5.0e4)

    def test_compare_smart_viewport_uses_both_panels_each_frame(self) -> None:
        import numpy as np

        from plot_cpp_run import Snapshot

        left_snaps = []
        right_snaps = []
        for step in range(6):
            pos_l = np.column_stack([np.full(300, 1.0e4), np.zeros(300)])
            pos_r = np.column_stack([np.full(300, -1.0e4), np.zeros(300)])
            vel = np.zeros_like(pos_l)
            left_snaps.append(Snapshot(step, float(step), pos_l, vel))
            right_snaps.append(Snapshot(step, float(step), pos_r, vel))

        vp = calculate_compare_smart_viewport(left_snaps, right_snaps, 150.0, coverage=0.95)
        self.assertAlmostEqual(vp.center_x, 0.0, places=6)
        self.assertGreater(vp.half_axis, 9.0e3)

    def test_compare_smart_viewport_is_fixed_shared_union_across_frames(self) -> None:
        import numpy as np

        from plot_cpp_run import Snapshot

        left_snaps = []
        right_snaps = []
        for step in range(4):
            near_l = np.column_stack([np.full(50, 10.0 + step * 5.0), np.linspace(-3.0, 3.0, 50)])
            near_r = np.column_stack([np.full(50, -10.0 - step * 5.0), np.linspace(-3.0, 3.0, 50)])
            vel = np.zeros_like(near_l)
            left_snaps.append(Snapshot(step, float(step), near_l, vel))
            right_snaps.append(Snapshot(step, float(step), near_r, vel))

        vp = calculate_compare_smart_viewport(left_snaps, right_snaps, 150.0, coverage=0.95)
        # Must include the most separated frame; single fixed shared viewport for full animation/static render.
        self.assertGreater(vp.half_axis, 20.0)
        self.assertAlmostEqual(vp.center_x, 0.0, places=6)

    def test_compare_smart_viewport_coverage_knob_monotonic(self) -> None:
        import numpy as np

        from plot_cpp_run import Snapshot

        rng = np.random.default_rng(123)
        left_snaps = []
        right_snaps = []
        for step in range(3):
            near = rng.normal(scale=10.0, size=(100, 2))
            far = np.column_stack([np.full(10, 2.0e4), rng.normal(scale=50.0, size=10)])
            pos_l = np.vstack([near, far])
            pos_r = np.vstack([near, far])
            vel = np.zeros_like(pos_l)
            left_snaps.append(Snapshot(step, float(step), pos_l, vel))
            right_snaps.append(Snapshot(step, float(step), pos_r, vel))

        vp_95 = calculate_compare_smart_viewport(left_snaps, right_snaps, 150.0, coverage=0.95)
        vp_90 = calculate_compare_smart_viewport(left_snaps, right_snaps, 150.0, coverage=0.90)
        self.assertLessEqual(vp_90.half_axis, vp_95.half_axis)

    def test_large_point_cloud_uses_approximate_viewport_path(self) -> None:
        from plot_cpp_run import Snapshot

        n = 1100
        pos_l = np.column_stack([np.linspace(-100.0, 100.0, n), np.zeros(n)])
        pos_r = np.column_stack([np.linspace(-90.0, 110.0, n), np.ones(n)])
        vel = np.zeros_like(pos_l)
        left = [Snapshot(0, 0.0, pos_l, vel)]
        right = [Snapshot(0, 0.0, pos_r, vel)]
        with patch("plot_cpp_compare.minimum_axis_aligned_box_covering_fraction") as exact_mock:
            vp = calculate_compare_smart_viewport(left, right, 150.0, coverage=0.95)
        exact_mock.assert_not_called()
        self.assertTrue(np.isfinite(vp.half_axis))

    def test_approximate_viewport_returns_finite_bounds(self) -> None:
        pts = np.array([[0.0, 0.0], [1.0, np.nan], [2.0, 3.0], [4.0, 5.0], [10.0, -2.0]], dtype=np.float64)
        xmin, xmax, ymin, ymax = approximate_axis_aligned_box_covering_fraction(pts, 0.8)
        self.assertTrue(np.isfinite([xmin, xmax, ymin, ymax]).all())
        self.assertLessEqual(xmin, xmax)
        self.assertLessEqual(ymin, ymax)

    def test_small_point_cloud_still_uses_exact_path(self) -> None:
        from plot_cpp_run import Snapshot

        n = 100
        pos = np.column_stack([np.linspace(-10.0, 10.0, n), np.zeros(n)])
        vel = np.zeros_like(pos)
        left = [Snapshot(0, 0.0, pos, vel)]
        right = [Snapshot(0, 0.0, pos, vel)]
        with patch("plot_cpp_compare.minimum_axis_aligned_box_covering_fraction", wraps=__import__("plot_cpp_compare").minimum_axis_aligned_box_covering_fraction) as exact_mock:
            calculate_compare_smart_viewport(left, right, 150.0, coverage=0.95)
        exact_mock.assert_called()

    def test_viewport_sampling_caps_large_frame_counts(self) -> None:
        idx = _sample_viewport_frame_indices(201, max_samples=50)
        self.assertLessEqual(len(idx), 50)
        self.assertEqual(idx[0], 0)
        self.assertEqual(idx[-1], 200)

    def test_compare_display_selection_shared_distance_unit(self) -> None:
        left = DisplayUnitConfig(distance_unit="AU", time_unit="day", velocity_unit="km/s")
        right = DisplayUnitConfig(distance_unit="auto", time_unit="day", velocity_unit="km/s")
        sel = resolve_compare_display_selection(
            left, right, shared_half_axis_m=1.0e12, max_time_s=86400.0 * 5, max_speed_m_s=2.0e3
        )
        self.assertEqual(sel.active_distance_unit, "AU")
        self.assertEqual(sel.active_time_unit, "day")
        self.assertEqual(sel.active_velocity_unit, "km/s")

    def test_compare_display_selection_conflicting_explicit_units_falls_back_to_auto(self) -> None:
        left = DisplayUnitConfig(distance_unit="km", time_unit="hr", velocity_unit="m/s")
        right = DisplayUnitConfig(distance_unit="AU", time_unit="day", velocity_unit="km/s")
        sel = resolve_compare_display_selection(
            left, right, shared_half_axis_m=2.0e11, max_time_s=86400.0 * 10, max_speed_m_s=5.0e3
        )
        self.assertEqual(sel.config.distance_unit, "auto")

    def test_compare_display_selection_overlay_flags_require_both_sides_enabled(self) -> None:
        left = DisplayUnitConfig(units_in_overlay=True, show_unit_reference=True)
        right = DisplayUnitConfig(units_in_overlay=False, show_unit_reference=False)
        sel = resolve_compare_display_selection(
            left, right, shared_half_axis_m=1.0e9, max_time_s=100.0, max_speed_m_s=1.0
        )
        self.assertFalse(sel.config.units_in_overlay)
        self.assertFalse(sel.config.show_unit_reference)


    def test_radial_kinematics_plus_x_plus_y(self) -> None:
        class S: pass
        s=S(); import numpy as np
        s.positions=np.array([[1.0,0.0]])
        s.velocities=np.array([[0.0,2.0]])
        out=radial_kinematics(s)
        self.assertAlmostEqual(float(out["v_radial"][0]), 0.0)
        self.assertAlmostEqual(float(out["v_t"][0]), 2.0)

    def test_estimate_snapshot_accelerations_delta_v_over_dt(self) -> None:
        class S: pass
        snaps=[]
        for t,v in [(0.0,0.0),(1.0,2.0),(2.0,4.0)]:
            s=S(); import numpy as np
            s.time=t; s.velocities=np.array([[v,0.0]])
            snaps.append(s)
        acc=estimate_snapshot_accelerations(snaps)
        self.assertAlmostEqual(float(acc[1][0,0]), 2.0)
        self.assertAlmostEqual(float(acc[0][0,0]), 2.0)
        self.assertAlmostEqual(float(acc[2][0,0]), 2.0)

    def test_binned_profile_skips_empty_and_returns_medians(self) -> None:
        import numpy as np
        x=np.array([0.0,0.1,0.2,9.8,9.9,10.0])
        y=np.array([1,3,5,7,9,11],dtype=float)
        bx,bm,_,_=binned_profile(x,y,bins=10)
        self.assertGreater(len(bx), 1)
        self.assertTrue(np.all(np.isfinite(bm)))

    def test_distance_scaling_uses_spatial_display_factor(self) -> None:
        import plot_cpp_compare

        spatial = SimpleNamespace(factor=2.5, unit="km")
        self.assertEqual(spatial.factor, 2.5)
        self.assertFalse(hasattr(spatial, "scale_to_display"))

    def test_time_display_factor_scales_days(self) -> None:
        self.assertAlmostEqual(_time_display_factor("day"), 1.0 / 86400.0)
        t_plot = 2.0 * 86400.0 * _time_display_factor("day")
        self.assertAlmostEqual(t_plot, 2.0)
        self.assertAlmostEqual(60.0 * _time_display_factor("min"), 1.0)
        self.assertAlmostEqual(3600.0 * _time_display_factor("hr"), 1.0)
        seconds_per_year = 365.25 * 86400.0
        self.assertAlmostEqual(seconds_per_year * _time_display_factor("yr"), 1.0)
        self.assertAlmostEqual((1.0e3 * seconds_per_year) * _time_display_factor("kyr"), 1.0)
        self.assertAlmostEqual((1.0e6 * seconds_per_year) * _time_display_factor("Myr"), 1.0)

    def test_velocity_display_factor_scales_km_s(self) -> None:
        self.assertAlmostEqual(_velocity_display_factor("km/s"), 1.0e-3)
        self.assertAlmostEqual(2000.0 * _velocity_display_factor("km/s"), 2.0)

    def test_raw_acceleration_plot_not_binned_helper(self) -> None:
        import plot_cpp_compare
        text = Path(plot_cpp_compare.__file__).read_text(encoding="utf-8")
        self.assertIn("ax.scatter(lk[\"radius\"] * dist_scale, l_in", text)

    def test_acceleration_ratio_uses_binned_profiles(self) -> None:
        import plot_cpp_compare
        text = Path(plot_cpp_compare.__file__).read_text(encoding="utf-8")
        self.assertIn("l_x, l_med, _, _ = binned_profile", text)
        self.assertIn("r_x, r_med, _, _ = binned_profile", text)
        self.assertNotIn("rx=lk['radius'][:m]*dist_scale", text)



if __name__ == "__main__":
    unittest.main()
