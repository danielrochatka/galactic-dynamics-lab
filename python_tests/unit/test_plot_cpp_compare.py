from __future__ import annotations

import json
import sys
import tempfile
import unittest
import importlib.util
from pathlib import Path
from unittest.mock import patch

# Repo root (parent of python_tests/)
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from plot_cpp_compare import (
    resolve_compare_display_selection,
    _resolve_side_run_dir,
    _mode_aware_compare_name,
    calculate_compare_smart_viewport,
    matched_steps_strict,
    render_compare,
)
from display_units import DisplayUnitConfig


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

    def test_resolve_side_run_dir_cpp_sim_relative_manifest(self) -> None:
        """Manifest stores paths relative to cpp_sim; compare_parent is .../outputs/RUN."""
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            run = root / "cpp_sim" / "outputs" / "RUN1"
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

    def test_mode_aware_compare_name_supports_legacy_run_info(self) -> None:
        left = {"simulation_mode": "galaxy", "physics_package": "TPFCore", "tpf_dynamics_mode": "direct_tpf"}
        right = {"simulation_mode": "galaxy", "physics_package": "Newtonian"}
        name = _mode_aware_compare_name("initial", left, right, ext="png")
        self.assertEqual(
            name,
            "galaxy_compare__tpfcore_direct_tpf_vs_newtonian__compare__initial_side_by_side.png",
        )

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


if __name__ == "__main__":
    unittest.main()
