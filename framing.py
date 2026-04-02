"""
Geometric square viewports for 2D scatter / animation (display layer only).

Single rule: axis-aligned quantile bounds on x and y, square envelope, margin.
SI meters in / SI meters out; multiply by display_units.SpatialDisplay.factor at draw time.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SquareViewport:
    """Square window in data (SI) space: max(|x - cx|, |y - cy|) <= half_axis."""

    center_x: float
    center_y: float
    half_axis: float

    @staticmethod
    def centered_origin(half_axis: float) -> SquareViewport:
        """Legacy behavior: square centered on the origin (e.g. compare scripts, fixed radius)."""
        h = float(half_axis)
        return SquareViewport(0.0, 0.0, h)


def compute_square_viewport(
    xy: np.ndarray,
    *,
    trim_fraction: float = 0.01,
    margin: float = 1.15,
    min_half_axis: float = 1.0,
    fallback_half_axis: float = 150.0,
) -> SquareViewport:
    """
    Robust square framing from point cloud (x, y).

    - x_lo, x_hi = quantile(x, p), quantile(x, 1-p) with p = trim_fraction (p=0 => min/max).
    - Same for y. Center = midpoint; half_axis = 0.5 * max(span_x, span_y) * margin.
    - Degenerate spans fall back to min/max on all finite points, then to fallback_half_axis.
    """
    fb = float(max(fallback_half_axis, min_half_axis, 1.0))
    if xy is None or getattr(xy, "size", 0) == 0 or len(xy) == 0:
        return SquareViewport(0.0, 0.0, fb)

    x = np.asarray(xy[:, 0], dtype=np.float64)
    y = np.asarray(xy[:, 1], dtype=np.float64)
    m = np.isfinite(x) & np.isfinite(y)
    if not np.any(m):
        return SquareViewport(0.0, 0.0, fb)

    xf = x[m]
    yf = y[m]
    p = float(trim_fraction)
    if p <= 0.0:
        x_lo, x_hi = float(np.min(xf)), float(np.max(xf))
        y_lo, y_hi = float(np.min(yf)), float(np.max(yf))
    else:
        p = min(p, 0.49)
        x_lo = float(np.quantile(xf, p))
        x_hi = float(np.quantile(xf, 1.0 - p))
        y_lo = float(np.quantile(yf, p))
        y_hi = float(np.quantile(yf, 1.0 - p))

    if x_lo > x_hi:
        x_lo, x_hi = x_hi, x_lo
    if y_lo > y_hi:
        y_lo, y_hi = y_hi, y_lo

    wx = x_hi - x_lo
    wy = y_hi - y_lo
    if not np.isfinite(wx) or wx <= 0.0:
        x_lo, x_hi = float(np.min(xf)), float(np.max(xf))
        wx = x_hi - x_lo
    if not np.isfinite(wy) or wy <= 0.0:
        y_lo, y_hi = float(np.min(yf)), float(np.max(yf))
        wy = y_hi - y_lo

    cx = 0.5 * (x_lo + x_hi)
    cy = 0.5 * (y_lo + y_hi)
    span = max(wx, wy)
    if not np.isfinite(span) or span <= 0.0:
        return SquareViewport(cx, cy, fb)

    half = 0.5 * float(span) * float(margin)
    half = max(half, 0.5 * float(min_half_axis), 1e-30)
    return SquareViewport(cx, cy, float(half))


def stack_particle_xy(
    snapshots: list,
    extra_xy: np.ndarray | None = None,
) -> np.ndarray:
    """Vertically stack snapshot.positions (each n×2); optionally append extra rows (e.g. fixed BH)."""
    if not snapshots:
        base = np.zeros((0, 2), dtype=np.float64)
    else:
        base = np.vstack([np.asarray(s.positions, dtype=np.float64) for s in snapshots])
    if extra_xy is not None and getattr(extra_xy, "size", 0) > 0:
        ex = np.asarray(extra_xy, dtype=np.float64).reshape(-1, 2)
        if base.size == 0:
            return ex
        return np.vstack([base, ex])
    return base


def global_viewport_from_snapshots(
    snapshots: list,
    *,
    extra_xy: np.ndarray | None = None,
    trim_fraction: float = 0.01,
    margin: float = 1.15,
    min_half_axis: float = 1.0,
    fallback_half_axis: float = 150.0,
) -> SquareViewport:
    """One viewport from all particles in all given snapshots (+ optional extra points)."""
    xy = stack_particle_xy(snapshots, extra_xy=extra_xy)
    return compute_square_viewport(
        xy,
        trim_fraction=trim_fraction,
        margin=margin,
        min_half_axis=min_half_axis,
        fallback_half_axis=fallback_half_axis,
    )


class DynamicViewportSmoother:
    """
    Exponential smoothing on viewport center and log(half_axis).
    First sample passes through unchanged (after margin / trim already applied upstream).
    """

    def __init__(self, alpha: float = 0.12, min_half_axis: float = 1.0) -> None:
        self.alpha = float(alpha)
        self.min_half_axis = float(min_half_axis)
        self._cx: float | None = None
        self._cy: float | None = None
        self._log_h: float | None = None

    def step(self, raw: SquareViewport) -> SquareViewport:
        a = self.alpha
        rx, ry, rh = float(raw.center_x), float(raw.center_y), float(raw.half_axis)
        rh = max(rh, self.min_half_axis, 1e-300)
        if self._cx is None:
            self._cx, self._cy = rx, ry
            self._log_h = np.log(rh)
            return SquareViewport(rx, ry, rh)

        self._cx = a * rx + (1.0 - a) * self._cx
        self._cy = a * ry + (1.0 - a) * self._cy
        lh = np.log(rh)
        self._log_h = a * lh + (1.0 - a) * self._log_h
        h = float(np.exp(self._log_h))
        h = max(h, self.min_half_axis)
        return SquareViewport(self._cx, self._cy, h)


def viewport_window_snapshots(snapshots: list, frame_index: int, window: int) -> list:
    """Inclusive window [frame_index - window + 1, frame_index] clipped to list bounds."""
    if window < 1:
        window = 1
    lo = max(0, int(frame_index) - int(window) + 1)
    return snapshots[lo : int(frame_index) + 1]
