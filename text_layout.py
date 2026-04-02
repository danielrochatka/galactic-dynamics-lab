"""Helpers to keep plot text readable and inside static image bounds."""

from __future__ import annotations

from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.text import Text


def _fit_text_to_pixel_width(
    fig: Figure,
    text_artist: Text,
    max_width_px: float,
    *,
    min_fontsize: float = 6.0,
) -> None:
    """Shrink text artist font size until it fits in `max_width_px`."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bbox = text_artist.get_window_extent(renderer=renderer)
    fs = float(text_artist.get_fontsize())
    while bbox.width > max_width_px and fs > min_fontsize:
        fs = max(min_fontsize, fs - 0.5)
        text_artist.set_fontsize(fs)
        fig.canvas.draw()
        bbox = text_artist.get_window_extent(renderer=renderer)


def set_fitted_title(
    ax: Axes,
    title: str,
    *,
    color: str | None = None,
    fontsize: float = 12.0,
    min_fontsize: float = 6.0,
    margin_px: float = 12.0,
) -> Text:
    """Set an axis title and scale font size down to fit plot width."""
    title_kwargs = {"fontsize": fontsize}
    if color is not None:
        title_kwargs["color"] = color
    title_artist = ax.set_title(title, **title_kwargs)
    fig = ax.figure
    fig.canvas.draw()
    ax_bbox = ax.get_window_extent(renderer=fig.canvas.get_renderer())
    _fit_text_to_pixel_width(
        fig,
        title_artist,
        max_width_px=max(1.0, ax_bbox.width - margin_px),
        min_fontsize=min_fontsize,
    )
    return title_artist


def add_fitted_footer(
    fig: Figure,
    text: str,
    *,
    x: float = 0.5,
    y: float = 0.02,
    ha: str = "center",
    va: str = "bottom",
    fontsize: float = 8.0,
    min_fontsize: float = 5.0,
    color: str | None = None,
    margin_px: float = 12.0,
) -> Text:
    """Add figure-level footer text and shrink it to stay inside the canvas."""
    text_kwargs = {"ha": ha, "va": va, "fontsize": fontsize}
    if color is not None:
        text_kwargs["color"] = color
    txt = fig.text(x, y, text, **text_kwargs)
    fig.canvas.draw()
    fig_w_px = fig.get_figwidth() * fig.dpi
    _fit_text_to_pixel_width(
        fig,
        txt,
        max_width_px=max(1.0, fig_w_px - margin_px),
        min_fontsize=min_fontsize,
    )
    return txt
