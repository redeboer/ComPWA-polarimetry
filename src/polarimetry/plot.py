"""Helper functions for `matplotlib`."""

from __future__ import annotations

import shutil
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt

if TYPE_CHECKING:
    from matplotlib.artist import Artist
    from matplotlib.axes import Axes
    from matplotlib.contour import QuadContourSet


def add_watermark(
    ax: Axes,
    x: float = 0.03,
    y: float = 0.03,
    fontsize: int | None = None,
    **kwargs,
) -> None:
    text = "LHCb\n" + R"$1.7\mathrm{~fb}^{-1}$"
    ax.text(x, y, text, size=fontsize, transform=ax.transAxes, **kwargs)


def get_contour_line(contour_set: QuadContourSet) -> Artist:
    (line_collection, *_), _ = contour_set.legend_elements()
    return line_collection


def use_mpl_latex_fonts(reset_mpl: bool = True) -> None:
    # cspell:ignore dejavusans fontset mathtext usetex
    if not _is_latex_allowed():
        return
    if reset_mpl:
        _wake_up_matplotlib()
    plt.rc("font", family="serif", serif="Helvetica")
    plt.rc("mathtext", fontset="dejavusans")
    plt.rc("text", usetex=True)


def _is_latex_allowed() -> bool:
    return shutil.which("latex") is not None


def _wake_up_matplotlib() -> None:
    """Somehow `plt.rc` does not work if a figure hasn't been created before..."""
    plt.figure()
    plt.close()
