"""Shared figure style and saving, so every figure follows the project's figure conventions (README).

Every figure is saved as PNG (200 dpi) in figures/<experiment>/ (no PDF copies; PNGs are enough).
Colours: one fixed order for all figures (a colour-blind-checked reference palette); always take them from COLORS
in this order, never pick hues by hand. Text stays black/grey, never in the line colour.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # draw to files only; no window needed
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent

COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK = "#0b0b0b"  # main text and the "total" line
INK_2 = "#52514e"  # secondary text
MUTED = "#898781"  # reference lines, streamlines
GRID = "#e1e0d9"  # hairline gridlines
AXIS = "#c3c2b7"  # axis lines


def use_style():
    """Quiet, print-friendly defaults: thin lines, hairline grid, no top/right frame."""
    plt.rcParams.update({
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10, "legend.fontsize": 9,
        "axes.edgecolor": AXIS, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": INK_2, "ytick.color": INK_2,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": "-",
        "lines.linewidth": 1.8, "lines.solid_capstyle": "round", "lines.dash_capstyle": "round",
        "legend.frameon": False,
    })


def plot_branch(ax, x, y, stable, color, label=None, lw=2.0, **kwargs):
    """Draw a branch of equilibria: SOLID where stable, DASHED where unstable (figure conventions).

    x, y, stable are equal-length sequences; the line is cut into pieces where the stability changes, and each
    piece is extended by one point so that the pieces join up.
    """
    x, y, stable = list(x), list(y), list(stable)
    start = 0
    first = True
    for i in range(1, len(x) + 1):
        if i == len(x) or stable[i] != stable[start]:
            end = min(i + 1, len(x))
            ax.plot(x[start:end], y[start:end], color=color, lw=lw, ls="-" if stable[start] else (0, (4, 2.5)),
                    label=label if first else None, **kwargs)
            first = False
            start = i


def ship_outline(state, s):
    """Corners of a simple top-view hull shape (units of L, earth axes: X upstream, Y to port), for drawing only.

    Rectangle with a pointed bow; the width is the real beam B/L.
    """
    import numpy as np  # local import keeps this module light for scripts that only need the style

    b = 0.5 * s.hull.B / s.hull.L
    body = np.array([[-0.5, -b], [0.36, -b], [0.5, 0.0], [0.36, b], [-0.5, b], [-0.5, -b]])
    psi = state[1]
    e_x = np.array([np.cos(psi), np.sin(psi)])
    e_y = np.array([-np.sin(psi), np.cos(psi)])
    gamma = state[0]
    O = s.L_H_L * np.array([-np.cos(gamma), np.sin(gamma)]) - s.x_A_L * e_x  # midship (see model_rigid.positions)
    return O + body[:, :1] * e_x + body[:, 1:] * e_y


def draw_mooring(ax, state, s, color, alpha=1.0, lw=1.4, fill=True):
    """Top view of buoy, hawser and ship at one instant (earth axes: X upstream to the right, Y to port = up).

    The current flows from right to left; the ship trails to the left of the buoy, bow pointing at the buoy.
    """
    import numpy as np

    hull = ship_outline(state, s)
    A = s.L_H_L * np.array([-np.cos(state[0]), np.sin(state[0])])
    ax.plot([0, A[0]], [0, A[1]], color=INK_2, lw=1.0, alpha=alpha)  # hawser
    if fill:
        ax.fill(hull[:, 0], hull[:, 1], color=color, alpha=0.18 * alpha, lw=0)
    ax.plot(hull[:, 0], hull[:, 1], color=color, lw=lw, alpha=alpha)
    ax.plot(0, 0, "o", ms=7, color=INK, zorder=6)  # buoy


def current_arrow(ax, x, y, length=0.3):
    """Arrow showing the current direction (right to left in our top views)."""
    ax.annotate("", xy=(x - length, y), xytext=(x, y),
                arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.4))
    ax.text(x - 0.5 * length, y + 0.04, "current", ha="center", va="bottom", color=INK_2, fontsize=8.5)


def save(fig, experiment, name):
    """Save `fig` as figures/<experiment>/<name>.png, then close it."""
    folder = ROOT / "figures" / experiment
    folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / f"{name}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"saved figures/{experiment}/{name}.png")
