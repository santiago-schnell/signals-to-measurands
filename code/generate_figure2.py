#!/usr/bin/env python3
"""Generate the signal-to-measurand chain (Figure 2)."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
OUT.mkdir(parents=True, exist_ok=True)

DARK = "#343a40"
TEXT = "#252a31"
MUTED = "#58616c"


BOXPAD = 0.010  # boxstyle padding: the visible border is this far outside (x, y, w, h)


def rounded_box(ax, x, y, w, h, facecolor, edgecolor=DARK, lw=2.0):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad={BOXPAD},rounding_size=0.020",
        linewidth=lw, edgecolor=edgecolor, facecolor=facecolor,
        zorder=2,
    )
    ax.add_patch(patch)
    return patch


def gap_arrow(ax, x0, y0, x1, y1):
    """Draw a clearly visible arrow entirely within the gap between boxes."""
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1),
        arrowstyle="-|>", mutation_scale=23,
        linewidth=2.2, color=DARK,
        shrinkA=0, shrinkB=0, zorder=4,
        connectionstyle="arc3,rad=0",
    ))


# Print-size typography: see generate_figure1.py.  The binding width is the
# SAVED BOUNDING BOX (~13.2 in here), not FIG_W, and the manuscript sets
# \linewidth = 6.524 in, so the scale factor is 0.493 and every source font
# must be at least 8/0.493 = 16.2 pt to clear the journal's 8 pt floor.  The
# smallest label below is 16.5 pt, i.e. 8.1 pt printed.  Box widths were
# enlarged to accommodate the larger labels.
FIG_W, FIG_H = 13.4, 7.2
SUB = 16.5       # -> 8.6 pt printed
BOXTITLE = 16.5  # -> 8.6 pt printed (bold)
CAPTION = 17.0   # -> 8.9 pt printed
HEAD2 = 20.0     # -> 10.4 pt printed (bold)
HEAD1 = 24.0     # -> 12.5 pt printed (bold)


def build():
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H), dpi=300)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(
        0.5, 0.945, "From a measured signal to a reported biochemical result",
        ha="center", va="center", fontsize=HEAD1, fontweight="bold", color=TEXT
    )

    # Upper measurement chain. The horizontal gaps are deliberately wide enough
    # to contain the arrows; no connector enters a box.  The outermost box edges
    # sit at xs[0] - BOXPAD and xs[-1] + w + BOXPAD, so MARGIN must exceed BOXPAD
    # (plus half the border linewidth) or the first and last borders fall outside
    # the axes and are clipped.
    MARGIN = 0.016
    w, h, y = 0.200, 0.285, 0.560
    gap = (1.0 - 2 * MARGIN - 4 * w) / 3.0
    xs = [MARGIN + i * (w + gap) for i in range(4)]
    assert xs[0] - BOXPAD > 0.0 and xs[-1] + w + BOXPAD < 1.0
    fills = ["#e8f0f7", "#f7f1e6", "#e8f4ec", "#f6e8e8"]
    edges = [DARK, DARK, DARK, "#b44f55"]
    titles = [
        "Material and assay",
        "Measured signal",
        "Measurement model",
        "Reported result",
    ]
    subtitles = [
        "system,\nconditions, and\nreactive fraction",
        "raw observations\nand calibration\ndata",
        "biochemical,\nobservation, and\nstatistical models",
        "estimate with\nan uncertainty\nstatement",
    ]

    for x, fc, ec, title, subtitle in zip(xs, fills, edges, titles, subtitles):
        rounded_box(ax, x, y, w, h, fc, ec, lw=2.4 if ec != DARK else 2.0)
        ax.text(x + w / 2, y + 0.212, title, ha="center", va="center",
                fontsize=BOXTITLE, fontweight="bold", color=TEXT, zorder=5)
        ax.text(x + w / 2, y + 0.098, subtitle, ha="center", va="center",
                fontsize=SUB, color=MUTED, linespacing=1.16, zorder=5)

    for i in range(3):
        gap_arrow(ax, xs[i] + w + BOXPAD, y + h / 2,
                  xs[i + 1] - BOXPAD, y + h / 2)

    ax.text(
        0.5, 0.525,
        "The protocol specifies what is done; the models specify how the result is inferred.",
        ha="center", va="center", fontsize=CAPTION, fontstyle="italic", color=DARK
    )

    ax.text(
        0.5, 0.435, "Three checks probe progressively more of the measurement chain",
        ha="center", va="center", fontsize=HEAD2, fontweight="bold", color=TEXT
    )

    # Lower reproducibility checks. Arrows sit in separate gaps and are larger
    # than in the previous version, so they remain visible at print size.
    y2, h2, w2 = 0.160, 0.200, 0.250
    x2s = [0.020, 0.375, 0.730]
    fills2 = ["#eaf1f7", "#dbe8f3", "#c8dbec"]
    titles2 = ["Repeatability", "Transfer replication",
               "Independent-method\nreproducibility"]
    subs2 = ["repeat under closely\nmatched local conditions",
             "same specified procedure\nin another laboratory",
             "different route to\nthe same measurand"]

    for x, fc, title, subtitle in zip(x2s, fills2, titles2, subs2):
        rounded_box(ax, x, y2, w2, h2, fc)
        ax.text(x + w2 / 2, y2 + 0.133, title, ha="center", va="center",
                fontsize=BOXTITLE, fontweight="bold", color=TEXT,
                linespacing=1.04, zorder=5)
        ax.text(x + w2 / 2, y2 + 0.050, subtitle, ha="center", va="center",
                fontsize=SUB, color=MUTED, linespacing=1.12, zorder=5)

    for i in range(2):
        gap_arrow(ax, x2s[i] + w2 + BOXPAD, y2 + h2 / 2,
                  x2s[i + 1] - BOXPAD, y2 + h2 / 2)

    ax.text(
        0.5, 0.090,
        "greater variation in laboratory, method, and assumptions  $\\longrightarrow$  stronger invariance claim",
        ha="center", va="center", fontsize=CAPTION, color=DARK
    )

    fig.subplots_adjust(left=0.012, right=0.988, bottom=0.035, top=0.985)
    return fig


def main() -> None:
    fig = build()
    for ext in ("pdf", "png", "eps"):
        kwargs = {"bbox_inches": "tight", "pad_inches": 0.08, "facecolor": "white"}
        if ext == "png":
            kwargs["dpi"] = 300
        fig.savefig(OUT / f"figure2-measurement-chain.{ext}", **kwargs)
    plt.close(fig)


if __name__ == "__main__":
    main()
