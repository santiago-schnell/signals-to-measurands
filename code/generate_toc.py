#!/usr/bin/env python3
"""
Graphical abstract (ACS table-of-contents graphic) for:
"From Signals to Measurands: A Measurement-Science Roadmap for
Reproducible Analytical Biochemistry"

Rendered at the ACS TOC size, 3.25 x 1.75 in, exact canvas. The title and
tagline font sizes are auto-fitted to the canvas width and each box is sized to
its text, so nothing is clipped at the edges (the prior version clipped the
title to "rom signals ...").

Outputs: toc-graphic.{eps,pdf,png}
"""
from __future__ import annotations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "dejavusans",
})
DPI = 600
DARK = "#2e333a"
MUTED = "#5f656e"
LINE = "#3b4148"
BOX_BLUE = "#e9f1f7"
BOX_GREEN = "#e9f6ee"
BOX_PINK = "#fbecec"
PINK_EDGE = "#b5505a"

FIG_W, FIG_H = 3.25, 1.75


def make_toc():
    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_axis_off()
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    figpx = FIG_W * fig.dpi

    def wfrac(s, fs, weight="normal"):
        t = fig.text(0.5, 0.5, s, fontsize=fs, fontweight=weight)
        w = t.get_window_extent(renderer=rend).width / figpx
        t.remove()
        return w

    def fit_fs(s, target, start, weight="normal", lo=6.0):
        fs = start
        while fs > lo and wfrac(s, fs, weight) > target:
            fs -= 0.2
        return fs

    def line_wfrac(label, fs, weight):
        return max(wfrac(ln, fs, weight) for ln in label.split("\n"))

    # ---- title (auto-fit so it never clips) ----
    title = "From signals to reproducible measurands"
    t_fs = fit_fs(title, 0.92, 10.5, "bold")
    ax.text(0.5, 0.865, title, ha="center", va="center", fontsize=t_fs,
            fontweight="bold", color=DARK)

    # ---- three-box chain, each box sized to its text, centered as a row ----
    boxes = [
        ("measured\nsignal", BOX_BLUE, DARK, 9.0, "bold"),
        ("calibration\n& inversion", BOX_GREEN, DARK, 9.0, "bold"),
        (r"$q \pm U(q)$", BOX_PINK, PINK_EDGE, 10.0, "bold"),
    ]
    hpad = 0.022
    widths = [line_wfrac(lbl, fs, wt) + 2 * hpad for lbl, _, _, fs, wt in boxes]
    gap = 0.052
    total = sum(widths) + gap * (len(boxes) - 1)
    x = (1.0 - total) / 2.0
    cy, h = 0.505, 0.30
    centers = []
    for (lbl, fc, ec, fs, wt), bw in zip(boxes, widths):
        ax.add_patch(FancyBboxPatch(
            (x, cy - h / 2), bw, h, boxstyle="round,pad=0,rounding_size=0.03",
            mutation_aspect=FIG_H / FIG_W, linewidth=1.9 if ec == PINK_EDGE else 1.7,
            edgecolor=ec, facecolor=fc, transform=ax.transAxes, clip_on=False, zorder=3))
        ax.text(x + bw / 2, cy, lbl, ha="center", va="center", fontsize=fs,
                fontweight=wt, color=DARK, linespacing=1.05, zorder=4)
        centers.append((x, x + bw))
        x += bw + gap

    for (l, r0), (l2, r1) in zip(centers[:-1], centers[1:]):
        ax.add_patch(FancyArrowPatch((r0 + 0.006, cy), (l2 - 0.006, cy),
                                     arrowstyle="-|>", mutation_scale=11,
                                     linewidth=1.9, color=LINE, transform=ax.transAxes,
                                     clip_on=False, zorder=5, shrinkA=0, shrinkB=0))

    # ---- tagline (auto-fit) ----
    tag = "identifiability  \u2022  uncertainty  \u2022  independent routes"
    g_fs = fit_fs(tag, 0.94, 8.2)
    ax.text(0.5, 0.145, tag, ha="center", va="center", fontsize=g_fs, color=MUTED)

    return fig


fig = make_toc()
for ext in ("eps", "pdf", "png"):
    fig.savefig(f"toc-graphic.{ext}", dpi=DPI, facecolor="white")
plt.close(fig)
print("TOC graphic written.")
