#!/usr/bin/env python3
"""
Figure 1 (rebuilt) for:
"From Signals to Measurands: A Measurement-Science Roadmap for
Reproducible Analytical Biochemistry"

Signal-to-measurand chain. Fixes the publication-quality issues in the prior
version: real gaps between boxes so the connecting arrows are fully visible,
no dotted leader lines crossing into boxes, stage descriptors anchored cleanly
above each stage, and an explicit bracket-and-arrow linking the chain to the
three reproducibility checks, which are drawn as a stringency ladder beneath it.

Outputs: figure1-measurement-chain.{eps,pdf,png}
"""
from __future__ import annotations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

# ---------- Styling (matches generate_figures.py) ----------
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "dejavusans",
    "axes.linewidth": 1.15,
})
DPI = 600
DARK = "#2e333a"
MUTED = "#6b727c"
LINE = "#3b4148"
# stage fills (calm, distinct)
BOX_BLUE = "#e9f1f7"
BOX_YELLOW = "#f8f3e6"
BOX_GREEN = "#e9f6ee"
BOX_PURPLE = "#f2ecf8"
BOX_PALE_BLUE = "#eaf2fb"
BOX_PINK = "#fbecec"
PINK_EDGE = "#b5505a"
# check-ladder tints (light -> deep)
CHK1 = "#eef3f8"
CHK2 = "#d9e6f1"
CHK3 = "#bed3e8"


def rbox(ax, x, y, w, h, *, fc, ec=DARK, lw=1.7, r=0.018):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
        linewidth=lw, edgecolor=ec, facecolor=fc,
        mutation_aspect=0.55, transform=ax.transAxes, clip_on=False, zorder=3))


def arrow(ax, x0, y0, x1, y1, *, lw=2.2, ms=17, color=LINE):
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=ms,
        linewidth=lw, color=color, transform=ax.transAxes,
        clip_on=False, zorder=5,
        shrinkA=0, shrinkB=0))


def make_figure1():
    fig = plt.figure(figsize=(13.4, 7.9), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_axis_off()

    ax.text(0.5, 0.940, "Reproducibility is engineered along the whole measurement chain",
            ha="center", va="center", fontsize=21, fontweight="bold", color=DARK)

    # ----- chain of stages -----
    stages = [
        ("Material\nstate",            "identity, purity,\nactive fraction",   BOX_BLUE),
        ("Raw\nsignal",                "blanks, standards,\ntraceability",     BOX_YELLOW),
        ("Calibration\n& correction",  "assumptions,\nregime of validity",     BOX_GREEN),
        ("Biochemical\nmodel",         "rate law\nor isotherm",                BOX_PURPLE),
        ("Statistical\ninverse problem","noise model,\nidentifiability",       BOX_PALE_BLUE),
        ("Measurand\n$q \\pm U(q)$",   "compatibility\nstatement",             BOX_PINK),
    ]
    n = len(stages)
    w, h, gap = 0.1280, 0.150, 0.0344
    x0 = 0.030
    ytop, ybot = 0.760, 0.760 - h
    cy = (ytop + ybot) / 2
    xs = [x0 + i * (w + gap) for i in range(n)]

    for i, (name, desc, fc) in enumerate(stages):
        x = xs[i]
        last = (i == n - 1)
        rbox(ax, x, ybot, w, h, fc=fc,
             ec=PINK_EDGE if last else DARK, lw=2.4 if last else 1.7)
        ax.text(x + w / 2, cy, name, ha="center", va="center",
                fontsize=12.5, fontweight="bold", color=DARK, linespacing=1.06,
                zorder=6)
        # descriptor anchored just above its box (no leader line)
        ax.text(x + w / 2, ytop + 0.014, desc, ha="center", va="bottom",
                fontsize=9.0, color=MUTED, linespacing=1.05, zorder=6)

    for i in range(n - 1):
        arrow(ax, xs[i] + w + 0.004, cy, xs[i + 1] - 0.004, cy)

    # ----- bracket grouping the chain as one procedure -----
    xL, xR = xs[0], xs[-1] + w
    yb = 0.598
    ax.plot([xL, xR], [yb, yb], color=LINE, lw=1.7, zorder=2)
    ax.plot([xL, xL], [yb, ybot], color=LINE, lw=1.7, zorder=2)
    ax.plot([xR, xR], [yb, ybot], color=LINE, lw=1.7, zorder=2)
    ax.text(0.5, 0.555,
            "A single pass through this chain is one measurement of the measurand.",
            ha="center", va="center", fontsize=12.5, color=DARK, style="italic")

    # ----- three reproducibility checks as a stringency ladder -----
    ax.text(0.5, 0.475, "Three reproducibility checks apply to the whole chain",
            ha="center", va="center", fontsize=14.5, fontweight="bold", color=DARK)

    checks = [
        ("Repeatability", "same operator, instrument, materials, and calibration \u2014 repeated", CHK1),
        ("Transfer replication", "the same specified procedure carried to another laboratory", CHK2),
        ("Independent-method reproducibility", "a different measurement route to the same measurand", CHK3),
    ]
    # size the boxes to their content (narrow), then centre them under the chain
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    figpx = fig.get_figwidth() * fig.dpi

    def wfrac(s, fs, weight="normal"):
        t = fig.text(0.5, 0.5, s, fontsize=fs, fontweight=weight)
        wpx = t.get_window_extent(renderer=rend).width
        t.remove()
        return wpx / figpx

    maxw = max(max(wfrac(nm, 13, "bold"), wfrac(ds, 10.5)) for nm, ds, _ in checks)
    lpad, rpad = 0.026, 0.026
    bw, bh = maxw + lpad + rpad, 0.094
    bx = (1.0 - bw) / 2.0
    centers = [0.388, 0.266, 0.144]
    for (name, desc, fc), c in zip(checks, centers):
        rbox(ax, bx, c - bh / 2, bw, bh, fc=fc, lw=1.6, r=0.012)
        ax.text(bx + lpad, c + 0.020, name, ha="left", va="center",
                fontsize=13, fontweight="bold", color=DARK, zorder=6)
        ax.text(bx + lpad, c - 0.021, desc, ha="left", va="center",
                fontsize=10.5, color="#3f444b", zorder=6)

    # stringency arrow placed just left of the centered boxes
    gx = bx - 0.045
    arrow(ax, gx, centers[0] + bh / 2 - 0.004, gx, centers[-1] - bh / 2 + 0.004,
          lw=2.0, ms=15, color=DARK)
    ax.text(gx - 0.034, (centers[0] + centers[-1]) / 2,
            "increasing variation,\nstronger claim", ha="center", va="center",
            rotation=90, fontsize=10, color=DARK, linespacing=1.05)

    return fig


fig = make_figure1()
for ext in ("eps", "pdf", "png"):
    fig.savefig(f"figure1-measurement-chain.{ext}",
                facecolor="white", bbox_inches="tight", pad_inches=0.14,
                **({"dpi": DPI} if ext == "png" else {}))
plt.close(fig)
print("Figure 1 (rebuilt) written.")
