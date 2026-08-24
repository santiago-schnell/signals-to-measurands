#!/usr/bin/env python3
"""
Figure 3 (rebuilt) for:
"From Signals to Measurands: A Measurement-Science Roadmap for
Reproducible Analytical Biochemistry"

Panel A: ligand depletion shifts the apparent binding curve.  When the tracer is
not negligible the curve is NOT a hyperbola in [R]_T -- it is the smaller root of
the mass-balance quadratic -- but it still reaches half tracer occupancy at the
displaced point [R]_{T,1/2} = K_D + [L]_T/2.  A midpoint-based or simple
hyperbolic reading therefore assigns an apparent constant displaced by [L]_T/2.
Panel B: normalized magnitude of the local sensitivity of the ideal
low-depletion bound fraction to ln K_D, exposing slope-amplitude and
saturation confounding away from [R]_T ~ K_D.

Fixes: interior gridlines removed; the "apparent midpoint" text and formula
moved into clear top-left space with a short leader to the true-K_D crossing;
midpoint markers and drop-lines make the rightward shift legible; the legend is
uncluttered and clear of the drop-lines; Panel B's three regime labels are off the
curve.

Outputs: figure3-ligand-binding.{eps,pdf,png}
"""
from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

# Print-size typography: see generate_figure1.py.  The binding width is the
# SAVED BOUNDING BOX (~13.25 in here), not FIG_W, and the manuscript sets
# \linewidth = 6.524 in, so the scale factor is 0.492 and every source font
# must be at least 8/0.492 = 16.2 pt.  The smallest label below is 17.0 pt,
# i.e. 8.4 pt printed.
FIG_W, FIG_H = 13.33, 6.8
ANNOT = 17.0                   # -> 8.9 pt printed
TICKS = 17.5                   # -> 9.2 pt printed
BASE = 18.0                    # -> 9.5 pt printed
AXLAB = 19.5                   # -> 10.2 pt printed
PANEL = 20.0                   # -> 10.5 pt printed (bold)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "dejavusans",
    "font.size": BASE,
    "axes.titlesize": BASE,
    "axes.labelsize": AXLAB,
    "xtick.labelsize": TICKS,
    "ytick.labelsize": TICKS,
    "axes.linewidth": 1.15,
    "xtick.major.width": 1.05,
    "ytick.major.width": 1.05,
})
DPI = 600
BLUE = "#1f77b4"
ORANGE = "#ff7f0e"
GREEN = "#2ca02c"
DARK = "#2e333a"
SHADE = "#e6f1f7"
FAINT = "#c7ccd2"


def bound_fraction_of_tracer(u, LT):
    """Fraction of labeled tracer bound in a depletion-aware titration.
    u = [R]_T/K_D (titrant, varying), LT = [L]_T/K_D (total tracer, fixed), K_D=1.
    RL = 1/2[(u+LT+1) - sqrt((u+LT+1)^2 - 4 u LT)]; returns RL/LT.
    Midpoint (RL=LT/2) at u = 1 + LT/2, i.e. [R]_{T,1/2}=K_D+[L]_T/2."""
    u = np.asarray(u, dtype=float)
    b = u + LT + 1.0
    disc = np.clip(b * b - 4.0 * u * LT, 0.0, None)
    RL = 0.5 * (b - np.sqrt(disc))
    return RL / LT


def make_figure3():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(FIG_W, FIG_H), dpi=DPI,
                                   gridspec_kw={"width_ratios": [1.0, 1.0]})

    u = np.logspace(-2.2, 2.2, 1400)

    # ---------- Panel A: depletion bias ----------
    ideal = u / (1.0 + u)
    dep1 = bound_fraction_of_tracer(u, 1.0)
    dep4 = bound_fraction_of_tracer(u, 4.0)

    ax1.axhline(0.5, color=FAINT, linewidth=1.2, zorder=1)
    ax1.plot(u, ideal, color=BLUE, linewidth=3.2, zorder=3, label=r"$[\mathrm{L}]_{\mathrm{T}} \ll K_D$")
    ax1.plot(u, dep1, color=ORANGE, linewidth=2.9, zorder=3, label=r"$[\mathrm{L}]_{\mathrm{T}} = K_D$")
    ax1.plot(u, dep4, color=GREEN, linewidth=2.9, linestyle=(0, (5, 2)), zorder=3,
             label=r"$[\mathrm{L}]_{\mathrm{T}} = 4\,K_D$")

    for xm, col in [(1.0, BLUE), (1.5, ORANGE), (3.0, GREEN)]:
        ax1.plot([xm, xm], [0.0, 0.5], color=col, linestyle="--", linewidth=1.3, zorder=2)
        ax1.plot([xm], [0.5], "o", color=col, ms=7, mec="white", mew=1.0, zorder=6)

    ax1.text(0.0115, 1.06, "apparent midpoint\nshifts right:\n"
             r"$[\mathrm{R}]_{\mathrm{T},1/2} = K_D + [\mathrm{L}]_{\mathrm{T}}/2$", transform=ax1.transData,
             ha="left", va="top", fontsize=ANNOT, color=DARK, linespacing=1.18)
    ax1.add_patch(FancyArrowPatch((0.7, 0.80), (1.0, 0.515), arrowstyle="-|>",
                                  mutation_scale=13, linewidth=1.4, color=DARK,
                                  zorder=5, shrinkA=0, shrinkB=0))

    ax1.set_xscale("log")
    ax1.set_xlim(9e-3, 1.5e2)
    ax1.set_ylim(-0.02, 1.12)
    ax1.set_xlabel(r"scaled titrant, $[\mathrm{R}]_{\mathrm{T}}/K_D$")
    ax1.set_ylabel("fraction of tracer bound")
    ax1.set_title("A  Depletion shifts the midpoint",
                  fontsize=PANEL, loc="left", fontweight="bold", pad=8)
    ax1.legend(frameon=False, loc="lower right", fontsize=ANNOT, handlelength=1.8,
               borderaxespad=0.9)

    # ---------- Panel B: normalized local sensitivity to ln K_D ----------
    # For f=u/(1+u), |partial f/partial ln K_D| = u/(1+u)^2 at fixed
    # [R]_T. Multiplication by four normalizes the maximum magnitude to one.
    sens = 4.0 * u / (1.0 + u) ** 2
    ax2.axvspan(0.3, 3.0, color=SHADE, zorder=0)
    ax2.axhline(0.0, color=FAINT, linewidth=1.0, zorder=1)
    #ax2.axvline(1.0, color=BLUE, linestyle="--", linewidth=1.4, zorder=2)
    # Draw the vertical line in two segments to leave a gap for the text
    ax2.vlines(x=[1.0, 1.0],           # The X-coordinates for both line segments
               ymin=[0.0, 0.36],       # Gap left for the annotation beneath the peak
               ymax=[0.06, 1.02],
               color=BLUE, 
               linestyle="--", 
               linewidth=1.4, 
               zorder=2)
    ax2.plot(u, sens, color=BLUE, linewidth=3.2, zorder=3)

    ax2.text(0.0115, 0.80, "slope-amplitude\nconfounding", ha="left", va="center",
             fontsize=ANNOT, color=DARK, linespacing=1.12)
    ax2.text(1.0, 0.20, "maximum local\ninformation near\n" r"$[\mathrm{R}]_{\mathrm{T}} \approx K_D$",
             ha="center", va="center", fontsize=ANNOT, color=DARK, linespacing=1.12)
    ax2.text(140.0, 0.80, "saturation\nconfounding", ha="right", va="center",
             fontsize=ANNOT, color=DARK, linespacing=1.12)

    ax2.set_xscale("log")
    ax2.set_xlim(9e-3, 1.5e2)
    ax2.set_ylim(-0.03, 1.07)
    ax2.set_xlabel(r"scaled titrant, $[\mathrm{R}]_{\mathrm{T}}/K_D$")
    ax2.set_ylabel(r"normalized $|\partial f/\partial\ln K_D|$")
    ax2.set_title("B  Normalized sensitivity to $K_D$",
                  fontsize=PANEL, loc="left", fontweight="bold", pad=8)

    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.135, top=0.90, wspace=0.215)
    return fig


def main() -> None:
    fig = make_figure3()
    for ext in ("eps", "pdf", "png"):
        fig.savefig(f"figure3-ligand-binding.{ext}",
                    facecolor="white", bbox_inches="tight", pad_inches=0.12,
                    **({"dpi": DPI} if ext == "png" else {}))
    plt.close(fig)
    print("Figure 3 (rebuilt) written.")


if __name__ == "__main__":
    main()
