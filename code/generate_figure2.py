#!/usr/bin/env python3
"""
Figure 2 (rebuilt) for:
"From Signals to Measurands: A Measurement-Science Roadmap for
Reproducible Analytical Biochemistry"

Fixes the publication-quality issues in the prior version: interior gridlines
removed (they crossed the curves and added clutter); panel letters masked so no
curve runs through them; Panel B uses color-coded labels placed directly on the
two sensitivity curves rather than a legend the blue line crossed; and the
regime annotations are moved into clear space away from the curves. The
meaningful shaded regime bands and the K_M reference line are kept.

Outputs: figure2-mm-progress-curves.{eps,pdf,png}
"""
from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "dejavusans",
    "font.size": 16,
    "axes.titlesize": 16,
    "axes.labelsize": 17,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
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
ZERO = "#9aa0a8"


def mm_progress_solution(tau, s0):
    """Solve ds/dtau = -s/(1+s), s(0)=s0 (s=[S]/K_M, tau=Vt/K_M).
    Implicit solution s + log s = s0 + log s0 - tau; Newton in z=log s."""
    tau = np.asarray(tau, dtype=float)
    c = s0 + np.log(s0) - tau
    z = np.full_like(tau, np.log(max(s0, 1e-12)))
    for _ in range(100):
        ez = np.exp(z)
        step = (ez + z - c) / (ez + 1.0)
        z -= step
        if np.nanmax(np.abs(step)) < 1e-12:
            break
    return np.exp(z)


def make_figure2():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.33, 6.8), dpi=DPI,
                                   gridspec_kw={"width_ratios": [1.0, 1.0]})

    # ---------- Panel A: progress trajectories ----------
    tau = np.linspace(0, 12, 600)
    curves = [(0.2, BLUE, r"$[\mathrm{S}]_0/K_M = 0.2$"),
              (3.0, ORANGE, r"$[\mathrm{S}]_0/K_M = 3$"),
              (10.0, GREEN, r"$[\mathrm{S}]_0/K_M = 10$")]
    ax1.axhspan(0.8, 1.2, color=SHADE, zorder=0)
    for s0, c, lab in curves:
        ax1.plot(tau, mm_progress_solution(tau, s0), color=c, linewidth=3.0,
                 label=lab, zorder=3)
    ax1.axhline(1.0, color=DARK, linestyle="--", linewidth=1.5, zorder=2)
    ax1.text(11.8, 1.40, r"$[\mathrm{S}] = K_M$", ha="right", va="bottom",
             fontsize=12.5, color=DARK)
    ax1.set_xlim(0, 12)
    ax1.set_ylim(0, 10.5)
    ax1.set_xlabel(r"scaled time, $Vt/K_M$")
    ax1.set_ylabel(r"substrate, $[\mathrm{S}]/K_M$")
    ax1.set_title("A  Trajectories traverse different regimes",
                  fontsize=15.5, loc="left", fontweight="bold", pad=8)
    ax1.legend(frameon=False, loc="upper right", fontsize=13.5, handlelength=1.7,
               borderaxespad=0.8)

    # ---------- Panel B: local sensitivities ----------
    s = np.logspace(-2, 2, 800)
    for a, b in [(1e-2, 1e-1), (0.7, 1.5), (10, 1e2)]:
        ax2.axvspan(a, b, color=SHADE, zorder=0)
    ax2.axhline(0.0, color=ZERO, linewidth=1.0, zorder=1)
    ax2.plot(s, np.ones_like(s), color=BLUE, linewidth=3.0, zorder=3)
    ax2.plot(s, -1.0 / (1.0 + s), color=ORANGE, linewidth=3.0, zorder=3)
    ax2.set_xscale("log")
    ax2.set_xlim(1e-2, 1e2)
    ax2.set_ylim(-1.18, 1.18)
    ax2.set_xlabel(r"instantaneous $[\mathrm{S}]/K_M$")
    ax2.set_ylabel("log-sensitivity")
    ax2.set_title("B  Sensitivities expose parameter confounding",
                  fontsize=15.5, loc="left", fontweight="bold", pad=8)

    # color-coded labels placed directly on the curves
    ax2.text(0.40, 1.035, r"$\partial\ln v\,/\,\partial\ln V$", color=BLUE,
             fontsize=14.5, ha="center", va="bottom")
    ax2.text(7.0, -0.34, r"$\partial\ln v\,/\,\partial\ln K_M$", color=ORANGE,
             fontsize=14.5, ha="center", va="center")

    # regime annotations, placed clear of both curves
    ax2.text(0.033, 0.50, "low $[\\mathrm{S}]$:\nonly $V/K_M$", fontsize=12.5, ha="center",
             va="center", color=DARK, linespacing=1.12)
    ax2.text(1.02, 0.60, "near $K_M$:\njoint information", fontsize=12.5,
             ha="center", va="center", color=DARK, linespacing=1.12)
    ax2.text(33.0, 0.50, "high $[\\mathrm{S}]$:\nmostly $V$", fontsize=12.5, ha="center",
             va="center", color=DARK, linespacing=1.12)

    fig.subplots_adjust(left=0.065, right=0.985, bottom=0.135, top=0.90, wspace=0.20)
    return fig


fig = make_figure2()
for ext in ("eps", "pdf", "png"):
    fig.savefig(f"figure2-mm-progress-curves.{ext}",
                facecolor="white", bbox_inches="tight", pad_inches=0.12,
                **({"dpi": DPI} if ext == "png" else {}))
plt.close(fig)
print("Figure 2 (rebuilt) written.")
