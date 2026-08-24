#!/usr/bin/env python3
"""
Worked exemplar for:
"From Signals to Measurands: A Measurement-Science Roadmap for
Reproducible Analytical Biochemistry"

A complete Michaelis-Menten signal-to-measurand record.

The same enzyme preparation and the same forward-model fit quality
(reduced chi-square ~ 1, structureless residuals) are analysed under two
designs:

  Design I  (uninformative): [S]_0 = 0.2 K_M, trajectory stays [S] << K_M.
  Design II (informative):   [S]_0 = 4 K_M, trajectory traverses [S] ~ K_M.

The initial substrate concentration is a separately calibrated input on the same
Beer-Lambert scale as the progress-curve observations; it is not estimated from
the first noisy progress-curve datum.

Forward adequacy holds in BOTH. Separate recovery of every parameter does not: a profile
deviance analysis shows K_M is practically non-identifiable in Design I (no finite
upper bound) but well constrained in Design II, while the specificity-
determining combination V/K_M is tightly determined in both.

The uncertainty budget then shows that the curve-fit (Type A) standard error
on k_cat is a small fraction of the combined standard uncertainty, which is
dominated by the Type B uncertainty in initial active-site concentration; and that
the molar absorption coefficient cancels in k_cat/K_M but not in k_cat or K_M.

Outputs: figure1-mm-case-study.{eps,pdf,png} and a printed numerical summary.
"""
from __future__ import annotations

import numpy as np
from pathlib import Path
from scipy.special import lambertw
from scipy.optimize import least_squares, minimize_scalar
from scipy.stats import chi2

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- Manuscript styling (matches generate_figures.py) ----------
# ---------------------------------------------------------------------------
# Print-size typography.
#
# ACS Measurement Science Au requires that no text in a figure fall below
# 8 pt at the size at which the figure is printed.  The effective printed
# size of a label is
#
#       effective_pt = source_pt x (printed_width / FIG_W)
#
# With FIG_W = 13.6 in and a printed width of 7.0 in (full ACS page width)
# the scale factor is 0.515, so every source font must be at least
# 8 / 0.515 = 15.6 pt.  ANNOT below is the floor used for all in-panel
# annotations and legends; the remaining sizes are set above it.
# ---------------------------------------------------------------------------
# The relevant width is the SAVED BOUNDING BOX, not FIG_W: bbox_inches="tight"
# trims and pads, so figure1 saves at about 14.77 in even though the canvas is
# 13.6 in.  The manuscript sets \linewidth = 6.524 in (letterpaper, 1 in
# margins), so the scale factor is 6.524/14.77 = 0.442 and every source font
# must be at least 8/0.442 = 18.1 pt to clear the journal's 8 pt floor.
FIG_W, FIG_H = 13.6, 10.2      # source canvas (inches); saved bbox ~14.77 in
PRINT_W = 6.524                # manuscript \linewidth (inches)
ANNOT = 18.5                   # -> 8.2 pt printed
TICKS = 18.5                   # -> 8.2 pt printed
BASE = 19.0                    # -> 8.4 pt printed
AXLAB = 20.5                   # -> 9.1 pt printed
PANEL = 21.0                   # -> 9.3 pt printed (bold)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "dejavusans",
    "font.size": BASE,
    "axes.titlesize": BASE,
    "axes.labelsize": AXLAB,
    "xtick.labelsize": TICKS,
    "ytick.labelsize": TICKS,
    "legend.fontsize": ANNOT,
    "axes.linewidth": 1.15,
    "xtick.major.width": 1.05,
    "ytick.major.width": 1.05,
})
DPI = 600
BLUE = "#1f77b4"      # Design II (informative)
ORANGE = "#ff7f0e"    # Design I  (uninformative)
GRID = "#d9d9d9"
SHADE = "#e6f1f7"
DARK = "#2e333a"
GREEN = "#2ca02c"
ZERO = "#9aa0a8"

# ---------- Ground truth and instrument ----------
KM_TRUE = 50.0e-6           # M
KCAT_TRUE = 50.0            # 1/s
E_TOTAL = 5.0e-9            # M  (nominal total enzyme)
ACTIVE_FRACTION = 0.80      # active-site titration result
E0_TRUE = ACTIVE_FRACTION * E_TOTAL     # 4.0 nM
V_TRUE = KCAT_TRUE * E0_TRUE            # 2.0e-7 M/s
EPS_EXT = 6220.0            # 1/(M cm)  NADH at 340 nm
PATH = 1.0                 # cm
SIGMA_A = 0.003            # AU, additive Gaussian on absorbance
SIGMA_S = SIGMA_A / (EPS_EXT * PATH)          # concentration noise (homoscedastic)

# Type B relative standard uncertainties (illustrative, realistic)
U_REL_EPS = 0.02           # molar absorption coefficient
U_REL_E0 = 0.08          # initial active-site concentration (active-site titration)

RNG = np.random.default_rng(20240530)


def substrate(t, V, KM, S0):
    """Exact MM progress curve via the Lambert-W (Schnell-Mendoza) form:
    [S](t) = K_M * W0( ([S]_0/K_M) exp(([S]_0 - V t)/K_M) )."""
    arg = (S0 / KM) * np.exp((S0 - V * t) / KM)
    return KM * np.real(lambertw(arg, k=0))


def simulate(S0, t):
    """Return (t, S_obs) for a design: noisy substrate concentrations from
    absorbance, S_obs = A_obs/(eps*path), A_obs = eps*path*S_true + N(0,sigma_A).

    Calibration convention. Every concentration in this exemplar is placed on
    one Beer-Lambert scale by division by eps*path. The initial substrate
    concentration is a separately calibrated input on that same scale,
    [S]_0 = A0_cal/(eps*path); it is not estimated from the first noisy
    progress-curve observation. A relative change in eps therefore rescales the
    observations and [S]_0 together, so V and K_M each scale as eps^-1 and
    V/K_M is invariant. The `epsilon_sensitivity` check verifies this covariance
    numerically. Repeatability uncertainty in A0_cal and uncertainty in the 1 cm
    path length are neglected in this illustrative budget."""
    S_true = substrate(t, V_TRUE, KM_TRUE, S0)
    A_obs = EPS_EXT * PATH * S_true + RNG.normal(0.0, SIGMA_A, size=t.shape)
    return A_obs / (EPS_EXT * PATH)


# V and K_M have different units and poorly matched numerical scales in the
# chosen parameterization, so a direct fit can stop short of the best numerical
# minimum. Optimizing log V and log K_M improves scaling, enforces
# positivity without bounds, and makes the tolerances scale-free.  Design I is
# the case that matters: its likelihood is nearly flat along V/K_M, so a direct
# fit can sit ~0.05% above the minimum while the point estimate is ~13% away.
LOG_LB = np.log((1e-9, 1e-7))
LOG_UB = np.log((1e-3, 1e-2))
MULTISTART = ((0.3, 0.2), (0.6, 0.5), (1.0, 1.0), (1.6, 2.0), (3.0, 5.0))


def resid_full(logp, t, S_obs, S0):
    V, KM = np.exp(logp)
    return (substrate(t, V, KM, S0) - S_obs) / SIGMA_S


def fit_full(t, S_obs, S0, p0=(V_TRUE, KM_TRUE)):
    """Multistart (V, KM) fit in log parameters; returns popt, pcov, SSR.

    The returned solution is the lowest residual sum of squares found by the
    stated bounded multistart procedure. It is used as the numerical reference
    SSR* for the profile deviance; no claim of analytically proved global
    optimality is made."""
    best = None
    for fV, fK in MULTISTART:
        start = np.log(np.clip((p0[0] * fV, p0[1] * fK),
                               np.exp(LOG_LB), np.exp(LOG_UB)))
        sol = least_squares(resid_full, x0=start, args=(t, S_obs, S0),
                            bounds=(LOG_LB, LOG_UB), method="trf",
                            xtol=1e-15, ftol=1e-15, gtol=1e-15)
        val = float(np.sum(sol.fun**2))
        if best is None or val < best[0]:
            best = (val, sol)
    SSRw, sol = best
    SSR = SSRw * SIGMA_S**2
    popt = np.exp(sol.x)
    # The Jacobian is with respect to log parameters; convert to natural
    # parameters by the chain rule, dtheta = theta * dlog(theta).
    J = sol.jac / popt          # columns rescaled: J_nat = J_log / theta
    JTJ = J.T @ J
    try:
        pcov = np.linalg.inv(JTJ)
    except np.linalg.LinAlgError:
        pcov = np.full((2, 2), np.nan)
    return popt, pcov, SSR


def _bounded_log_profile(objective, lower, upper):
    """Deterministic one-dimensional minimization in a bounded log parameter.

    A scalar bounded minimizer is substantially faster and more robust than a
    generic least-squares iteration for the 18,000 conditional fits in the
    sampling-schedule study. The tolerance is much tighter than required by the
    0.11-decade profile grid, while the finite iteration cap prevents a single
    nearly flat replicate from stalling the validation workflow.
    """
    result = minimize_scalar(
        objective, bounds=(np.log(lower), np.log(upper)), method="bounded",
        options={"xatol": 1e-10, "maxiter": 250})
    if not result.success or not np.isfinite(result.fun):
        raise RuntimeError(f"Conditional profile minimization failed: {result.message}")
    return float(result.fun) * SIGMA_S**2


def ssr_at_fixed_KM(KM, t, S_obs, S0):
    """Profile over V at fixed KM (bounded scalar minimization of log V)."""
    def objective(logV):
        residual = (substrate(t, np.exp(logV), KM, S0) - S_obs) / SIGMA_S
        return float(np.dot(residual, residual))
    return _bounded_log_profile(objective, 1e-9, 1e-3)


def ssr_at_fixed_ratio(psi, t, S_obs, S0):
    """Profile over KM at fixed psi = V/KM (bounded minimization of log KM)."""
    def objective(logKM):
        KM = np.exp(logKM)
        residual = (substrate(t, psi * KM, KM, S0) - S_obs) / SIGMA_S
        return float(np.dot(residual, residual))
    return _bounded_log_profile(objective, 1e-7, 1e-2)


def profile_curve(grid, ssr_fn, t, S_obs, S0, ssr_star):
    """Profile deviance D_p = [SSR_min(theta) - SSR*] / sigma^2.

    Referenced to the best-fit numerical minimum SSR* identified by the
    multistart procedure -- not to the lowest point of the finite grid, which
    would depend on where the grid happens to fall."""
    vals = np.array([ssr_fn(g, t, S_obs, S0) for g in grid])
    return np.maximum(vals - ssr_star, 0.0) / SIGMA_S**2


def ci_from_profile(grid, T, thresh):
    inside = grid[T <= thresh]
    if inside.size == 0:
        return None, None
    return inside.min(), inside.max()


# Nominal 95% likelihood-ratio threshold for one profiled parameter.
THRESH = chi2.ppf(0.95, df=1)


def limiting_exponential_fit(t, S_obs, S0):
    """Best fit of the K_M -> infinity limiting model S0*exp(-psi*t).

    Returns (psi_hat, SSR). This supplies the asymptotic profile-deviance check
    needed to establish whether the upper profile region truly remains open,
    rather than merely reaching the end of a finite K_M grid.
    """
    sol = least_squares(
        lambda lp: (S0 * np.exp(-np.exp(lp[0]) * t) - S_obs) / SIGMA_S,
        x0=[np.log(V_TRUE / KM_TRUE)],
        bounds=([np.log(1e-7)], [np.log(1.0)]),
        method="trf", xtol=1e-15, ftol=1e-15, gtol=1e-15)
    return float(np.exp(sol.x[0])), float(np.sum(sol.fun**2)) * SIGMA_S**2


def main() -> None:
    global RNG
    RNG = np.random.default_rng(20240530)
    # ---------- Designs ----------
    designs = {
        "I":  dict(label="Design I  ($0.2\\,K_M$)", color=ORANGE,
                   S0=0.2 * KM_TRUE, t=np.linspace(0.0, 1000.0, 50)),
        "II": dict(label="Design II ($4\\,K_M$)", color=BLUE,
                   S0=4.0 * KM_TRUE, t=np.linspace(0.0, 1800.0, 60)),
    }

    KM_GRID = np.unique(np.concatenate([
        np.logspace(np.log10(1e-6), np.log10(5e-2), 200),   # broad: 1 uM .. 50 mM
        np.linspace(5e-6, 90e-6, 240),                      # dense across the profile wells
    ]))
    PSI_GRID = np.logspace(np.log10(1.0e-3), np.log10(1.6e-2), 160)   # around V/KM=4e-3

    results = {}
    for key, d in designs.items():
        t = d["t"]
        S_obs = simulate(d["S0"], t)
        popt, pcov, SSR = fit_full(t, S_obs, d["S0"])
        Vhat, KMhat = popt
        dof = t.size - 2
        redchi = (SSR / SIGMA_S**2) / dof

        # Type A relative SEs
        seV = np.sqrt(pcov[0, 0]); seKM = np.sqrt(pcov[1, 1])
        uA_V = seV / Vhat
        uA_KM = seKM / KMhat
        # delta method for V/KM
        g = np.array([1.0 / KMhat, -Vhat / KMhat**2])      # grad of V/KM wrt (V,KM)
        var_ratio = g @ pcov @ g
        ratio = Vhat / KMhat
        uA_ratio = np.sqrt(var_ratio) / ratio

        # profiles, referenced to the best-fit multistart minimum SSR*
        T_KM = profile_curve(KM_GRID, ssr_at_fixed_KM, t, S_obs, d["S0"], SSR)
        KM_lo, KM_hi = ci_from_profile(KM_GRID, T_KM, THRESH)
        T_psi = profile_curve(PSI_GRID, ssr_at_fixed_ratio, t, S_obs, d["S0"], SSR)
        psi_lo, psi_hi = ci_from_profile(PSI_GRID, T_psi, THRESH)
        psi_inf, SSR_inf = limiting_exponential_fit(t, S_obs, d["S0"])
        Dp_inf = max(SSR_inf - SSR, 0.0) / SIGMA_S**2

        results[key] = dict(t=t, S_obs=S_obs, Vhat=Vhat, KMhat=KMhat, redchi=redchi,
                            uA_V=uA_V, uA_KM=uA_KM, uA_ratio=uA_ratio, ratio=ratio,
                            T_KM=T_KM, KM_lo=KM_lo, KM_hi=KM_hi,
                            T_psi=T_psi, psi_lo=psi_lo, psi_hi=psi_hi,
                            psi_inf=psi_inf, Dp_inf=Dp_inf,
                            S0=d["S0"], color=d["color"], label=d["label"])

    # ---------- Simulation-and-recovery ----------
    NSIM = 400
    recov = {}
    for key, d in designs.items():
        t = d["t"]
        Vs, KMs = [], []
        for _ in range(NSIM):
            S_true = substrate(t, V_TRUE, KM_TRUE, d["S0"])
            A_obs = EPS_EXT * PATH * S_true + RNG.normal(0.0, SIGMA_A, size=t.shape)
            S_obs = A_obs / (EPS_EXT * PATH)
            popt_r, _, _ = fit_full(t, S_obs, d["S0"])
            Vs.append(popt_r[0]); KMs.append(popt_r[1])
        Vs = np.array(Vs); KMs = np.array(KMs)
        ratios = Vs / KMs
        at_cap = int(np.sum(KMs >= 9.0e-3))
        recov[key] = dict(
            KM_med=np.median(KMs) * 1e6,
            KM_lo=np.percentile(KMs, 2.5) * 1e6, KM_hi=np.percentile(KMs, 97.5) * 1e6,
            ratio_med=np.median(ratios), ratio_lo=np.percentile(ratios, 2.5),
            ratio_hi=np.percentile(ratios, 97.5), at_cap=at_cap)

    # ---------- Uncertainty budget (use Design II, where K_M is identified) ----------
    b = results["II"]
    uA_V, uA_KM, uA_ratio = b["uA_V"], b["uA_KM"], b["uA_ratio"]
    uc_KM = np.hypot(uA_KM, U_REL_EPS)                          # eps enters K_M
    uc_kcat = np.sqrt(uA_V**2 + U_REL_EPS**2 + U_REL_E0**2)   # eps + [E]_0 enters kcat
    uc_kcatKM = np.hypot(uA_ratio, U_REL_E0)                  # eps cancels in kcat/KM


    # ---------- Calibration-scale check (see `simulate` for the convention) -------
    def epsilon_sensitivity(key, rel=1e-3):
        """d ln(theta) / d ln(eps) under the common concentration scale.

        Rescaling eps by (1+rel) divides every absorbance-derived concentration by
        (1+rel). Because the separately calibrated [S]_0 input is on the same scale,
        it is divided too. The
        expected result is -1 for V, -1 for K_M, and 0 for V/K_M, which is exactly
        what Eqs. (S15)-(S17) of the Supporting Information assume."""
        r = results[key]
        t, S_obs, S0 = r["t"], r["S_obs"], r["S0"]
        out = {}
        for sign in (+1, -1):
            f = 1.0 + sign * rel
            popt, _, _ = fit_full(t, S_obs / f, S0 / f)
            out[sign] = popt
        dln = np.log(1.0 + rel) - np.log(1.0 - rel)
        dV = (np.log(out[+1][0]) - np.log(out[-1][0])) / dln
        dK = (np.log(out[+1][1]) - np.log(out[-1][1])) / dln
        return dV, dK, dV - dK


    EPS_SENS = epsilon_sensitivity("II")

    # ---------- Figure 1 ----------
    def mm_progress_solution(tau, s0):
        """Solve ds/dtau = -s/(1+s), s(0)=s0, with s=[S]/K_M and tau=Vt/K_M."""
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


    fig, axes = plt.subplots(2, 2, figsize=(FIG_W, FIG_H), dpi=DPI)
    axA, axB, axC, axD = axes.ravel()

    # Panel A: substrate-depletion progress curves and fits
    for key in ("II", "I"):
        r = results[key]
        tt = np.linspace(0, r["t"].max(), 400)
        axA.plot(tt / 60.0, substrate(tt, r["Vhat"], r["KMhat"], r["S0"]) * 1e6,
                 color=r["color"], linewidth=2.6, zorder=3)
        y_obs = r["S_obs"] * 1e6
        positive = y_obs > 0.0
        axA.plot(r["t"][positive] / 60.0, y_obs[positive], "o", color=r["color"],
                 markersize=4.2, markerfacecolor="white", markeredgewidth=1.0,
                 label=r["label"], zorder=4)
        if np.any(~positive):
            # A log axis cannot display nonpositive concentration estimates. Retain
            # them in the fit and show them explicitly at a declared plotting floor.
            plotting_floor = 0.030
            axA.plot(r["t"][~positive] / 60.0,
                     np.full(np.sum(~positive), plotting_floor), marker="v",
                     linestyle="none", color=r["color"], markersize=6.0,
                     markerfacecolor="white", markeredgewidth=1.1, zorder=5)
    axA.axhspan(0.6 * KM_TRUE * 1e6, 1.6 * KM_TRUE * 1e6, color=SHADE, zorder=0)
    axA.axhline(KM_TRUE * 1e6, color=DARK, linestyle="--", linewidth=1.5, zorder=2)
    axA.text(0.6, KM_TRUE * 1e6 * 1.22, r"$K_M$", ha="left", va="bottom",
             fontsize=ANNOT, color=DARK)
    axA.set_yscale("log")
    axA.set_xlabel("time (min)")
    axA.set_ylabel(r"substrate $[\mathrm{S}]$ ($\mu$M)")
    axA.set_ylim(0.022, 420)
    axA.set_title("A  Both designs fit the signal well", loc="left",
                  fontweight="bold", fontsize=PANEL, pad=8)
    axA.legend(frameon=False, loc="lower left", fontsize=ANNOT)
    axA.text(0.97, 0.205,
             "reduced $\\chi^2$:\nI {:.2f}   II {:.2f}".format(
                 results["I"]["redchi"], results["II"]["redchi"]),
             transform=axA.transAxes, ha="right", va="top", fontsize=ANNOT, color=DARK)

    # Panel B: profile likelihoods for K_M
    for key in ("II", "I"):
        r = results[key]
        axB.plot(KM_GRID * 1e6, r["T_KM"], color=r["color"], linewidth=2.8,
                 label=r["label"])
    axB.axhline(THRESH, color=DARK, linestyle=":", linewidth=1.8)
    axB.text(4.5e4, THRESH * 1.12, "95% threshold (3.84)", fontsize=ANNOT,
             color=DARK, ha="right", va="bottom")
    axB.axvline(KM_TRUE * 1e6, color=DARK, linestyle="--", linewidth=1.3)
    axB.set_xscale("log")
    axB.set_xlim(1.0, 5.0e4)
    axB.set_ylim(0, 16)
    axB.set_xlabel(r"$K_M$ ($\mu$M)")
    axB.set_ylabel(r"profile deviance $D_p(K_M)$")
    axB.set_title("B  Only Design II bounds $K_M$ above", loc="left",
                  fontweight="bold", fontsize=PANEL, pad=8)
    axB.legend(frameon=False, loc="upper right", fontsize=ANNOT)
    axB.text(4.5e4, 2.60,
             rf"Design I: $D_p(\infty)={results['I']['Dp_inf']:.2f}<3.84$",
             fontsize=ANNOT, color=ORANGE, ha="right", va="center")

    # Panel C: dimensionless substrate-depletion progress curves
    tau = np.linspace(0, 7, 600)
    curves = [(0.2, ORANGE, r"$[\mathrm{S}]_0/K_M = 0.2$  (Design I)"),
              (1.0, GREEN, r"$[\mathrm{S}]_0/K_M = 1$"),
              (4.0, BLUE, r"$[\mathrm{S}]_0/K_M = 4$  (Design II)")]
    axC.axhspan(0.85, 1.15, color=SHADE, zorder=0)
    for s0, color, label in curves:
        axC.plot(tau, mm_progress_solution(tau, s0), color=color, linewidth=2.8,
                 label=label, zorder=3)
    axC.axhline(1.0, color=DARK, linestyle="--", linewidth=1.4, zorder=2)
    # Positioned above the shaded band so it no longer touches the green curve.
    axC.text(6.85, 1.30, r"$[\mathrm{S}] = K_M$", ha="right", va="bottom",
             fontsize=ANNOT, color=DARK)
    axC.set_xlim(0, 7)
    axC.set_ylim(0, 4.4)
    axC.set_xlabel(r"scaled time, $Vt/K_M$")
    axC.set_ylabel(r"substrate $[\mathrm{S}]/K_M$")
    axC.set_title("C  Designs span different regimes",
                  fontsize=PANEL, loc="left", fontweight="bold", pad=8)
    axC.legend(frameon=False, loc="upper right", fontsize=ANNOT, handlelength=1.7)

    # Panel D: local sensitivities
    s_grid = np.logspace(-2, 2, 800)
    for a, bnd in [(1e-2, 1e-1), (0.7, 1.5), (10, 1e2)]:
        axD.axvspan(a, bnd, color=SHADE, zorder=0)
    axD.axhline(0.0, color=ZERO, linewidth=1.0, zorder=1)
    axD.plot(s_grid, np.ones_like(s_grid), color=BLUE, linewidth=2.8, zorder=3)
    axD.plot(s_grid, -1.0 / (1.0 + s_grid), color=ORANGE, linewidth=2.8, zorder=3)
    axD.set_xscale("log")
    axD.set_xlim(1e-2, 1e2)
    axD.set_ylim(-1.20, 1.28)
    axD.set_xlabel(r"instantaneous $[\mathrm{S}]/K_M$")
    axD.set_ylabel("log-sensitivity")
    axD.set_title("D  Sensitivities explain the confounding",
                  fontsize=PANEL, loc="left", fontweight="bold", pad=8)
    axD.text(0.30, 1.115, r"$\partial\ln v/\partial\ln V$", color=BLUE,
             fontsize=ANNOT, ha="center", va="bottom")
    axD.text(6.0, -0.62, r"$\partial\ln v/\partial\ln K_M$", color=ORANGE,
             fontsize=ANNOT, ha="center", va="center")
    axD.text(0.030, 0.46, "low $[\\mathrm{S}]$:\nonly $V/K_M$", fontsize=ANNOT,
             ha="center", va="center", color=DARK, linespacing=1.10)
    axD.text(1.02, 0.56, "near $K_M$:\njoint information", fontsize=ANNOT,
             ha="center", va="center", color=DARK, linespacing=1.10)
    axD.text(34.0, 0.46, "high $[\\mathrm{S}]$:\nmostly $V$", fontsize=ANNOT,
             ha="center", va="center", color=DARK, linespacing=1.10)

    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.08, top=0.955,
                        wspace=0.26, hspace=0.34)
    for ext in ("eps", "pdf", "png"):
        fig.savefig(f"figure1-mm-case-study.{ext}", facecolor="white",
                    bbox_inches="tight", pad_inches=0.10,
                    **({"dpi": DPI} if ext == "png" else {}))
    plt.close(fig)

    # ---------- Numerical summary ----------
    def pct(x): return f"{100*x:.1f}%"

    print("="*70)
    print("GROUND TRUTH")
    print(f"  K_M = {KM_TRUE*1e6:.1f} uM | k_cat = {KCAT_TRUE:.0f} /s | "
          f"[E]_0 = {E0_TRUE*1e9:.1f} nM (phi={ACTIVE_FRACTION}) | "
          f"V = {V_TRUE*1e6:.3f} uM/s")
    print(f"  eps = {EPS_EXT:.0f} /M/cm | path = {PATH} cm | sigma_A = {SIGMA_A} AU "
          f"-> sigma_S = {SIGMA_S*1e6:.3f} uM")
    print(f"  V/K_M (true) = {V_TRUE/KM_TRUE*1e3:.3f} e-3 /s ; "
          f"k_cat/K_M (true) = {KCAT_TRUE/KM_TRUE:.3e} /M/s")
    print("="*70)
    for key in ("I", "II"):
        r = results[key]
        print(f"\n{r['label'].split(' ')[0]} {key}:  reduced chi^2 = {r['redchi']:.3f}  "
              f"(forward adequacy {'OK' if 0.7<r['redchi']<1.4 else 'check'})")
        print(f"  V_hat   = {r['Vhat']*1e6:.4f} uM/s   (Type A SE {pct(r['uA_V'])})")
        print(f"  KM_hat  = {r['KMhat']*1e6:.2f} uM     (Type A SE {pct(r['uA_KM'])})")
        if r["KM_hi"] is not None and r["KM_hi"] >= 4.0e-2:
            print(f"  K_M 95% profile CI: [{r['KM_lo']*1e6:.1f} uM, no finite upper bound]"
                  f"  -> NON-IDENTIFIABLE")
        else:
            print(f"  K_M 95% profile CI: [{r['KM_lo']*1e6:.1f}, {r['KM_hi']*1e6:.1f}] uM")
        print(f"  V/K_M 95% profile CI: [{r['psi_lo']*1e3:.2f}, {r['psi_hi']*1e3:.2f}] e-3 /s "
              f"(Type A SE {pct(r['uA_ratio'])})")
        if key == "I":
            print(f"  limiting exponential: D_p(infinity) = {r['Dp_inf']:.4f} "
                  f"(< {THRESH:.4f}; upper profile region remains open)")
    print("="*70)
    print("SIMULATION-AND-RECOVERY (Nsim = %d)" % NSIM)
    for key in ("I", "II"):
        rc = recov[key]
        print(f"  Design {key}:  K_M median {rc['KM_med']:.1f} uM, "
              f"95% spread [{rc['KM_lo']:.1f}, {rc['KM_hi']:.1f}] uM, "
              f"{rc['at_cap']}/{NSIM} runs at the 10 mM cap")
        print(f"             V/K_M median {rc['ratio_med']*1e3:.3f} e-3/s, "
              f"95% spread [{rc['ratio_lo']*1e3:.3f}, {rc['ratio_hi']*1e3:.3f}] e-3/s")
    print("="*70)
    print("UNCERTAINTY BUDGET (Design II; relative standard uncertainties)")
    print(f"  Type A (fit): u(V)={pct(uA_V)}, u(K_M)={pct(uA_KM)}, u(V/K_M)={pct(uA_ratio)}")
    print(f"  Type B: u(eps)/eps={pct(U_REL_EPS)}, u([E]_0)/[E]_0={pct(U_REL_E0)}")
    print(f"  u_c(K_M)/K_M       = {pct(uc_KM)}   [Type A (+) eps]")
    print(f"  u_c(k_cat)/k_cat   = {pct(uc_kcat)}   [Type A(V) (+) eps (+) [E]_0]")
    print(f"  u_c(k_cat/K_M)     = {pct(uc_kcatKM)}   [eps CANCELS; Type A(V/K_M) (+) [E]_0]")
    print(f"  Fit-only k_cat uncertainty would be {pct(uA_V)} -> understated "
          f"by ~{uc_kcat/uA_V:.1f}x")
    print("="*70)

    print("CALIBRATION-SCALE CHECK (separate [S]_0 input on common scale)")
    print(f"  d ln V     / d ln eps = {EPS_SENS[0]:+.4f}   (expected -1)")
    print(f"  d ln K_M   / d ln eps = {EPS_SENS[1]:+.4f}   (expected -1)")
    print(f"  d ln(V/K_M)/ d ln eps = {EPS_SENS[2]:+.4f}   (expected  0; eps cancels)")
    print("="*70)

    # ---------------------------------------------------------------------------
    # Deposit the generated observations and the headline results.
    #
    # The Perspective argues that a portable record must preserve the measured
    # signal, not merely the means of regenerating it.  Writing the observations
    # out makes this repository obey its own recommendation, and decouples the
    # long-term auditability of the exemplar from the NumPy generator and the
    # SciPy optimizer: a future reader can re-analyze these files even if the
    # random stream or the optimizer changes.
    # ---------------------------------------------------------------------------
    import csv
    import json

    DATA_DIR = Path("data")
    DATA_DIR.mkdir(exist_ok=True)

    for key, r in results.items():
        with open(DATA_DIR / f"design_{key}_observations.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["# Design", key])
            w.writerow(["# [S]_0 / M", f"{r['S0']:.6e}"])
            w.writerow(["# sigma_S / M", f"{SIGMA_S:.6e}"])
            w.writerow(["# seed", 20240530])
            w.writerow(["# note", "concentrations are A/(eps*path); values may be "
                                  "negative where additive absorbance noise exceeds "
                                  "the signal near depletion"])
            w.writerow(["t_s", "S_obs_M", "A_obs_AU"])
            for ti, si in zip(r["t"], r["S_obs"]):
                w.writerow([f"{ti:.6f}", f"{si:.9e}", f"{si * EPS_EXT * PATH:.9f}"])

    expected = {
        "_comment": "Headline values of the worked exemplar. Compare with a "
                    "relative tolerance of 1e-6 on a matching environment, or 1e-3 "
                    "across optimizer/BLAS versions. See tests/test_reproduction.py.",
        "seed": 20240530,
        "ground_truth": {"KM_M": KM_TRUE, "V_M_per_s": V_TRUE,
                         "kcat_per_s": KCAT_TRUE, "E0_M": E0_TRUE,
                         "sigma_S_M": SIGMA_S},
        "designs": {
            k: {"S0_M": r["S0"], "n_obs": int(r["t"].size),
                "t_max_s": float(r["t"][-1]), "reduced_chi2": float(r["redchi"]),
                "V_hat_M_per_s": float(r["Vhat"]), "KM_hat_M": float(r["KMhat"]),
                "uA_rel_V": float(r["uA_V"]), "uA_rel_KM": float(r["uA_KM"]),
                "uA_rel_ratio": float(r["uA_ratio"]),
                "KM_profile_lo_M": float(r["KM_lo"]),
                "KM_profile_hi_M": (None if r["KM_hi"] >= KM_GRID[-1] * 0.999
                                    else float(r["KM_hi"])),
                "psi_lo_per_s": float(r["psi_lo"]), "psi_hi_per_s": float(r["psi_hi"]),
                "Dp_infinity": float(r["Dp_inf"])}
            for k, r in results.items()},
        "recovery_nsim": NSIM,
        "recovery": {k: {kk: float(vv) for kk, vv in v.items()} for k, v in recov.items()},
        "budget_design_II": {"uc_rel_KM": float(uc_KM), "uc_rel_kcat": float(uc_kcat),
                             "uc_rel_kcat_over_KM": float(uc_kcatKM)},
        "epsilon_sensitivity": {"dlnV_dlneps": float(EPS_SENS[0]),
                                "dlnKM_dlneps": float(EPS_SENS[1]),
                                "dlnratio_dlneps": float(EPS_SENS[2])},
    }
    with open("current_results.json", "w") as fh:
        json.dump(expected, fh, indent=2, sort_keys=False)
        fh.write("\n")

    print("Deposited data/design_I_observations.csv, data/design_II_observations.csv,")
    print("and current_results.json")
    print("="*70)


if __name__ == "__main__":
    main()
