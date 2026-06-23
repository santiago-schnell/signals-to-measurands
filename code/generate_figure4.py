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

Forward adequacy holds in BOTH. Inverse security does not: a profile
likelihood shows K_M is practically non-identifiable in Design I (no finite
upper bound) but well constrained in Design II, while the specificity-
determining combination V/K_M is tightly determined in both.

The uncertainty budget then shows that the curve-fit (Type A) standard error
on k_cat is a small fraction of the combined standard uncertainty, which is
dominated by the Type B uncertainty in initial active-site concentration; and that
the molar absorption coefficient cancels in k_cat/K_M but not in k_cat or K_M.

Outputs: figure4-worked-exemplar.{eps,pdf,png} and a printed numerical summary.
"""
from __future__ import annotations

import numpy as np
from scipy.special import lambertw
from scipy.optimize import least_squares
from scipy.stats import chi2

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- Manuscript styling (matches generate_figures.py) ----------
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "dejavusans",
    "font.size": 16,
    "axes.titlesize": 16,
    "axes.labelsize": 17,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 12.5,
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
    absorbance, S_obs = A_obs/(eps*path), A_obs = eps*path*S_true + N(0,sigma_A)."""
    S_true = substrate(t, V_TRUE, KM_TRUE, S0)
    A_obs = EPS_EXT * PATH * S_true + RNG.normal(0.0, SIGMA_A, size=t.shape)
    return A_obs / (EPS_EXT * PATH)


def resid_full(p, t, S_obs, S0):
    V, KM = p
    return (substrate(t, V, KM, S0) - S_obs) / SIGMA_S


def fit_full(t, S_obs, S0, p0=(V_TRUE, KM_TRUE)):
    """Global (V, KM) fit; returns popt, pcov, SSR."""
    lb = (1e-9, 1e-7)
    ub = (1e-3, 1e-2)
    sol = least_squares(resid_full, x0=p0, args=(t, S_obs, S0),
                        bounds=(lb, ub), method="trf", xtol=1e-14, ftol=1e-14)
    SSR = float(np.sum(sol.fun**2)) * SIGMA_S**2
    # covariance from Gauss-Newton approximation J^T J (residuals are weighted)
    J = sol.jac
    JTJ = J.T @ J
    try:
        pcov = np.linalg.inv(JTJ)
    except np.linalg.LinAlgError:
        pcov = np.full((2, 2), np.nan)
    return sol.x, pcov, SSR


def ssr_at_fixed_KM(KM, t, S_obs, S0):
    """Profile over V at fixed KM."""
    sol = least_squares(lambda V: (substrate(t, V[0], KM, S0) - S_obs) / SIGMA_S,
                        x0=[V_TRUE], bounds=([1e-9], [1e-3]),
                        method="trf", xtol=1e-14, ftol=1e-14)
    return float(np.sum(sol.fun**2)) * SIGMA_S**2


def ssr_at_fixed_ratio(psi, t, S_obs, S0):
    """Profile over KM at fixed psi = V/KM (so V = psi*KM)."""
    sol = least_squares(
        lambda KM: (substrate(t, psi * KM[0], KM[0], S0) - S_obs) / SIGMA_S,
        x0=[KM_TRUE], bounds=([1e-7], [1e-2]),
        method="trf", xtol=1e-14, ftol=1e-14)
    return float(np.sum(sol.fun**2)) * SIGMA_S**2


def profile_curve(grid, ssr_fn, t, S_obs, S0):
    vals = np.array([ssr_fn(g, t, S_obs, S0) for g in grid])
    return (vals - vals.min()) / SIGMA_S**2  # delta(-2 log L)


def ci_from_profile(grid, T, thresh):
    inside = grid[T <= thresh]
    if inside.size == 0:
        return None, None
    return inside.min(), inside.max()


# ---------- Designs ----------
designs = {
    "I":  dict(label="Design I  ($[\\mathrm{S}]_0=0.2\\,K_M$)", color=ORANGE,
               S0=0.2 * KM_TRUE, t=np.linspace(0.0, 1000.0, 50)),
    "II": dict(label="Design II ($[\\mathrm{S}]_0=4\\,K_M$)", color=BLUE,
               S0=4.0 * KM_TRUE, t=np.linspace(0.0, 1800.0, 60)),
}

THRESH = chi2.ppf(0.95, df=1)   # 3.841
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

    # profiles
    T_KM = profile_curve(KM_GRID, ssr_at_fixed_KM, t, S_obs, d["S0"])
    KM_lo, KM_hi = ci_from_profile(KM_GRID, T_KM, THRESH)
    T_psi = profile_curve(PSI_GRID, ssr_at_fixed_ratio, t, S_obs, d["S0"])
    psi_lo, psi_hi = ci_from_profile(PSI_GRID, T_psi, THRESH)

    results[key] = dict(t=t, S_obs=S_obs, Vhat=Vhat, KMhat=KMhat, redchi=redchi,
                        uA_V=uA_V, uA_KM=uA_KM, uA_ratio=uA_ratio, ratio=ratio,
                        T_KM=T_KM, KM_lo=KM_lo, KM_hi=KM_hi,
                        T_psi=T_psi, psi_lo=psi_lo, psi_hi=psi_hi,
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
        sol = least_squares(resid_full, x0=(V_TRUE, KM_TRUE),
                            args=(t, S_obs, d["S0"]),
                            bounds=((1e-9, 1e-7), (1e-3, 1e-2)),
                            method="trf", xtol=1e-13, ftol=1e-13)
        Vs.append(sol.x[0]); KMs.append(sol.x[1])
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

# ---------- Figure 4 ----------
fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.0, 5.9), dpi=DPI,
                               gridspec_kw={"width_ratios": [1.05, 1.0]})

# Panel A: progress curves + fits
for key in ("II", "I"):
    r = results[key]
    tt = np.linspace(0, r["t"].max(), 400)
    axA.plot(tt / 60.0, substrate(tt, r["Vhat"], r["KMhat"], r["S0"]) * 1e6,
             color=r["color"], linewidth=2.6, zorder=3)
    axA.plot(r["t"] / 60.0, r["S_obs"] * 1e6, "o", color=r["color"],
             markersize=4.5, markerfacecolor="white", markeredgewidth=1.1,
             label=r["label"], zorder=4)
axA.axhspan(0.6 * KM_TRUE * 1e6, 1.6 * KM_TRUE * 1e6, color=SHADE, zorder=0)
axA.axhline(KM_TRUE * 1e6, color=DARK, linestyle="--", linewidth=1.6, zorder=2)
axA.text(axA.get_xlim()[1] * 0.985, KM_TRUE * 1e6 * 1.18, r"$K_M$",
         ha="right", va="bottom", fontsize=14, color=DARK)
axA.set_yscale("log")
axA.set_xlabel("time (min)")
axA.set_ylabel(r"substrate $[\mathrm{S}]$ ($\mu$M)")
axA.set_ylim(0.08, 320)
axA.set_title("A  Both designs fit the signal", loc="left", fontweight="bold",
              fontsize=15.5, pad=8)
axA.legend(frameon=False, loc="lower left", fontsize=12)
axA.text(0.97, 0.94,
         "reduced $\\chi^2$:\nI {:.2f}   II {:.2f}".format(
             results["I"]["redchi"], results["II"]["redchi"]),
         transform=axA.transAxes, ha="right", va="top", fontsize=12, color=DARK)

# Panel B: profile likelihoods for K_M
for key in ("II", "I"):
    r = results[key]
    axB.plot(KM_GRID * 1e6, r["T_KM"], color=r["color"], linewidth=2.8,
             label=r["label"])
axB.axhline(THRESH, color=DARK, linestyle=":", linewidth=1.8)
axB.text(4.2e4, THRESH * 1.06, "95% threshold (3.84)", fontsize=11.5,
         color=DARK, ha="right", va="bottom")
axB.axvline(KM_TRUE * 1e6, color=DARK, linestyle="--", linewidth=1.3)
axB.set_xscale("log")
axB.set_xlim(1.0, 5.0e4)
axB.set_ylim(0, 16)
axB.set_xlabel(r"$K_M$ ($\mu$M)")
axB.set_ylabel(r"profile $\Delta(-2\ln L)$")
axB.set_title("B  Only one design constrains $K_M$", loc="left",
              fontweight="bold", fontsize=15.5, pad=8)
axB.legend(frameon=False, loc="upper right", fontsize=12)
# label the flat tail directly (no arrow needed): K_M unbounded above in Design I
axB.text(1.5e3, 2.75, "Design I: no upper bound\n($K_M$ non-identifiable)",
         fontsize=12, color=ORANGE, ha="center", va="center")
axB.text(0.45, 0.60,
"$V/K_M$ estimable in both:\nI [{:.1f}, {:.1f}]   II [{:.1f}, {:.1f}] $\\times10^{{-3}}$ s$^{{-1}}$".format(
             results["I"]["psi_lo"] * 1e3, results["I"]["psi_hi"] * 1e3,
             results["II"]["psi_lo"] * 1e3, results["II"]["psi_hi"] * 1e3),
         transform=axB.transAxes, fontsize=11, color=DARK,
         va="top", ha="left")

fig.subplots_adjust(left=0.07, right=0.985, bottom=0.135, top=0.90, wspace=0.24)
for ext in ("eps", "pdf", "png"):
    fig.savefig(f"figure4-worked-exemplar.{ext}",
                dpi=DPI, facecolor="white")
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
