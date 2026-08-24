#!/usr/bin/env python3
"""
Sampling-schedule validation for:
"From Signals to Measurands: A Measurement-Science Roadmap for
Reproducible Analytical Biochemistry"

Produces Supporting Information Table S4.

WHY THIS IS A SEPARATE SCRIPT
-----------------------------
Designs I and II differ not only in the initial substrate concentration but
also in the observation window and the number of observations (50 points over
1000 s versus 60 points over 1800 s), because the two designs deplete
substrate on very different timescales.  A reader could reasonably ask whether
Design I fails because of its substrate regime or merely because it has less
data.

This study answers that.  For each row it generates NSIM synthetic data sets
under the indicated design and schedule, profiles K_M over
[1 uM, 10 mM], and records

  (i)  how often the nominal 95% profile region is still open at the 10 mM
       limit of the profiled range, and
  (ii) the median width of the region on a logarithmic scale.

Note the careful wording: reaching the grid limit without closure shows that
no upper bound was found WITHIN THE PROFILED RANGE.  The stronger asymptotic
statement -- that D_p approaches a finite limit below 3.841 as K_M -> infinity
at fixed V/K_M -- is made in the Supporting Information for the single worked
Design I data set, where the profile is evaluated to 50 mM.  It does not
follow from this table alone.

The computation is roughly 18,000 one-dimensional profile fits. It is
kept separate from `run_all.py` so that the principal figure workflow remains
short; on typical current hardware it takes tens of seconds to a few minutes.  Run it with

    python validate_schedules.py

from the repository root (or anywhere). It writes a machine-readable record to `figures/schedule_validation.json` and also prints the table values.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import generate_figure1 as core  # side-effect free: execution is under main()

REPO = HERE.parent
OUT = REPO / "figures" / "schedule_validation.json"


NSIM = 100
KM_GRID = np.logspace(np.log10(1e-6), np.log10(1e-2), 36)
GRID_STEP = np.log10(KM_GRID[1] / KM_GRID[0])

SCHEDULES = [
    (0.2, 1000.0, 50, "Design I,  as reported (n=50, 1000 s)"),
    (0.2, 1800.0, 60, "Design I,  Design II schedule (n=60, 1800 s)"),
    (0.2, 1800.0, 200, "Design I,  n=200 over 1800 s"),
    (0.2, 5000.0, 300, "Design I,  n=300 over 5000 s"),
    (4.0, 1800.0, 60, "Design II, as reported (n=60, 1800 s)"),
]


def schedule_study(core, S0, tmax, n, rng, label, row, total):
    t = np.linspace(0.0, tmax, n)
    S_true = core.substrate(t, core.V_TRUE, core.KM_TRUE, S0)
    still_open, log_widths = 0, []
    t0 = time.time()
    for i in range(NSIM):
        S_obs = S_true + rng.normal(0.0, core.SIGMA_S, size=t.size)
        _, _, ssr_star = core.fit_full(t, S_obs, S0)
        T = core.profile_curve(KM_GRID, core.ssr_at_fixed_KM, t, S_obs, S0,
                               ssr_star)
        lo, hi = core.ci_from_profile(KM_GRID, T, core.THRESH)
        if lo is None:
            continue
        if hi >= KM_GRID[-1] * 0.999:
            still_open += 1
            log_widths.append(np.log10(KM_GRID[-1] / lo))
        else:
            log_widths.append(np.log10(hi / lo))
        if (i + 1) % 20 == 0:
            done = (i + 1) / NSIM
            eta = (time.time() - t0) * (1 - done) / done
            print(f"    [{row}/{total}] {label[:34]:34s} "
                  f"{i + 1:3d}/{NSIM} replicates  (~{eta:4.0f} s left)",
                  flush=True)
    return still_open / NSIM, float(np.median(log_widths))


def main() -> None:
    print("Loaded the side-effect-free model and estimator from generate_figure1.py.",
          flush=True)
    print(f"Sampling-schedule validation: {len(SCHEDULES)} schedules x {NSIM} "
          f"replicates x {KM_GRID.size} profile points.")
    print("Expect tens of seconds to a few minutes. Progress is reported per schedule.\n",
          flush=True)

    rng = np.random.default_rng(20240531)
    rows = []
    t0 = time.time()
    for k, (fac, tmax, n, label) in enumerate(SCHEDULES, start=1):
        frac, width = schedule_study(core, fac * core.KM_TRUE, tmax, n, rng,
                                     label, k, len(SCHEDULES))
        rows.append((label, frac, width))

    print()
    print("=" * 70)
    print(f"SAMPLING-SCHEDULE ROBUSTNESS (Nsim = {NSIM} per row)")
    print("  'Open at limit' = the nominal 95% profile region for K_M is still")
    print("  open at the 10 mM limit of the profiled range. Width is the median")
    print(f"  span of that region in decades; the grid resolves {GRID_STEP:.2f}"
          " decades.")
    for label, frac, width in rows:
        w = f"< {GRID_STEP:.2f}" if width < GRID_STEP else f"{width:.2f}"
        print(f"  {label:44s} open at limit in {frac * 100:5.1f}% of "
              f"replicates; median width {w} decades")
    print("=" * 70)
    elapsed = time.time() - t0
    OUT.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "seed": 20240531,
        "nsim_per_row": NSIM,
        "profile_KM_min_M": float(KM_GRID[0]),
        "profile_KM_max_M": float(KM_GRID[-1]),
        "grid_step_decades": float(GRID_STEP),
        "rows": [
            {"label": label, "open_fraction": float(frac),
             "median_width_decades": float(width)}
            for label, frac, width in rows
        ],
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"Completed in {elapsed:.0f} s.")
    print("These are the values in Supporting Information Table S4.")
    print(f"Machine-readable results written to {OUT}.")


if __name__ == "__main__":
    main()
