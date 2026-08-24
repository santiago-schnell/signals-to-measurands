#!/usr/bin/env python3
"""
run_all.py  --  Reproduce the figures and the worked-exemplar numbers.

    "From Signals to Measurands: A Measurement-Science Roadmap for
     Reproducible Analytical Biochemistry"  (S. Schnell)

You do NOT need to understand this file to use it. Just run, from a terminal
opened in this folder:

    python run_all.py

(On macOS/Linux you may need to type "python3" instead of "python".)

What it does, in plain terms:
  1. Finds the four figure scripts in the "code/" folder.
  2. Runs each one, sending its output into the "figures/" folder, and
     streaming its progress to your screen as it goes.
  3. Saves a plain-text record of the numbers in
     "figures/numerical_summary.txt" so you can check them against
     Supporting Information Tables S1-S3 and S5.
  4. Writes the generated observations to "figures/data/" and the headline
     results to "figures/current_results.json".

No external data are required. The worked example creates its synthetic
observations from fixed seed 20240530, and the same observations are deposited
as CSV files for direct audit.

Expect roughly two minutes in total; Figure 1 does the bulk of the work
(400 simulation-and-recovery fits per design).

NOT INCLUDED HERE
-----------------
Supporting Information Table S4, the sampling-schedule comparison, is a
separate and slower study (about 18,000 profile fits). Run it on its own with

    python code/validate_schedules.py

after this script has finished.
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

# ----------------------------------------------------------------------
# Locations. Everything is relative to THIS file, so the repository works
# no matter where you put it on your computer.
# ----------------------------------------------------------------------
REPO = Path(__file__).resolve().parent
CODE = REPO / "code"
FIGURES = REPO / "figures"

# The scripts to run, in order, with a friendly description of each.
SCRIPTS = [
    ("generate_toc.py",      "Graphical abstract (table-of-contents graphic)"),
    ("generate_figure1.py",  "Figure 1  - Michaelis-Menten case study (prints the numbers)"),
    ("generate_figure2.py",  "Figure 2  - the signal-to-measurand chain"),
    ("generate_figure3.py",  "Figure 3  - ligand-receptor binding"),
]


def main() -> int:
    FIGURES.mkdir(exist_ok=True)

    # Quick, friendly check that the required libraries are installed.
    missing = []
    for pkg in ("numpy", "scipy", "matplotlib"):
        try:
            __import__(pkg)
        except ImportError:
            missing.append(pkg)
    if missing:
        print("ERROR: the following required package(s) are not installed: "
              + ", ".join(missing))
        print("Fix this by running:  pip install -r requirements.txt")
        return 1

    summary_lines = []
    print("Reproducing all figures. Output will appear in the 'figures' folder.\n")

    for filename, description in SCRIPTS:
        script_path = CODE / filename
        if not script_path.exists():
            print(f"  [SKIP] {filename} was not found in the code/ folder.")
            continue

        print(f"  Running {filename:22s}  ({description}) ...", flush=True)
        started = time.time()
        # We run each script with its working directory set to "figures/",
        # so every image it writes lands there automatically.  The output is
        # echoed line by line as it arrives rather than captured in silence,
        # so a long-running script never looks as though it has hung.
        proc = subprocess.Popen(
            [sys.executable, "-u", str(script_path)],
            cwd=str(FIGURES),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        captured = []
        for line in proc.stdout:
            captured.append(line)
            print("      " + line.rstrip(), flush=True)
        proc.wait()

        if proc.returncode != 0:
            print(f"    -> FAILED after {time.time() - started:.0f} s.")
            return 1
        print(f"    -> done in {time.time() - started:.0f} s.\n", flush=True)

        # Figure 1 prints the full numerical summary; keep it for the record.
        summary_lines.append(f"{'=' * 70}\n# OUTPUT OF {filename}\n{'=' * 70}\n")
        summary_lines.append("".join(captured))

    # Save the captured numbers so they can be compared with the paper.
    summary_path = FIGURES / "numerical_summary.txt"
    summary_path.write_text("".join(summary_lines))

    print("\nDone.")
    print(f"  Figures written to:  {FIGURES}")
    print(f"  Numbers written to:  {summary_path}")
    print("\nTo confirm the results reproduced, compare the values in")
    print("figures/numerical_summary.txt with Supporting Information")
    print("Tables S1-S3 and S5, or run:  python -m pytest tests")
    print("\nFor Supporting Information Table S4 (separate seed 20240531; slower):")
    print("    python code/validate_schedules.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
