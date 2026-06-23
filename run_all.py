#!/usr/bin/env python3
"""
run_all.py  --  Reproduce every figure and number in the paper.

    "From Signals to Measurands: A Measurement-Science Roadmap for
     Reproducible Analytical Biochemistry"  (S. Schnell)

You do NOT need to understand this file to use it. Just run, from a terminal
opened in this folder:

    python run_all.py

(On macOS/Linux you may need to type "python3" instead of "python".)

What it does, in plain terms:
  1. Finds the five figure scripts in the "code/" folder.
  2. Runs each one, sending its output into the "figures/" folder.
  3. Saves a plain-text record of the numbers in
     "figures/numerical_summary.txt" so you can check them against the
     tables in the paper's Supporting Information.

Nothing is downloaded and no external data are needed: the worked example
(Figure 4) creates its own data from a fixed random seed, so everyone who
runs this gets exactly the same numbers.
"""
from __future__ import annotations

import subprocess
import sys
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
    ("generate_figure1.py",  "Figure 1  - the signal-to-measurand chain"),
    ("generate_figure2.py",  "Figure 2  - Michaelis-Menten progress curves"),
    ("generate_figure3.py",  "Figure 3  - ligand-binding affinity"),
    ("generate_figure4.py",  "Figure 4  - worked exemplar (this one prints the numbers)"),
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

        print(f"  Running {filename:22s}  ({description}) ...")
        # We run each script with its working directory set to "figures/",
        # so every image it writes lands there automatically.
        result = subprocess.run(
            [sys.executable, str(script_path)],
            cwd=str(FIGURES),
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            print(f"    -> FAILED. Error message:\n{result.stderr}")
            return 1

        # Figure 4 prints the full numerical summary; keep it for the record.
        summary_lines.append(f"{'=' * 70}\n# OUTPUT OF {filename}\n{'=' * 70}\n")
        summary_lines.append(result.stdout)

    # Save the captured numbers so they can be compared with the paper.
    summary_path = FIGURES / "numerical_summary.txt"
    summary_path.write_text("".join(summary_lines))

    print("\nDone.")
    print(f"  Figures written to:  {FIGURES}")
    print(f"  Numbers written to:  {summary_path}")
    print("\nTo confirm the results reproduced, compare the values in")
    print("figures/numerical_summary.txt with Supporting Information Tables S1-S4.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
