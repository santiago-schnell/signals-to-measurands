# Reproducing the figures and results of *From Signals to Measurands*

This repository contains everything needed to regenerate **all figures and all
reported numbers** in the paper:

> **From Signals to Measurands: A Measurement-Science Roadmap for Reproducible
> Analytical Biochemistry**
> Santiago Schnell, manuscript prepared for *ACS Measurement Science Au* (2026).

It is written so that **you do not need any prior experience with Python**. If
you can open a terminal and copy–paste a few commands, you can reproduce the
results. The whole process takes about ten minutes, most of which is the
one-time installation of Python.

There is **no data to download**. The worked example in Figure 4 generates its
own synthetic data from a fixed random seed, so every run follows the same deterministic data-generation and analysis
path.

---

## What is in this repository

```
signals-to-measurands/
├── README.md            <- you are reading this
├── run_all.py           <- the single command that reproduces everything
├── requirements.txt     <- the Python packages you need
├── environment.yml      <- an alternative for Anaconda/conda users
├── LICENSE              <- terms of use (MIT)
├── CITATION.cff         <- how to cite this work
├── code/                <- the five scripts that draw the figures
│   ├── generate_toc.py
│   ├── generate_figure1.py
│   ├── generate_figure2.py
│   ├── generate_figure3.py
│   └── generate_figure4.py
└── figures/             <- the generated figures and the numerical summary
```

---

## Quick start (if you already have Python)

From a terminal opened in this folder:

```bash
pip install -r requirements.txt
python run_all.py
```

The figures and a file called `numerical_summary.txt` will appear in the
`figures/` folder. If `python` is not found, try `python3`.

---

## Step-by-step (no Python experience needed)

### 1. Install Python

Download Python 3.10 or newer from <https://www.python.org/downloads/> and
install it.

- **On Windows:** during installation, tick the box **"Add Python to PATH"**.
  This one checkbox saves a lot of trouble later.

To check it worked, open a terminal (see step 3) and type:

```bash
python --version
```

You should see something like `Python 3.12.3`. On macOS/Linux you may need to
type `python3 --version` instead.

### 2. Get this repository

Either:
- On the GitHub page, click the green **Code** button, choose **Download ZIP**,
  and unzip it; **or**
- if you know `git`: `git clone <repository-url>`.

### 3. Open a terminal *in this folder*

- **Windows:** open the unzipped folder in File Explorer, click in the address
  bar, type `cmd`, and press Enter.
- **macOS:** right-click the folder in Finder → *New Terminal at Folder*.
- **Linux:** right-click inside the folder → *Open in Terminal*.

You should now be "inside" the folder that contains `run_all.py`.

### 4. (Recommended) Create a virtual environment

A virtual environment is just a private box that keeps this project's packages
separate from the rest of your computer. It is optional but tidy.

**Windows:**
```bash
py -m venv .venv
.venv\Scripts\activate
```

**macOS / Linux:**
```bash
python3 -m venv .venv
source .venv/bin/activate
```

After activating, your prompt will usually show `(.venv)` at the start.

### 5. Install the required packages

```bash
pip install -r requirements.txt
```

(If you are using Anaconda instead, you can run
`conda env create -f environment.yml` then `conda activate signals-to-measurands`.)

### 6. Run it

```bash
python run_all.py
```

That's it. You will see each figure being produced, and at the end a message
telling you where the output went.

---

## What you should see

After running, the `figures/` folder contains, for each graphic, three file
formats — `.pdf` and `.eps` (vector, for the manuscript) and `.png` (for quick
viewing):

- `toc-graphic.*` — the graphical abstract
- `figure1-measurement-chain.*` — Figure 1
- `figure2-mm-progress-curves.*` — Figure 2
- `figure3-binding-affinity.*` — Figure 3
- `figure4-worked-exemplar.*` — Figure 4
- `numerical_summary.txt` — every number the analysis produces

### Confirming the results reproduced

Open `figures/numerical_summary.txt` and compare it with **Supporting
Information Tables S1–S4**. A few values to spot-check:

| Quantity | Expected value |
| --- | --- |
| Design I, reduced χ² | 1.02 (a good fit) |
| Design I, *K*<sub>M</sub> 95% profile interval | lower bound ≈ 13 µM, **no finite upper bound** (non-identifiable) |
| Design II, reduced χ² | 1.13 (also a good fit) |
| Design II, *K*<sub>M</sub> estimate / 95% interval | 49.5 µM / [48.4, 50.5] µM |
| *V*/*K*<sub>M</sub> 95% interval (Design II) | [3.97, 4.11] × 10⁻³ s⁻¹ |
| Combined uncertainty of *k*<sub>cat</sub> | 8.3% (vs. 0.4% from the curve fit alone) |

If those match, you have reproduced the paper.

---

## Which script makes which figure

| Script | Produces | Needs SciPy? |
| --- | --- | --- |
| `generate_toc.py` | Graphical abstract | no |
| `generate_figure1.py` | Figure 1 (measurement chain) | no |
| `generate_figure2.py` | Figure 2 (progress curves) | no |
| `generate_figure3.py` | Figure 3 (binding affinity) | no |
| `generate_figure4.py` | Figure 4 + the numerical summary | **yes** |

You can also run any one of them individually, for example:

```bash
cd figures
python ../code/generate_figure4.py
```

(`run_all.py` simply does this for all five and collects the output for you.)

---

## Why the results are exactly reproducible

The worked example draws its synthetic measurements from NumPy's random
generator seeded with the fixed value `20240530` (see the top of
`code/generate_figure4.py`). The fixed seed, recorded package versions, and common analysis path make the
run auditable and repeatable. Small numerical differences can occur with other
package or platform versions; use the supplied environment for exact comparison
with the deposited numerical summary.

---


## Methodological references

The progress-curve model in `generate_figure4.py` uses the Schnell--Mendoza
closed-form solution of the integrated Michaelis--Menten equation:

- Schnell, S.; Mendoza, C. *J. Theor. Biol.* **1997**, *187*, 207--212.
  DOI: 10.1006/jtbi.1997.0425.
- Goudar, C. T.; Sonnad, J. R.; Duggleby, R. G. *Biochim. Biophys. Acta*
  **1999**, *1429*, 377--383. DOI: 10.1016/S0167-4838(98)00247-7.

---

## Software versions used

The deposited figures were produced with:

- Python 3.12
- NumPy 2.4.4
- SciPy 1.17.1
- matplotlib 3.10.8

No proprietary fonts are required; the scripts use matplotlib's built-in
DejaVu Sans.

---

## Troubleshooting

- **`python: command not found`** — try `python3` instead. On Windows, try
  `py`. If none work, Python is not installed or (on Windows) was installed
  without "Add Python to PATH"; reinstall and tick that box.
- **`pip: command not found`** — use `python -m pip install -r requirements.txt`
  (or `python3 -m pip ...`).
- **A package "is not installed" message from `run_all.py`** — run
  `pip install -r requirements.txt` first.
- **Permission errors when installing** — use a virtual environment (step 4);
  it avoids needing administrator rights.
- **Nothing appears** — make sure the terminal is open *inside* the folder that
  contains `run_all.py` (step 3).

---

## How to cite

If you use this code or build on it, please cite both the article and this
repository. GitHub will offer a "Cite this repository" button generated from
`CITATION.cff`. The article reference is:

> Schnell, S. From Signals to Measurands: A Measurement-Science Roadmap for
> Reproducible Analytical Biochemistry. *ACS Meas. Sci. Au* **2026**.

---

## License

The code in this repository is released under the MIT License (see `LICENSE`),
which permits reuse with attribution.

---

## Contact

Santiago Schnell — santiago.schnell@dartmouth.edu
ORCID: <https://orcid.org/0000-0002-9477-3914>
