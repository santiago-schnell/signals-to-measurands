"""Tolerance-based reproduction tests for the worked exemplar.

Run from the repository root after ``python run_all.py``::

    python -m pytest tests

The generated results in ``figures/current_results.json`` are compared with an
immutable reference record committed under ``tests/reference_results.json``.
The separate sampling-schedule study has its own generated and reference records.

A fixed random seed fixes the synthetic observations, but optimizer results can
vary slightly with SciPy, BLAS/LAPACK, and platform floating-point behavior.
The tests therefore use tolerances for numerical values while checking the
scientific conclusions exactly.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
CURRENT = REPO / "figures" / "current_results.json"
REFERENCE = REPO / "tests" / "reference_results.json"
SCHEDULE_CURRENT = REPO / "figures" / "schedule_validation.json"
SCHEDULE_REFERENCE = REPO / "tests" / "reference_schedules.json"
RTOL = 1e-3
ATOL = 1e-12


def _load(path: Path) -> dict:
    with path.open() as fh:
        return json.load(fh)


@pytest.fixture(scope="module")
def results() -> dict:
    if not CURRENT.exists():
        pytest.skip("current_results.json not found; run `python run_all.py` first")
    return _load(CURRENT)


@pytest.fixture(scope="module")
def reference() -> dict:
    if not REFERENCE.exists():
        pytest.fail("immutable tests/reference_results.json is missing")
    return _load(REFERENCE)


def _close(a: float, b: float, rtol: float = RTOL, atol: float = ATOL) -> bool:
    return abs(a - b) <= atol + rtol * abs(b)


def _compare_leaf(actual, expected, path="root") -> None:
    """Recursively compare the committed scientific reference record."""
    if isinstance(expected, dict):
        assert isinstance(actual, dict), path
        for key, value in expected.items():
            assert key in actual, f"missing {path}.{key}"
            _compare_leaf(actual[key], value, f"{path}.{key}")
    elif isinstance(expected, list):
        assert isinstance(actual, list) and len(actual) == len(expected), path
        for i, (a, e) in enumerate(zip(actual, expected)):
            _compare_leaf(a, e, f"{path}[{i}]")
    elif isinstance(expected, float):
        assert isinstance(actual, (int, float)), path
        assert _close(float(actual), expected), f"{path}: {actual} != {expected}"
    else:
        assert actual == expected, f"{path}: {actual!r} != {expected!r}"


# --- immutable-reference comparison -----------------------------------------

def test_generated_results_match_reference(results, reference):
    _compare_leaf(results, reference)


# --- ground truth and noise model -------------------------------------------

def test_ground_truth(results):
    g = results["ground_truth"]
    assert _close(g["KM_M"], 50e-6)
    assert _close(g["V_M_per_s"], 2.0e-7)
    assert _close(g["kcat_per_s"], 50.0)
    assert _close(g["E0_M"], 4.0e-9)
    assert _close(g["sigma_S_M"], 0.003 / 6220.0)


# --- forward adequacy and identifiability -----------------------------------

@pytest.mark.parametrize("design", ["I", "II"])
def test_reduced_chi2_near_one(results, design):
    assert 0.7 < results["designs"][design]["reduced_chi2"] < 1.4


def test_design_I_has_no_finite_upper_bound(results):
    d = results["designs"]["I"]
    assert d["KM_profile_hi_M"] is None
    # The limiting exponential supplies the K_M -> infinity profile asymptote.
    assert _close(d["Dp_infinity"], 1.6773, rtol=0.01)
    assert d["Dp_infinity"] < 3.841


def test_design_I_has_a_finite_lower_bound(results):
    lo = results["designs"]["I"]["KM_profile_lo_M"]
    assert lo is not None and _close(lo, 12.8e-6, rtol=0.05)


def test_design_II_interval_is_bounded_and_contains_truth(results):
    d = results["designs"]["II"]
    assert _close(d["KM_profile_lo_M"], 48.4e-6, rtol=0.02)
    assert _close(d["KM_profile_hi_M"], 50.5e-6, rtol=0.02)
    assert d["KM_profile_lo_M"] < 50e-6 < d["KM_profile_hi_M"]


def test_both_designs_constrain_the_ratio(results):
    for design in ("I", "II"):
        d = results["designs"][design]
        assert d["psi_lo_per_s"] < 4.0e-3 < d["psi_hi_per_s"]


def test_design_II_point_estimates(results):
    d = results["designs"]["II"]
    assert _close(d["V_hat_M_per_s"], 1.9971816e-07)
    assert _close(d["KM_hat_M"], 4.9517765e-05)


# --- simulation-and-recovery and uncertainty budget ------------------------

def test_recovery_separates_the_designs(results):
    r = results["recovery"]
    assert r["I"]["at_cap"] > 20
    assert r["II"]["at_cap"] == 0
    assert r["II"]["KM_hi"] - r["II"]["KM_lo"] < 5.0
    assert r["I"]["KM_hi"] > 1000.0


def test_budget_dominated_by_active_site_concentration(results):
    b = results["budget_design_II"]
    a = results["designs"]["II"]
    assert _close(b["uc_rel_KM"], 0.024, rtol=0.05)
    assert _close(b["uc_rel_kcat"], 0.083, rtol=0.05)
    assert _close(b["uc_rel_kcat_over_KM"], 0.081, rtol=0.05)
    assert b["uc_rel_kcat"] > 15 * a["uA_rel_V"]


def test_epsilon_cancels_from_the_specificity_constant(results):
    e = results["epsilon_sensitivity"]
    assert abs(e["dlnV_dlneps"] + 1.0) < 1e-3
    assert abs(e["dlnKM_dlneps"] + 1.0) < 1e-3
    assert abs(e["dlnratio_dlneps"]) < 1e-3


# --- deposited observations --------------------------------------------------

@pytest.mark.parametrize("design,n", [("I", 50), ("II", 60)])
def test_observations_were_deposited(design, n):
    path = REPO / "figures" / "data" / f"design_{design}_observations.csv"
    if not path.exists():
        pytest.skip("observations not found; run `python run_all.py` first")
    rows = [ln for ln in path.read_text().splitlines()
            if ln and not ln.startswith("#")]
    assert rows[0].startswith("t_s,")
    assert len(rows) - 1 == n


def test_nonpositive_design_I_observations_are_preserved():
    path = REPO / "figures" / "data" / "design_I_observations.csv"
    if not path.exists():
        pytest.skip("observations not found; run `python run_all.py` first")
    values = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#") or line.startswith("t_s,"):
            continue
        values.append(float(line.split(",")[1]))
    assert sum(value <= 0.0 for value in values) == 3


# --- separate sampling-schedule validation ----------------------------------

def test_sampling_schedule_results_match_reference():
    if not SCHEDULE_CURRENT.exists():
        pytest.skip("schedule_validation.json not found; run `python code/validate_schedules.py`")
    if not SCHEDULE_REFERENCE.exists():
        pytest.fail("immutable tests/reference_schedules.json is missing")
    _compare_leaf(_load(SCHEDULE_CURRENT), _load(SCHEDULE_REFERENCE))
