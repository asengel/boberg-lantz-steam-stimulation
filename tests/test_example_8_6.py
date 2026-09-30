"""Validation tests for Green & Willhite Example 8.6."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from boberg_lantz import ModelInputs, run_model


def test_example_8_6_heated_radius():
    _, geometry = run_model(ModelInputs())
    # Published worked value: approximately 51.2 ft.
    assert abs(geometry["r_h"] - 51.2) < 0.1


def test_example_8_6_hot_oil_viscosity():
    results, _ = run_model(ModelInputs())
    # Published worked value: approximately 1.22 cp.
    assert abs(results.loc[0, "mu_oh (cp)"] - 1.22) < 0.03


def test_example_8_6_surface_oil_rate():
    results, _ = run_model(ModelInputs())
    # Published hand-calculation value: approximately 22.75 STB/D.
    assert abs(results.loc[0, "qo (STB/D)"] - 22.75) < 0.05


def test_example_8_6_reservoir_oil_rate():
    results, _ = run_model(ModelInputs())
    # Published hand-calculation value: approximately 25.71 RB/D.
    assert abs(results.loc[0, "qoh (RB/D)"] - 25.71) < 0.05


def test_example_8_6_first_heat_fraction():
    results, _ = run_model(ModelInputs())
    # Published worked value: approximately 0.038.
    assert abs(results.loc[0, "delta"] - 0.038) < 0.002
