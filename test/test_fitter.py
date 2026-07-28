"""Unit / end-to-end tests for :mod:`lyapower.fitter`."""

import numpy as np
import pytest

from lyapower import fitter

from conftest import FLUX_FILE


# --- Arinyo non-linear model terms -----------------------------------------


def test_D0_matches_formula():
    k, mu = 0.5, 0.6
    params = dict(k_nl=6.4, a_nl=0.4, k_p=13.0, a_p=1.6,
                  k_v0=1.2, a_v0=1.2, k_v1=1.5, a_v1=1.5)
    got = fitter.D0(k, mu, **params)
    expected = np.exp(
        (k / params["k_nl"]) ** params["a_nl"]
        - (k / params["k_p"]) ** params["a_p"]
        - ((k * mu) / (params["k_v0"] * (1 + (k / params["k_v1"])) ** params["a_v1"]))
        ** params["a_v0"]
    )
    assert got == pytest.approx(expected)


def test_D0_tends_to_one_at_small_k():
    # every term vanishes as k -> 0, so D0 -> 1
    val = fitter.D0(1e-12, 0.5, 6.4, 0.4, 13.0, 1.6, 1.2, 1.2, 1.5, 1.5)
    assert val == pytest.approx(1.0, abs=1e-3)


def test_D1_matches_formula():
    k, mu, plin = 0.5, 0.6, 20.0
    params = dict(q_1=0.8, q_2=0.2, k_v=1.2, a_v=0.5, b_v=1.6, k_p=13.0)
    got = fitter.D1(k, mu, linear_power_spectrum=plin, **params)
    delta2 = (1 / (2 * np.pi**2)) * k**3 * plin
    expected = np.exp(
        (params["q_1"] * delta2 + params["q_2"] * delta2**2)
        * (1 - ((k / params["k_v"]) ** params["a_v"]) * mu ** params["b_v"])
        - (k / params["k_p"]) ** 2
    )
    assert got == pytest.approx(expected)


# --- helpers ---------------------------------------------------------------


def test_rebin_matter_power_interpolates_and_nans_oob():
    k = np.array([0.1, 0.5, 1.0])
    p = np.array([10.0, 5.0, 2.0])
    # a wavenumber inside the range is interpolated
    assert fitter.rebin_matter_power(p, k, np.array([0.5]))[0] == pytest.approx(5.0)
    # out of range -> NaN
    assert np.isnan(fitter.rebin_matter_power(p, k, np.array([10.0]))[0])


def test_mask_data_sets_near_zero_to_nan():
    a = np.array([1.0, 2.0, 0.0, 4.0])
    b = np.array([10.0, 20.0, 30.0, 40.0])
    fitter.mask_data([0], a, b)  # mask where a is ~0
    assert np.isnan(a[2]) and np.isnan(b[2])
    # the other entries are untouched
    assert not np.isnan(a[0]) and not np.isnan(b[3])


# --- reader (end to end on the bundled flux file) --------------------------


def test_read_pfkmu_end_to_end():
    power = fitter.read_pfkmu(FLUX_FILE)
    assert power.k_array.shape[0] == 2
    assert len(power.power_array) == power.k_array.shape[1]
