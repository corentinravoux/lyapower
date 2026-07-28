"""Unit / end-to-end tests for :mod:`lyapower.power_spectra` (gimlet I/O)."""

import numpy as np
import pytest

from lyapower import power_spectra

from conftest import MATTER_FILE


# --- gimlet readers (end to end on the bundled data) -----------------------


def test_matter_from_gimlet(matter_spectrum):
    assert matter_spectrum.k_array.ndim == 1
    assert len(matter_spectrum.k_array) == 60
    assert np.all(matter_spectrum.power_array > 0)
    assert matter_spectrum.h_normalized is True


def test_flux_from_gimlet(flux_spectrum):
    # 3D flux P(k, mu): k_array is stacked [k, mu]
    assert flux_spectrum.k_array.shape[0] == 2
    assert flux_spectrum.k_array.shape[1] == len(flux_spectrum.power_array)
    assert flux_spectrum.power_array.ndim == 1


def test_init_from_ascii(tmp_path):
    # a plain two-column file (first row skipped as an edge/header row)
    path = tmp_path / "pk.txt"
    k = np.linspace(0.1, 1.0, 11)
    p = 1.0 / k
    np.savetxt(str(path), np.column_stack([k, p]))
    spec = power_spectra.PowerSpectrum.init_from_ascii_file(str(path))
    np.testing.assert_allclose(spec.k_array, k[1:])
    np.testing.assert_allclose(spec.power_array, p[1:])


# --- transforms ------------------------------------------------------------


def test_center_wavenumbers_1d(matter_spectrum):
    n = len(matter_spectrum.k_array)
    matter_spectrum.center_wavenumbers_1d()
    assert len(matter_spectrum.k_array) == n
    assert matter_spectrum.edge_stored is False


def test_rebin_arrays(matter_spectrum):
    matter_spectrum.rebin_arrays(10, operation="mean")
    assert len(matter_spectrum.k_array) == 10
    assert len(matter_spectrum.power_array) == 10


def test_cut_extremum_flux(flux_spectrum):
    n0 = len(flux_spectrum.power_array)
    flux_spectrum.cut_extremum(0.2, 3.0)
    assert len(flux_spectrum.power_array) < n0
    assert flux_spectrum.k_array[0].min() >= 0.2
    assert flux_spectrum.k_array[0].max() <= 3.0


def test_change_k_normalization(matter_spectrum):
    k0 = matter_spectrum.k_array.copy()
    # h/Mpc -> 1/Mpc divides by h, and the reverse multiplies by h
    matter_spectrum.change_k_normalization(wanted_h_normalized=False, h=0.7)
    assert matter_spectrum.h_normalized is False
    np.testing.assert_allclose(matter_spectrum.k_array, k0 / 0.7)
    matter_spectrum.change_k_normalization(wanted_h_normalized=True, h=0.7)
    np.testing.assert_allclose(matter_spectrum.k_array, k0)
