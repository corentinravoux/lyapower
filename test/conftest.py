"""Shared fixtures for the lyapower test suite.

Tests exercise the pure model/utility functions plus the gimlet readers against
the two files shipped in ``data_example`` (a matter P3D and a flux P(k, mu)),
so no simulation outputs or network access are needed. Matplotlib uses the Agg
backend.
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pytest

TEST_DIR = Path(__file__).resolve().parent
PACKAGE_ROOT = TEST_DIR.parent
DATA_DIR = PACKAGE_ROOT / "data_example"
MATTER_FILE = str(DATA_DIR / "plt00216rhom_ps3d.txt")
FLUX_FILE = str(DATA_DIR / "plt00216z_flux_pkmu.txt")

for _p in (str(PACKAGE_ROOT), str(TEST_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)


@pytest.fixture
def matter_spectrum():
    """A :class:`MatterPowerSpectrum` read from the bundled gimlet matter file."""
    from lyapower import power_spectra

    return power_spectra.MatterPowerSpectrum.init_from_gimlet(
        MATTER_FILE, specie="matter"
    )


@pytest.fixture
def flux_spectrum():
    """A :class:`FluxPowerSpectrum` read from the bundled gimlet flux P(k,mu) file."""
    from lyapower import power_spectra

    return power_spectra.FluxPowerSpectrum.init_3D_from_gimlet(
        FLUX_FILE, type_file="txt", kmu=True
    )
