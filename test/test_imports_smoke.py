"""Smoke tests: every shipped module byte-compiles and imports."""

import importlib
import py_compile

import pytest

from conftest import PACKAGE_ROOT

MODULES = ["utils_fitter", "power_spectra", "fitter", "CLASS"]

ALL_PY = sorted((PACKAGE_ROOT / "lyapower").glob("*.py")) \
    + sorted((PACKAGE_ROOT / "launchers").glob("*.py"))


@pytest.mark.parametrize("path", ALL_PY, ids=[p.name for p in ALL_PY])
def test_module_compiles(path):
    """Every shipped .py file byte-compiles (no syntax errors)."""
    py_compile.compile(str(path), doraise=True)


@pytest.mark.parametrize("module", MODULES)
def test_module_imports(module):
    """Each package module imports (classy is imported lazily where needed)."""
    importlib.import_module(f"lyapower.{module}")
