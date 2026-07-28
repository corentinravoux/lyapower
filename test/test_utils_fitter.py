"""Unit tests for :mod:`lyapower.utils_fitter`."""

import numpy as np
import pytest

from lyapower import utils_fitter as utils


def test_return_key():
    assert utils.return_key({"a": 1}, "a", 99) == 1
    assert utils.return_key({"a": 1}, "b", 99) == 99


def test_latex_float():
    assert utils.latex_float(100.0) == "1 \\times 10^{2}"
    assert "times" not in utils.latex_float(1.5)


def test_bin_ndarray_docstring_example():
    m = np.arange(0, 100, 1).reshape((10, 10))
    n = utils.bin_ndarray(m, new_shape=(5, 5), operation="sum")
    expected = np.array(
        [
            [22, 30, 38, 46, 54],
            [102, 110, 118, 126, 134],
            [182, 190, 198, 206, 214],
            [262, 270, 278, 286, 294],
            [342, 350, 358, 366, 374],
        ]
    )
    np.testing.assert_array_equal(n, expected)


def test_bin_ndarray_errors():
    m = np.ones((4, 4))
    with pytest.raises(ValueError):
        utils.bin_ndarray(m, (2, 2), operation="bogus")
    with pytest.raises(ValueError):
        utils.bin_ndarray(m, (2, 2, 2), operation="mean")


# --- error estimators ------------------------------------------------------


def test_error_estimator_uncorrelated():
    power = np.array([2.0, 4.0])
    err = utils.error_estimator(
        power, model="uncorrelated", bin_count=np.array([4.0, 16.0]), epsilon=0.1
    )
    # power * (1/sqrt(N) + epsilon)
    np.testing.assert_allclose(err, power * (1 / np.sqrt([4.0, 16.0]) + 0.1))


def test_error_estimator_constant():
    power = np.array([2.0, 4.0])
    err = utils.error_estimator(power, model="constant", epsilon=0.2)
    np.testing.assert_allclose(err, power * 0.2)


def test_error_estimator_computed_epsilon():
    power = np.array([2.0, 4.0])
    err = utils.error_estimator(
        power, model="computed_epsilon", bin_count=np.array([0.5, 0.5]), epsilon=0.1
    )
    np.testing.assert_allclose(err, np.array([0.5, 0.5]) + 0.1 * power)


def test_error_estimator_unknown_raises():
    with pytest.raises(KeyError):
        utils.error_estimator(np.array([1.0]), model="does_not_exist")


def test_error_estimator_missing_inputs():
    with pytest.raises(KeyError):
        utils.error_estimator(np.array([1.0]), model="uncorrelated")  # no bin_count
