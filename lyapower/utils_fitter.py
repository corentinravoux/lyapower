#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Date : 17/05/2019

Author: Corentin Ravoux

Description :
"""


#############################################################################
#############################################################################
########################## MODULE IMPORTATION ###############################
#############################################################################
#############################################################################


import logging
import multiprocessing as mp
import time
from functools import partial

import h5py
import numpy as np

#############################################################################
#############################################################################
############################### FUNCTIONS ###################################
#############################################################################
#############################################################################

lambdaLy = 1215.673123130217


def create_log(name="Python_Report", log_level="info"):
    log = Logger(name=name, log_level=log_level)
    log.setup_logging()
    return log


def create_report_log(name="Python_Report", log_level="info"):
    log = Logger(name=name, log_level=log_level)
    log.setup_report_logging()
    return log


_logging_handler = None


class Logger(object):
    def __init__(self, name="Python_Report", log_level="info"):
        self.name = name
        self.log_level = log_level

    def setup_logging(self):
        """
        Taken from https://nbodykit.readthedocs.io/
        Turn on logging, with the specified level.
        Parameters
        ----------
        log_level : 'info', 'debug', 'warning'
                the logging level to set; logging below this level is ignored
        """

        # This gives:
        #
        # [ 000000.43 ]   0: 06-28 14:49  measurestats	INFO	 Nproc = [2, 1, 1]
        # [ 000000.43 ]   0: 06-28 14:49  measurestats	INFO	 Rmax = 120

        levels = {
            "info": logging.INFO,
            "debug": logging.DEBUG,
            "warning": logging.WARNING,
        }

        logger = logging.getLogger()
        t0 = time.time()

        class Formatter(logging.Formatter):
            def format(self, record):
                s1 = "[ %09.2f ]: " % (time.time() - t0)
                return s1 + logging.Formatter.format(self, record)

        fmt = Formatter(
            fmt="%(asctime)s %(name)-15s %(levelname)-8s %(message)s",
            datefmt="%m-%d %H:%M ",
        )

        global _logging_handler
        if _logging_handler is None:
            _logging_handler = logging.StreamHandler()
            logger.addHandler(_logging_handler)

        _logging_handler.setFormatter(fmt)
        logger.setLevel(levels[self.log_level])

    def setup_report_logging(self):
        levels = {
            "info": logging.INFO,
            "debug": logging.DEBUG,
            "warning": logging.WARNING,
        }
        logging.basicConfig(
            filename=self.name,
            filemode="w",
            level=levels[self.log_level],
            format="%(asctime)s :: %(levelname)s :: %(message)s",
        )

    @staticmethod
    def add(line, level="info"):
        if level == "info":
            logging.info(line)
        if level == "warning":
            logging.warning(line)
        if level == "debug":
            logging.debug(line)

    @staticmethod
    def close():
        logging.shutdown()


def latex_float(float_input, decimals_input="{0:.2g}"):
    """
    example use:
    import matplotlib.pyplot as plt
    plt.figure(),plt.clf()
    plt.plot(np.array([1,2.]),'ko-',label="$P_0="+latex_float(7.63e-5)+'$'),
    plt.legend()
    """
    float_str = decimals_input.format(float_input)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str


def return_key(dictionary, string, default_value):
    return dictionary[string] if string in dictionary.keys() else default_value


def error_estimator(power, model="uncorrelated", **kwargs):
    if model == "uncorrelated":
        return error_estimator_uncorrelated(power, **kwargs)
    elif model == "constant":
        return error_estimator_constant(power, **kwargs)
    elif model == "computed":
        return error_estimator_computed(power, **kwargs)
    elif model == "computed_epsilon":
        return error_estimator_computed_epsilon(power, **kwargs)
    else:
        raise KeyError("model of error estimator not available")


def error_estimator_uncorrelated(power, **kwargs):
    epsilon = return_key(kwargs, "epsilon", None)
    bin_count = return_key(kwargs, "bin_count", None)
    if (bin_count is None) | (epsilon is None):
        raise KeyError("Need bin_count and epsilon")
    return power * ((1 / np.sqrt(bin_count)) + epsilon)


def error_estimator_constant(power, **kwargs):
    epsilon = return_key(kwargs, "epsilon", None)
    if epsilon is None:
        raise KeyError("Need bin_count and epsilon")
    return power * epsilon


def error_estimator_computed(power, **kwargs):
    bin_count = return_key(kwargs, "bin_count", None)
    if bin_count is None:
        raise KeyError("Need bin_count")
    return bin_count


def error_estimator_computed_epsilon(power, **kwargs):
    bin_count = return_key(kwargs, "bin_count", None)
    epsilon = return_key(kwargs, "epsilon", None)
    if bin_count is None:
        raise KeyError("Need bin_count")
    return bin_count + epsilon * power


def bin_ndarray(ndarray, new_shape, operation="mean"):
    """
    From : https://stackoverflow.com/questions/8090229/resize-with-averaging-or-rebin-a-numpy-2d-array/29042041
    Bins an ndarray in all axes based on the target shape, by summing or
    averaging.
    Number of output dimensions must match number of input dimensions.
    Example
    -------
    >>> m = np.arange(0,100,1).reshape((10,10))
    >>> n = bin_ndarray(m, new_shape=(5,5), operation='sum')
    >>> print(n)
    [[ 22  30  38  46  54]
    [102 110 118 126 134]
    [182 190 198 206 214]
    [262 270 278 286 294]
    [342 350 358 366 374]]
    """
    if not operation.lower() in ["sum", "mean", "average", "avg", "gauss"]:
        raise ValueError("Operation {} not supported.".format(operation))
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape, new_shape))
    compression_pairs = [(d, c // d) for d, c in zip(new_shape, ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    ndarray = ndarray.reshape(flattened)
    for i in range(len(new_shape)):
        if operation.lower() == "sum":
            ndarray = ndarray.sum(-1 * (i + 1))
        elif operation.lower() in ["mean", "average", "avg"]:
            ndarray = ndarray.mean(-1 * (i + 1))
        elif operation.lower() in ["gauss"]:
            if i != 0:
                raise KeyError("gaussian mean is not available for dim higher than 1")
            from scipy import signal

            newndarray = np.zeros(new_shape)
            gaussian_weights = signal.gaussian(
                int(ndarray.shape[1]), int(ndarray.shape[1]) / 4
            )
            for j in range(len(ndarray)):
                newndarray[j] = np.average(ndarray[j], axis=0, weights=gaussian_weights)
            ndarray = newndarray
    return ndarray


def rebin_slice(
    sim_name,
    new_shape_slice,
    index_rescaling,
    index,
    operation="mean",
    first_field="derived_fields",
    second_field="tau_red",
    transform=None,
):
    print("Treating ", index)
    sim = h5py.File(sim_name)
    full_slice = sim[first_field][second_field][
        index * index_rescaling : (index + 1) * index_rescaling, :, :
    ]
    if transform is not None:
        full_slice = transform(full_slice)
    return bin_ndarray(full_slice, new_shape_slice, operation=operation)


def rebin_simulation(
    sim_name,
    index_rescaling,
    operation="mean",
    number_worker=1,
    first_field="derived_fields",
    second_field="tau_red",
    transform_name=None,
):
    sim = h5py.File(sim_name)
    shape_sim = sim["domain"].attrs["shape"]
    shape_rebinned_field = (shape_sim / index_rescaling).astype(int)
    rebinned_field = np.zeros(shape_rebinned_field)
    new_shape_slice = (1, shape_rebinned_field[1], shape_rebinned_field[2])
    if transform_name == "exp":
        transform = lambda x: np.exp(-x)
    else:
        transform = None
    func = partial(
        rebin_slice,
        sim_name,
        new_shape_slice,
        index_rescaling,
        operation=operation,
        first_field=first_field,
        second_field=second_field,
        transform=transform,
    )
    if number_worker == 1:
        for i in range(shape_rebinned_field[0]):
            print("Treating ", i)
            rebinned_field[i, :, :] = func(i)
    else:

        with mp.Pool(number_worker) as p:
            results = p.map(func, np.arange(shape_rebinned_field[0]))
        rebinned_field = np.array(results)
    return rebinned_field
