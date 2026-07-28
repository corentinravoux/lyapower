"""Fitting routines for the 3D Lyman-alpha forest flux power spectrum.

This module implements the flux power spectrum model of
Arinyo-i-Prats et al. (2015), including the non-linear correction terms
D0 and D1 (and a BAO-damped D1 variant), builds ``Pf_model`` from a
linear matter power spectrum (obtained via CLASS or cosmoprimo, or read
directly from Nyx/gimlet outputs), sets up ``iminuit``-based cost
functions and Minuit fits, drives the end-to-end fitting pipeline
(``prepare_data``, ``fitter_k_mu``), reads gimlet power-spectrum output
files, and produces diagnostic plots and LaTeX summaries of fit results.

Conventions: wavenumbers ``k`` are in h/Mpc, ``mu = k_par / k``; ``b``
and ``beta`` are the linear bias and RSD parameter; ``non_linear_model``
selects between the "0" (D0), "1" (D1), "1_bao"/"1_BAO" (D1 + BAO
damping) and None (linear) flux power spectrum models.
"""

import os
import pickle

import numpy as np
import scipy.interpolate
from iminuit import Minuit
from scipy import integrate

from lyapower import CLASS, power_spectra, utils_fitter


def read_pfkmu_hdf5(
    filename, field_name, power_weighted=False, error_estimator=None, **kwargs
):
    """Read a 3D flux power spectrum P(k, mu) from an HDF5 gimlet output file.

    Args:
        filename (str): Path to the HDF5 file produced by gimlet.
        field_name (str): Name of the HDF5 field/dataset to read.
        power_weighted (bool): Whether the stored k/mu values are already
            power-weighted bin centers. Defaults to False.
        error_estimator (str or None): Name of the error estimator to use
            when building the error array. Defaults to None.
        **kwargs: Additional keyword arguments forwarded to
            `power_spectra.FluxPowerSpectrum.init_3D_from_gimlet`.

    Returns:
        power_spectra.FluxPowerSpectrum: The loaded 3D flux power spectrum.
    """
    power = power_spectra.FluxPowerSpectrum.init_3D_from_gimlet(
        filename,
        "hdf5",
        kmu=True,
        power_weighted=power_weighted,
        error_estimator=error_estimator,
        field_name=field_name,
        **kwargs,
    )
    return power


def read_pfkmu(filename, power_weighted=False, error_estimator=None, **kwargs):
    """Read a 3D flux power spectrum P(k, mu) from a text gimlet output file.

    Args:
        filename (str): Path to the text file produced by gimlet.
        power_weighted (bool): Whether the stored k/mu values are already
            power-weighted bin centers. Defaults to False.
        error_estimator (str or None): Name of the error estimator to use
            when building the error array. Defaults to None.
        **kwargs: Additional keyword arguments forwarded to
            `power_spectra.FluxPowerSpectrum.init_3D_from_gimlet`.

    Returns:
        power_spectra.FluxPowerSpectrum: The loaded 3D flux power spectrum,
        parametrized in (k, mu).
    """
    power = power_spectra.FluxPowerSpectrum.init_3D_from_gimlet(
        filename,
        "txt",
        kmu=True,
        power_weighted=power_weighted,
        error_estimator=error_estimator,
        **kwargs,
    )
    return power


def read_pfkperpkpar(filename, power_weighted=False, error_estimator=None, **kwargs):
    """Read a 3D flux power spectrum P(k_perp, k_par) from a text gimlet file.

    Args:
        filename (str): Path to the text file produced by gimlet.
        power_weighted (bool): Whether the stored k values are already
            power-weighted bin centers. Defaults to False.
        error_estimator (str or None): Name of the error estimator to use
            when building the error array. Defaults to None.
        **kwargs: Additional keyword arguments forwarded to
            `power_spectra.FluxPowerSpectrum.init_3D_from_gimlet`.

    Returns:
        power_spectra.FluxPowerSpectrum: The loaded 3D flux power spectrum,
        parametrized in (k_perp, k_par).
    """
    power = power_spectra.FluxPowerSpectrum.init_3D_from_gimlet(
        filename,
        "txt",
        kmu=False,
        power_weighted=power_weighted,
        error_estimator=error_estimator,
        **kwargs,
    )
    return power


def read_p1d(filename, power_weighted=False, error_estimator=None, **kwargs):
    """Read a 1D flux power spectrum P(k) from a gimlet output file.

    Args:
        filename (str): Path to the file produced by gimlet.
        power_weighted (bool): Whether the stored k values are already
            power-weighted bin centers. Defaults to False.
        error_estimator (str or None): Name of the error estimator to use
            when building the error array. Defaults to None.
        **kwargs: Additional keyword arguments forwarded to
            `power_spectra.FluxPowerSpectrum.init_1D_from_gimlet`.

    Returns:
        power_spectra.FluxPowerSpectrum: The loaded 1D flux power spectrum.
    """
    power = power_spectra.FluxPowerSpectrum.init_1D_from_gimlet(
        filename,
        power_weighted=power_weighted,
        error_estimator=error_estimator,
        **kwargs,
    )
    return power


def read_pk(filename, power_weighted=False, error_estimator=None, **kwargs):
    """Read a matter power spectrum P(k) from a gimlet output file.

    Args:
        filename (str): Path to the file produced by gimlet.
        power_weighted (bool): Whether the stored k values are already
            power-weighted bin centers. Defaults to False.
        error_estimator (str or None): Name of the error estimator to use
            when building the error array. Defaults to None.
        **kwargs: Additional keyword arguments forwarded to
            `power_spectra.MatterPowerSpectrum.init_from_gimlet`.

    Returns:
        power_spectra.MatterPowerSpectrum: The loaded matter power spectrum.
    """
    power = power_spectra.MatterPowerSpectrum.init_from_gimlet(
        filename,
        power_weighted=power_weighted,
        error_estimator=error_estimator,
        **kwargs,
    )
    return power


def rebin_matter_power(power_m, k_m, k_f):
    """Interpolate a matter power spectrum onto a new set of wavenumbers.

    Args:
        power_m (numpy.ndarray): Matter power spectrum values sampled at `k_m`.
        k_m (numpy.ndarray): Wavenumbers at which `power_m` is sampled.
        k_f (numpy.ndarray): Target wavenumbers (typically the flux power
            spectrum k-grid) to interpolate onto.

    Returns:
        numpy.ndarray: `power_m` linearly interpolated onto `k_f`; values
        outside the range of `k_m` are set to NaN.
    """
    power = scipy.interpolate.interp1d(
        k_m, power_m, bounds_error=False, fill_value=np.nan
    )
    power_m_rebin = power(k_f)
    return power_m_rebin


def mask_data(indexes, *args):
    """Mask near-zero entries in-place across a set of arrays.

    For each array indexed in `indexes`, flags array elements whose
    absolute value is smaller than 1e-10 times the array's mean, then
    sets those positions to NaN in every array in `args` (in-place).

    Args:
        indexes (iterable of int): Indices into `args` used to build the
            mask (elements are flagged as invalid if any of these arrays
            is near zero at that position).
        *args (numpy.ndarray): Arrays to be masked in-place; all must
            share the same shape as `args[0]`.

    Returns:
        None
    """
    mask = np.full(args[0].shape, False)
    for i in indexes:
        mask |= np.abs(args[i]) < np.abs(10**-10 * np.mean(args[i]))
    for i in range(len(args)):
        args[i][mask] = np.nan


### models & cost


def D0(k, mu, k_nl, a_nl, k_p, a_p, k_v0, a_v0, k_v1, a_v1):
    """Arinyo-i-Prats et al. (2015) "D0" non-linear correction term.

    Args:
        k (numpy.ndarray or float): Wavenumber, in h/Mpc.
        mu (numpy.ndarray or float): Cosine of the angle to the line of sight.
        k_nl (float): Non-linear growth scale.
        a_nl (float): Non-linear growth exponent.
        k_p (float): Pressure (smoothing) scale.
        a_p (float): Pressure exponent.
        k_v0 (float): Velocity non-linear scale normalization.
        a_v0 (float): Velocity non-linear exponent.
        k_v1 (float): Velocity scale used in the k-dependence of `k_v0`.
        a_v1 (float): Exponent of the k-dependence of `k_v0`.

    Returns:
        numpy.ndarray or float: The multiplicative non-linear correction
        D0(k, mu).
    """
    return np.exp(
        (k / k_nl) ** a_nl
        - (k / k_p) ** a_p
        - ((k * mu) / (k_v0 * (1 + (k / k_v1)) ** a_v1)) ** a_v0
    )


def D1(k, mu, q_1, q_2, k_v, a_v, b_v, k_p, linear_power_spectrum):
    """Arinyo-i-Prats et al. (2015) "D1" non-linear correction term.

    Args:
        k (numpy.ndarray or float): Wavenumber, in h/Mpc.
        mu (numpy.ndarray or float): Cosine of the angle to the line of sight.
        q_1 (float): Linear coefficient of the non-linear growth term.
        q_2 (float): Quadratic coefficient of the non-linear growth term.
        k_v (float): Velocity non-linear scale.
        a_v (float): Velocity non-linear exponent (k-dependence).
        b_v (float): Velocity non-linear exponent (mu-dependence).
        k_p (float): Pressure (smoothing) scale.
        linear_power_spectrum (numpy.ndarray or float): Linear matter
            power spectrum evaluated at `k`, used to build the
            dimensionless variance Delta^2(k).

    Returns:
        numpy.ndarray or float: The multiplicative non-linear correction
        D1(k, mu).
    """
    Delta_square = (1 / (2 * np.pi**2)) * k**3 * linear_power_spectrum
    non_linear_term = np.exp(
        (q_1 * Delta_square + q_2 * Delta_square**2)
        * (1 - ((k / k_v) ** a_v) * mu**b_v)
        - (k / k_p) ** 2
    )
    return non_linear_term


def Pl_class(k_array, settings, z, name="class"):
    """Compute the linear matter power spectrum at redshift `z` using CLASS.

    Args:
        k_array (numpy.ndarray): Wavenumbers (h/Mpc) spanning the desired
            range; only its min/max and length are used to set up the
            CLASS k-grid.
        settings (dict): CLASS configuration/settings dictionary passed
            to `CLASS.MyClass`.
        z (float): Redshift at which to evaluate the power spectrum.
        name (str): Output file base name used by `write_pk_tk`.
            Defaults to "class".

    Returns:
        power_spectra.MatterPowerSpectrum: The linear matter power
        spectrum computed by CLASS, h-normalized.
    """
    my_class = CLASS.MyClass(os.getcwd(), settings)
    kmin, kmax, nb_points = (
        np.log10(np.min(k_array)),
        np.log10(np.max(k_array)),
        2 * len(k_array),
    )
    (Power, _, _) = my_class.write_pk_tk(
        z, name, kmin=kmin, kmax=kmax, nb_points=nb_points, output=False
    )
    h_normalized = True
    power = power_spectra.MatterPowerSpectrum(
        k_array=Power[:, 0],
        power_array=Power[:, 1],
        dimension="1D",
        specie="matter",
        h_normalized=h_normalized,
    )
    return power


def Pl_cosmoprimo(k_array, settings, z):
    """Compute the (wiggle and no-wiggle) linear matter power spectra via cosmoprimo.

    Args:
        k_array (numpy.ndarray): Wavenumbers (h/Mpc) at which to evaluate
            the power spectra.
        settings (dict): Cosmoprimo/CLASS configuration dictionary passed
            to `CLASS.CosmoprimoInterface`.
        z (float): Redshift at which to evaluate the power spectra.

    Returns:
        tuple: `(power, power_no_bao)`, both
        `power_spectra.MatterPowerSpectrum` instances (h-normalized):
        the full (BAO wiggles included) linear power spectrum and its
        no-wiggle counterpart, both evaluated on `k_array`.
    """
    cosmoprimo = CLASS.CosmoprimoInterface(os.getcwd(), settings)
    Power_no_bao = cosmoprimo.Pl_class_cosmoprimo_no_wiggle(k_array, z)
    Power_bao = cosmoprimo.Pl_class_cosmoprimo(k_array, z)

    h_normalized = True
    power = power_spectra.MatterPowerSpectrum(
        k_array=k_array,
        power_array=Power_bao,
        dimension="1D",
        specie="matter",
        h_normalized=h_normalized,
    )
    power_no_bao = power_spectra.MatterPowerSpectrum(
        k_array=k_array,
        power_array=Power_no_bao,
        dimension="1D",
        specie="matter",
        h_normalized=h_normalized,
    )

    return power, power_no_bao


def Pm_normalized(pm_file, class_dict, z_simu, z_init, name="pmnorm"):
    """Build a matter power spectrum from simulation output, rescaled by CLASS growth.

    Reads a raw simulation matter power spectrum from `pm_file` and
    rescales it by the ratio of CLASS linear power spectra at `z_simu`
    and `z_init` (times an Omega_m/(Omega_m - Omega_b) coefficient), to
    correct for baryon growth suppression / normalize onto a fiducial
    linear growth.

    Args:
        pm_file (str): Path to the gimlet matter power spectrum file.
        class_dict (dict): CLASS configuration dictionary.
        z_simu (float): Simulation output redshift.
        z_init (float): Simulation initial-conditions redshift.
        name (str): Output file base name used by `write_pk_tk`.
            Defaults to "pmnorm".

    Returns:
        power_spectra.MatterPowerSpectrum: The rescaled, h-normalized
        matter power spectrum, sampled on the same k-grid as `pm_file`.
    """
    power_m = read_pk(pm_file)
    k_array = power_m.k_array
    my_class = CLASS.MyClass(os.getcwd(), class_dict)
    kmin, kmax, nb_points = (
        np.log10(np.min(k_array)),
        np.log10(np.max(k_array)),
        4 * len(k_array),
    )
    (Power, _, _) = my_class.write_pk_tk(
        z_simu,
        name,
        kmin=kmin,
        kmax=kmax,
        nb_points=nb_points,
        output=False,
        verbose=False,
    )
    (Power_init, _, _) = my_class.write_pk_tk(
        z_init,
        name,
        kmin=kmin,
        kmax=kmax,
        nb_points=nb_points,
        output=False,
        verbose=False,
    )
    interp_power = scipy.interpolate.interp1d(
        Power[:, 0], Power[:, 1], bounds_error=False, fill_value=np.nan
    )
    interp_power_init = scipy.interpolate.interp1d(
        Power_init[:, 0], Power_init[:, 1], bounds_error=False, fill_value=np.nan
    )
    coeff = my_class.model.Omega_m() / (
        my_class.model.Omega_m() - my_class.model.Omega_b()
    )
    power = (
        coeff
        * power_m.power_array
        * (interp_power(k_array) / interp_power_init(k_array))
    )
    h_normalized = True
    power = power_spectra.MatterPowerSpectrum(
        k_array=k_array,
        power_array=power,
        dimension="1D",
        specie="matter",
        h_normalized=h_normalized,
    )
    return power


def Pm(pm_file, name="pm"):
    """Read a raw matter power spectrum from a gimlet file and set its name.

    Args:
        pm_file (str): Path to the gimlet matter power spectrum file.
        name (str): Name to assign to the returned power spectrum object.
            Defaults to "pm".

    Returns:
        power_spectra.MatterPowerSpectrum: The matter power spectrum read
        from `pm_file`.
    """
    power_m = read_pk(pm_file)
    power_m.name = name
    return power_m


def Pf_model(
    linear_power_spectrum,
    non_linear_model="0",
    linear_power_spectrum_no_bao=None,
    integrate_model=True,
    N_mu_integration=1000,
    mu_max=1.0,
):
    """Build the flux power spectrum model function P_f(k, mu).

    Selects and returns a callable model of the flux power spectrum,
    P_f = b^2 (1 + beta*mu^2)^2 * P_linear * D(k, mu), where D is either
    the Arinyo-i-Prats "D0" or "D1" non-linear correction (optionally
    with BAO damping via `linear_power_spectrum_no_bao`), or omitted
    entirely for a linear model (`non_linear_model=None`). If
    `integrate_model` is True, the returned model integrates the
    analytic expression over each mu bin (using Simpson's rule) instead
    of evaluating it at bin-center mu values.

    Args:
        linear_power_spectrum (numpy.ndarray): Linear matter power
            spectrum evaluated on the flux power spectrum k-grid.
        non_linear_model (str or None): Which non-linear correction to
            use: "0" (D0), "1" (D1), "1_bao"/"1_BAO" (D1 + BAO damping),
            or None for a purely linear model. Defaults to "0".
        linear_power_spectrum_no_bao (numpy.ndarray or None): No-wiggle
            linear matter power spectrum, required for the BAO-damped
            variant. Defaults to None.
        integrate_model (bool): If True, integrate the model over each
            mu bin instead of evaluating at bin centers. Defaults to
            True.
        N_mu_integration (int): Number of mu sub-samples used per bin
            when `integrate_model` is True. Defaults to 1000.
        mu_max (float): Maximum mu value used to compute mu bin widths.
            Defaults to 1.0.

    Returns:
        callable: A model function `model(x, b, beta, ...)`, where
        `x = (k, mu)` and the remaining positional arguments are the
        non-linear model's nuisance parameters, returning the predicted
        flux power spectrum.
    """
    if integrate_model:
        if non_linear_model == "0":

            def modelD0(x, b, beta, k_nl, a_nl, k_p, a_p, k_v0, a_v0, k_v1, a_v1):
                """Flux power spectrum model with D0 non-linear correction, integrated over mu.

                Args:
                    x (tuple): `(k, mu)` arrays of wavenumbers and mu
                        bin-start values.
                    b (float): Linear bias.
                    beta (float): RSD parameter.
                    k_nl, a_nl, k_p, a_p, k_v0, a_v0, k_v1, a_v1 (float):
                        D0 non-linear correction parameters (see `D0`).

                Returns:
                    numpy.ndarray: The flux power spectrum evaluated at
                    `k`, integrated over each mu bin `[mu, mu + dmu]` via
                    Simpson's rule.
                """
                k, mu = x[0], x[1]
                mu_next_bin = mu + power_spectra.PowerSpectrum.compute_dmu(
                    mu, mu_max=mu_max
                )
                kmu = np.array(
                    [
                        np.transpose(np.tile(k, (N_mu_integration, 1))),
                        np.array(
                            [
                                np.linspace(mu[i], mu_next_bin[i], N_mu_integration)
                                for i in range(len(k))
                            ]
                        ),
                    ]
                )
                linear_power_spectrum_repeat = np.transpose(
                    np.tile(linear_power_spectrum, (N_mu_integration, 1))
                )

                def integrand(kmu):
                    """Integrand b^2 (1 + beta*mu^2)^2 P_lin D0(k, mu) for the mu integration.

                    Args:
                        kmu (numpy.ndarray): Array of shape `(2, ...)`
                            holding the k and mu sample grids over which
                            to evaluate the integrand.

                    Returns:
                        numpy.ndarray: The integrand values on the `kmu`
                        grid.
                    """
                    return (
                        b**2
                        * (1 + beta * kmu[1] ** 2) ** 2
                        * linear_power_spectrum_repeat
                        * D0(
                            kmu[0], kmu[1], k_nl, a_nl, k_p, a_p, k_v0, a_v0, k_v1, a_v1
                        )
                    )

                integrand_kmu = integrand(kmu)
                Pf_integrated = np.array(
                    [
                        integrate.simpson(integrand_kmu[i], x=kmu[1][i])
                        / (mu_next_bin[i] - mu[i])
                        for i in range(len(k))
                    ]
                )

                return Pf_integrated

            return modelD0

        elif non_linear_model == "1":

            def modelD1(x, b, beta, q_1, q_2, k_v, a_v, b_v, k_p):
                """Flux power spectrum model with D1 non-linear correction, integrated over mu.

                Args:
                    x (tuple): `(k, mu)` arrays of wavenumbers and mu
                        bin-start values.
                    b (float): Linear bias.
                    beta (float): RSD parameter.
                    q_1, q_2, k_v, a_v, b_v, k_p (float): D1 non-linear
                        correction parameters (see `D1`).

                Returns:
                    numpy.ndarray: The flux power spectrum evaluated at
                    `k`, integrated over each mu bin `[mu, mu + dmu]` via
                    Simpson's rule.
                """
                k, mu = x[0], x[1]
                mu_next_bin = mu + power_spectra.PowerSpectrum.compute_dmu(
                    mu, mu_max=mu_max
                )
                kmu = np.array(
                    [
                        np.transpose(np.tile(k, (N_mu_integration, 1))),
                        np.array(
                            [
                                np.linspace(mu[i], mu_next_bin[i], N_mu_integration)
                                for i in range(len(k))
                            ]
                        ),
                    ]
                )
                linear_power_spectrum_repeat = np.transpose(
                    np.tile(linear_power_spectrum, (N_mu_integration, 1))
                )

                def integrand(kmu):
                    """Integrand b^2 (1 + beta*mu^2)^2 P_lin D1(k, mu) for the mu integration.

                    Args:
                        kmu (numpy.ndarray): Array of shape `(2, ...)`
                            holding the k and mu sample grids over which
                            to evaluate the integrand.

                    Returns:
                        numpy.ndarray: The integrand values on the `kmu`
                        grid.
                    """
                    return (
                        b**2
                        * (1 + beta * kmu[1] ** 2) ** 2
                        * linear_power_spectrum_repeat
                        * D1(
                            kmu[0],
                            kmu[1],
                            q_1,
                            q_2,
                            k_v,
                            a_v,
                            b_v,
                            k_p,
                            linear_power_spectrum_repeat,
                        )
                    )

                integrand_kmu = integrand(kmu)
                Pf_integrated = np.array(
                    [
                        integrate.simpson(integrand_kmu[i], x=kmu[1][i])
                        / (mu_next_bin[i] - mu[i])
                        for i in range(len(k))
                    ]
                )

                return Pf_integrated

            return modelD1

        elif non_linear_model == "1_bao":  # LUCAS' MODEL
            wiggle_power_spectrum = linear_power_spectrum - linear_power_spectrum_no_bao

            def modelD1BAO(x, b, beta, q_1, q_2, k_v, a_v, b_v, k_p, S_p, S_t):
                """Flux power spectrum model with BAO-damped D1 correction, integrated over mu.

                Args:
                    x (tuple): `(k, mu)` arrays of wavenumbers and mu
                        bin-start values.
                    b (float): Linear bias.
                    beta (float): RSD parameter.
                    q_1, q_2, k_v, a_v, b_v, k_p (float): D1 non-linear
                        correction parameters (see `D1`).
                    S_p (float): BAO damping scale along the line of sight.
                    S_t (float): BAO damping scale transverse to the line
                        of sight.

                Returns:
                    numpy.ndarray: The flux power spectrum evaluated at
                    `k`, using the BAO-wiggle-damped linear power
                    spectrum, integrated over each mu bin via Simpson's
                    rule.
                """
                k, mu = x[0], x[1]

                Snl = (S_p * mu) ** 2 + (S_t**2 * (1 - mu**2))
                damped_linear_power_spectrum = (
                    linear_power_spectrum_no_bao
                    + wiggle_power_spectrum * np.exp((-((k * Snl) ** 2)) / 2)
                )

                mu_next_bin = mu + power_spectra.PowerSpectrum.compute_dmu(
                    mu, mu_max=mu_max
                )
                kmu = np.array(
                    [
                        np.transpose(np.tile(k, (N_mu_integration, 1))),
                        np.array(
                            [
                                np.linspace(mu[i], mu_next_bin[i], N_mu_integration)
                                for i in range(len(k))
                            ]
                        ),
                    ]
                )
                linear_power_spectrum_repeat = np.transpose(
                    np.tile(linear_power_spectrum, (N_mu_integration, 1))
                )
                damped_linear_power_spectrum_repeat = np.transpose(
                    np.tile(damped_linear_power_spectrum, (N_mu_integration, 1))
                )

                def integrand(kmu):
                    """Integrand b^2 (1 + beta*mu^2)^2 P_damped D1(k, mu) for the mu integration.

                    Args:
                        kmu (numpy.ndarray): Array of shape `(2, ...)`
                            holding the k and mu sample grids over which
                            to evaluate the integrand.

                    Returns:
                        numpy.ndarray: The integrand values on the `kmu`
                        grid, using the BAO-wiggle-damped linear power
                        spectrum.
                    """
                    return (
                        b**2
                        * (1 + beta * kmu[1] ** 2) ** 2
                        * damped_linear_power_spectrum_repeat
                        * D1(
                            kmu[0],
                            kmu[1],
                            q_1,
                            q_2,
                            k_v,
                            a_v,
                            b_v,
                            k_p,
                            linear_power_spectrum_repeat,
                        )
                    )

                integrand_kmu = integrand(kmu)

                Pf_integrated = np.array(
                    [
                        integrate.simpson(integrand_kmu[i], x=kmu[1][i])
                        / (mu_next_bin[i] - mu[i])
                        for i in range(len(k))
                    ]
                )

                return Pf_integrated

            return modelD1BAO

        elif non_linear_model == None:

            def modellinear(x, b, beta):
                """Linear flux power spectrum model, integrated over mu.

                Args:
                    x (tuple): `(k, mu)` arrays of wavenumbers and mu
                        bin-start values.
                    b (float): Linear bias.
                    beta (float): RSD parameter.

                Returns:
                    numpy.ndarray: The flux power spectrum evaluated at
                    `k`, integrated over each mu bin `[mu, mu + dmu]` via
                    Simpson's rule.
                """
                k, mu = x[0], x[1]
                mu_next_bin = mu + power_spectra.PowerSpectrum.compute_dmu(
                    mu, mu_max=mu_max
                )
                kmu = np.array(
                    [
                        np.transpose(np.tile(k, (N_mu_integration, 1))),
                        np.array(
                            [
                                np.linspace(mu[i], mu_next_bin[i], N_mu_integration)
                                for i in range(len(k))
                            ]
                        ),
                    ]
                )
                linear_power_spectrum_repeat = np.transpose(
                    np.tile(linear_power_spectrum, (N_mu_integration, 1))
                )

                def integrand(kmu):
                    """Integrand b^2 (1 + beta*mu^2)^2 P_lin for the mu integration.

                    Args:
                        kmu (numpy.ndarray): Array of shape `(2, ...)`
                            holding the k and mu sample grids over which
                            to evaluate the integrand.

                    Returns:
                        numpy.ndarray: The integrand values on the `kmu`
                        grid.
                    """
                    return (
                        b**2
                        * (1 + beta * kmu[1] ** 2) ** 2
                        * linear_power_spectrum_repeat
                    )

                integrand_kmu = integrand(kmu)
                Pf_integrated = np.array(
                    [
                        integrate.simpson(integrand_kmu[i], x=kmu[1][i])
                        / (mu_next_bin[i] - mu[i])
                        for i in range(len(k))
                    ]
                )
                return Pf_integrated

            return modellinear
    else:
        if non_linear_model == "0":

            def modelD0(x, b, beta, k_nl, a_nl, k_p, a_p, k_v0, a_v0, k_v1, a_v1):
                """Flux power spectrum model with D0 non-linear correction, evaluated at mu.

                Args:
                    x (tuple): `(k, mu)` arrays of wavenumbers and mu
                        values.
                    b (float): Linear bias.
                    beta (float): RSD parameter.
                    k_nl, a_nl, k_p, a_p, k_v0, a_v0, k_v1, a_v1 (float):
                        D0 non-linear correction parameters (see `D0`).

                Returns:
                    numpy.ndarray: The flux power spectrum evaluated
                    pointwise at `(k, mu)`.
                """
                k, mu = x[0], x[1]
                Pf = (
                    b**2
                    * (1 + beta * mu**2) ** 2
                    * linear_power_spectrum
                    * D0(k, mu, k_nl, a_nl, k_p, a_p, k_v0, a_v0, k_v1, a_v1)
                )
                return Pf

            return modelD0
        elif non_linear_model == "1":

            def modelD1(x, b, beta, q_1, q_2, k_v, a_v, b_v, k_p):
                """Flux power spectrum model with D1 non-linear correction, evaluated at mu.

                Args:
                    x (tuple): `(k, mu)` arrays of wavenumbers and mu
                        values.
                    b (float): Linear bias.
                    beta (float): RSD parameter.
                    q_1, q_2, k_v, a_v, b_v, k_p (float): D1 non-linear
                        correction parameters (see `D1`).

                Returns:
                    numpy.ndarray: The flux power spectrum evaluated
                    pointwise at `(k, mu)`.
                """
                k, mu = x[0], x[1]
                Pf = (
                    b**2
                    * (1 + beta * mu**2) ** 2
                    * linear_power_spectrum
                    * D1(k, mu, q_1, q_2, k_v, a_v, b_v, k_p, linear_power_spectrum)
                )
                return Pf

            return modelD1
        elif non_linear_model == "1_BAO":

            def modelD1BAO(x, b, beta, q_1, q_2, k_v, a_v, b_v, k_p, S_p, S_t):
                """Flux power spectrum model with BAO-damped D1 correction, evaluated at mu.

                Args:
                    x (tuple): `(k, mu)` arrays of wavenumbers and mu
                        values.
                    b (float): Linear bias.
                    beta (float): RSD parameter.
                    q_1, q_2, k_v, a_v, b_v, k_p (float): D1 non-linear
                        correction parameters (see `D1`).
                    S_p (float): BAO damping scale along the line of sight.
                    S_t (float): BAO damping scale transverse to the line
                        of sight.

                Returns:
                    numpy.ndarray: The flux power spectrum evaluated
                    pointwise at `(k, mu)`, using the BAO-wiggle-damped
                    linear power spectrum.
                """
                k, mu = x[0], x[1]

                Snl = (S_p * mu) ** 2 + (S_t**2 * (1 - mu**2))
                damped_linear_power_spectrum = (
                    linear_power_spectrum_no_bao
                    + wiggle_power_spectrum * np.exp((-((k * Snl) ** 2)) / 2)
                )
                Pf = (
                    b**2
                    * (1 + beta * mu**2) ** 2
                    * damped_linear_power_spectrum
                    * D1(k, mu, q_1, q_2, k_v, a_v, b_v, k_p, linear_power_spectrum)
                )
                return Pf

        elif non_linear_model == None:

            def modellinear(x, b, beta):
                """Linear flux power spectrum model, evaluated at mu.

                Args:
                    x (tuple): `(k, mu)` arrays; only `mu = x[1]` is used.
                    b (float): Linear bias.
                    beta (float): RSD parameter.

                Returns:
                    numpy.ndarray: The flux power spectrum evaluated
                    pointwise at `mu`.
                """
                mu = x[1]
                Pf = b**2 * (1 + beta * mu**2) ** 2 * linear_power_spectrum
                return Pf

            return modellinear


def custom_least_squares(model, data_x, data_y, data_yerr, non_linear_model="0"):
    """Build a plain least-squares iminuit cost function for a given model.

    Args:
        model (callable): Flux power spectrum model function, as
            returned by `Pf_model`.
        data_x (tuple): `(k, mu)` data coordinates.
        data_y (numpy.ndarray): Observed flux power spectrum values.
        data_yerr (numpy.ndarray): Uncertainties on `data_y`.
        non_linear_model (str or None): Selects which cost function
            signature to return: "0" (D0), "1" (D1), "1_bao" (D1 +
            BAO), or None (linear). Defaults to "0".

    Returns:
        callable: A cost function `cost(*params)` computing
        `sum(((data_y - model(data_x, *params)) / data_yerr) ** 2)`,
        ignoring NaNs, with a signature matching the selected model.
    """
    def costD0(b, beta, k_nl, a_nl, k_p, a_p, k_v0, a_v0, k_v1, a_v1):
        """Least-squares cost for the D0 model.

        Args:
            b, beta, k_nl, a_nl, k_p, a_p, k_v0, a_v0, k_v1, a_v1 (float):
                D0 model parameters (see `Pf_model` / `D0`).

        Returns:
            float: `sum(((data_y - model(data_x, ...)) / data_yerr) ** 2)`,
            ignoring NaNs.
        """
        ym = model(data_x, b, beta, k_nl, a_nl, k_p, a_p, k_v0, a_v0, k_v1, a_v1)
        z = (data_y - ym) / data_yerr
        return np.nansum(z**2)

    def costD1(b, beta, q_1, q_2, k_v, a_v, b_v, k_p):
        """Least-squares cost for the D1 model.

        Args:
            b, beta, q_1, q_2, k_v, a_v, b_v, k_p (float): D1 model
                parameters (see `Pf_model` / `D1`).

        Returns:
            float: `sum(((data_y - model(data_x, ...)) / data_yerr) ** 2)`,
            ignoring NaNs.
        """
        ym = model(data_x, b, beta, q_1, q_2, k_v, a_v, b_v, k_p)
        z = (data_y - ym) / data_yerr
        return np.nansum(z**2)

    def costD1BAO(b, beta, q_1, q_2, k_v, a_v, b_v, k_p, S_p, S_t):
        """Least-squares cost for the BAO-damped D1 model.

        Args:
            b, beta, q_1, q_2, k_v, a_v, b_v, k_p, S_p, S_t (float):
                BAO-damped D1 model parameters (see `Pf_model` / `D1`).

        Returns:
            float: `sum(((data_y - model(data_x, ...)) / data_yerr) ** 2)`,
            ignoring NaNs.
        """
        ym = model(data_x, b, beta, q_1, q_2, k_v, a_v, b_v, k_p, S_p, S_t)
        z = (data_y - ym) / data_yerr
        return np.nansum(z**2)

    def costlinear(b, beta):
        """Least-squares cost for the linear model.

        Args:
            b (float): Linear bias.
            beta (float): RSD parameter.

        Returns:
            float: `sum(((data_y - model(data_x, b, beta)) / data_yerr) ** 2)`,
            ignoring NaNs.
        """
        ym = model(data_x, b, beta)
        z = (data_y - ym) / data_yerr
        return np.nansum(z**2)

    if non_linear_model == "0":
        return costD0
    elif non_linear_model == "1":
        return costD1
    elif non_linear_model == "1_bao":
        return costD1BAO
    elif non_linear_model == None:
        return costlinear


def custom_least_squares_arinyo(model, data_x, data_y, data_yerr, non_linear_model="0"):
    """Build an Arinyo-style least-squares iminuit cost function for a model.

    Like `custom_least_squares`, but computes the residual as
    `((data_y**2 / ym) - ym) / data_yerr` instead of
    `(data_y - ym) / data_yerr`.

    Args:
        model (callable): Flux power spectrum model function, as
            returned by `Pf_model`.
        data_x (tuple): `(k, mu)` data coordinates.
        data_y (numpy.ndarray): Observed flux power spectrum values.
        data_yerr (numpy.ndarray): Uncertainties on `data_y`.
        non_linear_model (str or None): Selects which cost function
            signature to return: "0" (D0), "1" (D1), "1_bao" (D1 +
            BAO), or None (linear). Defaults to "0".

    Returns:
        callable: A cost function `cost(*params)` computing
        `sum((((data_y ** 2 / ym) - ym) / data_yerr) ** 2)` with
        `ym = model(data_x, *params)`, ignoring NaNs.
    """
    def costD0(b, beta, k_nl, a_nl, k_p, a_p, k_v0, a_v0, k_v1, a_v1):
        """Arinyo-style least-squares cost for the D0 model.

        Args:
            b, beta, k_nl, a_nl, k_p, a_p, k_v0, a_v0, k_v1, a_v1 (float):
                D0 model parameters (see `Pf_model` / `D0`).

        Returns:
            float: `sum((((data_y**2/ym) - ym) / data_yerr) ** 2)` with
            `ym = model(data_x, ...)`, ignoring NaNs.
        """
        ym = model(data_x, b, beta, k_nl, a_nl, k_p, a_p, k_v0, a_v0, k_v1, a_v1)
        z = ((data_y**2 / ym) - ym) / data_yerr
        return np.nansum(z**2)

    def costD1(b, beta, q_1, q_2, k_v, a_v, b_v, k_p):
        """Arinyo-style least-squares cost for the D1 model.

        Args:
            b, beta, q_1, q_2, k_v, a_v, b_v, k_p (float): D1 model
                parameters (see `Pf_model` / `D1`).

        Returns:
            float: `sum((((data_y**2/ym) - ym) / data_yerr) ** 2)` with
            `ym = model(data_x, ...)`, ignoring NaNs.
        """
        ym = model(data_x, b, beta, q_1, q_2, k_v, a_v, b_v, k_p)
        z = ((data_y**2 / ym) - ym) / data_yerr
        return np.nansum(z**2)

    def costD1BAO(b, beta, q_1, q_2, k_v, a_v, b_v, k_p, S_p, S_t):
        """Arinyo-style least-squares cost for the BAO-damped D1 model.

        Args:
            b, beta, q_1, q_2, k_v, a_v, b_v, k_p, S_p, S_t (float):
                BAO-damped D1 model parameters (see `Pf_model` / `D1`).

        Returns:
            float: `sum((((data_y**2/ym) - ym) / data_yerr) ** 2)` with
            `ym = model(data_x, ...)`, ignoring NaNs.
        """
        ym = model(data_x, b, beta, q_1, q_2, k_v, a_v, b_v, k_p, S_p, S_t)
        z = ((data_y**2 / ym) - ym) / data_yerr
        return np.nansum(z**2)

    def costlinear(b, beta):
        """Arinyo-style least-squares cost for the linear model.

        Args:
            b (float): Linear bias.
            beta (float): RSD parameter.

        Returns:
            float: `sum((((data_y**2/ym) - ym) / data_yerr) ** 2)` with
            `ym = model(data_x, b, beta)`, ignoring NaNs.
        """
        ym = model(data_x, b, beta)
        z = ((data_y**2 / ym) - ym) / data_yerr
        return np.nansum(z**2)

    if non_linear_model == "0":
        return costD0
    elif non_linear_model == "1":
        return costD1
    elif non_linear_model == "1_bao":
        return costD1BAO
    elif non_linear_model == None:
        return costlinear


def cost_function(model, data_x, data_y, data_yerr, cost_name, non_linear_model="0"):
    """Dispatch to the requested iminuit cost function builder.

    Args:
        model (callable): Flux power spectrum model function.
        data_x (tuple): `(k, mu)` data coordinates.
        data_y (numpy.ndarray): Observed flux power spectrum values.
        data_yerr (numpy.ndarray): Uncertainties on `data_y`.
        cost_name (str): Which cost function family to build: "least"
            (`custom_least_squares`) or "least_arinyo"
            (`custom_least_squares_arinyo`).
        non_linear_model (str or None): Non-linear model selector
            forwarded to the chosen cost-function builder. Defaults to
            "0".

    Returns:
        callable or None: The selected cost function, or None if
        `cost_name` is not recognized.
    """
    if cost_name == "least":
        return custom_least_squares(
            model, data_x, data_y, data_yerr, non_linear_model=non_linear_model
        )
    elif cost_name == "least_arinyo":
        return custom_least_squares_arinyo(
            model, data_x, data_y, data_yerr, non_linear_model=non_linear_model
        )


### Minuit stuff


def run_minuit(
    data_x,
    data_y,
    data_yerr,
    minuit_parameters,
    minuit_limits,
    power_l_rebin,
    non_linear_model="0",
    cost_name="least",
    ncall=100,
    fix_args=None,
    launch_minos=False,
    sigma_minos=None,
    var_minos=None,
    integrate_model=True,
    N_mu_integration=1000,
    mu_max=1.0,
    power_l_no_bao_rebin=None,
):
    """Build the flux power spectrum model and cost function, and run a Minuit fit.

    Constructs the model via `Pf_model`, builds the corresponding cost
    function via `cost_function`, initializes an `iminuit.Minuit`
    instance, applies parameter limits and fixed parameters, then runs
    MIGRAD (and optionally HESSE and MINOS).

    Args:
        data_x (tuple): `(k, mu)` data coordinates.
        data_y (numpy.ndarray): Observed flux power spectrum values.
        data_yerr (numpy.ndarray): Uncertainties on `data_y`.
        minuit_parameters (dict): Initial parameter values/settings
            passed to `Minuit(cost, **minuit_parameters)`.
        minuit_limits (list or None): Sequence of `(name, limits)`
            pairs applied via `minuit.limits[name] = limits`.
        power_l_rebin (numpy.ndarray): Linear power spectrum rebinned
            onto the flux power spectrum k-grid, passed to `Pf_model`.
        non_linear_model (str or None): Non-linear model selector ("0",
            "1", "1_bao", or None). Defaults to "0".
        cost_name (str): Cost function family ("least" or
            "least_arinyo"). Defaults to "least".
        ncall (int): Maximum number of function calls for MIGRAD (and
            MINOS). Defaults to 100.
        fix_args (list or None): Parameter names to fix (hold constant).
            Defaults to None.
        launch_minos (bool): Whether to run MINOS after HESSE. Defaults
            to False.
        sigma_minos (float or None): Sigma level passed to `run_minos`.
            Defaults to None.
        var_minos (str or list or None): Parameter(s) to run MINOS on.
            Defaults to None.
        integrate_model (bool): Whether `Pf_model` should integrate over
            mu bins. Defaults to True.
        N_mu_integration (int): Number of mu sub-samples per bin when
            integrating. Defaults to 1000.
        mu_max (float): Maximum mu used for bin-width computation.
            Defaults to 1.0.
        power_l_no_bao_rebin (numpy.ndarray or None): No-wiggle linear
            power spectrum rebinned onto the flux power spectrum
            k-grid, required for the BAO-damped model. Defaults to
            None.

    Returns:
        iminuit.Minuit: The fitted Minuit instance after MIGRAD (and
        optionally HESSE/MINOS).
    """
    model = Pf_model(
        power_l_rebin,
        non_linear_model=non_linear_model,
        linear_power_spectrum_no_bao=power_l_no_bao_rebin,
        integrate_model=integrate_model,
        N_mu_integration=N_mu_integration,
        mu_max=mu_max,
    )
    cost = cost_function(
        model, data_x, data_y, data_yerr, cost_name, non_linear_model=non_linear_model
    )
    minuit = Minuit(cost, **minuit_parameters)
    minuit.errordef = 1
    if minuit_limits is not None:
        for i in range(len(minuit_limits)):
            minuit.limits[minuit_limits[i][0]] = minuit_limits[i][1]
    if fix_args is not None:
        for i in range(len(fix_args)):
            minuit.fixed[fix_args[i]] = True
    print(run_migrad(minuit, ncall=ncall))
    run_hesse(minuit)
    if launch_minos:
        run_minos(minuit, sigma_minos, ncall=ncall, var_minos=var_minos)
    return minuit


def run_migrad(minuit, ncall=1000):
    """Run the iminuit MIGRAD minimizer.

    Args:
        minuit (iminuit.Minuit): The Minuit instance to minimize.
        ncall (int): Maximum number of function calls. Defaults to 1000.

    Returns:
        iminuit.Minuit: The Minuit instance after MIGRAD, as returned by
        `minuit.migrad`.
    """
    return minuit.migrad(ncall)


def run_minos(minuit, sigma, ncall=1000, var_minos=None):
    """Run the iminuit MINOS error analysis.

    Args:
        minuit (iminuit.Minuit): The Minuit instance to analyze.
        sigma (float): Confidence-level sigma passed to `minuit.minos`.
        ncall (int): Maximum number of function calls. Defaults to 1000.
        var_minos (str or list or None): Parameter(s) to run MINOS on;
            None runs it on all free parameters. Defaults to None.

    Returns:
        iminuit.Minuit: The Minuit instance after MINOS, as returned by
        `minuit.minos`.
    """
    return minuit.minos(var=var_minos, sigma=sigma, ncall=ncall)


def run_hesse(minuit):
    """Run the iminuit HESSE error analysis.

    Args:
        minuit (iminuit.Minuit): The Minuit instance to analyze.

    Returns:
        iminuit.Minuit: The Minuit instance after HESSE, as returned by
        `minuit.hesse`.
    """
    return minuit.hesse()


### Main


def prepare_data(
    pf_file,
    pk_file,
    power_weighted=False,
    class_dict=None,
    z_simu=None,
    z_init=None,
    kmax=None,
    kmin=None,
    name_pm_file=None,
    error_estimator=None,
    use_wavenumber_centers=True,
    use_mu_centers=False,
    **kwargs,
):
    """Load, center, cut, and rebin the flux and linear power spectra for fitting.

    Reads the 3D flux power spectrum from `pf_file`, optionally
    recenters its wavenumber/mu bin values, trims it to
    `[kmin, kmax]`, optionally rebins it, then builds the linear
    matter power spectrum (from CLASS, a normalized/raw simulation
    power spectrum, cosmoprimo, or a gimlet file depending on
    `pk_file`) and rebins it onto the flux power spectrum k-grid.

    Args:
        pf_file (str): Path to the gimlet 3D flux power spectrum file.
        pk_file (str): Source of the linear matter power spectrum:
            "class", "pmnorm", "pm", "cosmoprimo", or a path to a
            gimlet matter power spectrum file.
        power_weighted (bool): Whether the flux power spectrum bins are
            already power-weighted. Defaults to False.
        class_dict (dict or None): CLASS/cosmoprimo settings
            dictionary, required when `pk_file` is "class", "pmnorm",
            or "cosmoprimo". Defaults to None.
        z_simu (float or None): Simulation redshift, required when
            `pk_file` is "class", "pmnorm", or "cosmoprimo". Defaults
            to None.
        z_init (float or None): Initial-conditions redshift, required
            when `pk_file` is "pmnorm". Defaults to None.
        kmax (float or None): Maximum wavenumber kept. Defaults to
            None.
        kmin (float or None): Minimum wavenumber kept. Defaults to
            None.
        name_pm_file (str or None): Path to the raw simulation matter
            power spectrum file, required when `pk_file` is "pmnorm"
            or "pm". Defaults to None.
        error_estimator (str or None): Error estimator forwarded to
            `read_pfkmu`. Defaults to None.
        use_wavenumber_centers (bool): Whether to recenter wavenumber
            bins via `power_f.center_wavenumbers_2d`. Defaults to True.
        use_mu_centers (bool): Whether to recenter mu bins via
            `power_f.center_mu_2d`. Defaults to False.
        **kwargs: Additional keyword arguments; may include a "rebin"
            dict with keys "nb_bin", "loglin", "k_loglin" controlling
            `power_f.rebin_2d_arrays`, and are otherwise forwarded to
            `read_pfkmu`.

    Returns:
        tuple: `(power_f, power_l, power_l_rebin, data_x, data_y,
        data_yerr, power_l_no_bao, power_l_no_bao_rebin)` — the flux
        power spectrum object, the linear matter power spectrum
        object, the linear power spectrum rebinned onto the flux
        k-grid, the fit x/y/yerr data arrays, and (when applicable) the
        no-wiggle linear power spectrum object and its rebinned
        counterpart (both None otherwise).

    Raises:
        KeyError: If `power_f.error_array` is None, i.e. no valid
            `error_estimator` was supplied.
    """
    power_f = read_pfkmu(
        pf_file,
        power_weighted=power_weighted,
        error_estimator=error_estimator,
        **kwargs,
    )
    if use_wavenumber_centers:
        if power_weighted:
            print("Wavenumbers already power weighted, skipping centering")
        else:
            power_f.center_wavenumbers_2d()

    if use_mu_centers:
        if power_weighted:
            print("Mu already power weighted, skipping centering")
        else:
            print(
                "Mu centering not properly tested. I suggest to use integration over mu"
            )
            power_f.center_mu_2d()

    power_f.cut_extremum(kmin, kmax)

    rebin = utils_fitter.return_key(kwargs, "rebin", None)
    if rebin is not None:
        power_f.rebin_2d_arrays(
            rebin["nb_bin"],
            operation="mean",
            loglin=rebin["loglin"],
            k_loglin=rebin["k_loglin"],
        )
    power_l_no_bao = None
    if pk_file == "class":
        power_l = Pl_class(power_f.k_array[0], class_dict, z_simu, name="class")
    elif pk_file == "pmnorm":
        power_l = Pm_normalized(name_pm_file, class_dict, z_simu, z_init, name="pmnorm")
    elif pk_file == "pm":
        power_l = Pm(name_pm_file, name="pm")
    elif pk_file == "cosmoprimo":
        power_l, power_l_no_bao = Pl_cosmoprimo(power_f.k_array[0], class_dict, z_simu)
    else:
        power_l = read_pk(pk_file, power_weighted=power_weighted)
    if use_wavenumber_centers:
        power_l.edge_stored = True

    power_l_rebin = rebin_matter_power(
        power_l.power_array, power_l.k_array, power_f.k_array[0]
    )
    power_l_no_bao_rebin = None
    if power_l_no_bao is not None:
        power_l_no_bao_rebin = rebin_matter_power(
            power_l_no_bao.power_array, power_l_no_bao.k_array, power_f.k_array[0]
        )

    data_x = power_f.k_array
    data_y = power_f.power_array
    if power_f.error_array is None:
        raise KeyError("Choose an error_estimator")
    data_yerr = power_f.error_array

    return (
        power_f,
        power_l,
        power_l_rebin,
        data_x,
        data_y,
        data_yerr,
        power_l_no_bao,
        power_l_no_bao_rebin,
    )


def fitter_k_mu(
    pf_file,
    pk_file,
    minuit_parameters,
    minuit_limits,
    power_weighted=False,
    class_dict=None,
    z_simu=None,
    z_init=None,
    non_linear_model="0",
    cost_name="least",
    ncall=100,
    kmax=None,
    kmin=None,
    launch_minos=None,
    var_minos=None,
    sigma_minos=None,
    name_pm_file=None,
    error_estimator=None,
    fix_args=None,
    integrate_model=True,
    N_mu_integration=1000,
    mu_max=1.0,
    use_wavenumber_centers=True,
    use_mu_centers=False,
    **kwargs,
):
    """End-to-end flux power spectrum fit in (k, mu): load data and run Minuit.

    Validates the combination of `power_weighted`,
    `use_wavenumber_centers`, `use_mu_centers` and `integrate_model`
    options (printing warnings for inconsistent/untested combinations),
    calls `prepare_data` to load and rebin the flux and linear power
    spectra, then calls `run_minuit` to perform the fit.

    Args:
        pf_file (str): Path to the gimlet 3D flux power spectrum file.
        pk_file (str): Source of the linear matter power spectrum (see
            `prepare_data`).
        minuit_parameters (dict): Initial Minuit parameter
            values/settings.
        minuit_limits (list or None): Sequence of `(name, limits)`
            pairs for Minuit parameter limits.
        power_weighted (bool): Whether the flux power spectrum bins are
            already power-weighted. Defaults to False.
        class_dict (dict or None): CLASS/cosmoprimo settings
            dictionary. Defaults to None.
        z_simu (float or None): Simulation redshift. Defaults to None.
        z_init (float or None): Initial-conditions redshift. Defaults
            to None.
        non_linear_model (str or None): Non-linear model selector ("0",
            "1", "1_bao", or None). Defaults to "0".
        cost_name (str): Cost function family ("least" or
            "least_arinyo"). Defaults to "least".
        ncall (int): Maximum number of MIGRAD/MINOS function calls.
            Defaults to 100.
        kmax (float or None): Maximum wavenumber kept. Defaults to
            None.
        kmin (float or None): Minimum wavenumber kept. Defaults to
            None.
        launch_minos (bool or None): Whether to run MINOS after HESSE.
            Defaults to None.
        var_minos (str or list or None): Parameter(s) to run MINOS on.
            Defaults to None.
        sigma_minos (float or None): Sigma level passed to MINOS.
            Defaults to None.
        name_pm_file (str or None): Path to the raw simulation matter
            power spectrum file. Defaults to None.
        error_estimator (str or None): Error estimator forwarded to
            `prepare_data`. Defaults to None.
        fix_args (list or None): Parameter names to fix in Minuit.
            Defaults to None.
        integrate_model (bool): Whether `Pf_model` should integrate over
            mu bins. Defaults to True.
        N_mu_integration (int): Number of mu sub-samples per bin when
            integrating. Defaults to 1000.
        mu_max (float): Maximum mu used for bin-width computation.
            Defaults to 1.0.
        use_wavenumber_centers (bool): Whether to recenter wavenumber
            bins. Defaults to True.
        use_mu_centers (bool): Whether to recenter mu bins. Defaults to
            False.
        **kwargs: Additional keyword arguments forwarded to
            `prepare_data`.

    Returns:
        tuple: `(minuit, power_f, power_l, power_l_rebin,
        non_linear_model, power_l_no_bao, power_l_no_bao_rebin)` — the
        fitted Minuit instance, the flux and linear power spectrum
        objects, the rebinned linear power spectrum, the non-linear
        model name, and (when applicable) the no-wiggle linear power
        spectrum and its rebinned counterpart.
    """
    if power_weighted is False:
        if use_wavenumber_centers is False:
            print(
                "You are computing the model at the wavenumber edges, it will lead to wrong values, "
                "please choose the option use_wavenumber_centers, "
                "or use the power_weighted option"
            )
        if (use_mu_centers is False) & (integrate_model is False):
            print(
                "You are computing the model at the mu edges, it will lead to wrong values, "
                "please choose between the options integrate_model or use_mu_centers, "
                "or use the power_weighted option"
            )
    else:
        print("Power weigthed fit not tested.")

    if use_mu_centers & integrate_model:
        print(
            "You are integrating the model between centered mu values, "
            "please choose between the options integrate_model or use_mu_centers"
        )
    (
        power_f,
        power_l,
        power_l_rebin,
        data_x,
        data_y,
        data_yerr,
        power_l_no_bao,
        power_l_no_bao_rebin,
    ) = prepare_data(
        pf_file,
        pk_file,
        power_weighted=power_weighted,
        class_dict=class_dict,
        z_simu=z_simu,
        z_init=z_init,
        kmax=kmax,
        kmin=kmin,
        name_pm_file=name_pm_file,
        error_estimator=error_estimator,
        use_wavenumber_centers=use_wavenumber_centers,
        use_mu_centers=use_mu_centers,
        **kwargs,
    )
    minuit = run_minuit(
        data_x,
        data_y,
        data_yerr,
        minuit_parameters,
        minuit_limits,
        power_l_rebin,
        non_linear_model=non_linear_model,
        cost_name=cost_name,
        ncall=ncall,
        fix_args=fix_args,
        launch_minos=launch_minos,
        var_minos=var_minos,
        sigma_minos=sigma_minos,
        integrate_model=integrate_model,
        N_mu_integration=N_mu_integration,
        mu_max=mu_max,
        power_l_no_bao_rebin=power_l_no_bao_rebin,
    )
    return (
        minuit,
        power_f,
        power_l,
        power_l_rebin,
        non_linear_model,
        power_l_no_bao,
        power_l_no_bao_rebin,
    )


def compute_kna(minuit, power_l, eps, nloopmax=1000):
    """Iteratively solve for the non-linear wavenumber k_na from a fit.

    Solves the fixed-point equation
    `k_na = ((2*pi)^2 * k_v^a_v * ln(1+beta) / (q_1 * P_lin(k_na)))^(1/(3+a_v))`
    by simple fixed-point iteration starting from `k_na = 3`.

    Args:
        minuit (iminuit.Minuit): Fitted Minuit instance holding values
            for "a_v", "k_v", "beta", and "q_1".
        power_l (power_spectra.MatterPowerSpectrum): Linear matter
            power spectrum used to evaluate P_lin(k_na) via
            interpolation.
        eps (float): Convergence tolerance on successive iterates of
            k_na.
        nloopmax (int): Maximum number of iterations. Defaults to 1000.

    Returns:
        float: The converged value of k_na.

    Raises:
        Warning: If the iteration does not converge within `nloopmax`
            steps.
    """
    values = dict(minuit.values.items())
    av, kv, beta, q1 = values["a_v"], values["k_v"], values["beta"], values["q_1"]
    power_l_interp = scipy.interpolate.interp1d(
        power_l.k_array, power_l.power_array, bounds_error=False, fill_value=np.nan
    )
    kna0 = 3
    knai = kna0
    diff = np.inf
    n = 0
    while (diff > eps) & (n < nloopmax):
        knai2 = (
            ((2 * np.pi) ** 2 * kv**av * np.log(1 + beta)) / (q1 * power_l_interp(knai))
        ) ** (1 / (3 + av))
        diff = abs(knai2 - knai)
        knai = knai2
        n += 1
    if n == nloopmax:
        raise Warning("Maximal loop iteration reached")
    return knai


### Plots


def plot_pl(power_l):
    """Open a new plot and draw the 1D linear matter power spectrum.

    Args:
        power_l (power_spectra.MatterPowerSpectrum): Linear matter power
            spectrum to plot.

    Returns:
        None
    """
    power_l.open_plot()
    power_l.plot_1d_pk()


def plot_pf(power_f, mu_bin, legend):
    """Open a new plot and draw the 2D flux power spectrum for given mu bins.

    Args:
        power_f (power_spectra.FluxPowerSpectrum): Flux power spectrum
            to plot.
        mu_bin (list): Mu bin indices/values to plot.
        legend (bool or list): Legend labels/flag forwarded to
            `power_f.plot_2d_pk`.

    Returns:
        None
    """
    power_f.open_plot()
    power_f.plot_2d_pk(mu_bin, legend=legend)


def plot_pf_pm(power_f, power_m, mu_bin, legend):
    """Plot the flux-to-matter power spectrum ratio for given mu bins.

    Rebins `power_m` onto the flux power spectrum k-grid, divides
    `power_f`'s power and error arrays by the rebinned matter power
    spectrum (in-place), and plots the resulting ratio.

    Args:
        power_f (power_spectra.FluxPowerSpectrum): Flux power spectrum;
            its `power_array`/`error_array` are overwritten with the
            ratio.
        power_m (power_spectra.MatterPowerSpectrum): Matter power
            spectrum used as the denominator.
        mu_bin (list): Mu bin indices/values to plot.
        legend (bool or list): Legend labels/flag forwarded to
            `power_f.plot_2d_pk`.

    Returns:
        None
    """
    power_f.open_plot()
    power_m_rebin = rebin_matter_power(
        power_m.power_array, power_m.k_array, power_f.k_array[0]
    )
    power_f.power_array = power_f.power_array / power_m_rebin
    power_f.error_array = power_f.error_array / power_m_rebin
    power_f.plot_2d_pk(mu_bin, legend=legend, ps="x")


def plot_fit(
    minuit,
    power_f,
    power_l_rebin,
    non_linear_model,
    mu_bin,
    legend,
    name_out="fit_results",
    integrate_model=True,
    N_mu_integration=1000,
    power_l_no_bao_rebin=None,
    plot_no_bao_ratio=False,
    mu_max=1.0,
    **kwargs,
):
    """Plot fitted flux/matter power spectrum ratio against the data and save results.

    Rebuilds the model via `Pf_model` using the best-fit Minuit
    parameters, computes the data-to-linear and model-to-linear power
    spectrum ratios, plots both (data as points, model as lines) per mu
    bin, saves the figure as PDF and PNG, and writes the fit parameters
    to a text file.

    Args:
        minuit (iminuit.Minuit): Fitted Minuit instance providing
            best-fit parameter values.
        power_f (power_spectra.FluxPowerSpectrum): Flux power spectrum
            data.
        power_l_rebin (numpy.ndarray): Linear power spectrum rebinned
            onto the flux power spectrum k-grid.
        non_linear_model (str or None): Non-linear model selector used
            to rebuild the model via `Pf_model`.
        mu_bin (list): Mu bin indices/values to plot.
        legend (bool or list): Legend labels/flag forwarded to
            `power2.plot_2d_pk`/`power1.plot_2d_pk`.
        name_out (str): Base filename (without extension) for the saved
            plot and parameter files. Defaults to "fit_results".
        integrate_model (bool): Whether `Pf_model` should integrate over
            mu bins. Defaults to True.
        N_mu_integration (int): Number of mu sub-samples per bin when
            integrating. Defaults to 1000.
        power_l_no_bao_rebin (numpy.ndarray or None): No-wiggle linear
            power spectrum rebinned onto the flux k-grid, used for the
            BAO-damped model and/or as the ratio denominator when
            `plot_no_bao_ratio` is True. Defaults to None.
        plot_no_bao_ratio (bool): If True, use `power_l_no_bao_rebin` as
            the denominator for the plotted ratios instead of
            `power_l_rebin`. Defaults to False.
        mu_max (float): Maximum mu used for bin-width computation.
            Defaults to 1.0.
        **kwargs: Additional keyword arguments forwarded to the
            plotting calls (`open_plot`, `plot_2d_pk`).

    Returns:
        None
    """
    model = Pf_model(
        power_l_rebin,
        non_linear_model=non_linear_model,
        linear_power_spectrum_no_bao=power_l_no_bao_rebin,
        integrate_model=integrate_model,
        N_mu_integration=N_mu_integration,
        mu_max=mu_max,
    )
    minuit_params = []
    for i in range(len(minuit.parameters)):
        minuit_params.append(minuit.params[minuit.parameters[i]].value)

    if plot_no_bao_ratio:
        power_l_plot = power_l_no_bao_rebin
    else:
        power_l_plot = power_l_rebin

    pf_over_pm_data = power_f.power_array / power_l_plot
    if power_f.error_array is not None:
        error_array = power_f.error_array / power_l_plot
    else:
        error_array = None
    pf_over_pm_model = model(power_f.k_array, *minuit_params) / power_l_plot

    color = [f"C{i}" for i in range(len(mu_bin))]

    h_normalized = power_f.h_normalized
    power1 = power_spectra.FluxPowerSpectrum(
        k_array=power_f.k_array,
        power_array=pf_over_pm_model,
        dimension="3D",
        h_normalized=h_normalized,
    )
    power1.open_plot(**kwargs)

    power2 = power_spectra.FluxPowerSpectrum(
        k_array=power_f.k_array,
        power_array=pf_over_pm_data,
        error_array=error_array,
        dimension="3D",
        h_normalized=h_normalized,
    )
    power2.plot_2d_pk(
        mu_bin,
        color=color,
        ps="x",
        linestyle=["None" for i in range(len(mu_bin))],
        **kwargs,
    )

    power1.plot_2d_pk(mu_bin, color=color, legend=legend, **kwargs)

    power1.save_plot(f"{name_out}.pdf")
    power1.save_plot(f"{name_out}.png", format_out="png")
    power1.close_plot()

    with open(f"{name_out}_param.txt", "w") as f:
        print(minuit.params, file=f)


### Latex


def minuit_to_latex(minuit, name=""):
    """obsolete"""
    try:
        file = open(f"{name}_minuit_matrix.tex", "w")
        file.write(minuit.latex_matrix().__str__())
        file.close()
    except:
        print("no minuit matrix")

    file = open(f"{name}_minuit_params.tex", "w")
    file.write(minuit.latex_param().__str__())
    file.close()
