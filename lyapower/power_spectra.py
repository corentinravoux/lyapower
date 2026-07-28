#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 16:29:56 2019

@author: cravoux

Power-spectrum data objects and I/O for the ``lyapower`` package.

This module defines the container classes used to represent 1D and 3D
(or k, mu binned) power spectra measured from Nyx hydrodynamical
simulations post-processed by ``gimlet``: :class:`PowerSpectrum` (base
class), :class:`MatterPowerSpectrum` and :class:`FluxPowerSpectrum`.
It provides readers for gimlet / genpk / ascii power-spectrum files,
rebinning utilities (1D and 2D in k or k, mu), plotting helpers, and
the "splicing" routines (:func:`splice_1D`, :func:`splice_3D`) that
combine multiple-resolution simulation boxes into a single spliced
power spectrum, following Arinyo-i-Prats et al. 2015.

Conventions: wavenumber ``k`` is expressed in h/Mpc or 1/Mpc (see
:meth:`PowerSpectrum.change_k_normalization`), ``mu`` is
k_parallel / k, and power spectra are stored as P(k) (1D) or
P(k, mu) (3D, with ``k_array`` stacked as ``[k, mu]``).
"""

import math

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import cm
from scipy.interpolate import interp1d

from lyapower import utils_fitter as utils


class PowerSpectrum(object):
    """Base container for a measured power spectrum.

    Holds the wavenumber array (1D ``k`` or stacked ``[k, mu]`` for
    2D spectra), the power values, an optional error array, and
    bookkeeping metadata (source file, simulation box size, and
    whether ``k`` is h-normalized). Provides file readers
    (:meth:`init_from_genpk_file`, :meth:`init_from_ascii_file`),
    rebinning, plotting and unit-conversion utilities shared by
    :class:`MatterPowerSpectrum` and :class:`FluxPowerSpectrum`.
    """

    def __init__(
        self,
        k_array=None,
        power_array=None,
        error_array=None,
        file_init=None,
        size_box=None,
        h_normalized=None,
    ):
        """Initialize a PowerSpectrum.

        Args:
            k_array: Wavenumber array (1D) or stacked ``[k, mu]`` array (2D).
            power_array: Power spectrum values, same length as ``k_array``.
            error_array: Optional array of uncertainties on ``power_array``.
            file_init: Optional path of the file the spectrum was read from.
            size_box: Optional simulation box size (Mpc/h).
            h_normalized: Whether ``k_array`` is expressed in h/Mpc (True)
                or 1/Mpc (False).
        """
        self.k_array = k_array
        self.power_array = power_array
        self.file_init = file_init
        self.size_box = size_box
        self.h_normalized = h_normalized
        self.error_array = error_array
        self.edge_stored = True

    @classmethod
    def init_from_genpk_file(cls, name_file, size_box):
        """Load a GenPk format power spectum, plotting the DM and the neutrinos (if present)
        Does not plot baryons."""
        # Load DM P(k)
        matpow = np.loadtxt(name_file)
        if size_box is None:
            raise KeyError(
                "To load a GenPk output, the size of the box in Mpc.h-1 must be given"
            )
        scale = 2 * math.pi / size_box
        # Adjust Fourier convention to match CAMB.
        simk = matpow[1:, 0] * scale
        Pk = matpow[1:, 1] / scale**3 * (2 * math.pi) ** 3
        h_normalized = True
        return cls(
            k_array=simk,
            power_array=Pk,
            file_init=name_file,
            size_box=size_box,
            h_normalized=h_normalized,
        )

    @classmethod
    def init_from_ascii_file(cls, name_file):
        """Load a power spectrum from a plain two-column ascii file.

        The first row of the file is skipped (treated as a header/edge
        row); column 0 is used as ``k`` and column 1 as the power.

        Args:
            name_file: Path to the ascii file with columns [k, power].

        Returns:
            PowerSpectrum: New instance with ``h_normalized=True`` and
            ``size_box=None``.
        """
        file_power = np.loadtxt(name_file)
        k_array = file_power[1:, 0]
        power_array = file_power[1:, 1]
        h_normalized = True
        return cls(
            k_array=k_array,
            power_array=power_array,
            file_init=name_file,
            size_box=None,
            h_normalized=h_normalized,
        )

    def center_wavenumbers_1d(self):
        """Replace 1D bin-edge wavenumbers with bin-center wavenumbers.

        Recomputes ``self.k_array`` in place as the midpoints between
        consecutive stored values, extrapolating the last bin center
        from the final spacing. Sets ``self.edge_stored = False``.
        """
        k_centers = (self.k_array[1:] + self.k_array[:-1]) / 2
        last_value = (3 * self.k_array[-1] - self.k_array[-2]) / 2
        self.k_array = np.concatenate([k_centers, [last_value]])
        self.edge_stored = False

    def center_wavenumbers_2d(self):
        """Replace 2D bin-edge k values with bin-center k values, per mu bin.

        For each unique ``mu`` value in ``self.k_array[1]``, converts the
        corresponding ``k`` edges (``self.k_array[0]``) to bin centers
        in place, extrapolating the last center from the final spacing.
        Sets ``self.edge_stored = False``.
        """
        mus = np.unique(self.k_array[1])
        for mu in mus:
            mask = self.k_array[1] == mu
            k_edge = self.k_array[0][mask]
            k_centers = (k_edge[1:] + k_edge[:-1]) / 2
            last_value = (3 * k_edge[-1] - k_edge[-2]) / 2
            self.k_array[0][mask] = np.concatenate([k_centers, [last_value]])
        self.edge_stored = False

    @staticmethod
    def compute_dmu(mu, mu_max=1.0):
        """Compute the bin width in mu for each entry of a sorted mu array.

        Args:
            mu: 1D array of mu bin-edge values (assumed sorted per group).
            mu_max: Upper bound of the mu range, used to close the last bin
                and to fix up any negative widths caused by mu wrap-around
                (e.g. between successive k groups). Defaults to 1.0.

        Returns:
            numpy.ndarray: Array of mu bin widths, same length as ``mu``.
        """
        dmu = mu[1:] - mu[:-1]
        dmu = np.concatenate([dmu, [mu_max - mu[-1]]], axis=0)
        mask = dmu < 0
        dmu[mask] = mu_max - mu[mask]
        return dmu

    def center_mu_2d(self, mu_max=1.0):
        """Shift stored mu edges to mu bin centers, in place.

        Args:
            mu_max: Upper bound of the mu range, forwarded to
                :meth:`compute_dmu`. Defaults to 1.0.
        """
        dmu = PowerSpectrum.compute_dmu(self.k_array[1], mu_max=mu_max)
        self.k_array[1] = self.k_array[1] + dmu / 2
        self.edge_stored = False

    def rebin_arrays(self, nb_bin, operation="mean"):
        """Rebin a 1D power spectrum onto a new log-spaced k grid, in place.

        Builds ``nb_bin`` log-spaced k values spanning the current
        ``k_array`` range and, for each resulting interval, aggregates
        the power values that fall in it (falling back to the nearest
        original point when an interval is empty). Updates
        ``self.k_array`` and ``self.power_array`` in place.

        Args:
            nb_bin: Number of new k bins (i.e. length of the new k grid).
            operation: Aggregation to apply within each bin: ``"mean"``
                (also ``"average"``/``"avg"``) for a plain mean, or
                ``"gauss"`` for a Gaussian-weighted average. Defaults to
                ``"mean"``.
        """
        new_k = np.logspace(
            np.log10(min(self.k_array)), np.log10(max(self.k_array)), nb_bin
        )
        new_Pk = np.zeros(new_k.shape)
        for i in range(nb_bin - 1):
            mask = (self.k_array >= new_k[i]) & (self.k_array < new_k[i + 1])
            if operation.lower() in ["mean", "average", "avg"]:
                if len(self.power_array[mask]) == 0:
                    nearest_index = np.argmin(np.abs(self.k_array - new_k[i]))
                    new_Pk[i] = self.power_array[nearest_index]
                else:
                    new_Pk[i] = np.mean(self.power_array[mask])
            elif operation.lower() in ["gauss"]:
                from scipy import signal

                gaussian_weights = signal.gaussian(
                    int(len(self.power_array[mask])),
                    int(len(self.power_array[mask])) / 4,
                )
                if len(self.power_array[mask]) == 0:
                    nearest_index = np.argmin(np.abs(self.k_array - new_k[i]))
                    new_Pk[i] = self.power_array[nearest_index]
                else:
                    new_Pk[i] = np.average(
                        self.power_array[mask], axis=0, weights=gaussian_weights
                    )
        self.k_array = new_k
        self.power_array = new_Pk

    def rebin_2d_arrays(self, nb_bin, operation="mean", loglin=False, k_loglin=None):
        """Rebin a 2D (k, mu) power spectrum onto a new log-spaced k grid.

        For each mu value, rebins the k axis onto ``nb_bin`` log-spaced
        bins spanning the (optionally restricted) k range, aggregating
        power (and error, if present) within each bin. Updates
        ``self.k_array``, ``self.power_array`` and ``self.error_array``
        in place.

        Args:
            nb_bin: Number of new k bins along each mu slice.
            operation: Aggregation to apply within each bin: ``"mean"``
                (also ``"average"``/``"avg"``) or ``"gauss"`` for a
                Gaussian-weighted average. Defaults to ``"mean"``.
            loglin: If True, only rebin the k range above ``k_loglin``
                and leave lower-k points untouched (log-linear hybrid
                binning). Defaults to False.
            k_loglin: k threshold above which rebinning is applied when
                ``loglin`` is True. Required if ``loglin`` is True.
        """
        if loglin:
            mask_k_rebin = self.k_array[0] > k_loglin
            mu_old = self.k_array[1][mask_k_rebin]
            k_old = self.k_array[0][mask_k_rebin]
            power_old = self.power_array[mask_k_rebin]
            if self.error_array is not None:
                error_old = self.error_array[mask_k_rebin]
        else:
            mu_old = self.k_array[1]
            k_old = self.k_array[0]
            power_old = self.power_array
            if self.error_array is not None:
                error_old = self.error_array
        new_k = np.logspace(np.log10(min(k_old)), np.log10(max(k_old)), nb_bin)
        bin_centers = np.array(
            [0.5 * (new_k[i] + new_k[i + 1]) for i in range(len(new_k) - 1)]
        )
        mus = np.unique(mu_old)
        new_Pk = np.zeros(len(mus) * len(bin_centers))
        if self.error_array is not None:
            new_error_array = np.zeros(len(mus) * len(bin_centers))
        for j in range(len(mus)):
            for i in range(nb_bin - 1):
                mask = k_old >= new_k[i]
                mask &= k_old < new_k[i + 1]
                mask &= mu_old == mus[j]
                if operation.lower() in ["mean", "average", "avg"]:
                    if len(power_old[mask]) == 0:
                        nearest_index = np.argmin(np.abs(k_old[0] - new_k[i]))
                        new_Pk[i * len(mus) + j] = power_old[nearest_index]
                        if self.error_array is not None:
                            new_error_array[i * len(mus) + j] = error_old[nearest_index]
                    else:
                        new_Pk[i * len(mus) + j] = np.mean(power_old[mask])
                        if self.error_array is not None:
                            new_error_array[i * len(mus) + j] = np.mean(error_old[mask])
                elif operation.lower() in ["gauss"]:
                    from scipy import signal

                    gaussian_weights = signal.gaussian(
                        int(len(power_old[mask])), int(len(power_old[mask])) / 4
                    )
                    if len(power_old[mask]) == 0:
                        nearest_index = np.argmin(np.abs(k_old - new_k[i]))
                        new_Pk[i * len(mus) + j] = power_old[nearest_index]
                        if self.error_array is not None:
                            new_error_array[i * len(mus) + j] = error_old[nearest_index]
                    else:
                        new_Pk[i * len(mus) + j] = np.average(
                            power_old[mask], axis=0, weights=gaussian_weights
                        )
                        if self.error_array is not None:
                            new_error_array[i * len(mus) + j] = np.average(
                                error_old[mask], axis=0, weights=gaussian_weights
                            )

        new_2d_k = np.transpose(
            [
                [bin_centers[i], mus[j]]
                for i in range(len(bin_centers))
                for j in range(len(mus))
            ]
        )
        if loglin:
            self.k_array = np.array(
                [
                    np.concatenate([self.k_array[0][~mask_k_rebin], new_2d_k[0]]),
                    np.concatenate([self.k_array[1][~mask_k_rebin], new_2d_k[1]]),
                ]
            )
            self.power_array = np.concatenate([self.power_array[~mask_k_rebin], new_Pk])
            if self.error_array is not None:
                self.error_array = np.concatenate(
                    [self.error_array[~mask_k_rebin], new_error_array]
                )
        else:
            self.k_array = new_2d_k
            self.power_array = new_Pk
            if self.error_array is not None:
                self.error_array = new_error_array

    def cut_extremum(self, kmin, kmax):
        """Restrict a 2D power spectrum to a k range, in place.

        Filters ``self.k_array``, ``self.power_array`` and (if present)
        ``self.error_array`` to keep only entries whose ``k_array[0]``
        lies within ``[kmin, kmax]``.

        Args:
            kmin: Lower k bound (inclusive), or None to skip the lower cut.
            kmax: Upper k bound (inclusive), or None to skip the upper cut.
        """
        mask = np.full(self.power_array.shape, True)
        if kmin is not None:
            mask &= self.k_array[0, :] >= kmin
        if kmax is not None:
            mask &= self.k_array[0, :] <= kmax
        self.k_array = np.transpose(np.transpose(self.k_array)[mask])
        self.power_array = self.power_array[mask]
        if self.error_array is not None:
            self.error_array = self.error_array[mask]

    def put_label(
        self,
        ax,
        xunit=True,
        yunit=True,
        y_label=r"$P$",
        x_label=r"$k$",
        labelsize_x=12,
        labelsize_y=12,
        fontsize=12,
    ):
        """Set axis labels (with units) on a matplotlib Axes.

        Args:
            ax: Matplotlib Axes to label.
            xunit: If True, append the k unit to the x label. Defaults to True.
            yunit: If True, append the power unit to the y label. Defaults to True.
            y_label: Base y-axis label (unit suffix appended if ``yunit``).
                Defaults to ``r"$P$"``.
            x_label: Base x-axis label (unit suffix appended if ``xunit``).
                Defaults to ``r"$k$"``.
            labelsize_x: Tick label font size for the x axis. Defaults to 12.
            labelsize_y: Tick label font size for the y axis. Defaults to 12.
            fontsize: Font size for the axis labels. Defaults to 12.
        """
        ylab, xlab = "", ""
        if yunit:
            if self.h_normalized:
                ylab = r" $[h^{-3}$" + r"$\cdot$" + "$\mathrm{Mpc}^3]$"
            else:
                ylab = r" $[\mathrm{Mpc}^3$]"
        if xunit:
            if self.h_normalized:
                xlab = r" $[h$" + r"$\cdot$" + "$\mathrm{Mpc}^{-1}]$"
            else:
                xlab = r" $[\mathrm{Mpc}^{-1}]$"
        ax.tick_params(axis="y", labelsize=labelsize_y)
        ax.set_ylabel(f"{y_label}{ylab}", fontsize=fontsize)
        ax.tick_params(axis="x", labelsize=labelsize_x)
        ax.set_xlabel(f"{x_label}{xlab}", fontsize=fontsize)

    def prepare_axes(self, kwargs):
        """Resolve the main and comparison Axes to plot on from kwargs.

        Reads the ``"ax"`` key (list of Axes) from ``kwargs``, defaulting
        to the current figure's axes. If none are found, uses
        ``plt.gca()`` for both; if one is found, it is used for both the
        main plot and the comparison (ratio) panel; if two or more are
        found, the first is the main plot and the second the comparison
        panel.

        Args:
            kwargs: Keyword-argument dict, inspected via
                ``utils.return_key(kwargs, "ax", ...)``.

        Returns:
            tuple: ``(ax_to_plot, ax_comparison)`` matplotlib Axes.
        """
        ax = utils.return_key(kwargs, "ax", plt.gcf().get_axes())
        if len(ax) == 0:
            ax_to_plot = plt.gca()
            ax_comparison = plt.gca()
        elif len(ax) == 1:
            ax_to_plot = ax[0]
            ax_comparison = ax[0]
        else:
            ax_to_plot = ax[0]
            ax_comparison = ax[1]
        return (ax_to_plot, ax_comparison)

    def plot_1d_pk(self, **kwargs):
        """Plot this 1D power spectrum, with an optional comparison ratio panel.

        Plots ``power_array`` vs ``k_array`` on the main axes and,
        if a ``"comparison"`` PowerSpectrum is given, plots the relative
        difference ``(comparison - self) / comparison`` (interpolated onto
        the comparison's k grid) on the comparison axes.

        Args:
            **kwargs: Options read via ``utils.return_key``:
                ``comparison`` (PowerSpectrum, default None): spectrum to
                    compare against.
                ``ax`` (list of Axes): passed to :meth:`prepare_axes`.
                ``color`` (default None): line/marker color.
                ``ps`` (default None): marker style.
                ``ls`` (default ``"-"``): line style.
                ``xscale``/``yscale`` (default ``"log"``): axis scales.
                ``x_min_lim``/``x_max_lim``/``y_min_lim``/``y_max_lim``:
                    main-axes limits (default None).
                ``x_min_lim_comparison`` etc.: comparison-axes limits.
                ``legend`` (default []): legend labels.
                ``legend_elements`` (default None): custom legend handles.
        """
        comparison = utils.return_key(kwargs, "comparison", None)
        (ax_to_plot, ax_comparison) = self.prepare_axes(kwargs)
        self.put_label(ax_to_plot)

        color = utils.return_key(kwargs, "color", None)

        if comparison is not None:
            power_array_comparison = interp1d(
                self.k_array, self.power_array, bounds_error=False, fill_value=np.NaN
            )(comparison.k_array)
            ax_comparison.plot(
                comparison.k_array,
                (comparison.power_array - power_array_comparison)
                / comparison.power_array,
                color=color,
            )
        ax_to_plot.plot(
            self.k_array,
            self.power_array,
            marker=utils.return_key(kwargs, "ps", None),
            linestyle=utils.return_key(kwargs, "ls", "-"),
            color=color,
        )

        xscale = utils.return_key(kwargs, "xscale", "log")
        yscale = utils.return_key(kwargs, "yscale", "log")

        ax_to_plot.set_xscale(xscale)
        ax_to_plot.set_yscale(yscale)

        x_min_lim = utils.return_key(kwargs, "x_min_lim", None)
        x_max_lim = utils.return_key(kwargs, "x_max_lim", None)
        y_min_lim = utils.return_key(kwargs, "y_min_lim", None)
        y_max_lim = utils.return_key(kwargs, "y_max_lim", None)

        ax_to_plot.set_xlim(left=x_min_lim, right=x_max_lim)
        ax_to_plot.set_ylim(bottom=y_min_lim, top=y_max_lim)
        ax_to_plot.legend(
            utils.return_key(kwargs, "legend", []),
            handles=utils.return_key(kwargs, "legend_elements", None),
        )

        if comparison is not None:
            x_min_lim_comparison = utils.return_key(
                kwargs, "x_min_lim_comparison", None
            )
            x_max_lim_comparison = utils.return_key(
                kwargs, "x_max_lim_comparison", None
            )
            y_min_lim_comparison = utils.return_key(
                kwargs, "y_min_lim_comparison", None
            )
            y_max_lim_comparison = utils.return_key(
                kwargs, "y_max_lim_comparison", None
            )

            ax_comparison.set_xlim(
                left=x_min_lim_comparison, right=x_max_lim_comparison
            )
            ax_comparison.set_ylim(
                bottom=y_min_lim_comparison, top=y_max_lim_comparison
            )
        plt.gcf().tight_layout()

    def plot_2d_pk(self, bin_edges, **kwargs):
        """Plot this 2D (k, mu) power spectrum, one line per mu bin.

        For each value in ``bin_edges`` (matched against
        ``self.k_array[1]``), plots power (with error bars, if
        ``error_array`` is set) vs k, optionally multiplied by
        ``k**3 / (2*pi**2)``. If a ``"comparison"`` PowerSpectrum is
        given, also plots the relative difference on a second (ratio)
        axes.

        Args:
            bin_edges: Sequence of mu values identifying which mu slices
                of ``self.k_array[1]`` to plot.
            **kwargs: Options read via ``utils.return_key`` (non-exhaustive):
                ``comparison`` (default None): spectrum to compare against.
                ``k_multiplication`` (default False): multiply power by
                    ``k**3 / (2*pi**2)``.
                ``ax``: passed to :meth:`prepare_axes`.
                ``color``: list of per-mu-bin colors (indexed by position
                    in ``bin_edges``).
                ``linestyle`` (default ``"-"`` for each bin): list of
                    per-mu-bin line styles.
                ``ps`` (default None): marker style.
                ``error_bar_comparison`` (default True): whether to draw
                    error bars on the comparison panel.
                ``x_unit``/``y_unit``/``x_label``/``y_label``/
                    ``labelsize_x``/``labelsize_y``/``fontsize`` and their
                    ``*_comparison`` counterparts: forwarded to
                    :meth:`put_label` for the main and comparison axes.
                ``xscale``/``yscale`` (default ``"log"``): axis scales.
                ``x_min_lim``/``x_max_lim``/``y_min_lim``/``y_max_lim`` and
                    ``*_comparison`` variants: axis limits.
                ``legend`` (default []), ``legend_elements`` (default None).
        """
        comparison = utils.return_key(kwargs, "comparison", None)
        k_multiplication = utils.return_key(kwargs, "k_multiplication", False)
        (ax_to_plot, ax_comparison) = self.prepare_axes(kwargs)
        self.put_label(
            ax_to_plot,
            xunit=utils.return_key(kwargs, "x_unit", True),
            yunit=utils.return_key(kwargs, "y_unit", True),
            x_label=utils.return_key(kwargs, "x_label", r"$k$"),
            y_label=utils.return_key(kwargs, "y_label", r"$P$"),
            labelsize_x=utils.return_key(kwargs, "labelsize_x", 12),
            labelsize_y=utils.return_key(kwargs, "labelsize_y", 12),
            fontsize=utils.return_key(kwargs, "fontsize", 12),
        )
        if comparison is not None:
            self.put_label(
                ax_comparison,
                xunit=utils.return_key(kwargs, "x_unit_comparison", True),
                yunit=utils.return_key(kwargs, "y_unit_comparison", True),
                x_label=utils.return_key(kwargs, "x_label_comparison", r"$k$"),
                y_label=utils.return_key(kwargs, "y_label_comparison", r"$P$"),
                labelsize_x=utils.return_key(kwargs, "labelsize_x_comparison", 12),
                labelsize_y=utils.return_key(kwargs, "labelsize_y_comparison", 12),
                fontsize=utils.return_key(kwargs, "fontsize_comparison", 12),
            )

        for i in range(len(bin_edges)):
            mask = (self.k_array[1] == bin_edges[i]) & (self.power_array != 0.0)
            c = kwargs["color"][i] if "color" in kwargs.keys() else None
            ls = utils.return_key(
                kwargs, "linestyle", ["-" for i in range(len(bin_edges))]
            )[i]
            if k_multiplication:
                factor_multiplication = self.k_array[0][mask] ** 3 / (2 * np.pi**2)
            else:
                factor_multiplication = 1
            if comparison is not None:
                error_bar_comparison = utils.return_key(
                    kwargs, "error_bar_comparison", True
                )
                mask_comparison = (comparison.k_array[1] == bin_edges[i]) & (
                    comparison.power_array != 0.0
                )
                power_array_comparison = interp1d(
                    self.k_array[0][mask],
                    self.power_array[mask],
                    bounds_error=False,
                    fill_value=np.NaN,
                )(comparison.k_array[0][mask_comparison])
                if (
                    (self.error_array is not None)
                    & (comparison.error_array is not None)
                    & error_bar_comparison
                ):
                    error_array_comparison = interp1d(
                        self.k_array[0][mask],
                        self.error_array[mask],
                        bounds_error=False,
                        fill_value=np.NaN,
                    )(comparison.k_array[0][mask_comparison])
                    ax_comparison.errorbar(
                        comparison.k_array[0][mask_comparison],
                        (
                            comparison.power_array[mask_comparison]
                            - power_array_comparison
                        )
                        / comparison.power_array[mask_comparison],
                        (
                            power_array_comparison
                            / comparison.power_array[mask_comparison]
                        )
                        * np.sqrt(
                            (
                                comparison.error_array[mask_comparison]
                                / comparison.power_array[mask_comparison]
                            )
                            ** 2
                            + (error_array_comparison / power_array_comparison) ** 2
                        ),
                        marker=utils.return_key(kwargs, "ps", None),
                        linestyle=ls,
                        color=c,
                    )
                    ax_comparison.plot(
                        [
                            np.min(comparison.k_array[0][mask_comparison]),
                            np.max(comparison.k_array[0][mask_comparison]),
                        ],
                        [0, 0],
                        "k-",
                        alpha=0.5,
                    )
                else:
                    ax_comparison.plot(
                        comparison.k_array[0][mask_comparison],
                        (
                            comparison.power_array[mask_comparison]
                            - power_array_comparison
                        )
                        / comparison.power_array[mask_comparison],
                        marker=utils.return_key(kwargs, "ps", None),
                        linestyle=ls,
                        color=c,
                    )
                    ax_comparison.plot(
                        [
                            np.min(comparison.k_array[0][mask_comparison]),
                            np.max(comparison.k_array[0][mask_comparison]),
                        ],
                        [0, 0],
                        "k-",
                        alpha=0.5,
                    )
            if self.error_array is not None:
                ax_to_plot.errorbar(
                    self.k_array[0][mask],
                    self.power_array[mask] * factor_multiplication,
                    self.error_array[mask] * factor_multiplication,
                    marker=utils.return_key(kwargs, "ps", None),
                    linestyle=ls,
                    color=c,
                )
            else:
                ax_to_plot.plot(
                    self.k_array[0][mask],
                    self.power_array[mask] * factor_multiplication,
                    marker=utils.return_key(kwargs, "ps", None),
                    linestyle=ls,
                    color=c,
                )
        xscale = utils.return_key(kwargs, "xscale", "log")
        yscale = utils.return_key(kwargs, "yscale", "log")

        ax_to_plot.set_xscale(xscale)
        ax_to_plot.set_yscale(yscale)

        x_min_lim = utils.return_key(kwargs, "x_min_lim", None)
        x_max_lim = utils.return_key(kwargs, "x_max_lim", None)
        y_min_lim = utils.return_key(kwargs, "y_min_lim", None)
        y_max_lim = utils.return_key(kwargs, "y_max_lim", None)

        ax_to_plot.set_xlim(left=x_min_lim, right=x_max_lim)
        ax_to_plot.set_ylim(bottom=y_min_lim, top=y_max_lim)

        if comparison is not None:
            x_min_lim_comparison = utils.return_key(
                kwargs, "x_min_lim_comparison", None
            )
            x_max_lim_comparison = utils.return_key(
                kwargs, "x_max_lim_comparison", None
            )
            y_min_lim_comparison = utils.return_key(
                kwargs, "y_min_lim_comparison", None
            )
            y_max_lim_comparison = utils.return_key(
                kwargs, "y_max_lim_comparison", None
            )

            ax_comparison.set_xlim(
                left=x_min_lim_comparison, right=x_max_lim_comparison
            )
            ax_comparison.set_ylim(
                bottom=y_min_lim_comparison, top=y_max_lim_comparison
            )

        ax_to_plot.legend(
            utils.return_key(kwargs, "legend", []),
            handles=utils.return_key(kwargs, "legend_elements", None),
            fontsize=utils.return_key(kwargs, "fontsize", 12),
        )
        plt.gcf().tight_layout()

    def plot_several_power_spectrum(self, Pks, k_space, name, legend):
        """Plot several power spectra sharing the same k grid to a new figure.

        Creates a new figure, log-log plots each spectrum in ``Pks``
        against ``k_space`` with a rainbow color cycle, and saves the
        result as ``"<name>matter_power_spectrum.pdf"``.

        Args:
            Pks: Sequence of power arrays, one per spectrum to plot.
            k_space: Shared wavenumber array (x axis) for all spectra.
            name: Filename prefix for the saved PDF.
            legend: Sequence of legend labels, one per spectrum.
        """
        plt.figure()
        color = cm.rainbow(np.linspace(0, 1, len(Pks)))
        for i in range(len(Pks)):
            plt.loglog(k_space, np.array(Pks[i]), "b", color=color[i])
        plt.grid()
        plt.legend(legend)
        plt.savefig(name + "matter_power_spectrum.pdf", format="pdf")

    def plot_comparison_spectra(
        self, list_spectra, label_list, diff_extremums=0.1, normalize=True
    ):
        """Create a two-panel figure comparing this spectrum to others.

        Builds a figure with a main panel (top, 3/4 height) showing all
        spectra and a ratio panel (bottom, 1/4 height) showing each
        spectrum's fractional difference to ``self`` (delegated to
        :meth:`add_comparison_spectra`), then applies axis labels,
        scales, legend and ratio-panel y-limits.

        Args:
            list_spectra: Sequence of PowerSpectrum instances to compare
                against ``self`` (the reference).
            label_list: Legend labels, one per plotted spectrum
                (including the reference, first).
            diff_extremums: Symmetric y-limit for the ratio panel.
                Defaults to 0.1.
            normalize: If True, plot ``k**3 * P(k) / 2 * pi**2`` instead
                of ``P(k)`` (as computed in :meth:`add_comparison_spectra`).
                Defaults to True.

        Returns:
            numpy.ndarray: The two-element array of matplotlib Axes
            ``[main_ax, ratio_ax]``.
        """
        fig, ax = plt.subplots(
            2, 1, gridspec_kw={"height_ratios": [3, 1]}, sharex=True, figsize=(8, 6)
        )  # note that height ratios can be used to scale the size of top vs bottom part
        self.add_comparison_spectra(list_spectra, ax, normalize=normalize)
        ax[0].set_title(r"...")
        if normalize:
            ax[0].set_ylabel(r"$\Delta_m^2$")
        else:
            ax[0].set_ylabel(r"$P_m$")
        if self.h_normalized:
            ax[1].set_xlabel("k (h Mpc-1)")
        else:
            ax[1].set_xlabel("k (Mpc-1)")
        ax[1].set_ylabel(r"$\Delta_m^2/\Delta_{m,ref}^2-1$")
        ax[0].set_xscale("log")
        ax[0].set_yscale("log")
        if len(label_list) <= 5:
            ax[0].legend(label_list)
        else:
            ax[0].legend(label_list, ncol=2)
        ax[1].set_ylim(-diff_extremums, diff_extremums)
        return ax

    def add_comparison_spectra(self, list_spectra, ax, normalize=True):
        """Overlay this spectrum and others on existing axes, with a ratio panel.

        Plots ``self`` and each spectrum in ``list_spectra`` on
        ``ax[0]``, and each spectrum's ratio to ``self`` (interpolated
        onto ``self``'s k grid, minus 1) on ``ax[1]``.

        Args:
            list_spectra: Sequence of PowerSpectrum instances to overlay
                and compare against ``self``.
            ax: Two-element sequence of matplotlib Axes,
                ``[main_ax, ratio_ax]``.
            normalize: If True, plot ``k**3 * P(k) / 2 * pi**2`` instead
                of ``P(k)``. Defaults to True.
        """
        kref = self.k_array
        if normalize:
            kpkref = (self.k_array**3 * self.power_array) / 2 * (np.pi) ** 2
        else:
            kpkref = self.power_array
        karr, kpkarr = [], []
        karr.append(kref)
        kpkarr.append(kpkref)
        for i in range(len(list_spectra)):
            karr.append(list_spectra[i].k_array)
            if normalize:
                kpkarr.append(
                    (list_spectra[i].k_array ** 3 * list_spectra[i].power_array)
                    / 2
                    * (np.pi) ** 2
                )
            else:
                kpkarr.append(list_spectra[i].power_array)
        # karr is your array of x values, i.e. a numpy array with shape (nlines,nvalues)
        # kpkarr is your array of y values same shape (nlines,nvalues)
        # larr is your array of labels (nlines)
        # kref,kpkref are the reference values (nvalues)
        for k, kpk in zip(karr, kpkarr):
            interp = interp1d(k, kpk, bounds_error=False)
            ax[0].plot(k, kpk)
            ax[1].plot(kref, (interp(kref) / kpkref) - 1)

    def save_plot(self, nameout, format_out="pdf", fig=None):
        """Save a matplotlib figure to disk.

        Args:
            nameout: Output file path.
            format_out: File format passed to ``Figure.savefig``.
                Defaults to ``"pdf"``.
            fig: Figure to save. Defaults to the current figure
                (``plt.gcf()``) if None.
        """
        if fig is None:
            fig = plt.gcf()
        fig.savefig(nameout, format=format_out)

    def close_plot(self, fig=None):
        """Close the current matplotlib figure.

        Args:
            fig: Unused (accepted for API symmetry with :meth:`save_plot`);
                the current figure is always resolved internally and
                ``plt.close()`` closes the current figure regardless.
        """
        if fig is None:
            fig = plt.gcf()
        plt.close()

    def open_plot(self, **kwargs):
        """Create a new matplotlib figure, optionally applying a style.

        Args:
            **kwargs: Options read via ``utils.return_key``:
                ``style`` (default None): matplotlib style name passed to
                    ``plt.style.use``.
                ``figsize`` (default ``(8, 6)``): figure size in inches.

        Returns:
            matplotlib.figure.Figure: The newly created figure.
        """
        style = utils.return_key(kwargs, "style", None)
        if style is not None:
            plt.style.use(style)
        figsize = utils.return_key(kwargs, "figsize", (8, 6))
        fig = plt.figure(figsize=figsize)
        return fig

    def open_subplot(self, x=2, y=1, figsize=(8, 6)):
        """Create a new figure with a grid of x-by-y subplots sharing the x axis.

        Args:
            x: Number of subplot rows. Defaults to 2.
            y: Number of subplot columns. Defaults to 1.
            figsize: Figure size in inches. Defaults to ``(8, 6)``.

        Returns:
            matplotlib.figure.Figure: The newly created figure.
        """
        fig, ax = plt.subplots(x, y, sharex=True, figsize=figsize)
        return fig

    def show_plot(self):
        """Display the current matplotlib figure (``plt.show()``)."""
        plt.show()

    def get_k_value(self, k):
        """Interpolate and print the power at a given k value (or values).

        Args:
            k: Wavenumber value(s) at which to evaluate the power
                spectrum. Must lie within the range of ``self.k_array``
                (``bounds_error=True``).

        Returns:
            float or numpy.ndarray: Interpolated power value(s).
        """
        interp = interp1d(self.k_array, self.power_array, bounds_error=True)
        print(interp(k))
        return interp(k)

    def change_k_normalization(self, wanted_h_normalized, h):
        """Convert ``k_array`` between h/Mpc and 1/Mpc units, in place.

        Args:
            wanted_h_normalized: Target normalization: True for h/Mpc,
                False for 1/Mpc.
            h: Dimensionless Hubble parameter used for the conversion
                (``k_array`` is divided by ``h`` to go from h/Mpc to
                1/Mpc, multiplied by ``h`` for the reverse).

        Returns:
            tuple: Empty tuple ``()`` in all cases (no-op sentinel).

        Raises:
            KeyError: If ``self.h_normalized`` is None (unknown current
                normalization).
        """
        if self.h_normalized is None:
            raise KeyError("The actual normalization of the k vector is not know")
        if self.h_normalized:
            if wanted_h_normalized:
                return ()
            else:
                self.k_array = self.k_array / h
                self.h_normalized = False
                return ()
        else:
            if wanted_h_normalized:
                self.k_array = self.k_array * h
                self.h_normalized = True
                return ()
            else:
                return ()


class MatterPowerSpectrum(PowerSpectrum):
    """Matter power spectrum for a given species, 1D or 3D.

    Extends :class:`PowerSpectrum` with a ``dimension`` (``"1D"`` or
    ``"3D"``) and a ``specie`` label (e.g. dark matter, baryons,
    neutrinos), and adds gimlet-format I/O.
    """

    def __init__(self, dimension, specie, **kwargs):
        """Initialize a MatterPowerSpectrum.

        Args:
            dimension: Spectrum dimensionality, either ``"1D"`` or ``"3D"``.
            specie: Label of the matter species this spectrum belongs to.
            **kwargs: Forwarded to :class:`PowerSpectrum.__init__`.

        Raises:
            KeyError: If ``dimension`` is not ``"1D"`` or ``"3D"``.
        """
        super(MatterPowerSpectrum, self).__init__(**kwargs)
        if dimension not in ["1D", "3D"]:
            raise KeyError(
                "Dimension of spectrum not available, please choose between 1D and 3D"
            )
        self.dimension = dimension
        self.specie = specie

    @classmethod
    def init_from_gimlet(
        cls,
        namefile,
        specie="unknown",
        power_weighted=False,
        error_estimator=None,
        **kwargs,
    ):
        """Pm(k) gimlet file contains
        - k: edge (higher) of the k bin considered
        - bincount: number of mode (pairs) computed in the bin
        - pwk: power weighted k
        - power: power of the bin"""
        f = np.loadtxt(namefile)
        k, bincount, pwk, power = f[:, 0], f[:, 1], f[:, 2], f[:, 3]
        if error_estimator is not None:
            error = utils.error_estimator(
                power, model=error_estimator, bin_count=bincount, **kwargs
            )
        else:
            error = None
        if power_weighted:
            k_array = pwk
        else:
            k_array = k
        dimension = "3D"
        h_normalized = True
        return cls(
            dimension,
            specie,
            k_array=k_array,
            power_array=power,
            error_array=error,
            file_init=namefile,
            size_box=None,
            h_normalized=h_normalized,
        )

    def write_to_gimlet(self, name_out, power_weighted=False):
        """Write this spectrum to a gimlet-format Pm(k) ascii file.

        Writes columns ``[k, bincount, pwk, power]`` where either ``k``
        or ``pwk`` (power-weighted k) is populated from ``self.k_array``
        depending on ``power_weighted``, and ``bincount`` is taken from
        ``self.error_array`` if present, else zeros.

        Args:
            name_out: Output file path.
            power_weighted: If True, store ``self.k_array`` in the
                power-weighted-k column (``pwk``) and zero the plain
                ``k`` column; if False, do the reverse. Defaults to False.
        """
        if power_weighted:
            pwk = self.k_array
            k = np.zeros(self.k_array.shape)
        else:
            pwk = np.zeros(self.k_array.shape)
            k = self.k_array
        if self.error_array is not None:
            bincount = self.error_array
        else:
            bincount = np.zeros(self.power_array.shape)
        power = self.power_array
        out = np.transpose(np.stack([k, bincount, pwk, power]))
        np.savetxt(name_out, out)


class FluxPowerSpectrum(PowerSpectrum):
    """Lyman-alpha flux power spectrum, 1D or 3D (k or k, mu binned).

    Extends :class:`PowerSpectrum` with a ``dimension`` (``"1D"`` or
    ``"3D"``) and adds gimlet-format readers/writers for both P(k) and
    P(k, mu) or P(k_perp, k_par) binnings, plus multi-file averaging.
    """

    def __init__(self, dimension, **kwargs):
        """Initialize a FluxPowerSpectrum.

        Args:
            dimension: Spectrum dimensionality, either ``"1D"`` or ``"3D"``.
            **kwargs: Forwarded to :class:`PowerSpectrum.__init__`.

        Raises:
            KeyError: If ``dimension`` is not ``"1D"`` or ``"3D"``.
        """
        super(FluxPowerSpectrum, self).__init__(**kwargs)
        if dimension not in ["1D", "3D"]:
            raise KeyError(
                "Dimension of spectrum not available, please choose between 1D and 3D"
            )
        self.dimension = dimension

    @classmethod
    def init_1D_from_gimlet(
        cls,
        namefile,
        power_weighted=False,
        error_estimator=None,
        error_stored=False,
        **kwargs,
    ):
        """Load a 1D flux power spectrum P(k) from a gimlet ascii file.

        The file is expected to have columns
        ``[k_edge, bincount, pwk_edge, power]``.

        Args:
            namefile: Path to the gimlet Pf(k) ascii file.
            power_weighted: If True, use the power-weighted k column
                (``pwk_edge``) as ``k_array``; otherwise use the plain
                k-edge column. Defaults to False.
            error_estimator: Name of the error model to pass to
                ``utils.error_estimator`` (e.g. ``"uncorrelated"``,
                ``"constant"``, ``"computed"``, ``"computed_epsilon"``).
                If None and ``error_stored`` is False, no error is set.
            error_stored: If True, use the ``bincount`` column directly
                as the error array instead of computing one. Defaults
                to False.
            **kwargs: Forwarded to ``utils.error_estimator`` (e.g.
                ``epsilon``).

        Returns:
            FluxPowerSpectrum: New 1D instance with ``h_normalized=True``.
        """
        pk_array = np.loadtxt(namefile)
        k_edge, bincount, pwk_edge, power = (
            pk_array[:, 0],
            pk_array[:, 1],
            pk_array[:, 2],
            pk_array[:, 3],
        )

        if power_weighted:
            k_array = pwk_edge
        else:
            k_array = k_edge
        if (error_estimator is not None) & (not (error_stored)):
            error = utils.error_estimator(
                power, model=error_estimator, bin_count=bincount, **kwargs
            )
        elif error_stored:
            error = bincount
        else:
            error = None
        dimension = "1D"
        h_normalized = True
        return cls(
            dimension,
            k_array=k_array,
            power_array=power,
            error_array=error,
            file_init=namefile,
            size_box=None,
            h_normalized=h_normalized,
        )

    @classmethod
    def init_3D_from_gimlet(
        cls,
        namefile,
        type_file,
        kmu=True,
        power_weighted=False,
        error_estimator=None,
        field_name=None,
        error_stored=False,
        **kwargs,
    ):
        """Load a 3D flux power spectrum P(k, mu) or P(k_perp, k_par) from gimlet.

        Reads a gimlet output (plain text or HDF5 dataset) with 6
        columns, parsed either as ``(k, mu, ...)`` via
        :meth:`init_kmu` or ``(k_perp, k_par, ...)`` via
        :meth:`init_kperpar`, and builds a 2-row ``k_array`` of
        ``[k1, k2]``.

        Args:
            namefile: Path to the gimlet Pf(k, mu) file.
            type_file: File format, either ``"txt"`` or ``"hdf5"``.
            kmu: If True, parse columns as (k, mu, bincount, pwk, pwmu,
                power) via :meth:`init_kmu`; if False, parse as
                (k_perp, k_par, bincount, pwkperp, pwkpar, power) via
                :meth:`init_kperpar`. Defaults to True.
            power_weighted: If True, use the power-weighted k1/k2 columns
                as ``k_array``; otherwise use the plain edge columns.
                Defaults to False.
            error_estimator: Name of the error model to pass to
                ``utils.error_estimator``. If None and ``error_stored``
                is False, no error is set.
            field_name: Dataset name to read within the HDF5 file when
                ``type_file == "hdf5"``.
            error_stored: If True, use the ``bincount`` column directly
                as the error array instead of computing one. Defaults
                to False.
            **kwargs: Forwarded to ``utils.error_estimator``.

        Returns:
            FluxPowerSpectrum: New 3D instance with ``h_normalized=True``.
        """
        if type_file == "txt":
            pk_array = np.loadtxt(namefile)
        elif type_file == "hdf5":
            file = h5py.File(namefile, "r")[field_name]
            pk_array = np.array(list(zip(*file))).transpose()
        if kmu:
            (
                k1_edge,
                k2_edge,
                bincount,
                pwk1,
                pwk2,
                power,
            ) = FluxPowerSpectrum.init_kmu(pk_array)
        else:
            (
                k1_edge,
                k2_edge,
                bincount,
                pwk1,
                pwk2,
                power,
            ) = FluxPowerSpectrum.init_kperpar(pk_array)
        if power_weighted:
            k1_array, k2_array = pwk1, pwk2
        else:
            k1_array, k2_array = k1_edge, k2_edge
        if (error_estimator is not None) & (not (error_stored)):
            error = utils.error_estimator(
                power, model=error_estimator, bin_count=bincount, **kwargs
            )
        elif error_stored:
            error = bincount
        else:
            error = None
        k_array = np.stack([k1_array, k2_array])
        dimension = "3D"
        h_normalized = True
        return cls(
            dimension,
            k_array=k_array,
            power_array=power,
            error_array=error,
            file_init=namefile,
            size_box=None,
            h_normalized=h_normalized,
        )

    @staticmethod
    def init_kmu(pk_array):
        """Pf(k,mu) gimlet file contains
        - k_edge: edge (lower) of the k bin considered
        - mu_edge: edge (lower) of the mu bin considered (mu positive)
        - bincount: number of mode (pairs) computed in the bin
        - pwk: power weighted k
        - pwmu: power weighted mu
        - power: power of the bin"""
        k_edge, mu_edge, bincount, pwk, pwmu, power = (
            pk_array[:, 0],
            pk_array[:, 1],
            pk_array[:, 2],
            pk_array[:, 3],
            pk_array[:, 4],
            pk_array[:, 5],
        )
        return (k_edge, mu_edge, bincount, pwk, pwmu, power)

    @staticmethod
    def init_kperpar(pk_array):
        """Pf(kperp,kpar) gimlet file contains
        - k_perp: edge (lower) of the k perp bin considered
        - k_par: edge (lower) of the k par bin considered
        - bincount: number of mode (pairs) computed in the bin
        - pwkperp: power weighted k perp
        - pwkpar: power weighted k par
        - power: power of the bin"""
        k_perp, k_par, bincount, pwkperp, pwkpar, power = (
            pk_array[:, 0],
            pk_array[:, 1],
            pk_array[:, 2],
            pk_array[:, 3],
            pk_array[:, 4],
            pk_array[:, 5],
        )
        return (k_perp, k_par, bincount, pwkperp, pwkpar, power)

    @staticmethod
    def compute_mean_gimlet(namefile, namemean, type_file, kmu=True, field_name=None):
        """Average several gimlet 3D power-spectrum files and write the result.

        Reads each file in ``namefile`` (plain text or HDF5), parses it
        with :meth:`init_kmu` or :meth:`init_kperpar`, sums bincount,
        the two power-weighted-k columns, and power across files, then
        divides the weighted-k and power sums by the number of files
        (bincount is left as a raw sum) and writes the result to
        ``namemean`` in gimlet ascii format
        ``[k1_edge, k2_edge, bincount, pwk1, pwk2, power]`` (k edges
        taken from the first file only).

        Args:
            namefile: Sequence of input gimlet file paths to average.
            namemean: Output ascii file path for the averaged spectrum.
            type_file: File format of the inputs, ``"txt"`` or ``"hdf5"``.
            kmu: If True, parse as (k, mu, ...) via :meth:`init_kmu`; if
                False, parse as (k_perp, k_par, ...) via
                :meth:`init_kperpar`. Defaults to True.
            field_name: Dataset name to read within each HDF5 file when
                ``type_file == "hdf5"``.
        """
        if type_file == "txt":
            pk_array = np.loadtxt(namefile[0])
        elif type_file == "hdf5":
            file = h5py.File(namefile[0], "r")[field_name]
            pk_array = np.array(list(zip(*file))).transpose()

        if kmu:
            (
                k1_edge_mean,
                k2_edge_mean,
                bincount_mean,
                pwk1_mean,
                pwk2_mean,
                power_mean,
            ) = FluxPowerSpectrum.init_kmu(pk_array)
        else:
            (
                k1_edge_mean,
                k2_edge_mean,
                bincount_mean,
                pwk1_mean,
                pwk2_mean,
                power_mean,
            ) = FluxPowerSpectrum.init_kperpar(pk_array)
        for j in range(1, len(namefile)):
            if type_file == "txt":
                pk_array = np.loadtxt(namefile[j])
            elif type_file == "hdf5":
                file = h5py.File(namefile[j], "r")[field_name]
                pk_array = np.array(list(zip(*file))).transpose()
            if kmu:
                indiv_value = FluxPowerSpectrum.init_kmu(pk_array)
            else:
                indiv_value = FluxPowerSpectrum.init_kperpar(pk_array)

            bincount, pwk1, pwk2, power = (
                indiv_value[2],
                indiv_value[3],
                indiv_value[4],
                indiv_value[5],
            )
            bincount_mean = bincount_mean + bincount
            pwk1_mean = pwk1_mean + pwk1
            pwk2_mean = pwk2_mean + pwk2
            power_mean = power_mean + power
        pwk1_mean = pwk1_mean / len(namefile)
        pwk2_mean = pwk2_mean / len(namefile)
        power_mean = power_mean / len(namefile)
        out = np.transpose(
            np.stack(
                [
                    k1_edge_mean,
                    k2_edge_mean,
                    bincount_mean,
                    pwk1_mean,
                    pwk2_mean,
                    power_mean,
                ]
            )
        )
        np.savetxt(namemean, out)

    def write_to_gimlet(self, name_out, power_weighted=False):
        """Write this 3D spectrum to a gimlet-format Pf(k, mu) ascii file.

        Writes columns
        ``[k1_edge, k2_edge, bincount, pwk1, pwk2, power]`` where either
        the plain edges or the power-weighted columns are populated
        from ``self.k_array`` depending on ``power_weighted``, and
        ``bincount`` is taken from ``self.error_array`` if present,
        else zeros.

        Args:
            name_out: Output file path.
            power_weighted: If True, store ``self.k_array`` rows in the
                power-weighted columns (``pwk1``, ``pwk2``) and zero the
                plain edge columns; if False, do the reverse. Defaults
                to False.
        """
        if power_weighted:
            pwk1 = self.k_array[0]
            pwk2 = self.k_array[1]
            k1_edge = np.zeros(self.k_array[0].shape)
            k2_edge = np.zeros(self.k_array[1].shape)
        else:
            pwk1 = np.zeros(self.k_array[0].shape)
            pwk2 = np.zeros(self.k_array[1].shape)
            k1_edge = self.k_array[0]
            k2_edge = self.k_array[1]
        if self.error_array is not None:
            bincount = self.error_array
        else:
            bincount = np.zeros(self.power_array.shape)
        power = self.power_array
        out = np.transpose(np.stack([k1_edge, k2_edge, bincount, pwk1, pwk2, power]))
        np.savetxt(name_out, out)


def init_spectrum(type_init, filename, boxsize=None):
    """
    The init_spectrum function takes in a type_init and filename, and returns a PowerSpectrum object.

    Args:
        type_init: Determine what type of file is being read in
        filename: Specify the file to read from
        boxsize: Specify the boxsize of the simulation

    Returns:
        A powerspectrum object

    Doc Author:
        Trelent
    """
    if type_init == "GENPK":
        spectrum = PowerSpectrum.init_from_genpk_file(filename, boxsize)
    elif type_init == "ASCII":
        spectrum = PowerSpectrum.init_from_ascii_file(filename)
    return spectrum


def launch_comparison_power_spectra(
    list_file,
    type_file,
    label_list,
    name_out,
    diff_extremums=0.1,
    rebin=None,
    rebin_method=None,
    flux_factor=None,
    normalize=True,
    size_box=None,
    wanted_normalization=None,
    h_normalization=None,
):
    """Load, optionally rebin/rescale, and plot several spectra against a reference.

    Loads ``list_file[0]`` as the reference spectrum and each remaining
    entry of ``list_file`` as a comparison spectrum (all via
    :func:`init_spectrum`), optionally rebins them (:meth:`PowerSpectrum.rebin_arrays`),
    rescales power by a per-file ``flux_factor``, and converts k
    normalization, then plots them with
    :meth:`PowerSpectrum.plot_comparison_spectra` and saves the figure.

    Args:
        list_file: Sequence of file paths; ``list_file[0]`` is the
            reference spectrum, the rest are compared against it.
        type_file: File type passed to :func:`init_spectrum`
            (``"GENPK"`` or ``"ASCII"``).
        label_list: Legend labels, one per spectrum (including the
            reference).
        name_out: Output path for the saved comparison figure.
        diff_extremums: Symmetric y-limit for the ratio panel. Defaults
            to 0.1.
        rebin: Number of bins to rebin each spectrum to, or None to skip
            rebinning. Defaults to None.
        rebin_method: Aggregation method forwarded to
            :meth:`PowerSpectrum.rebin_arrays` (e.g. ``"mean"``, ``"gauss"``).
        flux_factor: Sequence of per-file multiplicative factors applied
            to each spectrum's power array, or None to skip. Defaults
            to None.
        normalize: Forwarded to
            :meth:`PowerSpectrum.plot_comparison_spectra`. Defaults to
            True.
        size_box: Box size (Mpc/h) forwarded to :func:`init_spectrum`
            for GENPK files.
        wanted_normalization: Target h-normalization forwarded to
            :meth:`PowerSpectrum.change_k_normalization`, or None to
            skip conversion.
        h_normalization: Hubble parameter ``h`` forwarded to
            :meth:`PowerSpectrum.change_k_normalization`.
    """
    reference_spectrum = init_spectrum(type_file, list_file[0], boxsize=size_box)
    if rebin is not None:
        reference_spectrum.rebin_arrays(rebin, operation=rebin_method)
    if flux_factor is not None:
        reference_spectrum.power_array = reference_spectrum.power_array * flux_factor[0]
    if wanted_normalization is not None:
        reference_spectrum.change_k_normalization(wanted_normalization, h_normalization)
    list_spectra = []
    for i in range(1, len(list_file)):
        ps = init_spectrum(type_file, list_file[i], boxsize=size_box)
        if rebin is not None:
            ps.rebin_arrays(rebin, operation=rebin_method)
        if flux_factor is not None:
            ps.power_array = ps.power_array * flux_factor[i]
        if wanted_normalization is not None:
            ps.change_k_normalization(wanted_normalization, h_normalization)
        list_spectra.append(ps)
    reference_spectrum.plot_comparison_spectra(
        list_spectra, label_list, diff_extremums=diff_extremums, normalize=normalize
    )
    reference_spectrum.save_plot(name_out)
    reference_spectrum.close_plot()


def launch_comparison_power_spectra_different_ref(
    list_file,
    type_file,
    label_list,
    name_out,
    diff_extremums=0.1,
    rebin=None,
    rebin_method=None,
    flux_factor=None,
    normalize=True,
    size_box=None,
    wanted_normalization=None,
    h_normalization=None,
):
    """Compare several groups of spectra, each against its own reference.

    Like :func:`launch_comparison_power_spectra`, but operates on a
    list of groups (``list_file[j]``), each with its own reference
    (``list_file[j][0]``), rebin settings, flux factors and target
    normalization (all indexed by ``j``). All groups are overlaid on
    the same figure: the first group creates the comparison figure via
    :meth:`PowerSpectrum.plot_comparison_spectra`, subsequent groups are
    added via :meth:`PowerSpectrum.add_comparison_spectra`. The figure
    is saved using the last group's reference spectrum.

    Args:
        list_file: Sequence of groups; each group is a sequence of file
            paths whose first entry is that group's reference spectrum.
        type_file: Sequence of file types (one per group), passed to
            :func:`init_spectrum`.
        label_list: Legend labels forwarded to
            :meth:`PowerSpectrum.plot_comparison_spectra` for the first group.
        name_out: Output path for the saved comparison figure.
        diff_extremums: Symmetric y-limit for the ratio panel. Defaults
            to 0.1.
        rebin: Sequence of rebin bin counts (one per group, or None
            entries to skip rebinning that group).
        rebin_method: Aggregation method forwarded to
            :meth:`PowerSpectrum.rebin_arrays`.
        flux_factor: Sequence of per-group sequences of per-file power
            multipliers, or None entries to skip.
        normalize: Forwarded to the comparison plotting calls. Defaults
            to True.
        size_box: Box size (Mpc/h) forwarded to :func:`init_spectrum`
            for GENPK files.
        wanted_normalization: Sequence of target h-normalizations (one
            per group), or None to skip conversion.
        h_normalization: Hubble parameter ``h`` forwarded to
            :meth:`PowerSpectrum.change_k_normalization`.
    """
    for j in range(len(list_file)):
        reference_spectrum = init_spectrum(
            type_file[j], list_file[j][0], boxsize=size_box
        )
        if rebin[j] is not None:
            reference_spectrum.rebin_arrays(rebin[j], operation=rebin_method)
        if flux_factor[j] is not None:
            reference_spectrum.power_array = (
                reference_spectrum.power_array * flux_factor[j][0]
            )
        if wanted_normalization is not None:
            reference_spectrum.change_k_normalization(
                wanted_normalization[j], h_normalization
            )
        list_spectra = []
        for i in range(1, len(list_file[j])):
            ps = init_spectrum(type_file[j], list_file[j][i], boxsize=size_box)
            if rebin[j] is not None:
                ps.rebin_arrays(rebin[j], operation=rebin_method)
            if flux_factor[j] is not None:
                ps.power_array = ps.power_array * flux_factor[j][i]
            if wanted_normalization is not None:
                ps.change_k_normalization(wanted_normalization[j], h_normalization)
            list_spectra.append(ps)
        if j == 0:
            ax = reference_spectrum.plot_comparison_spectra(
                list_spectra,
                label_list,
                diff_extremums=diff_extremums,
                normalize=normalize,
            )
        else:
            reference_spectrum.add_comparison_spectra(
                list_spectra, ax, normalize=normalize
            )
    reference_spectrum.save_plot(name_out)
    reference_spectrum.close_plot()


def compute_k_extremums(
    power_Ll,
    power_Sl,
    power_Ss,
    tol=0.01,
):
    """
    The compute_k_extremums function takes in the power spectra of the long-long, short-long and short-short modes.
    It then computes a list of k_max values for each mu bin. The k_max value is defined as the maximum value of
    k where P(k) = P(Ss)(k). This function is used to compute an upper limit on our integration range when computing
    the covariance matrix.

    Args:
        power_Ll: Compute the k_max for each mu bin
        power_Sl: Compute the interpolation of power_sl
        power_Ss: Find the minimum and maximum k values for each mu bin
        tol: Compute the k_max value
        : Compute the upper limit on our integration range when computing

    Returns:
        A list of k_max values for each mu bin

    Doc Author:
        Trelent
    """
    mu_bins = np.unique(power_Ll.k_array[1])
    k_max = []
    for i in range(len(mu_bins)):
        mask_mu_Ll = power_Ll.k_array[1] == mu_bins[i]
        mask_mu_Sl = power_Sl.k_array[1] == mu_bins[i]
        mask_mu_Ss = power_Ss.k_array[1] == mu_bins[i]

        min_k = np.min(power_Ss.k_array[0][mask_mu_Ss])
        max_k = np.max(power_Ss.k_array[0][mask_mu_Ss])

        mask_k_Ll = (
            mask_mu_Ll & (power_Ll.k_array[0] >= min_k) & (power_Ll.k_array[0] <= max_k)
        )

        power_Sl_interp = interp1d(
            power_Sl.k_array[0][mask_mu_Sl], power_Sl.power_array[mask_mu_Sl]
        )
        mask_k_select = (
            np.abs(
                (
                    power_Ll.power_array[mask_k_Ll]
                    - power_Sl_interp(power_Ll.k_array[0][mask_k_Ll])
                )
                / power_Ll.power_array[mask_k_Ll]
            )
            < tol
        )
        k_max.append(np.max(power_Ll.k_array[0][mask_k_Ll][mask_k_select]))

    return k_max


def compute_k_extremums_1D(
    power_Ll,
    power_Sl,
    power_Ss,
    tol=0.01,
):
    """
    The compute_k_extremums_1D function computes the maximum k value for a given power spectrum.

    Args:
        power_Ll: Compute the k_max value
        power_Sl: Compute the maximum k value
        power_Ss: Compute the maximum k value
        tol: Determine the maximum k value
        : Compute the maximum k value for a given power spectrum

    Returns:
        The maximum k value

    Doc Author:
        Trelent
    """
    min_k = np.min(power_Ss.k_array)
    max_k = np.max(power_Ss.k_array)

    mask_k_Ll = (power_Ll.k_array >= min_k) & (power_Ll.k_array <= max_k)

    power_Sl_interp = interp1d(power_Sl.k_array, power_Sl.power_array)
    mask_k_select = (
        np.abs(
            (
                power_Ll.power_array[mask_k_Ll]
                - power_Sl_interp(power_Ll.k_array[mask_k_Ll])
            )
            / power_Ll.power_array[mask_k_Ll]
        )
        < tol
    )

    k_max = np.max(power_Ll.k_array[mask_k_Ll][mask_k_select])

    return k_max


def splice_1D(
    power_Ll,
    power_Sl,
    power_Ss,
    size_small,
    size_large,
    N_large,
    use_nyquist=False,
    tol=0.01,
):
    """L,S = Large or Small size
    l,s = large or small number of particles/resolution elements
    splice the Ll box, using resolved Sl box and splicing Ss box"""
    kmin_S = 2 * np.pi / size_small
    if use_nyquist:
        knyq_L = N_large * np.pi / size_large
        k_max = knyq_L / 4
    else:
        k_max = compute_k_extremums_1D(power_Ll, power_Sl, power_Ss, tol=tol)
    power = []
    k_array = []
    error = None
    if (power_Ll.error_array is not None) & (power_Sl.error_array is not None):
        error = []
    power_Ll_interp = interp1d(power_Ll.k_array, power_Ll.power_array)
    power_Ss_interp = interp1d(power_Ss.k_array, power_Ss.power_array)
    power_Sl_interp = interp1d(power_Sl.k_array, power_Sl.power_array)

    ## low k:  k <= kminS
    mask_k_Ll = power_Ll.k_array <= kmin_S
    power.append(
        power_Ll.power_array[mask_k_Ll]
        * (power_Sl_interp(kmin_S) / power_Ss_interp(kmin_S))
    )
    k_array.append(power_Ll.k_array[mask_k_Ll])
    if error is not None:
        error.append(power_Ll.error_array[mask_k_Ll])

    ## mid k:  kminS < k <= kNyqL / 4
    mask_k_Ll = (power_Ll.k_array > kmin_S) & (power_Ll.k_array <= k_max)
    k = power_Ll.k_array[mask_k_Ll]
    power.append(
        power_Ll.power_array[mask_k_Ll] * (power_Sl_interp(k) / power_Ss_interp(k))
    )
    k_array.append(k)
    if error is not None:
        error.append(power_Ll.error_array[mask_k_Ll])

    ## large k:  k > kNyqL / 4
    mask_k_Sl = power_Sl.k_array > k_max
    power.append(
        power_Sl.power_array[mask_k_Sl]
        * (power_Ll_interp(k_max) / power_Ss_interp(k_max))
    )
    k_array.append(power_Sl.k_array[mask_k_Sl])
    if error is not None:
        error.append(power_Sl.error_array[mask_k_Sl])
        error = np.concatenate(error, axis=0)
    power = np.concatenate(power, axis=0)
    k_array = np.concatenate(k_array, axis=0)

    power_spectrum = MatterPowerSpectrum(
        "1D",
        "matter",
        k_array=k_array,
        power_array=power,
        error_array=error,
        h_normalized=True,
    )
    return power_spectrum


def splice_3D(
    power_Ll,
    power_Sl,
    power_Ss,
    size_small,
    size_large,
    N_large,
    use_nyquist=False,
    impose_kmin_coeff=None,
    impose_kmin=None,
    impose_kmax=None,
    tol=0.01,
    power_pwk_Ll=None,
    power_pwk_Sl=None,
):
    """L,S = Large or Small size
    l,s = large or small number of particles/resolution elements
    splice the Ll box, using resolved Sl box and splicing Ss box"""
    if impose_kmin is not None:
        kmin_S = impose_kmin
    else:
        if impose_kmin_coeff is None:
            kmin_S = 2 * np.pi / size_small
        else:
            kmin_S = impose_kmin_coeff * 2 * np.pi / size_small
    if impose_kmax is not None:
        k_max = impose_kmax
    else:
        if use_nyquist:
            knyq_L = N_large * np.pi / size_large
            k_max = knyq_L / 4
        else:
            k_max_array = compute_k_extremums(power_Ll, power_Sl, power_Ss, tol=tol)
            k_max = np.min(k_max_array)
    mu_value = []
    power = []
    k_array = []
    error = None
    if power_pwk_Ll is not None:
        k_array_pwk = []
        mu_array_pwk = []

    if (power_Ll.error_array is not None) & (power_Sl.error_array is not None):
        error = []
    mu_bins = np.unique(power_Ll.k_array[1])
    for i in range(len(mu_bins)):
        power_mu, k_array_mu = [], []
        if error is not None:
            error_mu = []
        if power_pwk_Ll is not None:
            k_array_pwk_mu = []
            mu_array_pwk_mu = []            
        mask_mu_Ll = power_Ll.k_array[1] == mu_bins[i]
        mask_mu_Sl = power_Sl.k_array[1] == mu_bins[i]
        mask_mu_Ss = power_Ss.k_array[1] == mu_bins[i]
        mu_value.append(power_Ll.k_array[1][mask_mu_Ll].mean())
        power_Ll_interp = interp1d(
            power_Ll.k_array[0][mask_mu_Ll], power_Ll.power_array[mask_mu_Ll]
        )
        power_Ss_interp = interp1d(
            power_Ss.k_array[0][mask_mu_Ss], power_Ss.power_array[mask_mu_Ss]
        )
        power_Sl_interp = interp1d(
            power_Sl.k_array[0][mask_mu_Sl], power_Sl.power_array[mask_mu_Sl]
        )

        ## low k:  k <= kmin_S
        mask_k_Ll = power_Ll.k_array[0][mask_mu_Ll] <= kmin_S
        power_mu.append(
            power_Ll.power_array[mask_mu_Ll][mask_k_Ll]
            * (power_Sl_interp(kmin_S) / power_Ss_interp(kmin_S))
        )
        k_array_mu.append(power_Ll.k_array[0][mask_mu_Ll][mask_k_Ll])
        if power_pwk_Ll is not None:
            k_array_pwk_mu.append(power_pwk_Ll.k_array[0][mask_mu_Ll][mask_k_Ll])
            mu_array_pwk_mu.append(power_pwk_Ll.k_array[1][mask_mu_Ll][mask_k_Ll])

        if error is not None:
            error_mu.append(power_Ll.error_array[mask_mu_Ll][mask_k_Ll])

        ## mid k:  kmin_S < k <= kNyqL / 4
        mask_k_Ll = (power_Ll.k_array[0][mask_mu_Ll] > kmin_S) & (
            power_Ll.k_array[0][mask_mu_Ll] <= k_max
        )
        k = power_Ll.k_array[0][mask_mu_Ll][mask_k_Ll]
        power_mu.append(
            power_Ll.power_array[mask_mu_Ll][mask_k_Ll]
            * (power_Sl_interp(k) / power_Ss_interp(k))
        )
        k_array_mu.append(k)
        if power_pwk_Ll is not None:
            k_array_pwk_mu.append(power_pwk_Ll.k_array[0][mask_mu_Ll][mask_k_Ll])
            mu_array_pwk_mu.append(power_pwk_Ll.k_array[1][mask_mu_Ll][mask_k_Ll])

        if error is not None:
            error_mu.append(power_Ll.error_array[mask_mu_Ll][mask_k_Ll])

        ## large k:  k > kNyqL / 4
        mask_k_Sl = power_Sl.k_array[0][mask_mu_Sl] > k_max
        power_mu.append(
            power_Sl.power_array[mask_mu_Sl][mask_k_Sl]
            * (power_Ll_interp(k_max) / power_Ss_interp(k_max))
        )
        k_array_mu.append(power_Sl.k_array[0][mask_mu_Sl][mask_k_Sl])
        if power_pwk_Ll is not None:
            k_array_pwk_mu.append(power_pwk_Sl.k_array[0][mask_mu_Sl][mask_k_Sl])
            mu_array_pwk_mu.append(power_pwk_Sl.k_array[1][mask_mu_Sl][mask_k_Sl])

        if error is not None:
            error_mu.append(power_Sl.error_array[mask_mu_Sl][mask_k_Sl])
            error_mu = np.concatenate(error_mu, axis=0)
        power_mu = np.concatenate(power_mu, axis=0)
        k_array_mu = np.concatenate(k_array_mu, axis=0)
        power.append(power_mu)
        k_array.append(k_array_mu)
        if power_pwk_Ll is not None:
            k_array_pwk_mu = np.concatenate(k_array_pwk_mu, axis=0)
            mu_array_pwk_mu = np.concatenate(mu_array_pwk_mu, axis=0)
            k_array_pwk.append(k_array_pwk_mu)
            mu_array_pwk.append(mu_array_pwk_mu)
        if error is not None:
            error.append(error_mu)

    power_spliced = np.array(
        [power[i][j] for j in range(len(power[i])) for i in range(len(mu_value))]
    )
    if power_pwk_Ll is not None:
        k_spliced = np.transpose(
            np.array(
                [
                    [k_array_pwk[i][j], mu_array_pwk[i][j]]
                    for j in range(len(power[i]))
                    for i in range(len(mu_value))
                ]
            )
        )
    else:
        k_spliced = np.transpose(
            np.array(
                [
                    [k_array[i][j], mu_value[i]]
                    for j in range(len(power[i]))
                    for i in range(len(mu_value))
                ]
            )
        )
    error_spliced = None
    if error is not None:
        error_spliced = np.array(
            [error[i][j] for j in range(len(error[i])) for i in range(len(mu_value))]
        )
    power_spectrum = FluxPowerSpectrum(
        "3D",
        k_array=k_spliced,
        power_array=power_spliced,
        error_array=error_spliced,
        file_init=None,
        size_box=None,
        h_normalized=True,
    )
    return power_spectrum, k_max, kmin_S


def verif_slicing(power_verif, power_spliced, mu_bins, name_out, style=None):
    """Plot the fractional residual between a spliced spectrum and a reference.

    For each mu bin in ``mu_bins``, interpolates ``power_spliced`` onto
    ``power_verif``'s k grid and plots the fractional difference
    ``(power_verif - power_spliced) / power_verif`` vs k (semilog-x),
    shading a +/-5% band, marking k=8, and saving the figure.

    Args:
        power_verif: Reference (e.g. fully-resolved) FluxPowerSpectrum
            to validate the splicing against.
        power_spliced: Spliced FluxPowerSpectrum (e.g. output of
            :func:`splice_3D`) to compare to the reference.
        mu_bins: Sequence of (up to 4) mu bin values to plot, matched
            against ``k_array[1]`` of both spectra.
        name_out: Output file path for the saved figure.
        style: Optional matplotlib style name passed to
            ``plt.style.use``. Defaults to None.
    """
    if style is not None:
        plt.style.use(style)
    plt.figure(figsize=(9, 6))
    size = 17
    size_ticks = 14
    for i in range(len(mu_bins)):
        mu = mu_bins[i]
        mask = power_verif.k_array[1] == mu
        mask_comparison = power_spliced.k_array[1] == mu
        power_array_comparison = interp1d(
            power_spliced.k_array[0][mask_comparison],
            power_spliced.power_array[mask_comparison],
            bounds_error=False,
            fill_value=np.NaN,
        )(power_verif.k_array[0][mask])

        delta_P = (
            power_verif.power_array[mask] - power_array_comparison
        ) / power_verif.power_array[mask]
        plt.semilogx(power_verif.k_array[0][mask], delta_P)
    plt.xlabel(r"$k$" + r" $[h$" + r"$\cdot$" + r"$\mathrm{Mpc}^{-1}]$", fontsize=size)
    plt.ylabel(
        r"$\left[P_{\alpha}(\mathrm{true}) - P_{\alpha}(\mathrm{splice})\right] / P_{\alpha}(\mathrm{true})}$",
        fontsize=size,
    )
    plt.fill_between(power_verif.k_array[0][mask], -0.05, 0.05, alpha=0.1, color="k")
    plt.gca().margins(x=0)
    plt.gca().tick_params(axis="x", labelsize=size_ticks)
    plt.gca().tick_params(axis="y", labelsize=size_ticks)
    plt.ylim([-0.15, 0.15])
    plt.plot([8, 8], [-0.15, 0.15], "k-")
    legend = [
        r"0.0 < $|\mu|$ < 0.25",
        r"0.25 < $|\mu|$ < 0.5",
        r"0.5 < $|\mu|$ < 0.75",
        r"0.75 < $|\mu|$ < 1.0",
    ]
    plt.legend(legend, fontsize=size)
    plt.tight_layout()
    power_spliced.save_plot(name_out)
