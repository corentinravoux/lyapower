#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from lyapower import fitter


name_p3d_file = "/global/cfs/cdirs/desi/users/ravouxco/accel2/solene_home/6144_600/gimlet_nyx/plt00713/plt00713_fB14meandirection_flux_pkmu.txt"
redshift = 2.0


####### Linear power spectra param

class_settings = {
    "h": 0.675,
    "Omega_b": 0.0487,
    "Omega_cdm": 0.2613,
    "k_pivot": 0.05,
    "sigma8": 0.83,
    "n_s": 0.96,
    "output": "dTk,mPk",
    "P_k_max_h/Mpc": 200,
    "z_max_pk": 10,
}


####### Minuit parameters

non_linear_model = "1_bao"

non_linear_params = {
    "q_1": 0,
    "q_2": 0,
    "k_v": 1,
    "a_v": 1,
    "b_v": 1,
    "k_p": 10,
    " S_p": 6,
    "S_t": 6,
}

non_linear_limits = [
    ("q_1", (0, 10)),
    ("q_2", (-10, 10)),
    ("k_v", (10**-2, 10**3)),
    ("a_v", (0.1, 2)),
    ("b_v", (0.1, 2)),
    ("k_p", (10**-2, 10**3)),
    ("S_p", (0, 100)),
    ("S_t", (0, 100)),
]

linear_params = {"b": 0.1, "beta": 1}

linear_limits = [("b", (10**-3, 1)), ("beta", (0.1, 10))]


minuit_parameters = {}
minuit_limits = []

minuit_parameters.update(linear_params)
minuit_limits = minuit_limits + linear_limits

minuit_parameters.update(non_linear_params)
minuit_limits = minuit_limits + non_linear_limits


### Fitter params


kmin, kmax = None, 10
cost_name = "least"
ncall = 10000
epsilon = 0.05
rebin = {
    "nb_bin": 20,
    "loglin": True,
    "k_loglin": 0.09,
}
filename_pk = "class"
fix_args = None
error_estimator = "uncorrelated"
name_pm_file = None
z_init = None


####### Plot parameters

name_out = "p3d_fit"

mu_bins = [0.0, 0.25, 0.5, 0.75]

mu_bins_legend = [
    r"0.0 < $|\mu|$ < 0.25",
    r"0.25 < $|\mu|$ < 0.5",
    r"0.5 < $|\mu|$ < 0.75",
    r"0.75 < $|\mu|$ < 1.0",
]

plt_args = {
    "style": "~/2_Software/desi_ec/style.mplstyle",
    "figsize": (8, 6),
    "fontsize": 16,
    "labelsize_x": 14,
    "labelsize_y": 14,
    "y_label": r"$P_{\alpha}(k,\mu)/P_{\mathrm{m}}(k)$",
    "y_unit": False,
}


(
    minuit,
    power_f,
    power_m,
    linear_power_spectrum,
    non_linear_model,
) = fitter.fitter_k_mu(
    name_p3d_file,
    filename_pk,
    minuit_parameters,
    minuit_limits,
    class_dict=class_settings,
    z_simu=redshift,
    non_linear_model=non_linear_model,
    cost_name=cost_name,
    ncall=ncall,
    kmin=kmin,
    kmax=kmax,
    error_estimator=error_estimator,
    epsilon=epsilon,
    rebin=rebin,
    name_pm_file=name_pm_file,
    z_init=z_init,
    fix_args=fix_args,
)

fitter.plot_fit(
    minuit,
    power_f,
    power_m,
    linear_power_spectrum,
    non_linear_model,
    mu_bins,
    mu_bins_legend,
    name_out=name_out,
    **plt_args,
)
