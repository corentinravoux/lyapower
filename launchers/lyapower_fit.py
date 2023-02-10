#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from lyapower import fitter
import glob, os, pickle


models = [
    "1_class",
    "1_class_fixedq2",
    "1_pmnorm",
    "1_pmnorm_fixedq2",
    "0_class",
    "0_pmnorm",
]

main_path = "/global/cfs/cdirs/desi/users/ravouxco/accel2/solene_home"
path_splice = "/global/cfs/cdirs/desi/users/ravouxco/accel2/p3d_plots/splices"

flux_normalized = True
mu_binning = ""

directions = ["x", "y", "z", "meandirection"]

redshifts = [5.0, 4.0, 3.0, 2.5, 2.0]

grid1 = os.path.join(main_path, "1536_150/gimlet_nyx")
grid2 = os.path.join(main_path, "3072_150/gimlet_nyx")
grid3 = os.path.join(main_path, "3072_300/gimlet_nyx")
grid4 = os.path.join(main_path, "6144_150/gimlet_nyx")
grid5 = os.path.join(main_path, "6144_300/gimlet_nyx")
grid6 = os.path.join(main_path, "6144_600/gimlet_nyx")

gridA = os.path.join(path_splice, "splicing_A")
gridB = os.path.join(path_splice, "splicing_B")
gridD = os.path.join(path_splice, "splicing_D")
gridF = os.path.join(path_splice, "splicing_F")


dict_redshift_grid1 = {5.0: "366", 4.0: "402", 3.0: "461", 2.5: "510", 2.0: "580"}
dict_redshift_grid2 = {5.0: "484", 4.0: "579", 3.0: "723", 2.5: "839", 2.0: "1006"}
dict_redshift_grid3 = {5.0: "389", 4.0: "439", 3.0: "526", 2.5: "586", 2.0: "666"}
dict_redshift_grid4 = {5.0: "781", 4.0: "1032", 3.0: "1424"}
dict_redshift_grid5 = {5.0: "522", 4.0: "657", 3.0: "861", 2.5: "1023", 2.0: "1226"}
dict_redshift_grid6 = {5.0: "400", 4.0: "458", 3.0: "549", 2.5: "617", 2.0: "713"}

avail_redshift_grid1 = [5.0, 4.0, 3.0, 2.5, 2.0]
avail_redshift_grid2 = [5.0, 4.0, 3.0, 2.5, 2.0]
avail_redshift_grid3 = [5.0, 4.0, 3.0, 2.5, 2.0]
avail_redshift_grid4 = [5.0, 4.0, 3.0]
avail_redshift_grid5 = [5.0, 4.0, 3.0, 2.5, 2.0]
avail_redshift_grid6 = [5.0, 4.0, 3.0, 2.5, 2.0]
avail_redshift_gridA = [5.0, 4.0, 3.0, 2.5, 2.0]
avail_redshift_gridB = [5.0, 4.0, 3.0, 2.5, 2.0]
avail_redshift_gridD = [5.0, 4.0, 3.0]
avail_redshift_gridF = [5.0, 4.0, 3.0]

legend_grid1 = "1536_150mpc"
legend_grid2 = "3072_150mpc"
legend_grid3 = "3072_300mpc"
legend_grid4 = "6144_150mpc"
legend_grid5 = "6144_300mpc"
legend_grid6 = "6144_600mpc"
legend_gridA = "splicing_A"
legend_gridB = "splicing_B"
legend_gridD = "splicing_D"
legend_gridF = "splicing_F"

mu_bins = [0.0, 0.25, 0.5, 0.75]

mu_bins_legend = [
    r"0.0 < $|\mu|$ < 0.25",
    r"0.25 < $|\mu|$ < 0.5",
    r"0.5 < $|\mu|$ < 0.75",
    r"0.75 < $|\mu|$ < 1.0",
]

compute_kna = False
tol_kna = 0.0001

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
}


####### Minuit parameters


non_linear_params_0 = {
    "k_nl": 1,
    "a_nl": 1,
    "k_p": 1,
    "a_p": 1,
    "k_v0": 1,
    "a_v0": 1,
    "k_v1": 1,
    "a_v1": 1,
}
non_linear_limits_0 = [
    ("k_nl", (10**-1, 10**2)),
    ("a_nl", (0.001, 4)),
    ("k_p", (10**-1, 10**2)),
    ("a_p", (0.001, 4)),
    ("k_v0", (10**-1, 10**2)),
    ("a_v0", (0.001, 4)),
    ("k_v1", (10**-1, 10**2)),
    ("a_v1", (0.001, 4)),
]


non_linear_params_1 = {"q_1": 0, "q_2": 0, "k_v": 1, "a_v": 1, "b_v": 1, "k_p": 10}

# Old parameters used for JZ
# non_linear_limits_1 = [
#     ("q_1", (0, 2)),
#     ("q_2", (-2, 2)),
#     ("k_v", (10**-8, 10**3)),
#     ("a_v", (-2, 2)),
#     ("b_v", (0, 2)),
#     ("k_p", (10**-1, 10**3)),
# ]


# Good convergence
non_linear_limits_1 = [
    ("q_1", (0, 10)),
    ("q_2", (-10, 10)),
    ("k_v", (10**-2, 10**3)),
    ("a_v", (0.1, 2)),
    ("b_v", (0.1, 2)),
    ("k_p", (10**-2, 10**3)),
]


linear_params = {"b": 0.1, "beta": 1}

linear_limits = [("b", (10**-3, 1)), ("beta", (0.1, 10))]


### Out params

base_name_out = "p3d_fit"

plt_args = {
    "style": "~/2_Software/desi_ec/style.mplstyle",
    "figsize": (8, 6),
    "fontsize": 16,
    "labelsize_x": 14,
    "labelsize_y": 14,
    "y_label": r"$P_{\alpha}(k,\mu)/P_{\mathrm{m}}(k)$",
    "y_unit": False,
}


### Fitter fixed params


kmin, kmax = None, 10

cost_name = "least"
ncall = 10000

power_weighted = False

epsilon = 0.05
rebin = {
    "nb_bin": 20,
    "loglin": True,
    "k_loglin": 0.09,
}


for model in models:

    ### Fitter model params

    if "pmnorm" in model:
        grid_to_plot = [1, 2, 3, 4, 5]
        filename_pk = "pmnorm"
        pm_number = "002"
        redshift_pm = 198.0
        class_settings.update({"z_max_pk": redshift_pm})

    else:
        grid_to_plot = [1, 2, 3, 4, 5, 6, "A", "B", "D", "F"]
        filename_pk = "class"
        pm_number = None
        redshift_pm = None
        class_settings.update({"z_max_pk": 10})

    minuit_parameters = {}
    minuit_parameters.update(linear_params)

    minuit_limits = []
    minuit_limits = minuit_limits + linear_limits

    if "0" in model:
        non_linear_model = "0"
        non_linear_params = non_linear_params_0
        non_linear_limits = non_linear_limits_0
        minuit_parameters.update(non_linear_params)
        minuit_limits = minuit_limits + non_linear_limits
        label_model = f"nonlinear{non_linear_model}"
    elif "1" in model:
        non_linear_model = "1"
        non_linear_params = non_linear_params_1
        non_linear_limits = non_linear_limits_1
        minuit_parameters.update(non_linear_params)
        minuit_limits = minuit_limits + non_linear_limits
        label_model = f"nonlinear{non_linear_model}"
    else:
        label_model = "linear"

    if "fixedq2" in model:
        fix_args = ["q_2"]
        add = "fixed_q2"
    else:
        fix_args = None
        add = ""

    path_out = os.path.join(
        os.getcwd(), f"fit_plots_model_{label_model}_{filename_pk}{add}"
    )
    os.makedirs(path_out, exist_ok=True)
    for direction in directions:
        for redshift in redshifts:
            for i in range(len(grid_to_plot)):
                i_box = grid_to_plot[i]
                if redshift in eval(f"avail_redshift_grid{i_box}"):
                    if type(i_box) != int:
                        error_estimator = "computed_epsilon"
                        str_box = f"""{eval(f'legend_grid{i_box}')}_spliced_spectra_{direction}_z{redshift}.txt"""
                        str_box_m = f"""{eval(f'legend_grid{i_box}')}_spliced_matter_power_spectra.txt"""
                    else:
                        error_estimator = "uncorrelated"
                        str_box = f"""*/plt*{eval(f'dict_redshift_grid{i_box}')[redshift]}{'_fB14' if flux_normalized else ''}{direction}_flux_pkmu{mu_binning}.txt"""
                        str_box_m = f"""*/plt*{pm_number}rhom_ps3d{mu_binning}.txt"""

                    name_box = os.path.join(eval(f"grid{i_box}"), str_box)
                    name_box_pm = os.path.join(eval(f"grid{i_box}"), str_box_m)

                    try:
                        locals()[f"file_grid{i_box}"] = glob.glob(name_box)[0]
                        if filename_pk == "pmnorm":
                            locals()[f"file_pm_grid{i_box}"] = glob.glob(name_box_pm)[0]
                    except:
                        print(f"grid {i_box} not found")
                        print(f"Was trying to find {name_box}")
                        if filename_pk == "pmnorm":
                            print(f"and {name_box_pm}")

                    if filename_pk == "class":
                        name_pm_file = None
                        z_init = None
                    else:
                        name_pm_file = eval(f"file_pm_grid{grid_to_plot[i]}")
                        z_init = redshift_pm

                    name_p3d_file = eval(f"file_grid{grid_to_plot[i]}")
                    legend_grid = eval(f"legend_grid{grid_to_plot[i]}")
                    name_out = os.path.join(
                        path_out,
                        f"{base_name_out}_{legend_grid}_z{redshift}_{label_model}_{direction}_{filename_pk}{mu_binning}",
                    )

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
                    minuit_params, minuit_errors = [], []
                    for i in range(len(minuit.parameters)):
                        minuit_params.append(minuit.params[minuit.parameters[i]].value)
                        minuit_errors.append(minuit.params[minuit.parameters[i]].error)
                    reduced_chi2 = minuit.fval / (
                        len(power_f.power_array) - len(minuit.values)
                    )

                    pickle.dump(
                        (
                            minuit_params,
                            minuit_errors,
                            power_f,
                            power_m,
                            linear_power_spectrum,
                            non_linear_model,
                            reduced_chi2,
                        ),
                        open(f"{name_out}_fit_result.pickle", "wb"),
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

                    if compute_kna:
                        kna = fitter.compute_kna(minuit, power_m, tol_kna)
                        print(r"$k_{na} =$ ", kna)
