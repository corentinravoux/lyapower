#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from lyapower import fitter
import glob, os
from lyapower import power_spectra

### Grids


main_path = "/global/cfs/cdirs/desi/users/ravouxco/accel2/solene_home"

mu_binning = ""
directions = ["x", "y", "z", "meandirection"]
flux_normalized = True
redshifts = [2.0, 2.5, 3.0, 4.0, 5.0]


splice_pm = False
pm_number = "002"


plot_splice = False
plot_splice_pm = False
verification_splice = False


grid1 = os.path.join(main_path, "1536_150/gimlet_nyx")
grid2 = os.path.join(main_path, "3072_150/gimlet_nyx")
grid3 = os.path.join(main_path, "3072_300/gimlet_nyx")
grid4 = os.path.join(main_path, "6144_150/gimlet_nyx")
grid5 = os.path.join(main_path, "6144_300/gimlet_nyx")
grid6 = os.path.join(main_path, "6144_600/gimlet_nyx")


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

legend_grid1 = "1536_150mpc"
legend_grid2 = "3072_150mpc"
legend_grid3 = "3072_300mpc"
legend_grid4 = "6144_150mpc"
legend_grid5 = "6144_300mpc"
legend_grid6 = "6144_600mpc"


### Splicing options
# Hr_S = grid_to_plot[0]
# Lr_L = grid_to_plot[1]
# Lr_S = grid_to_plot[2]


name_splicing = "splicing_A"
size_small = 300
size_large = 600
N_small = 3072
N_large = 6144
grid_to_plot = [5, 6, 3]
grid_verification = None


name_splicing = "splicing_B"
size_small = 150
size_large = 300
N_small = 1536
N_large = 3072
grid_to_plot = [2, 3, 1]
grid_verification = 5

name_splicing = "splicing_D"
size_small = 150
size_large = 600
N_small = 1536
N_large = 6144
grid_to_plot = [4, 6, 1]
grid_verification = None


name_splicing = "splicing_F"
size_small = 150
size_large = 300
N_small = 1536
N_large = 6144
grid_to_plot = [4, 5, 2]
grid_verification = None


# Others - do not work, three different N

# name_splicing = "splicing_C"
# size_small = 150
# size_large = 600
# N_small = 1536
# N_large = 6144
# grid_to_plot = [2,6,1]
# grid_verification = None


# name_splicing = "splicing_E"
# size_small = 150
# size_large = 300
# N_small = 1536
# N_large = 6144
# grid_to_plot = [4,3,1]
# grid_verification = None


use_nyquist = False
impose_kmin_coeff = 5

power_weighted = False
error_estimator = "uncorrelated"
epsilon = 0.0

### Plot options

mu_bins = [[0.0], [0.25], [0.5], [0.75]]


kwargs = {
    "x_min_lim": 0.001,
    "x_max_lim": 100,
    "y_max_lim": 0.1,
    "y_min_lim": 10**-7,
    "linestyle": ["-", "--"],
    "k_multiplication": True,
    "y_label": r"$\Delta(k) = \frac{k^{3} P(k)}{2\pi^{2}}$",
    "y_unit": False,
    "y_label_comparison": r"$\frac{\Delta(Spliced) - \Delta(i)}{\Delta(Spliced)}$",
    "y_unit_comparison": False,
    "y_max_lim_comparison": 0.20,
    "y_min_lim_comparison": -0.20,
    "error_bar_comparison": False,
}
kwargs_pm = {
    "y_max_lim_comparison": 0.05,
    "y_min_lim_comparison": -0.05,
    "error_bar_comparison": False,
}

kwargs2 = {
    "k_multiplication": True,
    "y_label": r"$\Delta(k) = \frac{k^{3} P(k)}{2\pi^{2}}$",
    "y_unit": False,
}


str_plot = "splicing_plots"
str_splice = "splices"


if __name__ == "__main__":
    for direction in directions:
        str_box = (
            f"{'_fB14' if flux_normalized else ''}{direction}_flux_pkmu{mu_binning}.txt"
        )
        str_box_m = f"rhom_ps3d{mu_binning}.txt"

        for redshift in redshifts:
            if (
                (redshift in eval(f"avail_redshift_grid{grid_to_plot[0]}"))
                & (redshift in eval(f"avail_redshift_grid{grid_to_plot[1]}"))
                & (redshift in eval(f"avail_redshift_grid{grid_to_plot[2]}"))
            ):
                name_out_splice = (
                    f"{name_splicing}_z{redshift}_{direction}_spliced_spectra.pdf"
                )
                path_plots = os.path.join(os.getcwd(), str_plot)

                if plot_splice | plot_splice_pm | verification_splice:
                    os.makedirs(path_plots, exist_ok=True)
                    name_out_splice = os.path.join(path_plots, name_out_splice)

                path_splices = os.path.join(os.getcwd(), str_splice)
                dir_out_splice = os.path.join(path_splices, name_splicing)
                os.makedirs(dir_out_splice, exist_ok=True)

                for i_box in grid_to_plot:
                    try:
                        locals()[f"file_grid{i_box}"] = glob.glob(
                            os.path.join(
                                eval(f"grid{i_box}"),
                                f"*/plt*{eval(f'dict_redshift_grid{i_box}')[redshift]}{str_box}",
                            )
                        )[0]
                        if splice_pm:
                            locals()[f"file_pm_grid{i_box}"] = glob.glob(
                                os.path.join(
                                    eval(f"grid{i_box}"),
                                    f"*/plt*{pm_number}{str_box_m}",
                                )
                            )[0]
                    except:
                        print(f"grid {i_box} not found")
                        print(
                            f"""Was trying to find {os.path.join(eval(f"grid{i_box}"),f"*/plt*{eval(f'dict_redshift_grid{i_box}')[redshift]}{str_box}")}"""
                        )
                        print(
                            f"""and {os.path.join(eval(f"grid{i_box}"),f"*/plt*{pm_number}{str_box_m}")}"""
                        )

                if verification_splice:
                    try:
                        locals()[f"file_grid{grid_verification}"] = glob.glob(
                            os.path.join(
                                eval(f"grid{grid_verification}"),
                                f"*/plt*{eval(f'dict_redshift_grid{grid_verification}')[redshift]}{str_box}",
                            )
                        )[0]
                        if splice_pm:
                            locals()[f"file_pm_grid{grid_verification}"] = glob.glob(
                                os.path.join(
                                    eval(f"grid{grid_verification}"),
                                    f"*/plt*{pm_number}{str_box_m}",
                                )
                            )[0]
                    except:
                        print(f"grid {grid_verification} not found")
                        print(
                            f"""Was trying to find {os.path.join(eval(f"grid{grid_verification}"),f"*/plt*{eval(f'dict_redshift_grid{grid_verification}')[redshift]}{str_box}")}"""
                        )
                        print(
                            f"""and {os.path.join(eval(f"grid{grid_verification}"),f"*/plt*{pm_number}{str_box_m}")}"""
                        )

                grid_Hr_S = eval(f"file_grid{grid_to_plot[0]}")
                grid_Sr_L = eval(f"file_grid{grid_to_plot[1]}")
                grid_Sr_S = eval(f"file_grid{grid_to_plot[2]}")

                power_Sr_L = fitter.read_pfkmu(
                    grid_Sr_L,
                    power_weighted=power_weighted,
                    error_estimator=error_estimator,
                    epsilon=epsilon,
                )
                power_Hr_S = fitter.read_pfkmu(
                    grid_Hr_S,
                    power_weighted=power_weighted,
                    error_estimator=error_estimator,
                    epsilon=epsilon,
                )
                power_Sr_S = fitter.read_pfkmu(
                    grid_Sr_S,
                    power_weighted=power_weighted,
                    error_estimator=error_estimator,
                    epsilon=epsilon,
                )

                power_spliced, k_max, kmin_S = power_spectra.splice_3D(
                    power_Sr_L,
                    power_Hr_S,
                    power_Sr_S,
                    size_small,
                    size_large,
                    N_large,
                    use_nyquist=use_nyquist,
                    impose_kmin_coeff=impose_kmin_coeff,
                )
                print(
                    f"For direction {direction}, redshift {redshift}, kmax = {k_max}, kmin = {kmin_S}"
                )

                splice_file_out = os.path.join(
                    dir_out_splice,
                    f"{name_splicing}_spliced_spectra_{direction}_z{redshift}.txt",
                )
                power_spliced.write_to_gimlet(splice_file_out, power_weighted=False)

                if splice_pm:
                    grid_m_Hr_S = eval(f"file_pm_grid{grid_to_plot[0]}")
                    grid_m_Sr_L = eval(f"file_pm_grid{grid_to_plot[1]}")
                    grid_m_Sr_S = eval(f"file_pm_grid{grid_to_plot[2]}")

                    power_m_Sr_L = fitter.read_pk(
                        grid_m_Sr_L,
                        power_weighted=power_weighted,
                        error_estimator=error_estimator,
                        epsilon=epsilon,
                    )
                    power_m_Hr_S = fitter.read_pk(
                        grid_m_Hr_S,
                        power_weighted=power_weighted,
                        error_estimator=error_estimator,
                        epsilon=epsilon,
                    )
                    power_m_Sr_S = fitter.read_pk(
                        grid_m_Sr_S,
                        power_weighted=power_weighted,
                        error_estimator=error_estimator,
                        epsilon=epsilon,
                    )

                    power_m_spliced = power_spectra.splice_1D(
                        power_m_Sr_L,
                        power_m_Hr_S,
                        power_m_Sr_S,
                        size_small,
                        size_large,
                        N_large,
                        use_nyquist=use_nyquist,
                        impose_kmin_coeff=impose_kmin_coeff,
                    )

                    splice_file_pm_out = os.path.join(
                        dir_out_splice,
                        f"{name_splicing}_spliced_matter_power_spectra.txt",
                    )
                    power_m_spliced.write_to_gimlet(
                        splice_file_pm_out, power_weighted=power_weighted
                    )

                if plot_splice:
                    for i in range(len(mu_bins)):
                        mu = mu_bins[i]
                        name_out = f"{name_splicing}_direction_{direction}_z{redshift}_mu{mu[0]}.pdf"

                        name_out = os.path.join(os.getcwd(), str_plot, name_out)

                        kwargs["comparison"] = power_spliced
                        legend = []
                        for j in range(len(grid_to_plot)):
                            color = f"C{j}"
                            legend.append(eval(f"legend_grid{grid_to_plot[j]}"))
                            power_f = fitter.read_pfkmu(
                                eval(f"file_grid{grid_to_plot[j]}"),
                                power_weighted=power_weighted,
                                error_estimator=error_estimator,
                                epsilon=epsilon,
                            )

                            if j == 0:
                                power_f.open_subplot(figsize=(8, 8))
                            power_f.plot_2d_pk(mu, color=[color, color], **kwargs)

                        legend.append("Spliced spectra")
                        kwargs["comparison"] = None

                        power_spliced.plot_2d_pk(
                            mu, color=[f"C{j+1}", f"C{j+1}"], legend=legend, **kwargs
                        )
                        power_spliced.save_plot(name_out)

                    l = [
                        r"0.0 < $|\mu|$ < 0.25",
                        r"0.25 < $|\mu|$ < 0.5",
                        r"0.5 < $|\mu|$ < 0.75",
                        r"0.75 < $|\mu|$ < 1.0",
                    ]

                    power_spliced.open_plot()
                    power_spliced.plot_2d_pk(
                        [0.0, 0.25, 0.5, 0.75], legend=l, **kwargs2
                    )
                    power_spliced.save_plot(name_out_splice)

                if plot_splice_pm:
                    legend = []

                    kwargs_pm["comparison"] = power_m_spliced

                    for j in range(len(grid_to_plot)):
                        color = f"C{j}"
                        legend.append(eval(f"legend_grid{grid_to_plot[j]}"))
                        power_m = fitter.read_pk(
                            eval(f"file_pm_grid{grid_to_plot[j]}"),
                            power_weighted=power_weighted,
                            error_estimator=error_estimator,
                            epsilon=epsilon,
                        )

                        if j == 0:
                            power_m.open_subplot(figsize=(8, 8))
                        power_m.plot_1d_pk(color=color, **kwargs_pm)

                    legend.append("Spliced spectra")
                    kwargs_pm["comparison"] = None

                    name_out_pm = (
                        f"{name_splicing}_direction_{direction}_pm_z{redshift}.pdf"
                    )

                    name_out_pm = os.path.join(os.getcwd(), str_plot, name_out_pm)

                    power_m_spliced.plot_1d_pk(
                        color=f"C{j+1}", legend=legend, **kwargs_pm
                    )
                    power_m_spliced.save_plot(name_out_pm)

                if verification_splice:
                    for i in range(len(mu_bins)):
                        mu = mu_bins[i]
                        name_out = f"{name_splicing}_direction{direction}_z{redshift}_mu{mu[0]}_verification.pdf"

                        name_out = os.path.join(os.getcwd(), str_plot, name_out)

                        kwargs["comparison"] = power_spliced
                        legend = []
                        color = f"C{0}"
                        legend.append(eval(f"legend_grid{grid_verification}"))
                        power_f = fitter.read_pfkmu(
                            eval(f"file_grid{grid_verification}"),
                            power_weighted=power_weighted,
                            error_estimator=error_estimator,
                            epsilon=epsilon,
                        )

                        power_f.open_subplot(figsize=(8, 8))
                        power_f.plot_2d_pk(mu, color=["C0", "C0"], **kwargs)

                        legend.append("Spliced spectra")
                        kwargs["comparison"] = None

                        power_spliced.plot_2d_pk(
                            mu, color=["C1", "C1"], legend=legend, **kwargs
                        )
                        power_spliced.save_plot(name_out)

                    power_f = fitter.read_pfkmu(
                        eval(f"file_grid{grid_verification}"),
                        power_weighted=power_weighted,
                        error_estimator=error_estimator,
                        epsilon=epsilon,
                    )
                    name_out = f"{name_splicing}_direction{direction}_z{redshift}_verification.pdf"
                    name_out = os.path.join(os.getcwd(), str_plot, name_out)
                    power_spectra.verif_slicing(
                        power_f, power_spliced, mu_bins, name_out, style=None
                    )
