#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from lyapower import fitter
import glob, os

main_path = "/global/cfs/cdirs/desi/users/ravouxco/accel2/solene_home"
plot_pm = False

mu_binning = ""
direction = "meandirection"
flux_normalized = True
redshifts = [5.0, 4.0, 3.0, 2.5, 2.0]

grid_to_plot = [1, 2, 3, 4, 5, 6]


str_box = f"{'_fB14' if flux_normalized else ''}{direction}_flux_pkmu{mu_binning}.txt"
str_box_m = f"rhom_ps3d{mu_binning}.txt"

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


power_weighted = False
error_estimator = "uncorrelated"
epsilon = 0.0


mu_bins = [0.0, 0.25, 0.5, 0.75]
legend = [
    r"0.0 < $|\mu|$ < 0.25",
    r"0.25 < $|\mu|$ < 0.5",
    r"0.5 < $|\mu|$ < 0.75",
    r"0.75 < $|\mu|$ < 1.0",
]

kwargs = {
    "x_min_lim": 0.003,
    "x_max_lim": 100,
    "y_max_lim": None,
    "y_min_lim": None,
    "k_multiplication": True,
    "y_label": r"$\Delta^2_{\alpha}(k,\mu) = \frac{k^{3} P_{\alpha}(k,\mu)}{2\pi^{2}}$",
    "y_unit": False,
    "style": "~/2_Software/desi_ec/style.mplstyle",
    "figsize": (8, 6),
    "fontsize": 20,
    "labelsize_x": 18,
    "labelsize_y": 18,
}


y_min_lim_dict = {
    5.0: 5 * 10 ** (-6),
    4.0: 10 ** (-6),
    3.0: 4 * 10 ** (-7),
    2.5: 2 * 10 ** (-7),
    2.0: 10 ** (-7),
}
y_max_lim_dict = {5.0: 1, 4.0: 0.4, 3.0: 0.2, 2.5: 0.1, 2.0: 0.03}

kwargs_m = {
    "x_min_lim": 0.003,
    "x_max_lim": 100,
    "k_multiplication": True,
    "y_label": r"$\Delta(k) = \frac{k^{3} P(k)}{2\pi^{2}}$",
    "y_unit": False,
    "figsize": (8, 6),
    "fontsize": 20,
    "labelsize_x": 18,
    "labelsize_y": 18,
}

if __name__ == "__main__":

    os.makedirs("plots", exist_ok=True)

    for redshift in redshifts:
        for i in range(len(grid_to_plot)):
            if redshift in eval(f"avail_redshift_grid{grid_to_plot[i]}"):
                kwargs.update(
                    {
                        "y_max_lim": y_max_lim_dict[redshift],
                        "y_min_lim": y_min_lim_dict[redshift],
                    }
                )
                try:
                    i_box = grid_to_plot[i]
                    locals()[f"file_grid{i_box}"] = glob.glob(
                        os.path.join(
                            eval(f"grid{i_box}"),
                            f"*/plt*{eval(f'dict_redshift_grid{i_box}')[redshift]}{str_box}",
                        )
                    )[0]
                    if plot_pm:
                        locals()[f"file_grid{i_box}_m"] = glob.glob(
                            os.path.join(
                                eval(f"grid{i_box}"),
                                f"*/plt*{eval(f'dict_redshift_grid{i_box}')[redshift]}{str_box_m}",
                            )
                        )[0]
                except:
                    print(f"grid {i_box} not found")
                    print(
                        f"""Was trying to find {os.path.join(eval(f"grid{i_box}"),f"*/plt*{eval(f'dict_redshift_grid{i_box}')[redshift]}{str_box}")}"""
                    )
                    print(
                        f"""and {os.path.join(eval(f"grid{i_box}"),f"*/plt*{eval(f'dict_redshift_grid{i_box}')[redshift]}{str_box_m}")}"""
                    )

                print(f'Plotting grid {eval(f"file_grid{grid_to_plot[i]}")}')
                power_f = fitter.read_pfkmu(
                    eval(f"file_grid{grid_to_plot[i]}"),
                    power_weighted=power_weighted,
                    error_estimator=error_estimator,
                    epsilon=epsilon,
                )
                power_f.open_plot()
                power_f.plot_2d_pk(mu_bins, legend=legend, **kwargs)
                power_f.save_plot(
                    os.path.join(
                        os.getcwd(),
                        "plots",
                        eval(f"legend_grid{grid_to_plot[i]}")
                        + f"_redshift{redshift}_direction{direction}.pdf",
                    )
                )

                if plot_pm:
                    power_m = fitter.read_pk(
                        eval(f"file_grid{grid_to_plot[i]}_m"),
                        power_weighted=power_weighted,
                        error_estimator=error_estimator,
                        epsilon=epsilon,
                    )
                    power_m.open_plot()
                    power_m.plot_1d_pk(**kwargs_m)
                    power_m.save_plot(
                        os.path.join(
                            os.getcwd(),
                            "plots",
                            eval(f"legend_grid{grid_to_plot[i]}")
                            + f"_rhom_redshift{redshift}.pdf",
                        )
                    )
