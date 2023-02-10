from lyapower import power_spectra
import glob, os

main_path = "/global/cfs/cdirs/desi/users/ravouxco/accel2/solene_home"

mu_binning = ""
flux_normalized = True

str_box = f"{'_fB14' if flux_normalized else ''}*_flux_pkmu{mu_binning}.txt"

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


grid_to_compute = [1, 2, 3, 4, 5, 6]
redshifts = [5.0, 4.0, 3.0, 2.5, 2.0]


if __name__ == "__main__":

    for redshift in redshifts:
        for i in range(len(grid_to_compute)):
            if redshift in eval(f"avail_redshift_grid{grid_to_compute[i]}"):
                try:
                    i_box = grid_to_compute[i]
                    locals()[f"file_grid{i_box}"] = glob.glob(
                        os.path.join(
                            eval(f"grid{i_box}"),
                            f"*/plt*{eval(f'dict_redshift_grid{i_box}')[redshift]}{str_box}",
                        )
                    )
                except:
                    print(f"grid {i_box} not found")
                    print(
                        f"""Was trying to find {os.path.join(eval(f"grid{i_box}"),f"*/plt*{eval(f'dict_redshift_grid{i_box}')[redshift]}{str_box}")}"""
                    )

                files = eval(f"file_grid{grid_to_compute[i]}")
                name_mean = (
                    files[0].split("_flux_pkmu.txt")[0][:-1]
                    + "meandirection_flux_pkmu.txt"
                )
                power_spectra.FluxPowerSpectrum.compute_mean_gimlet(
                    files, name_mean, "txt", kmu=True
                )
