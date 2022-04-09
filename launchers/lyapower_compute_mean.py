from lyapower import power_spectra
import glob,os


mu_binning = ""
redshift = "2.0"


grid1 = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/2560_500mpc/"
grid2 = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/2560_250mpc/"
grid3 = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/2560_125mpc/"
grid4 = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/1280_500mpc/"
grid5 = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/1280_250mpc/"
grid6 = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/1280_125mpc/"
grid7 = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/1280_62.5mpc/"
grid8 = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/1280_62.5mpc_michael_treecool/"
grid9 = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/640_125mpc/"

legend_grid1 = "2560_500mpc"
legend_grid2 = "2560_250mpc"
legend_grid3 = "2560_125mpc"
legend_grid4 = "1280_500mpc"
legend_grid5 = "1280_250mpc"
legend_grid6 = "1280_125mpc"
legend_grid7 = "1280_62.5mpc"
legend_grid8 = "1280_62.5mpc_michael_treecool"
legend_grid9 = "640_125mpc"


str_box = f"{redshift}_gimlet_michael_rescaled_model0_fbar*_T0_box_gamma_box_*_flux_pkmu{mu_binning}.txt"

file_grid1 = glob.glob(os.path.join(grid1,f"*/plt*_z{str_box}"))
file_grid2 = glob.glob(os.path.join(grid2,f"*/plt*_z{str_box}"))
file_grid3 = glob.glob(os.path.join(grid3,f"*/plt*_z{str_box}"))
file_grid4 = glob.glob(os.path.join(grid4,f"*/plt*_z{str_box}"))
file_grid5 = glob.glob(os.path.join(grid5,f"*/plt*_z{str_box}"))
file_grid6 = glob.glob(os.path.join(grid6,f"*/plt*_z{str_box}"))
file_grid7 = glob.glob(os.path.join(grid7,f"*/plt*_z{str_box}"))
file_grid8 = glob.glob(os.path.join(grid8,f"*/plt*_z{str_box}"))
file_grid9 = glob.glob(os.path.join(grid9,f"*/plt*_z{str_box}"))



grid_to_compute = ["1","2","3","4","5","6","7","8","9"]



if __name__ == "__main__":
    import numpy as np
    os.makedirs("plots",exist_ok=True)

    for i in range(len(grid_to_compute)):
        files = eval(f"file_grid{grid_to_compute[i]}")
        name_mean = files[0].split("_flux_pkmu.txt")[0][:-1] + "mean_flux_pkmu.txt"
        power_spectra.FluxPowerSpectrum.compute_mean_gimlet(files,
                                                            name_mean,
                                                            "txt",
                                                            kmu=True)
