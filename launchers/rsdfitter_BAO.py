import glob,os
from lyapower import fitter, CLASS

### Grids

mu_binning = ""
direction = "y"
redshift = 2.0

str_box = f"{redshift}_gimlet_michael_rescaled_model0_fbar*_T0_box_gamma_box_{direction}_flux_pkmu{mu_binning}.txt"

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

file_grid1 = glob.glob(os.path.join(grid1,f"*/plt*_z{str_box}"))[0]
file_grid2 = glob.glob(os.path.join(grid2,f"*/plt*_z{str_box}"))[0]
file_grid3 = glob.glob(os.path.join(grid3,f"*/plt*_z{str_box}"))[0]
file_grid4 = glob.glob(os.path.join(grid4,f"*/plt*_z{str_box}"))[0]
file_grid5 = glob.glob(os.path.join(grid5,f"*/plt*_z{str_box}"))[0]
file_grid6 = glob.glob(os.path.join(grid6,f"*/plt*_z{str_box}"))[0]
file_grid7 = glob.glob(os.path.join(grid7,f"*/plt*_z{str_box}"))[0]
file_grid8 = glob.glob(os.path.join(grid8,f"*/plt*_z{str_box}"))[0]
file_grid9 = glob.glob(os.path.join(grid9,f"*/plt*_z{str_box}"))[0]

redshift_pm = 48.0

str_box_pm = f"{redshift_pm}_gimlet_michael_rho_matter_ps3d.txt"


file_pm_grid1 = glob.glob(os.path.join(grid1,f"*/plt*_z{str_box_pm}"))[0]
file_pm_grid2 = glob.glob(os.path.join(grid2,f"*/plt*_z{str_box_pm}"))[0]
file_pm_grid3 = glob.glob(os.path.join(grid3,f"*/plt*_z{str_box_pm}"))[0]
file_pm_grid4 = glob.glob(os.path.join(grid4,f"*/plt*_z{str_box_pm}"))[0]
file_pm_grid5 = glob.glob(os.path.join(grid5,f"*/plt*_z{str_box_pm}"))[0]
file_pm_grid6 = glob.glob(os.path.join(grid6,f"*/plt*_z{str_box_pm}"))[0]
file_pm_grid7 = glob.glob(os.path.join(grid7,f"*/plt*_z{str_box_pm}"))[0]
file_pm_grid8 = glob.glob(os.path.join(grid8,f"*/plt*_z{str_box_pm}"))[0]
file_pm_grid9 = glob.glob(os.path.join(grid9,f"*/plt*_z{str_box_pm}"))[0]

power_weighted = False
error_estimator="uncorrelated"
epsilon = 0.0



### Class options


class_settings = {'h' : 0.6732117,
        'Omega_b' : 0.04938682464547351,
        'Omega_cdm' : 0.2650127795050758,
        'k_pivot' : 0.05,
        'A_s' :  0.2100549e-08,
        'n_s' : 0.9660499,
}

monofonic_default = {
        "z_max_pk" : 199,
        "P_k_max_h/Mpc" : 222.88,
        "output" : 'dTk,vTk,mPk',
        "extra metric transfer functions" : 'yes',
        "gauge" : 'synchronous',
        "Omega_k" : 0,
        "Omega_fld" : 0,
        "Omega_scf" : 0,
        "N_eff" : 3.046,
        "N_ncdm" : 0,
        "P_k_ini type" : 'analytic_Pk',
        "alpha_s" : 0,
        "T_cmb" : 2.7255,
        "YHe" : 0.248,
        "reio_parametrization" : 'reio_none',
        "k_per_decade_for_pk" : 100,
        "k_per_decade_for_bao" : 100,
        "compute damping scale" : 'yes',
        "tol_perturb_integration" : 1e-08,
        "tol_background_integration" : 1e-09,
        "hyper_flat_approximation_nu" : 7000,
        "transfer_neglect_delta_k_S_t0" : 0.17,
        "transfer_neglect_delta_k_S_t1" : 0.05,
        "transfer_neglect_delta_k_S_t2" : 0.17,
        "transfer_neglect_delta_k_S_e" : 0.13,
        "delta_l_max" : 1000,
        "background_verbose" : 0,
        "thermodynamics_verbose" : 0,
        "perturbations_verbose" : 0,
        "transfer_verbose" : 0,
        "primordial_verbose" : 0,
        "spectra_verbose" : 0,
        "nonlinear_verbose" : 0,
        "lensing_verbose" : 0,
        "output_verbose" : 0,
        'format' : 'class',
        }


class_settings.update(monofonic_default)





import numpy as np

if __name__ == "__main__":

    from cosmoprimo import Cosmology,Fourier,PowerSpectrumBAOFilter


    power_1 = fitter.read_pfkmu(file_grid1,
                                 power_weighted=power_weighted,
                                 error_estimator=error_estimator,
                                 epsilon=epsilon)
    cosmoprimo = CLASS.CosmoprimoInterface(os.getcwd(),class_settings)

    mask_mu = power_1.k_array[1] == 0.0
    k_array = power_1.k_array[0][mask_mu]

    k_inf = k = np.logspace(-3,2,1000)

    # k_array = k_inf
    pk_no_wiggle = cosmoprimo.Pl_class_cosmoprimo_no_wiggle(k_array,redshift)
    pk_wiggle = cosmoprimo.Pl_class_cosmoprimo(k_array,redshift)

    import matplotlib.pyplot as plt

    plt.loglog(k_array,power_1.power_array[mask_mu])
    plt.loglog(k_array,pk_no_wiggle)


    # plt.loglog(k_array,power_1.power_array[mask_mu]/pk_no_wiggle)
    # plt.loglog(k_array,power_1.power_array[mask_mu]/pk_wiggle)

    # plt.loglog(k_array,pk_wiggle/pk_no_wiggle)
    #
    # plt.loglog(k_array,pk_no_wiggle)
    # plt.loglog(k_array,pk_wiggle)
