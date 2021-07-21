from rsd_fitter import fitter
import glob
import os

### Box params


mu_binning = ""
direction = "y"
redshift = 2.0

str_box = f"{redshift}_gimlet_michael_rescaled_model0_fbar*_T0_box_gamma_box_{direction}_flux_pkmu{mu_binning}.txt"
str_spliced_box = f"_spliced_spectra_z{redshift}.txt"

grid1 = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/2560_500mpc/"
grid2 = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/2560_250mpc/"
grid3 = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/2560_125mpc/"
grid4 = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/1280_500mpc/"
grid5 = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/1280_250mpc/"
grid6 = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/1280_125mpc/"
grid7 = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/1280_62.5mpc/"
grid8 = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/1280_62.5mpc_michael_treecool/"
grid9 = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/640_125mpc/"
gridA = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/splicing_A/"
gridB = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/splicing_B/"
gridC = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid/splicing_C/"

legend_grid1 = "2560_500mpc"
legend_grid2 = "2560_250mpc"
legend_grid3 = "2560_125mpc"
legend_grid4 = "1280_500mpc"
legend_grid5 = "1280_250mpc"
legend_grid6 = "1280_125mpc"
legend_grid7 = "1280_62.5mpc"
legend_grid8 = "1280_62.5mpc_michael_treecool"
legend_grid9 = "640_125mpc"
legend_gridA = "splicing_A"
legend_gridB = "splicing_B"
legend_gridC = "splicing_C"

file_grid1 = glob.glob(os.path.join(grid1,f"*/plt*_z{str_box}"))[0]
file_grid2 = glob.glob(os.path.join(grid2,f"*/plt*_z{str_box}"))[0]
file_grid3 = glob.glob(os.path.join(grid3,f"*/plt*_z{str_box}"))[0]
file_grid4 = glob.glob(os.path.join(grid4,f"*/plt*_z{str_box}"))[0]
file_grid5 = glob.glob(os.path.join(grid5,f"*/plt*_z{str_box}"))[0]
file_grid6 = glob.glob(os.path.join(grid6,f"*/plt*_z{str_box}"))[0]
file_grid7 = glob.glob(os.path.join(grid7,f"*/plt*_z{str_box}"))[0]
file_grid8 = glob.glob(os.path.join(grid8,f"*/plt*_z{str_box}"))[0]
file_grid9 = glob.glob(os.path.join(grid9,f"*/plt*_z{str_box}"))[0]
file_gridA = glob.glob(os.path.join(gridA,f"splicing_A*{str_spliced_box}"))[0]
file_gridB = glob.glob(os.path.join(gridB,f"splicing_B*{str_spliced_box}"))[0]
file_gridC = glob.glob(os.path.join(gridC,f"splicing_C*{str_spliced_box}"))[0]



redshift_pm = 48.0

str_box_pm = f"{redshift_pm}_gimlet_michael_rho_matter_ps3d.txt"
str_spliced_box_pm = "spliced_matter_power_spectra.txt"


file_pm_grid1 = glob.glob(os.path.join(grid1,f"*/plt*_z{str_box_pm}"))[0]
file_pm_grid2 = glob.glob(os.path.join(grid2,f"*/plt*_z{str_box_pm}"))[0]
file_pm_grid3 = glob.glob(os.path.join(grid3,f"*/plt*_z{str_box_pm}"))[0]
file_pm_grid4 = glob.glob(os.path.join(grid4,f"*/plt*_z{str_box_pm}"))[0]
file_pm_grid5 = glob.glob(os.path.join(grid5,f"*/plt*_z{str_box_pm}"))[0]
file_pm_grid6 = glob.glob(os.path.join(grid6,f"*/plt*_z{str_box_pm}"))[0]
file_pm_grid7 = glob.glob(os.path.join(grid7,f"*/plt*_z{str_box_pm}"))[0]
file_pm_grid8 = glob.glob(os.path.join(grid8,f"*/plt*_z{str_box_pm}"))[0]
file_pm_grid9 = glob.glob(os.path.join(grid9,f"*/plt*_z{str_box_pm}"))[0]
file_pm_gridA = glob.glob(os.path.join(gridA,f"splicing_A_{str_spliced_box_pm}"))[0]
file_pm_gridB = glob.glob(os.path.join(gridB,f"splicing_B_{str_spliced_box_pm}"))[0]
file_pm_gridC = glob.glob(os.path.join(gridC,f"splicing_C_{str_spliced_box_pm}"))[0]


grid_michael = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/Michael_grid/all_model_outputs_copy.hdf5"
grid_michael_field_name = 'fiducial/redshift_2.2/rescale_Fbar_fiducial/3d power kmu/default binning'
legend_grid_michael = "4096_80mpc"


filename_pk = "pmnorm" # "class"



### Fitter params


non_linear_model="1"
kmin,kmax = None , 8

cost_name="least"
var = 1
sigma = 1
ncall = 10000

var_minos = None
power_weighted = False
error_estimator="computed" #"uncorrelated"
# error_estimator="uncorrelated"
epsilon = 0.05
rebin = {"nb_bin" : 20,
         "loglin" : True,
         "k_loglin" : 0.09}

fix_args = None # ["q_2"]
eps = 0.0001

mu_bins_legend = [r"0.0 < $|\mu|$ < 0.25",r"0.25 < $|\mu|$ < 0.5",r"0.5 < $|\mu|$ < 0.75",r"0.75 < $|\mu|$ < 1.0"]
mu_bins = [0.0,0.25,0.5,0.75]



####### Linear power spectra param


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


####### Minuit parameters



non_linear_params_0 = {"k_nl" : 1 ,
                       "a_nl" : 1 ,
                       "k_p" : 1 ,
                       "a_p" : 1 ,
                       "k_v0" : 1 ,
                       "a_v0" : 1 ,
                       "k_v1" : 1 ,
                       "a_v1" : 1}
non_linear_errors_0 = {"limit_k_nl" : (10**-1,10**2) ,
                       "limit_a_nl" : (0,4) ,
                       "limit_k_p" : (10**-1,10**2) ,
                       "limit_a_p" : (0,4) ,
                       "limit_k_v0" : (10**-1,10**2) ,
                       "limit_a_v0" : (0,4) ,
                       "limit_k_v1" : (10**-1,10**2) ,
                       "limit_a_v1" : (0,4)}


non_linear_params_1 = {"q_1" : 0,
                       "q_2" : 0,
                       "k_v" : 1,
                       "a_v": 1,
                       "b_v" :1,
                       "k_p" : 10}

non_linear_errors_1 = {"limit_q_1" : (0,2) ,
                       "limit_q_2" : (-2,2) ,
                       "limit_k_v" : (10**-8,10**3) ,
                       "limit_a_v" : (-2,2) ,
                       "limit_b_v" : (0,2) ,
                       "limit_k_p" : (10**-1,10**3)}


linear_params = {"b" : 0.1 ,
                 "beta" : 1}
linear_errors = {"limit_b" : (10**-3,1),
                 "limit_beta" : (0.1,10)}


minuit_parameters = {}
minuit_parameters.update(linear_params)
minuit_parameters.update(linear_errors)


if(non_linear_model == "0"):
    non_linear_params = non_linear_params_0
    non_linear_errors = non_linear_errors_0
    minuit_parameters.update(non_linear_params)
    minuit_parameters.update(non_linear_errors)
    label_model = f"nonlinear{non_linear_model}"
elif(non_linear_model == "1"):
    non_linear_params = non_linear_params_1
    non_linear_errors = non_linear_errors_1
    minuit_parameters.update(non_linear_params)
    minuit_parameters.update(non_linear_errors)
    label_model = f"nonlinear{non_linear_model}"
else:
    label_model = "linear"


### Out params

grid_to_plot = ["1","2","3","4","5","6","7","8","9"]
grid_to_plot = ["2","3","6"]
grid_to_plot = ["A","B","C"]
grid_to_plot = ["2"]
grid_to_plot = ["B"]

base_name_out = "p3d_fit"

if __name__ == "__main__":

    os.makedirs(os.path.join(os.getcwd(),"fit_plots"),exist_ok=True)


    for i in range(len(grid_to_plot)):
        if(filename_pk == "class"):
            name_pm_file = None
            z_init = None
        else:
            name_pm_file = eval(f"file_pm_grid{grid_to_plot[i]}")
            z_init = redshift_pm

        name_p3d_file = eval(f"file_grid{grid_to_plot[i]}")
        legend_grid  = eval(f"legend_grid{grid_to_plot[i]}")
        class_settings.update({"z_pk" : f"99, {redshift}, 0.0"})
        name_out = os.path.join(os.getcwd(),
                                "fit_plots",
                                f"{base_name_out}_{legend_grid}_z{redshift}_{label_model}_{direction}{mu_binning}")

        (minuit,power_f,power_m,model) = fitter.fitter_k_mu(name_p3d_file,
                                                            filename_pk,
                                                            minuit_parameters,
                                                            sigma,
                                                            class_dict=class_settings,
                                                            z_simu=redshift,
                                                            non_linear_model=non_linear_model,
                                                            cost_name=cost_name,
                                                            ncall=ncall,
                                                            kmin=kmin,
                                                            kmax=kmax,
                                                            var_minos=var_minos,
                                                            error_estimator=error_estimator,
                                                            epsilon=epsilon,
                                                            rebin=rebin,
                                                            name_pm_file=name_pm_file,
                                                            z_init=z_init,
                                                            fix_args=fix_args)
        fitter.plot_fit(minuit,
                        power_f,
                        power_m,
                        model,
                        mu_bins,
                        mu_bins_legend,
                        name_out=name_out)
        kna = fitter.compute_kna(minuit,power_m,eps)
        print(kna)
