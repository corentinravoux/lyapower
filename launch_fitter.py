import rsd_fitter



filename = "./plt00225z_flux_pkmu.txt"





filename_pk = "class"


class_settings = {'h' : 0.6732117,
        'T_cmb' : 2.7255,
        'omega_b' : 0.022383,
        'N_ur' : 3.046,
        'omega_cdm' : 0.1201074794045102,
        'Omega_dcdmdr' : 0.0,
        'N_ncdm' : 0 , 
        'Omega_k' : 0.,
        'Omega_fld' : 0,
        'Omega_scf' : 0,
        'YHe' : 'BBN',
        'recombination' : 'RECFAST',
        'reio_parametrization' : 'reio_camb',
        'z_reio' : 11.357,
        'reionization_exponent' : 1.5,
        'reionization_width' : 0.5,
        'helium_fullreio_redshift' : 3.5,
        'helium_fullreio_width' : 0.5,
        'annihilation' : 0.,
        'decay' : 0.,
        'output' : 'tCl,pCl,lCl,nCl,sCl,mPk,dTk,vTk',
        'modes' : 's',
        'lensing' : 'yes',
        'ic' : 'ad',
        'gauge' : 'synchronous',
        'P_k_ini type' : 'analytic_Pk',
        'k_pivot' : 0.05,
        'A_s' :  2.100549e-9,
        'n_s' : 0.9660499,
        'alpha_s' : 0.,
        'l_max_scalars' : 2500,
        'headers' : 'yes',
        'format' : 'class',
        'write background' : 'no',
        'write thermodynamics' : 'no',
        'write primordial' : 'no',
        'write parameters' : 'yeap',
        'extra metric transfer functions' : 'yes',
        'input_verbose' : 0,
        'background_verbose' : 0,
        'thermodynamics_verbose' : 0,
        'perturbations_verbose' : 0,
        'transfer_verbose' : 1,
        'primordial_verbose' : 0,
        'spectra_verbose' : 0,
        'nonlinear_verbose' : 0,
        'lensing_verbose' : 0,
        'output_verbose' : 0,
        'root' : './',
        'P_k_max_h/Mpc' : 1000,
        'z_pk' : 3.810644783
        }

z_simu = 3.810644783



non_linear_model="1"
cost_name="least"
var = 1
sigma = 1
ncall = 10000
kmin = 0.5
kmax = 10
var_minos = None

power_weighted = False

# ## "0" model 

# non_linear_params = {"k_nl" : 1 ,"a_nl" : 1 ,"k_p" : 1 , "a_p" : 1 , "k_v0" : 1 , "a_v0" : 1 , "k_v1" : 1 , "a_v1" : 1}
# non_linear_errors = {"limit_k_nl" : (10**-2,10**2) ,"limit_a_nl" : (0,4) ,"limit_k_p" : (10**-2,10**2) , "limit_a_p" : (0,4) , "limit_k_v0" : (10**-2,10**2) , "limit_a_v0" : (0,4) , "limit_k_v1" : (10**-2,10**2) , "limit_a_v1" : (0,4)}


# ## "1" model 

non_linear_params = {"q_1" : 0,"q_2" : 0,"k_v" : 1,"a_v": 1,"b_v" :1,"k_p" : 1}
non_linear_errors = {"limit_q_1" : (-2,2) ,"limit_q_2" : (-2,2) ,"limit_k_v" : (10**-3,10**3) , "limit_a_v" : (0,4) , "limit_b_v" : (0,4) , "limit_k_p" : (10**-3,10**3)}



linear_params = {"b" : 0.1 , "beta" : 1}
linear_errors = {"error_b" : 0.1 , "error_beta" : 1}


minuit_parameters = {}
minuit_parameters.update(linear_params)
minuit_parameters.update(linear_errors)
minuit_parameters.update(non_linear_params)
minuit_parameters.update(non_linear_errors)

if __name__ == "__main__":

    (minuit,power_f,power_m,model) = rsd_fitter.fitter_k_mu(filename,filename_pk,minuit_parameters,sigma,class_dict=class_settings,z_simu=z_simu,non_linear_model=non_linear_model,cost_name=cost_name,ncall=ncall,kmin=kmin,kmax=kmax,var_minos=var_minos)

    # rsd_fitter.plot_pm(power_m)
    # rsd_fitter.plot_pf(power_f)
    # rsd_fitter.plot_pf_pm(power_f,power_m)
    rsd_fitter.plot_fit(minuit,power_f,power_m,model)