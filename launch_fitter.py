import rsd_fitter



filename = "./plt00216z_flux_pkmu.txt"
filename_pk = "./plt00216rhom_ps3d.txt"

non_linear_model=None
cost_name="least"
var = 1
sigma = 1
ncall = 1000



## "0" model 

non_linear_params = {"k_nl" : 1 ,"a_nl" : 1 ,"k_p" : 1 , "a_p" : 1 , "k_v0" : 1 , "a_v0" : 1 , "k_v1" : 1 , "a_v1" : 1}
non_linear_errors = {"limit_k_nl" : (10**-2,10**2) ,"limit_a_nl" : (0,4) ,"limit_k_p" : (10**-2,10**2) , "limit_a_p" : (0,4) , "limit_k_v0" : (10**-2,10**2) , "limit_a_v0" : (0,4) , "limit_k_v1" : (10**-2,10**2) , "limit_a_v1" : (0,4)}


# ## "1" model 

# non_linear_params = {"q_1" : 1,"q_2" : 1,"k_v" : 1,"a_v": 1,"b_v" :1,"k_p" : 1}
# non_linear_errors = {"error_q_1" : 1,"error_q_2" : 1,"error_k_v" : 1,"error_a_v": 1,"error_b_v" :1,"error_k_p" : 1}



linear_params = {"b" : 0.1 , "beta" : 1}
linear_errors = {"error_b" : 0.1 , "error_beta" : 1}


minuit_parameters = {}
minuit_parameters.update(linear_params)
minuit_parameters.update(linear_errors)
# minuit_parameters.update(non_linear_params)
# minuit_parameters.update(non_linear_errors)

if __name__ == "__main__":

    (minuit,power_f,power_m) = rsd_fitter.fitter_k_mu(filename,filename_pk,minuit_parameters,var,sigma,non_linear_model=non_linear_model,cost_name=cost_name,ncall=ncall)

    
    rsd_fitter.plot(power_f,power_m)
