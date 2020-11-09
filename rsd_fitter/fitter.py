import numpy as np
import os
from iminuit import Minuit
import scipy.interpolate
from rsd_fitter import power_spectra
from rsd_fitter import CLASS



def read_pfkmu(filename,power_weighted=False,error_estimator=None,**kwargs):
    power = power_spectra.FluxPowerSpectrum.init_from_gimlet(filename,kmu=True,power_weighted=power_weighted,error_estimator=error_estimator,**kwargs)
    return(power)

def read_pfkperpkpar(filename,power_weighted=False,error_estimator=None,**kwargs):
    power = power_spectra.FluxPowerSpectrum.init_from_gimlet(filename,kmu=False,power_weighted=power_weighted,error_estimator=error_estimator,**kwargs)
    return(power)

def read_pk(filename,power_weighted=False,error_estimator=None,**kwargs):
    power = power_spectra.MatterPowerSpectrum.init_from_gimlet(filename,power_weighted=power_weighted,error_estimator=error_estimator,**kwargs)
    return(power)



def rebin_matter_power(power_m,k_m,k_f):
    power = scipy.interpolate.interp1d(k_m,power_m,bounds_error=False,fill_value=np.nan)
    power_m_rebin = power(k_f)
    return(power_m_rebin)




def mask_data(indexes,*args):
    mask = np.full(args[0].shape,False)
    for i in indexes:
        mask |= np.abs(args[i]) < np.abs(10**-10 * np.mean(args[i]))
    for i in range(len(args)):
        args[i][mask] = np.nan




### models & cost

def D0(k,mu,k_nl,a_nl,k_p,a_p,k_v0,a_v0,k_v1,a_v1):
    return(np.exp((k/k_nl)**a_nl - (k/k_p)**a_p - ((k*mu)/(k_v0 * (1 + (k/k_v1))**a_v1))**a_v0))

def D1(k,mu,q_1,q_2,k_v,a_v,b_v,k_p,linear_power_spectrum):
    Delta = (1/(2*np.pi**2))*k**3 * linear_power_spectrum
    non_linear_term = np.exp((q_1*Delta**2 + q_2*Delta**4)*(1-((k/k_v)**a_v)*mu**b_v) - (k/k_p)**2)
    return(non_linear_term)


def Pl(k):
    """ Put to one since the fitted data is Pf/Pl """
    return(1)


def Pl_class(k_array,settings,z,name="class"):
    my_class =CLASS.MyClass(os.getcwd(),settings)
    kmin,kmax,nb_points = np.log10(np.min(k_array)),np.log10(np.max(k_array)),2*len(k_array)
    (Power,Transfer,sigma_8) = my_class.write_pk_tk(z,name,kmin=kmin,kmax=kmax,nb_points=nb_points)
    h_normalized = True
    power = power_spectra.MatterPowerSpectrum(k_array=Power[:,0],power_array=Power[:,1],dimension="1D",specie="matter",h_normalized=h_normalized)
    return(power)


def Pf_model(linear_power_spectrum=None,non_linear_model="0"):

    def modelD0(x,b,beta,k_nl,a_nl,k_p,a_p,k_v0,a_v0,k_v1,a_v1):
        k,mu = x[0],x[1]
        Pf = b**2 * (1 + beta * mu**2)**2 * Pl(k) * D0(k,mu,k_nl,a_nl,k_p,a_p,k_v0,a_v0,k_v1,a_v1)
        return(Pf)


    def modelD1(x,b,beta,q_1,q_2,k_v,a_v,b_v,k_p):
        k,mu = x[0],x[1]
        Pf = b**2 * (1 + beta * mu**2)**2 * Pl(k) * D1(k,mu,q_1,q_2,k_v,a_v,b_v,k_p,linear_power_spectrum)
        return(Pf)

    def modellinear(x,b,beta):
        k,mu = x[0],x[1]
        Pf = b**2 * (1 + beta * mu**2)**2 * Pl(k)
        return(Pf)


    if(non_linear_model == "0"):
        return(modelD0)
    elif(non_linear_model == "1"):
        return(modelD1)
    elif(non_linear_model == None):
        return(modellinear)


def custom_least_squares(model,data_x,data_y,data_yerr,non_linear_model="0"):

    def costD0(b,beta,k_nl,a_nl,k_p,a_p,k_v0,a_v0,k_v1,a_v1):
        ym = model(data_x,b,beta,k_nl,a_nl,k_p,a_p,k_v0,a_v0,k_v1,a_v1)
        z = (data_y - ym) / data_yerr
        return np.nansum(z ** 2)

    def costD1(b,beta,q_1,q_2,k_v,a_v,b_v,k_p):
        ym = model(data_x,b,beta,q_1,q_2,k_v,a_v,b_v,k_p)
        z = (data_y - ym) / data_yerr
        return np.nansum(z ** 2)

    def costlinear(b,beta):
        ym = model(data_x,b,beta)
        z = (data_y - ym) / data_yerr
        return np.nansum(z ** 2)

    if(non_linear_model == "0"):
        return(costD0)
    elif(non_linear_model == "1"):
        return(costD1)
    elif(non_linear_model == None):
        return(costlinear)

def cost_function(model,data_x,data_y,data_yerr,cost_name,non_linear_model="0"):
    if(cost_name =="least"):
        return(custom_least_squares(model,data_x,data_y,data_yerr,non_linear_model=non_linear_model))



### Minuit stuff


def run_minuit(data_x,data_y,data_yerr,minuit_parameters,sigma,linear_power_spectrum,non_linear_model="0",cost_name="least",ncall=100,var_minos=None):
    model = Pf_model(linear_power_spectrum=linear_power_spectrum,non_linear_model=non_linear_model)
    cost = cost_function(model,data_x,data_y,data_yerr,cost_name,non_linear_model=non_linear_model)
    minuit_parameters.update({"errordef":1, "pedantic":False, "print_level": 0})
    minuit = Minuit(cost,**minuit_parameters)
    print(run_migrad(minuit,ncall=ncall))
    run_hesse(minuit)
    run_minos(minuit,sigma,ncall=ncall,var_minos=var_minos)
    return(minuit,model)

def run_migrad(minuit,ncall=1000):
    return(minuit.migrad(ncall,resume=True))

def run_minos(minuit,sigma,ncall=1000,var_minos=None):
    if(minuit.valid):
        return(minuit.minos(var=var_minos,sigma=sigma,ncall=ncall))

def run_hesse(minuit):
    return(minuit.hesse())




### Main




def fitter_k_mu(pf_file,pk_file,minuit_parameters,sigma,power_weighted=False,class_dict=None,z_simu=None,non_linear_model="0",cost_name="least",ncall=100,kmax=None,kmin=None,var_minos=None,error_estimator=None,**kwargs):
    power_f = read_pfkmu(pf_file,power_weighted=power_weighted,error_estimator=error_estimator,**kwargs)
    power_f.cut_extremum(kmin,kmax)
    if(pk_file == "class"):
        power_l = Pl_class(power_f.k_array[0],class_dict,z_simu,name="class")
    else:
        power_l = read_pk(pk_file,power_weighted=power_weighted)
    power_l_rebin = rebin_matter_power(power_l.power_array,power_l.k_array,power_f.k_array[0])
    # mask_data([3],k_edge, mu_edge, bincount, pwk, pwmu, power)
    data_x = power_f.k_array
    data_y = power_f.power_array / power_l_rebin
    if(power_f.error_array is not None):
        data_yerr = power_f.error_array
    else: data_yerr = 0.001 * data_y
    minuit,model = run_minuit(data_x,data_y,data_yerr,minuit_parameters,sigma,power_l_rebin,non_linear_model=non_linear_model,cost_name=cost_name,ncall=ncall,var_minos=var_minos)
    return(minuit,power_f,power_l,model)


### Plots




def plot_pl(power_l):
    power_l.open_plot()
    power_l.plot_1d_pk()

def plot_pf(power_f,mu_bin,legend):
    power_f.open_plot()
    power_f.plot_2d_pk(mu_bin,legend=legend)


def plot_pf_pm(power_f,power_m,mu_bin,legend):
    power_f.open_plot()
    power_m_rebin = rebin_matter_power(power_m.power_array,power_m.k_array,power_f.k_array[0])
    power_f.power_array = power_f.power_array / power_m_rebin
    power_f.plot_2d_pk(mu_bin,legend=legend,ps="x")


def plot_fit(minuit,power_f,power_l,model,mu_bin,legend,name_out="fit_results"):
    minuit_params = []
    for i in range(len(minuit.params)):
        minuit_params.append(minuit.params[i].value)
    power_l_rebin = rebin_matter_power(power_l.power_array,power_l.k_array,power_f.k_array[0])
    pf_over_pm_data = power_f.power_array / power_l_rebin
    if(power_f.error_array is not None): error_array = power_f.error_array/ power_l_rebin
    else: error_array = None
    pf_over_pm_model = model(power_f.k_array,*minuit_params)

    color = [f"C{i}" for i in range(len(mu_bin))]

    h_normalized = power_f.h_normalized
    power1 = power_spectra.FluxPowerSpectrum(k_array=power_f.k_array,power_array= pf_over_pm_model,dimension="3D",h_normalized=h_normalized)
    power1.plot_2d_pk(mu_bin,color=color)


    h_normalized = power_l.h_normalized
    power2 = power_spectra.FluxPowerSpectrum(k_array=power_f.k_array,power_array= pf_over_pm_data,error_array=error_array,dimension="3D",h_normalized=h_normalized)
    power2.plot_2d_pk(mu_bin,legend=legend,ps="x",ls="None",color=color,y_label=r"$P_f/P_l$")
    power2.save_fig(name_out)

    minuit_to_latex(minuit,name=name_out)



### Latex


def minuit_to_latex(minuit,name=""):
<<<<<<< HEAD
    file = open(f"{name}_minuit_matrix.tex",'w')
    file.write(minuit.latex_matrix().__str__())
    file.close()
=======
    file = open(f"{name}_minuit_results.tex",'w')
    file.write(minuit.latex_matrix().__str__())
    file.close()

>>>>>>> 15797c5e56295cd9677fe9d0940a39014c5a82b7

    file = open(f"{name}_minuit_params.tex",'w')
    file.write(minuit.latex_param().__str__())
    file.close()







###
