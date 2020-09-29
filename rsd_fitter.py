import numpy as np
from iminuit import Minuit
import scipy.interpolate
import power_spectra




def read_pfkmu(filename,power_weighted=False):
    power = power_spectra.FluxPowerSpectrum.init_from_gimlet(filename,kmu=True,power_weighted=power_weighted)
    return(power)

def read_pfkperpkpar(filename,power_weighted=False):
    power = power_spectra.FluxPowerSpectrum.init_from_gimlet(filename,kmu=False,power_weighted=power_weighted)
    return(power)

def read_pk(filename,power_weighted=False):
    power = power_spectra.MatterPowerSpectrum.init_from_gimlet(filename,power_weighted=power_weighted)
    return(power)



def rebin_matter_power(power_m,k_m,k_f):
    print(power_m.shape,k_m.shape)
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

def D1(k,mu,q_1,q_2,k_v,a_v,b_v,k_p,matter_power_spectrum):
    Delta = (1/(2*np.pi**2))*k**3 * matter_power_spectrum
    non_linear_term = np.exp((q_1*Delta**2 + q_2*Delta**4)*(1-((k/k_v)**a_v)*mu**b_v) - (k/k_p)**2)
    return(non_linear_term)


def Pl(k):
    """ Put to one since the fitted data is Pf/Pl """
    return(1)


def Pf_model(matter_power_spectrum,non_linear_model="0"):
    
    def modelD0(x,b,beta,k_nl,a_nl,k_p,a_p,k_v0,a_v0,k_v1,a_v1):
        k,mu = x[0],x[1]
        Pf = b**2 * (1 + beta * mu**2)**2 * Pl(k) * D0(k,mu,k_nl,a_nl,k_p,a_p,k_v0,a_v0,k_v1,a_v1)
        return(Pf)


    def modelD1(x,b,beta,q_1,q_2,k_v,a_v,b_v,k_p):
        k,mu = x[0],x[1]
        Pf = b**2 * (1 + beta * mu**2)**2 * Pl(k) * D1(k,mu,q_1,q_2,k_v,a_v,b_v,k_p,matter_power_spectrum)
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
        z = (data_y - ym) / data_yerr ** 2
        return np.nansum(z ** 2)

    def costD1(b,beta,q_1,q_2,k_v,a_v,b_v,k_p):
        ym = model(data_x,b,beta,q_1,q_2,k_v,a_v,b_v,k_p)
        z = (data_y - ym) / data_yerr ** 2
        return np.nansum(z ** 2)

    def costlinear(b,beta):
        ym = model(data_x,b,beta)
        z = (data_y - ym) / data_yerr ** 2
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


def run_minuit(data_x,data_y,data_yerr,minuit_parameters,var,sigma,matter_power_spectrum,non_linear_model="0",cost_name="least",ncall=100):
    model = Pf_model(matter_power_spectrum,non_linear_model=non_linear_model)
    cost = cost_function(model,data_x,data_y,data_yerr,cost_name,non_linear_model=non_linear_model)
    minuit_parameters.update({"errordef":1, "pedantic":False, "print_level": 0})
    minuit = Minuit(cost,**minuit_parameters)
    run_migrad(minuit,ncall=ncall)
    # run_minos(minuit,var,sigma,ncall=ncall)
    # run_hesse(minuit)
    return(minuit)

def run_migrad(minuit,ncall=100):
    minuit.migrad(ncall,resume=True)

def run_minos(minuit,var,sigma,ncall=100):
    return(minuit.minos(var,sigma,ncall))

def run_hesse(minuit):
    minuit.hesse()




### Main




def fitter_k_mu(pf_file,pk_file,minuit_parameters,var,sigma,power_weighted=False,non_linear_model="0",cost_name="least",ncall=100):
    power_f = read_pfkmu(pf_file,power_weighted=power_weighted)
    power_m = read_pk(pk_file,power_weighted=power_weighted)
    power_m_rebin = rebin_matter_power(power_m.power_array,power_m.k_array,power_f.k_array[0])
    # mask_data([3],k_edge, mu_edge, bincount, pwk, pwmu, power)
    data_x = power_f.k_array
    data_y = power_f.power_array / power_m_rebin
    data_yerr = 0.001 * data_y
    minuit = run_minuit(data_x,data_y,data_yerr,minuit_parameters,var,sigma,power_m,non_linear_model=non_linear_model,cost_name=cost_name,ncall=ncall)
    return(minuit,power_f,power_m)






def plot(power_f,power_m):
    power_m.plot_1d_pk()
    power_m.close_plot()
    power_f.plot_2d_pk([0.25,0.5,0.75])







