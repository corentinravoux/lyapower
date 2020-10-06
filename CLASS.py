#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 14:12:40 2019

@author: cravoux
"""
import numpy as np
from scipy.interpolate import interp1d



header_tk_class = """ Transfer functions T_i(k)  at redshift z={} in CLASS convention for adiabatic (AD) mode (normalized to initial curvature=1)
            for k= {} to {} h/Mpc,
            number of wavenumbers equal to {}
            d_i   stands for (delta rho_i/rho_i)(k,z) with above normalization 
            d_tot stands for (delta rho_tot/rho_tot)(k,z) with rho_Lambda NOT included in rho_tot
            (note that this differs from the transfer function output from CAMB/CMBFAST, which gives the same
             quantities divided by -k^2 with k in Mpc^-1; use format=camb to match CAMB)
            \n"""
header_tk_camb =  """ Transfer functions T_i(k) at redshift z={} in CAMB convention for adiabatic (AD) mode (normalized to initial curvature=1)
            for k= {} to {} h/Mpc,
            number of wavenumbers equal to {}
            k vector equidistant in log space. Suitable for 2LPT if no header. Suitable for cosmicic
            (note that this differs from the transfer function output from CLASS, which gives the same
            quantities multiplied by -k^2 with k in Mpc^-1; use format=class to match CLASS)
            \n"""
header_pk_class =  """ Power spectrum in the class format (in h3.Mpc-3) at redshift z={}
            for k= {} to {} h/Mpc,
            number of wavenumbers equal to {}
            k (in h.Mpc-1) & Pk (in h3.Mpc-3)
            \n"""
header_pk_camb =  """ Power spectrum in the camb format (in h3.Mpc-3) at redshift z={}
            for k= {} to {} h/Mpc,
            number of wavenumbers equal to {}
            k vector equidistant in log space. Suitable for 2LPT if no header and same size than Pk. Suitable for Cosmicic.
            k (in h.Mpc-1) & Pk (in h3.Mpc-3)
            \n"""






class MyClass(object):
    def __init__(self,pwd,settings,hierarchy="EQUI",class_version="public"):
        if(class_version == "public"):
            from classy import Class
        elif(class_version == "pk1D"):
            from classy_pk1D import Class

        self.pwd = pwd
        self.settings = settings
        self.model = Class()
        self.model.set(settings)

        self.define_output_format()

        self.hierarchy = hierarchy
        self.class_exe_path = "/local/home/cravoux/Software/class/class_public/class"
        


    def define_output_format(self):
        if("format" not in self.settings.keys()): 
            self.output_format = "class"
        else:
            self.output_format = self.settings["format"]
            

    def write_pk_tk(self,z,name,kmin=-4,kmax= 3,nb_points=2000):
        """ return Power Spectra with k in h.Mpc-1 and P in h3.Mpc-3 in the convention of CAMB and CLASS (the one needed for MP-Gadget).
            return Transfer function in the convention of CLASS (needed for MP-Genic) or CAMB (needed for 2LPT) !!!
            note : . CLASS format for MP-GENIC, CAMB format for 2LPT, CAMB format for COSMICIC
                   . for 2LPT & COSMICIC: no header! k vector need to be equidistant in log-k space and length of pk must be the same than tk
        """        
        self.model.compute()
        Power,sigma_8 = self.write_pk(name,z,kmin=kmin,kmax=kmax,nb_points=nb_points)
        Transfer = self.write_tk(name,z)
        return(Power,Transfer,sigma_8)
    

    def write_pk_tk_grid_variation(self,z,name,varying_dictionary,derived_params=None,kmin=-4,kmax= 3,nb_points=2000):
        """ return Power Spectra with k in h.Mpc-1 and P in h3.Mpc-3 in the convention of CAMB and CLASS (the one needed for MP-Gadget).
            return Transfer function in the convention of CLASS (needed for MP-Genic) or CAMB (needed for 2LPT) !!!
            note : . CLASS format for MP-GENIC, CAMB format for 2LPT, CAMB format for COSMICIC
                   . for 2LPT & COSMICIC: no header! k vector need to be equidistant in log-k space and length of pk must be the same than tk
        """        
        dict_deriv = {}
        if(derived_params is not None):
            for i in range(len(derived_params)):
                dict_deriv[derived_params[i]] = []
        list_key = list(varying_dictionary.keys())
        for i in range(len(varying_dictionary[list_key[0]])):
            dict_to_set = {list_key[j] : varying_dictionary[list_key[j]][i] for j in range(len(list_key))}
            self.model.set(dict_to_set)
            self.model.compute()
            Power,sigma_8 = self.write_pk(name,z,kmin=kmin,kmax=kmax,nb_points=nb_points)
            Transfer = self.write_tk(name,z)
            if(derived_params is not None):
                par = self.model.get_current_derived_parameters(derived_params)
                print(par.keys())
                for key in par.keys():
                    dict_deriv[key].append(par[key])
        return(Power,Transfer,sigma_8,dict_deriv)

    def write_pk_tk_neutrino_mass(self,z,name,kmin=-4,kmax= 3,nb_points=2000,nb_neutrinos=3,neutrino_mass = None,method = "cbnu"):
        """ return Power Spectra with k in h.Mpc-1 and P in h3.Mpc-3 in the convention of CAMB and CLASS (the one needed for MP-Gadget).
            return Transfer function in the convention of CLASS (needed for MP-Genic) or CAMB (needed for 2LPT) !!!
        """
        m_ncdm = ""
        for j in range(nb_neutrinos):
            m_ncdm =m_ncdm + "0.0,"
        m_ncdm = m_ncdm[0:-1]
        self.model.set({'N_ncdm' : nb_neutrinos,'m_ncdm' : m_ncdm})
        self.model.compute()
        if (neutrino_mass is not None):
            omega_cdm_ref = self.model.Omega0_cdm()
            omega_nu_ref = self.model.Omega_nu
            m_ncdm = ""
            for j in range(nb_neutrinos):
                m_ncdm =m_ncdm +  str(neutrino_mass/nb_neutrinos) + ","
            m_ncdm = m_ncdm[0:-1]
            if(method == "cbnu"):
                self.model.set({'N_ncdm' : nb_neutrinos,'m_ncdm':m_ncdm,'omega_cdm':str((omega_cdm_ref - self.get_omeganu_from_mass(neutrino_mass) + omega_nu_ref)*self.model.h()**2)})
            if(method == "cb"):
                self.model.set({'N_ncdm' : nb_neutrinos,'m_ncdm':m_ncdm})
            self.model.compute()
            
        Power,sigma_8 = self.write_pk(name,z,kmin=kmin,kmax=kmax,nb_points=nb_points)
        Transfer = self.write_tk(name,z)
        return(Power,Transfer,sigma_8)



    def write_pk(self,name,z,kmin=-4,kmax=3,nb_points=2000,header_output=True):
        if(self.output_format == "class"):
            (k_space,Pk,sigma_8) = self.compute_power_spectrum(z,kmin=kmin,kmax=kmax,nb_points=nb_points) # k in h/Mpc
            header = header_pk_class.format(z,kmin,kmax,nb_points)
        elif(self.output_format == 'camb'): 
            tr = self.model.get_transfer(z,output_format='camb')
            k = np.array(tr['k (h/Mpc)']) 
            kpower=np.logspace(*(np.log10(k[[0,-1]]))*(1-1e-5),len(k)) #k in h/Mpc
            (k_space,Pk,sigma_8) = self.compute_power_spectrum(z,k_array=kpower)
            header = header_pk_camb.format(z,np.min(k_space),np.max(k_space),len(k_space))
        Power = np.stack([k_space,Pk],axis=1)
        if(header_output):
            np.savetxt("pk_{}_z{}.dat".format(name,z),Power,header = header)     
        else :
            np.savetxt("pk_{}_z{}.dat".format(name,z),Power)     
        return(Power,sigma_8)
        


    
    def write_tk(self,name,z,header_output=True):
        if(self.output_format == "class"):
            Transfer_dict = self.model.get_transfer(z,output_format='class')
            Transfer = np.stack([Tk for key,Tk in Transfer_dict.items()],axis=1)
            header = header_tk_class.format(z,np.min(Transfer[:,0]),np.max(Transfer[:,0]),len(Transfer[:,0]))
        elif(self.output_format == 'camb'): 
            Transfer_dict = self.model.get_transfer(z,output_format='camb')
            Transfer = np.stack([Tk for key,Tk in Transfer_dict.items()],axis=1)
            Transfer = self.interp_to_equidistant_log_space_camb_format(Transfer)
            header = header_tk_camb.format(z,np.min(Transfer[:,0]),np.max(Transfer[:,0]),len(Transfer[:,0]))
        list_key  = [key for key,Tk in Transfer_dict.items()]
        head = ''
        for i in range(len(list_key)):
            head = head + "            " + str(i+1) + ":" + list_key[i]
        header = header + head
        if(header_output):        
            np.savetxt("tk_{}_z{}.dat".format(name,z),Transfer,header = header,delimiter ="    ")
        else :
            np.savetxt("tk_{}_z{}.dat".format(name,z),Transfer,delimiter ="    ")
        return(Transfer)



    def compute_power_spectrum(self,z,kmin=-4,kmax=3,nb_points=2000,k_array=None,verbose=True):
        if(k_array is None):k_space = np.logspace(kmin,kmax,num=nb_points)
        else : k_space = k_array
        if(verbose):
            print("Omega matter = " + str(self.model.Omega_m()))
            print("Omega lambda = " + str(self.model.Omega_Lambda()))
            print("Omega baryon = " + str(self.model.Omega_b()))
            print("Omega dm = " + str(self.model.Omega0_cdm()))
            print("Omega k = " + str(self.model.Omega0_k()))
            print("Omega sum = " + str(self.model.Omega0_cdm() + self.model.Omega_b() +self.model.Omega_nu))
            print("Omega nu = " + str(self.model.Omega_nu))
            print("Omega rad = " + str(self.model.Omega_g()))
            sigma_8 = self.model.sigma(8/self.model.h(),z)
            print("sigma 8 at (z={}) = {}".format(z,sigma_8))
        Pk = []
        h=self.model.h()
        for k in k_space :
            Pk.append(self.model.pk(k*h,z)*h**3)
        return(k_space,Pk,sigma_8)
        
    
    def interp_to_equidistant_log_space_camb_format(self,Transfer):
        k_space = Transfer[:,0] 
        k_log_space=np.logspace(*(np.log10(k_space[[0,-1]]))*(1-1e-5),len(k_space)) #k in h/Mpc
        for i in range(1,Transfer.shape[-1]):
            interp = interp1d(k_space,Transfer[:,i], bounds_error=True)
            Transfer[:,i]=interp(k_log_space)
        Transfer[:,0]= k_log_space
        return(Transfer)



    def sigmaR_all_species_at(self,R,z):
        self.model.compute()
        return self.model.sigma(R/self.model.h(),z)

    
    def free_structure(self):
        self.model.struct_cleanup()

    def close_class(self):
        self.model.empty()
    