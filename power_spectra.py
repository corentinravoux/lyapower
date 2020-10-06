#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 16:29:56 2019

@author: cravoux
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import math
from matplotlib.pyplot import cm




class PowerSpectrum(object):
    
    
    def __init__(self,k_array=None,power_array=None,file_init=None,size_box=None):
        self.k_array = k_array
        self.power_array = power_array
        self.file_init = file_init
        self.size_box = size_box
        self.h_normalized = None
        

    @classmethod
    def init_from_genpk_file(cls,name_file,size_box):
        """Load a GenPk format power spectum, plotting the DM and the neutrinos (if present)
        Does not plot baryons."""
        #Load DM P(k)
        matpow=np.loadtxt(name_file)
        if(size_box is None): raise KeyError("To load a GenPk output, the size of the box in Mpc.h-1 must be given")
        scale=2*math.pi/size_box
        #Adjust Fourier convention to match CAMB.
        simk=matpow[1:,0]*scale
        Pk=matpow[1:,1]/scale**3*(2*math.pi)**3
        return(cls(k_array=simk,power_array=Pk,file_init=name_file,size_box=size_box))


    @classmethod
    def init_from_ascii_file(cls,name_file):
        file_power=np.loadtxt(name_file)
        k_array = file_power[1:,0]
        power_array=file_power[1:,1]
        power = cls(k_array=k_array,power_array=power_array,file_init=name_file,size_box=None)
        power.h_normalized = True
        return(power)


    
    

    def rebin_arrays(self,nb_bin,operation="mean"):
        if nb_bin is None : return()
        new_k = np.logspace(np.log10(min(self.k_array)),np.log10(max(self.k_array)),nb_bin)
        new_Pk = np.zeros(new_k.shape)
        for i in range(nb_bin-1):
            mask = (self.k_array >= new_k[i]) & (self.k_array < new_k[i+1])
            if operation.lower() in ["mean", "average", "avg"]:
                if(len(self.power_array[mask])==0):
                    nearest_index = np.argmin(np.abs(self.k_array - new_k[i]))
                    new_Pk[i] = self.power_array[nearest_index]
                else :
                    new_Pk[i] = np.mean(self.power_array[mask])
            elif operation.lower() in ["gauss"]:
                from scipy import signal
                gaussian_weights = signal.gaussian(int(len(self.power_array[mask])),int(len(self.power_array[mask]))/4)
                if(len(self.power_array[mask])==0):
                    nearest_index = np.argmin(np.abs(self.k_array - new_k[i]))
                    new_Pk[i] = self.power_array[nearest_index]
                else :
                    new_Pk[i] = np.average(self.power_array[mask],axis=0,weights=gaussian_weights)
        self.k_array=new_k
        self.power_array=new_Pk



    def plot_1d_pk(self,color=None, ls="-",flux_factor = None,ps =None):
        plt.ylabel("P(k) /(h-3 Mpc3)")
        if(self.h_normalized): plt.xlabel("k (h Mpc-1)")
        else: plt.xlabel("k (Mpc-1)")
        plt.title("Power spectrum")
        plt.loglog(self.k_array,self.power_array,marker = ps, linestyle=ls, color=color)        


    def plot_2d_pk(self,bin_edges,color=None, ls="-",flux_factor = None,ps =None):
        bin_start = 0.0
        for i in range(len(bin_edges)):
            mask = (self.k_array[1] < bin_edges[i]) & (self.k_array[1] >= bin_start)
            plt.loglog(self.k_array[0][mask],self.power_array[mask],marker = ps, linestyle=ls, color=color)
            bin_start = bin_edges[i]
        plt.ylabel("P(k) /(h-3 Mpc3)")
        if(self.h_normalized): plt.xlabel("k (h Mpc-1)")
        else: plt.xlabel("k (Mpc-1)")
        plt.title("Power spectrum")




    def plot_several_power_spectrum(self,Pks,k_space,name,legend):
        plt.figure()
        color = cm.rainbow(np.linspace(0, 1, len(Pks)))
        for i in range(len(Pks)):
            plt.loglog(k_space,np.array(Pks[i]),'b',color=color[i])
        plt.grid()
        plt.legend(legend)
        plt.savefig(name +"matter_power_spectrum.pdf",format="pdf")



    def plot_comparison_spectra(self,list_spectra,label_list,diff_extremums=0.1,normalize=True):
        fig, ax = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, sharex=True,figsize=(8,6))  #note that height ratios can be used to scale the size of top vs bottom part
        self.add_comparison_spectra(list_spectra,ax,normalize=normalize)
        ax[0].set_title(r'...')
        if(normalize): ax[0].set_ylabel(r'$\Delta_m^2$')
        else: ax[0].set_ylabel(r'$P_m$')
        if(self.h_normalized): ax[1].set_xlabel("k (h Mpc-1)")
        else: ax[1].set_xlabel("k (Mpc-1)")
        ax[1].set_ylabel(r'$\Delta_m^2/\Delta_{m,ref}^2-1$')
        ax[0].set_xscale('log')
        ax[0].set_yscale('log')
        if(len(label_list)<=5):ax[0].legend(label_list)
        else:ax[0].legend(label_list,ncol=2)
        ax[1].set_ylim(-diff_extremums,diff_extremums)
        return(ax)
        
        
    def add_comparison_spectra(self,list_spectra,ax,normalize=True):
        kref = self.k_array
        if(normalize): kpkref = ( self.k_array**3 * self.power_array )/2*(np.pi)**2
        else: kpkref = self.power_array
        karr,kpkarr = [], []
        karr.append(kref)
        kpkarr.append(kpkref)
        for i in range(len(list_spectra)):
            karr.append(list_spectra[i].k_array)
            if(normalize): kpkarr.append((list_spectra[i].k_array**3 * list_spectra[i].power_array)/2*(np.pi)**2)
            else: kpkarr.append(list_spectra[i].power_array)
        #karr is your array of x values, i.e. a numpy array with shape (nlines,nvalues)
        #kpkarr is your array of y values same shape (nlines,nvalues)
        #larr is your array of labels (nlines)
        #kref,kpkref are the reference values (nvalues)
        for k, kpk in zip(karr, kpkarr):
            interp = interp1d(k, kpk, bounds_error=False)
            ax[0].plot(k, kpk)
            ax[1].plot(kref,(interp(kref)/kpkref)-1)



    def save_fig(self,name_out):
        fig =plt.gcf()
        fig.savefig("{}.pdf".format(name_out),format="pdf")

    def close_fig(self):
        fig =plt.gcf()
        fig.close()


    def save_plot(self,nameout,format_out = "pdf"):
        plt.savefig(nameout,format=format_out)

    def close_plot(self):
        plt.close()
        
    def get_k_value(self,k):
        interp = interp1d(self.k_array, self.power_array, bounds_error=True)
        print(interp(k))
        return(interp(k))


    def change_k_normalization(self,wanted_h_normalized,h):
        if(self.h_normalized is None):
            raise KeyError("The actual normalization of the k vector is not know")
        if(self.h_normalized):
            if(wanted_h_normalized):
                return()
            else :
                self.k_array = self.k_array/h
                self.h_normalized = False
                return()
        else:
            if(wanted_h_normalized):
                self.k_array = self.k_array * h
                self.h_normalized = True
                return()
            else:
                return()
                


class MatterPowerSpectrum(PowerSpectrum):
    
    def __init__(self,dimension,specie,**kwargs):
        super(MatterPowerSpectrum,self).__init__(**kwargs)
        if dimension not in ["1D","3D"]: raise KeyError("Dimension of spectrum not available, please choose between 1D and 3D")
        self.dimension = dimension
        self.specie = specie


    @classmethod
    def init_from_gimlet(cls,namefile,specie="unknown",power_weighted=False):
        """ Pm(k) gimlet file contains 
         - k: edge (higher) of the k bin considered
         - bincount: number of mode (pairs) computed in the bin 
         - pwk: power weighted k
         - power: power of the bin"""
        f = np.loadtxt(namefile)
        k, pwk,power = f[:,0],f[:,2],f[:,3]
        if (power_weighted) : k_array = pwk
        else : k_array = k
        dimension = "3D"
        return(cls(dimension,specie,k_array=k_array,power_array=power,file_init=namefile,size_box=None))




class FluxPowerSpectrum(PowerSpectrum):
    
    
    def __init__(self,dimension,**kwargs):
        super(FluxPowerSpectrum,self).__init__(**kwargs)
        if dimension not in ["1D","3D"]: raise KeyError("Dimension of spectrum not available, please choose between 1D and 3D")
        self.dimension = dimension


    @classmethod
    def init_from_gimlet(cls,namefile,kmu=True,power_weighted=False):
        if(kmu): (k1_edge, k2_edge, pwk1, pwk2, power) = cls.init_kmu(namefile)
        else: (k1_edge, k2_edge, pwk1, pwk2, power) = cls.init_kperpar(namefile)
        if (power_weighted) : k1_array , k2_array = pwk1,pwk2
        else : k1_array , k2_array = k1_edge,k2_edge
        k_array = np.stack([k1_array, k2_array])
        dimension = "3D"
        return(cls(dimension,k_array=k_array,power_array=power,file_init=namefile,size_box=None))        
        
        
    @classmethod
    def init_kmu(cls,namefile):
        """ Pf(k,mu) gimlet file contains 
         - k_edge: edge (higher) of the k bin considered
         - mu_edge: edge (lower) of the mu bin considered (mu positive)
         - bincount: number of mode (pairs) computed in the bin 
         - pwk: power weighted k
         - pwmu: power weighted mu
         - power: power of the bin"""
        f = np.loadtxt(namefile)
        k_edge, mu_edge, pwk, pwmu,power = f[:,0],f[:,1],f[:,3],f[:,4],f[:,5]
        return(k_edge, mu_edge, pwk, pwmu, power)



    @classmethod
    def init_kperpar(cls,namefile):
        """ Pf(kperp,kpar) gimlet file contains 
         - k_perp: edge (higher) of the k perp bin considered
         - k_par: edge (higher) of the k par bin considered
         - bincount: number of mode (pairs) computed in the bin 
         - pwkperp: power weighted k perp
         - pwkpar: power weighted k par
         - power: power of the bin"""
        f = np.loadtxt(namefile)
        k_perp, k_par, pwkperp, pwkpar,power = f[:,0],f[:,1],f[:,3],f[:,4],f[:,5]
        return(k_perp, k_par, pwkperp, pwkpar,power)






def init_spectrum(type_init,filename,boxsize=None):
    if(type_init == "GENPK"):
        spectrum = PowerSpectrum.init_from_genpk_file(filename, boxsize)
    elif(type_init == "ASCII"):
        spectrum = PowerSpectrum.init_from_ascii_file(filename)
    return(spectrum)



def launch_comparison_power_spectra(list_file,type_file,label_list,name_out,diff_extremums=0.1,rebin=None,rebin_method=None,flux_factor=None,normalize=True,size_box=None,wanted_normalization=None,h_normalization=None):
    reference_spectrum = init_spectrum(type_file,list_file[0],boxsize=size_box)
    if(rebin is not None): reference_spectrum.rebin_arrays(rebin,operation=rebin_method)
    if(flux_factor is not None) : reference_spectrum.power_array = reference_spectrum.power_array * flux_factor[0]
    if(wanted_normalization is not None) : reference_spectrum.change_k_normalization(wanted_normalization,h_normalization)
    list_spectra = []
    for i in range(1,len(list_file)):
        ps = init_spectrum(type_file,list_file[i],boxsize=size_box)
        if(rebin is not None): ps.rebin_arrays(rebin,operation=rebin_method)
        if(flux_factor is not None) : ps.power_array = ps.power_array * flux_factor[i]
        if(wanted_normalization is not None) : ps.change_k_normalization(wanted_normalization,h_normalization)
        list_spectra.append(ps)
    reference_spectrum.plot_comparison_spectra(list_spectra,label_list,diff_extremums=diff_extremums,normalize=normalize)
    reference_spectrum.save_fig(name_out)
    reference_spectrum.close_plot()


def launch_comparison_power_spectra_different_ref(list_file,type_file,label_list,name_out,diff_extremums=0.1,rebin=None,rebin_method=None,flux_factor=None,normalize=True,size_box=None,wanted_normalization=None,h_normalization=None):
    for j in range(len(list_file)):
        reference_spectrum = init_spectrum(type_file[j],list_file[j][0],boxsize=size_box)
        if(rebin[j] is not None): reference_spectrum.rebin_arrays(rebin[j],operation=rebin_method)
        if(flux_factor[j] is not None) : reference_spectrum.power_array = reference_spectrum.power_array * flux_factor[j][0]
        if(wanted_normalization is not None) : reference_spectrum.change_k_normalization(wanted_normalization[j],h_normalization)
        list_spectra = []
        for i in range(1,len(list_file[j])):
            ps = init_spectrum(type_file[j],list_file[j][i],boxsize=size_box)
            if(rebin[j] is not None): ps.rebin_arrays(rebin[j],operation=rebin_method)
            if(flux_factor[j] is not None) : ps.power_array = ps.power_array * flux_factor[j][i]
            if(wanted_normalization is not None) : ps.change_k_normalization(wanted_normalization[j],h_normalization)
            list_spectra.append(ps)
        if(j==0): ax = reference_spectrum.plot_comparison_spectra(list_spectra,label_list,diff_extremums=diff_extremums,normalize=normalize)
        else:reference_spectrum.add_comparison_spectra(list_spectra,ax,normalize=normalize)
    reference_spectrum.save_fig(name_out)
    reference_spectrum.close_plot()










