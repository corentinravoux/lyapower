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
from lyapower import utils_fitter as utils
import h5py







class PowerSpectrum(object):


    def __init__(self,k_array=None,power_array=None,error_array=None,file_init=None,size_box=None,h_normalized=None):
        self.k_array = k_array
        self.power_array = power_array
        self.file_init = file_init
        self.size_box = size_box
        self.h_normalized = h_normalized
        self.error_array = error_array

        # if(k_array is not None):
        #     print(mask[~mask])
        #     mask = self.k_array[0] != np.nan
        #     self.k_array = np.transpose(np.transpose(self.k_array)[mask])
        #     self.power_array = self.power_array[mask]



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
        h_normalized = True
        return(cls(k_array=simk,power_array=Pk,file_init=name_file,size_box=size_box,h_normalized=h_normalized))


    @classmethod
    def init_from_ascii_file(cls,name_file):
        file_power=np.loadtxt(name_file)
        k_array = file_power[1:,0]
        power_array=file_power[1:,1]
        h_normalized = True
        return(cls(k_array=k_array,power_array=power_array,file_init=name_file,size_box=None,h_normalized = h_normalized))





    def rebin_arrays(self,nb_bin,operation="mean"):
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



    def rebin_2d_arrays(self,
                        nb_bin,
                        operation="mean",
                        loglin=False,
                        k_loglin=None):
        if(loglin):
            mask_k_rebin = self.k_array[0] > k_loglin
            mu_old = self.k_array[1][mask_k_rebin]
            k_old = self.k_array[0][mask_k_rebin]
            power_old = self.power_array[mask_k_rebin]
            if(self.error_array is not None):
                error_old = self.error_array[mask_k_rebin]
        else:
            mu_old = self.k_array[1]
            k_old = self.k_array[0]
            power_old = self.power_array
            if(self.error_array is not None):
                error_old = self.error_array
        new_k = np.logspace(np.log10(min(k_old)),np.log10(max(k_old)),nb_bin)
        bin_centers= np.array([0.5 * (new_k[i] + new_k[i+1]) for i in range(len(new_k)-1)])
        mus = np.unique(mu_old)
        new_Pk = np.zeros(len(mus) * len(bin_centers))
        if(self.error_array is not None):
            new_error_array = np.zeros(len(mus) * len(bin_centers))
        for j in range(len(mus)):
            for i in range(nb_bin-1):
                mask = (k_old >= new_k[i])
                mask &= (k_old < new_k[i+1])
                mask &= (mu_old == mus[j])
                if operation.lower() in ["mean", "average", "avg"]:
                    if(len(power_old[mask])==0):
                        nearest_index = np.argmin(np.abs(k_old[0] - new_k[i]))
                        new_Pk[i*len(mus) + j] = power_old[nearest_index]
                        if(self.error_array is not None):
                            new_error_array[i*len(mus) + j] = error_old[nearest_index]
                    else :
                        new_Pk[i*len(mus) + j] = np.mean(power_old[mask])
                        if(self.error_array is not None):
                            new_error_array[i*len(mus) + j] = np.mean(error_old[mask])
                elif operation.lower() in ["gauss"]:
                    from scipy import signal
                    gaussian_weights = signal.gaussian(int(len(power_old[mask])),int(len(power_old[mask]))/4)
                    if(len(power_old[mask])==0):
                        nearest_index = np.argmin(np.abs(k_old - new_k[i]))
                        new_Pk[i*len(mus) + j] = power_old[nearest_index]
                        if(self.error_array is not None):
                            new_error_array[i*len(mus) + j] = error_old[nearest_index]
                    else :
                        new_Pk[i*len(mus) + j] = np.average(power_old[mask],axis=0,weights=gaussian_weights)
                        if(self.error_array is not None):
                            new_error_array[i*len(mus) + j] = np.average(error_old[mask],axis=0,weights=gaussian_weights)


        new_2d_k = np.transpose([[bin_centers[i],mus[j]] for i in range(len(bin_centers)) for j in range(len(mus))])
        if(loglin):
            self.k_array = np.array([np.concatenate([self.k_array[0][~mask_k_rebin],new_2d_k[0]]),
                                     np.concatenate([self.k_array[1][~mask_k_rebin],new_2d_k[1]])])
            self.power_array = np.concatenate([self.power_array[~mask_k_rebin],new_Pk])
            if(self.error_array is not None):
                self.error_array = np.concatenate([self.error_array[~mask_k_rebin],new_error_array])
        else:
            self.k_array=new_2d_k
            self.power_array=new_Pk
            if(self.error_array is not None):
                self.error_array = new_error_array


    def cut_extremum(self,kmin,kmax):
        mask = np.full(self.power_array.shape,True)
        if(kmin is not None):
            mask &= self.k_array[0,:] >= kmin
        if(kmax is not None):
            mask &= self.k_array[0,:] <= kmax
        self.k_array = np.transpose(np.transpose(self.k_array)[mask])
        self.power_array = self.power_array[mask]
        if(self.error_array is not None):self.error_array = self.error_array[mask]


    def put_label(self,ax,xunit=True,yunit=True,y_label ="P",x_label="k",labelsize_x=12,labelsize_y=12,fontsize=12):
        ylab,xlab = "",""
        if(yunit):
            if(self.h_normalized):
                ylab = r" [$\rm{h}^{-3}.\rm{Mpc}^3$]"
            else:
                ylab = r" [$\rm{Mpc}^3$]"
        if(xunit):
            if(self.h_normalized):
                xlab = r" [$\rm{h}.\rm{Mpc}^{-1}$]"
            else:
                xlab = r" [$\rm{Mpc}^{-1}$]"
        ax.tick_params(axis='y', labelsize=labelsize_y)
        ax.set_ylabel(f"{y_label}{ylab}", fontsize=fontsize)
        ax.tick_params(axis='y', labelsize=labelsize_x)
        ax.set_xlabel(f"{x_label}{xlab}", fontsize=fontsize)



    def prepare_axes(self,kwargs):
        ax = utils.return_key(kwargs,"ax", plt.gcf().get_axes())
        if(len(ax) == 0):
            ax_to_plot = plt.gca()
            ax_comparison = plt.gca()
        elif(len(ax) == 1):
            ax_to_plot = ax[0]
            ax_comparison = ax[0]
        else:
            ax_to_plot = ax[0]
            ax_comparison = ax[1]
        return(ax_to_plot,ax_comparison)



    def plot_1d_pk(self,**kwargs):
        style = utils.return_key(kwargs,"style",None)
        if style is not None:
            plt.style.use(style)
        comparison = utils.return_key(kwargs,"comparison",None)
        (ax_to_plot,ax_comparison) = self.prepare_axes(kwargs)
        self.put_label(ax_to_plot)

        color = utils.return_key(kwargs,"color",None)


        if(comparison is not None):
            power_array_comparison = interp1d(self.k_array,
                                              self.power_array,
                                              bounds_error=False,
                                              fill_value=np.NaN)(comparison.k_array)
            ax_comparison.plot(comparison.k_array,
                               (comparison.power_array - power_array_comparison)/comparison.power_array,
                               color=color)
        ax_to_plot.set_title("Power spectrum")
        ax_to_plot.plot(self.k_array,
                        self.power_array,
                        marker = utils.return_key(kwargs,"ps",None),
                        linestyle= utils.return_key(kwargs,"ls","-"),
                        color=color)

        xscale = utils.return_key(kwargs,"xscale","log")
        yscale = utils.return_key(kwargs,"yscale","log")

        ax_to_plot.set_xscale(xscale)
        ax_to_plot.set_yscale(yscale)


        x_min_lim = utils.return_key(kwargs,"x_min_lim",None)
        x_max_lim = utils.return_key(kwargs,"x_max_lim",None)
        y_min_lim = utils.return_key(kwargs,"y_min_lim",None)
        y_max_lim = utils.return_key(kwargs,"y_max_lim",None)

        ax_to_plot.set_xlim(left=x_min_lim,right=x_max_lim)
        ax_to_plot.set_ylim(bottom=y_min_lim,top=y_max_lim)
        ax_to_plot.legend(utils.return_key(kwargs,"legend",[]),handles=utils.return_key(kwargs,"legend_elements",None))


        if(comparison is not None):
            x_min_lim_comparison = utils.return_key(kwargs,"x_min_lim_comparison",None)
            x_max_lim_comparison = utils.return_key(kwargs,"x_max_lim_comparison",None)
            y_min_lim_comparison = utils.return_key(kwargs,"y_min_lim_comparison",None)
            y_max_lim_comparison = utils.return_key(kwargs,"y_max_lim_comparison",None)

            ax_comparison.set_xlim(left=x_min_lim_comparison,right=x_max_lim_comparison)
            ax_comparison.set_ylim(bottom=y_min_lim_comparison,top=y_max_lim_comparison)
        plt.gcf().tight_layout()





    def plot_2d_pk(self,bin_edges,**kwargs):
        style = utils.return_key(kwargs,"style",None)
        if style is not None:
            plt.style.use(style)
        comparison = utils.return_key(kwargs,"comparison",None)
        k_multiplication = utils.return_key(kwargs,"k_multiplication",False)
        (ax_to_plot,ax_comparison) = self.prepare_axes(kwargs)
        self.put_label(ax_to_plot,
                       xunit = utils.return_key(kwargs,"x_unit",True),
                       yunit = utils.return_key(kwargs,"y_unit",True),
                       x_label = utils.return_key(kwargs,"x_label","k"),
                       y_label = utils.return_key(kwargs,"y_label","P"),
                       labelsize_x = utils.return_key(kwargs,"labelsize_x",12),
                       labelsize_y = utils.return_key(kwargs,"labelsize_y",12),
                       fontsize = utils.return_key(kwargs,"fontsize",12))
        if(comparison is not None):
            self.put_label(ax_comparison,
                           xunit = utils.return_key(kwargs,"x_unit_comparison",True),
                           yunit = utils.return_key(kwargs,"y_unit_comparison",True),
                           x_label = utils.return_key(kwargs,"x_label_comparison","k"),
                           y_label = utils.return_key(kwargs,"y_label_comparison","Pk"),
                           labelsize_x = utils.return_key(kwargs,"labelsize_x_comparison",12),
                           labelsize_y = utils.return_key(kwargs,"labelsize_y_comparison",12),
                           fontsize = utils.return_key(kwargs,"fontsize_comparison",12))

        for i in range(len(bin_edges)):
            mask = self.k_array[1] == bin_edges[i]
            c = kwargs["color"][i] if "color" in kwargs.keys() else None
            ls = utils.return_key(kwargs,"linestyle",["-" for i in range(len(bin_edges))])[i]
            if(k_multiplication): factor_multiplication = self.k_array[0][mask]**3 / (2 * np.pi**2)
            else: factor_multiplication = 1
            if(comparison is not None):
                error_bar_comparison = utils.return_key(kwargs,"error_bar_comparison",True)
                mask_comparison = comparison.k_array[1] == bin_edges[i]
                power_array_comparison = interp1d(self.k_array[0][mask],
                                                  self.power_array[mask],
                                                  bounds_error=False,
                                                  fill_value=np.NaN)(comparison.k_array[0][mask_comparison])
                if((self.error_array is not None)&(comparison.error_array is not None)&error_bar_comparison):
                    error_array_comparison = interp1d(self.k_array[0][mask],
                                                      self.error_array[mask],
                                                      bounds_error=False,
                                                      fill_value=np.NaN)(comparison.k_array[0][mask_comparison])
                    ax_comparison.errorbar(comparison.k_array[0][mask_comparison],
                                           (comparison.power_array[mask_comparison] - power_array_comparison)/comparison.power_array[mask_comparison],
                                           (power_array_comparison/comparison.power_array[mask_comparison]) * np.sqrt((comparison.error_array[mask_comparison]/comparison.power_array[mask_comparison])**2 +  (error_array_comparison/power_array_comparison)**2),
                                           marker = utils.return_key(kwargs,"ps",None),
                                           linestyle= ls,
                                           color=c)
                    ax_comparison.plot([np.min(comparison.k_array[0][mask_comparison]),np.max(comparison.k_array[0][mask_comparison])],
                                       [0,0],"k-",alpha=0.5)
                else:
                    ax_comparison.plot(comparison.k_array[0][mask_comparison],
                                       (comparison.power_array[mask_comparison] - power_array_comparison)/comparison.power_array[mask_comparison],
                                       marker = utils.return_key(kwargs,"ps",None),
                                       linestyle= ls,
                                       color=c)
                    ax_comparison.plot([np.min(comparison.k_array[0][mask_comparison]),np.max(comparison.k_array[0][mask_comparison])],
                                       [0,0],"k-",alpha=0.5)
            if(self.error_array is not None):
                ax_to_plot.errorbar(self.k_array[0][mask],
                                    self.power_array[mask]*factor_multiplication,
                                    self.error_array[mask]*factor_multiplication,
                                    marker = utils.return_key(kwargs,"ps",None),
                                    linestyle= ls,
                                    color=c)
            else:
                ax_to_plot.plot(self.k_array[0][mask],
                                self.power_array[mask]*factor_multiplication,
                                marker = utils.return_key(kwargs,"ps",None),
                                linestyle= ls,
                                color=c)
        xscale = utils.return_key(kwargs,"xscale","log")
        yscale = utils.return_key(kwargs,"yscale","log")

        ax_to_plot.set_xscale(xscale)
        ax_to_plot.set_yscale(yscale)

        x_min_lim = utils.return_key(kwargs,"x_min_lim",None)
        x_max_lim = utils.return_key(kwargs,"x_max_lim",None)
        y_min_lim = utils.return_key(kwargs,"y_min_lim",None)
        y_max_lim = utils.return_key(kwargs,"y_max_lim",None)

        ax_to_plot.set_xlim(left=x_min_lim,right=x_max_lim)
        ax_to_plot.set_ylim(bottom=y_min_lim,top=y_max_lim)

        if(comparison is not None):
            x_min_lim_comparison = utils.return_key(kwargs,"x_min_lim_comparison",None)
            x_max_lim_comparison = utils.return_key(kwargs,"x_max_lim_comparison",None)
            y_min_lim_comparison = utils.return_key(kwargs,"y_min_lim_comparison",None)
            y_max_lim_comparison = utils.return_key(kwargs,"y_max_lim_comparison",None)

            ax_comparison.set_xlim(left=x_min_lim_comparison,right=x_max_lim_comparison)
            ax_comparison.set_ylim(bottom=y_min_lim_comparison,top=y_max_lim_comparison)



        ax_to_plot.legend(utils.return_key(kwargs,"legend",[]),handles=utils.return_key(kwargs,"legend_elements",None))
        ax_to_plot.set_title(utils.return_key(kwargs,"title","Power spectrum"))
        plt.gcf().tight_layout()

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




    def save_plot(self,nameout,format_out = "pdf",fig=None):
        if(fig is None):
            fig = plt.gcf()
        fig.savefig(nameout,format=format_out)

    def close_plot(self,fig=None):
        if(fig is None):
            fig = plt.gcf()
        plt.close()

    def open_plot(self):
        fig = plt.figure()
        return(fig)

    def open_subplot(self,x=2, y=1,figsize=(8,6)):
        fig, ax = plt.subplots(x,y, sharex=True,figsize=figsize)
        return(fig)


    def show_plot(self):
        plt.show()



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
    def init_from_gimlet(cls,namefile,specie="unknown",power_weighted=False,error_estimator=None,**kwargs):
        """ Pm(k) gimlet file contains
         - k: edge (higher) of the k bin considered
         - bincount: number of mode (pairs) computed in the bin
         - pwk: power weighted k
         - power: power of the bin"""
        f = np.loadtxt(namefile)
        k, bincount, pwk,power = f[:,0],f[:,1],f[:,2],f[:,3]
        if(error_estimator is not None):
            error = utils.error_estimator(power,model=error_estimator,bin_count=bincount,**kwargs)
        else: error = None
        if (power_weighted) : k_array = pwk
        else : k_array = k
        dimension = "3D"
        h_normalized  = True
        return(cls(dimension,specie,k_array=k_array,power_array=power,error_array=error,file_init=namefile,size_box=None,h_normalized=h_normalized))


    def write_to_gimlet(self,name_out,power_weighted=False):
        if(power_weighted):
            pwk  = self.k_array
            k = np.zeros(self.k_array.shape)
        else:
            pwk  = np.zeros(self.k_array.shape)
            k = self.k_array
        if(self.error_array is not None):
            bincount = self.error_array
        else:
            bincount = np.zeros(self.power_array.shape)
        power = self.power_array
        out =  np.transpose(np.stack([k, bincount, pwk,power]))
        np.savetxt(name_out,out)



class FluxPowerSpectrum(PowerSpectrum):


    def __init__(self,dimension,**kwargs):
        super(FluxPowerSpectrum,self).__init__(**kwargs)
        if dimension not in ["1D","3D"]: raise KeyError("Dimension of spectrum not available, please choose between 1D and 3D")
        self.dimension = dimension


    @classmethod
    def init_from_gimlet(cls,
                         namefile,
                         type_file,
                         kmu=True,
                         power_weighted=False,
                         error_estimator=None,
                         field_name=None,
                         error_stored = False,
                         **kwargs):
        if(type_file == "txt"):
            pk_array = np.loadtxt(namefile)
        elif(type_file == "hdf5"):
            file = h5py.File(namefile,"r")[field_name]
            pk_array = np.array(list(zip(*file))).transpose()
        if(kmu): (k1_edge, k2_edge, bincount, pwk1, pwk2, power) = cls.init_kmu(pk_array)
        else: (k1_edge, k2_edge, bincount, pwk1, pwk2, power) = cls.init_kperpar(pk_array)
        if (power_weighted) : k1_array , k2_array = pwk1,pwk2
        else : k1_array , k2_array = k1_edge,k2_edge
        if(error_estimator is not None)&(not(error_stored)):
            error = utils.error_estimator(power,model=error_estimator,bin_count=bincount,**kwargs)
        elif(error_stored):
            error = bincount
        else: error = None
        k_array = np.stack([k1_array, k2_array])
        dimension = "3D"
        h_normalized  = True
        return(cls(dimension,
                   k_array=k_array,
                   power_array=power,
                   error_array=error,
                   file_init=namefile,
                   size_box=None,
                   h_normalized=h_normalized))


    @classmethod
    def init_kmu(cls,pk_array):
        """ Pf(k,mu) gimlet file contains
         - k_edge: edge (higher) of the k bin considered
         - mu_edge: edge (lower) of the mu bin considered (mu positive)
         - bincount: number of mode (pairs) computed in the bin
         - pwk: power weighted k
         - pwmu: power weighted mu
         - power: power of the bin"""
        k_edge, mu_edge, bincount, pwk, pwmu,power = pk_array[:,0],pk_array[:,1],pk_array[:,2],pk_array[:,3],pk_array[:,4],pk_array[:,5]
        return(k_edge, mu_edge, bincount , pwk, pwmu, power)



    @classmethod
    def init_kperpar(cls,pk_array):
        """ Pf(kperp,kpar) gimlet file contains
         - k_perp: edge (higher) of the k perp bin considered
         - k_par: edge (higher) of the k par bin considered
         - bincount: number of mode (pairs) computed in the bin
         - pwkperp: power weighted k perp
         - pwkpar: power weighted k par
         - power: power of the bin"""
        k_perp, k_par, bincount, pwkperp, pwkpar,power = pk_array[:,0],pk_array[:,1],pk_array[:,2],pk_array[:,3],pk_array[:,4],pk_array[:,5]
        return(k_perp, k_par, bincount, pwkperp, pwkpar,power)



    def write_to_gimlet(self,name_out,power_weighted=False):
        if(power_weighted):
            pwk1  = self.k_array[0]
            pwk2  = self.k_array[1]
            k1_edge = np.zeros(self.k_array[0].shape)
            k2_edge = np.zeros(self.k_array[1].shape)
        else:
            pwk1  = np.zeros(self.k_array[0].shape)
            pwk2  = np.zeros(self.k_array[1].shape)
            k1_edge = self.k_array[0]
            k2_edge = self.k_array[1]
        if(self.error_array is not None):
            bincount = self.error_array
        else:
            bincount = np.zeros(self.power_array.shape)
        power = self.power_array
        out =  np.transpose(np.stack([k1_edge, k2_edge, bincount, pwk1, pwk2, power]))
        np.savetxt(name_out,out)



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
    reference_spectrum.save_plot(name_out)
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
    reference_spectrum.save_plot(name_out)
    reference_spectrum.close_plot()




def compute_k_extremums(power_Ll,power_Sl,power_Ss,size_small,size_large,N_small,N_large):

    mu_bins = np.unique(power_Ll.k_array[1])
    k_max = []
    for i in range(len(mu_bins)):
        mask_mu_Ll = power_Ll.k_array[1]  == mu_bins[i]
        mask_mu_Sl = power_Sl.k_array[1]  == mu_bins[i]
        mask_mu_Ss = power_Ss.k_array[1]  == mu_bins[i]

        min_k = np.min(power_Ss.k_array[0][mask_mu_Ss])
        max_k = np.max(power_Ss.k_array[0][mask_mu_Ss])

        mask_k_Ll = mask_mu_Ll & (power_Ll.k_array[0] >=min_k) & (power_Ll.k_array[0] <= max_k)


        power_Sl_interp = interp1d(power_Sl.k_array[0][mask_mu_Sl],power_Sl.power_array[mask_mu_Sl])
        mask_k_select = np.abs((power_Ll.power_array[mask_k_Ll]-power_Sl_interp(power_Ll.k_array[0][mask_k_Ll]))/power_Ll.power_array[mask_k_Ll]) < 0.01

        k_max.append(np.max(power_Ll.k_array[0][mask_k_Ll][mask_k_select]))

    return(k_max)


def compute_k_extremums_1D(power_Ll,power_Sl,power_Ss):


    min_k = np.min(power_Ss.k_array)
    max_k = np.max(power_Ss.k_array)

    mask_k_Ll = (power_Ll.k_array >=min_k) & (power_Ll.k_array <= max_k)


    power_Sl_interp = interp1d(power_Sl.k_array,power_Sl.power_array)
    mask_k_select = np.abs((power_Ll.power_array[mask_k_Ll]-power_Sl_interp(power_Ll.k_array[mask_k_Ll]))/power_Ll.power_array[mask_k_Ll]) < 0.01

    k_max = np.max(power_Ll.k_array[mask_k_Ll][mask_k_select])

    return(k_max)


def splice_1D(power_Ll,power_Sl,power_Ss,size_small,size_large,N_small,N_large,use_nyquist=False):
    """ L,S = Large or Small size
        l,s = large or small number of particles/resolution elements
        splice the Ll box, using resolved Sl box and splicing Ss box """
    kmin_S = 2*np.pi / size_small
    if(use_nyquist):
        knyq_L = N_large*np.pi / size_large
        k_max = knyq_L/4
    else:
        k_max = compute_k_extremums_1D(power_Ll,power_Sl,power_Ss)
    power = []
    k_array = []
    error = None
    if((power_Ll.error_array is not None)&(power_Sl.error_array is not None)):
        error = []
    power_Ll_interp = interp1d(power_Ll.k_array,power_Ll.power_array)
    power_Ss_interp = interp1d(power_Ss.k_array,power_Ss.power_array)
    power_Sl_interp = interp1d(power_Sl.k_array,power_Sl.power_array)

    ## low k:  k <= kminS
    mask_k_Ll = power_Ll.k_array <= kmin_S
    power.append(power_Ll.power_array[mask_k_Ll] * (power_Sl_interp(kmin_S)/power_Ss_interp(kmin_S)))
    k_array.append(power_Ll.k_array[mask_k_Ll])
    if(error is not None):
        error.append(power_Ll.error_array[mask_k_Ll])

    ## mid k:  kminS < k <= kNyqL / 4
    mask_k_Ll = (power_Ll.k_array > kmin_S) &  (power_Ll.k_array <= k_max)
    k = power_Ll.k_array[mask_k_Ll]
    power.append(power_Ll.power_array[mask_k_Ll] * (power_Sl_interp(k)/power_Ss_interp(k)))
    k_array.append(k)
    if(error is not None):
        error.append(power_Ll.error_array[mask_k_Ll])

    ## large k:  k > kNyqL / 4
    mask_k_Sl = power_Sl.k_array > k_max
    power.append(power_Sl.power_array[mask_k_Sl] * (power_Ll_interp(k_max)/power_Ss_interp(k_max)))
    k_array.append(power_Sl.k_array[mask_k_Sl])
    if(error is not None):
        error.append(power_Sl.error_array[mask_k_Sl])
        error = np.concatenate(error,axis=0)
    power = np.concatenate(power,axis=0)
    k_array = np.concatenate(k_array,axis=0)



    power_spectrum = MatterPowerSpectrum("1D",
                                         "matter",
                                         k_array=k_array,
                                         power_array=power,
                                         error_array=error,
                                         h_normalized=True)
    return(power_spectrum)






def splice_3D(power_Ll,power_Sl,power_Ss,size_small,size_large,N_small,N_large,use_nyquist=False):
    """ L,S = Large or Small size
        l,s = large or small number of particles/resolution elements
        splice the Ll box, using resolved Sl box and splicing Ss box """
    kmin_S = 2*np.pi / size_small
    if(use_nyquist):
        knyq_L = N_large*np.pi / size_large
        k_max = knyq_L/4
    else:
        k_max_array = compute_k_extremums(power_Ll,power_Sl,power_Ss,size_small,size_large,N_small,N_large)
        k_max = np.min(k_max_array)
    mu_value = []
    power = []
    k_array = []
    error = None
    if((power_Ll.error_array is not None)&(power_Sl.error_array is not None)):
        error = []
    mu_bins = np.unique(power_Ll.k_array[1])
    for i in range(len(mu_bins)):
        power_mu, k_array_mu = [],[]
        if(error is not None):
            error_mu = []
        mask_mu_Ll = power_Ll.k_array[1]  == mu_bins[i]
        mask_mu_Sl = power_Sl.k_array[1]  == mu_bins[i]
        mask_mu_Ss = power_Ss.k_array[1]  == mu_bins[i]
        mu_value.append(power_Ll.k_array[1][mask_mu_Ll].mean())
        power_Ll_interp = interp1d(power_Ll.k_array[0][mask_mu_Ll],power_Ll.power_array[mask_mu_Ll])
        power_Ss_interp = interp1d(power_Ss.k_array[0][mask_mu_Ss],power_Ss.power_array[mask_mu_Ss])
        power_Sl_interp = interp1d(power_Sl.k_array[0][mask_mu_Sl],power_Sl.power_array[mask_mu_Sl])

        ## low k:  k <= kminS
        mask_k_Ll = power_Ll.k_array[0][mask_mu_Ll] <= kmin_S
        power_mu.append(power_Ll.power_array[mask_mu_Ll][mask_k_Ll] * (power_Sl_interp(kmin_S)/power_Ss_interp(kmin_S)))
        k_array_mu.append(power_Ll.k_array[0][mask_mu_Ll][mask_k_Ll])
        if(error is not None):
            error_mu.append(power_Ll.error_array[mask_mu_Ll][mask_k_Ll])

        ## mid k:  kminS < k <= kNyqL / 4
        mask_k_Ll = (power_Ll.k_array[0][mask_mu_Ll] > kmin_S) &  (power_Ll.k_array[0][mask_mu_Ll] <= k_max)
        k = power_Ll.k_array[0][mask_mu_Ll][mask_k_Ll]
        power_mu.append(power_Ll.power_array[mask_mu_Ll][mask_k_Ll] * (power_Sl_interp(k)/power_Ss_interp(k)))
        k_array_mu.append(k)
        if(error is not None):
            error_mu.append(power_Ll.error_array[mask_mu_Ll][mask_k_Ll])

        ## large k:  k > kNyqL / 4
        mask_k_Sl = power_Sl.k_array[0][mask_mu_Sl] > k_max
        power_mu.append(power_Sl.power_array[mask_mu_Sl][mask_k_Sl] * (power_Ll_interp(k_max)/power_Ss_interp(k_max)))
        k_array_mu.append(power_Sl.k_array[0][mask_mu_Sl][mask_k_Sl])
        if(error is not None):
            error_mu.append(power_Sl.error_array[mask_mu_Sl][mask_k_Sl])
            error_mu = np.concatenate(error_mu,axis=0)
        power_mu = np.concatenate(power_mu,axis=0)
        k_array_mu = np.concatenate(k_array_mu,axis=0)

        power.append(power_mu)
        k_array.append(k_array_mu)
        if(error is not None):
            error.append(error_mu)


    power_spliced = np.array([power[i][j]
                              for j in range(len(power[i]))
                              for i in range(len(mu_value))])
    k_spliced = np.transpose(np.array([[k_array[i][j],mu_value[i]]
                                        for j in range(len(power[i]))
                                        for i in range(len(mu_value))]))
    error_spliced = None
    if(error is not None):
        error_spliced = np.array([error[i][j]
                                  for j in range(len(error[i]))
                                  for i in range(len(mu_value))])
    power_spectrum = FluxPowerSpectrum("3D",
                                       k_array=k_spliced,
                                       power_array=power_spliced,
                                       error_array=error_spliced,
                                       file_init=None,
                                       size_box=None,
                                       h_normalized=True)
    return(power_spectrum)
