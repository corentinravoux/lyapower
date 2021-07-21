from rsd_fitter import fitter
import glob,os
from rsd_fitter import power_spectra

### Grids

mu_binning = ""
direction = "y"
redshift = "2.0"

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


### Splicing options
name_splicing = "splicing_A"
size_small = 125
size_large = 250
N_small = 1280
N_large = 2560
grid_to_plot = ["2","3","6"]
grid_Ll = file_grid2
grid_Sl = file_grid3
grid_Ss = file_grid6
grid_m_Ll = file_pm_grid2
grid_m_Sl = file_pm_grid3
grid_m_Ss = file_pm_grid6


# name_splicing = "splicing_B"
# size_small = 250
# size_large = 500
# N_small = 1280
# N_large = 2560
# grid_to_plot = ["1","2","5"]
# grid_Ll = file_grid1
# grid_Sl = file_grid2
# grid_Ss = file_grid5
# grid_m_Ll = file_pm_grid1
# grid_m_Sl = file_pm_grid2
# grid_m_Ss = file_pm_grid5

# name_splicing = "splicing_C"
# size_small = 125
# size_large = 500
# N_small = 640
# N_large = 2560
# grid_to_plot = ["1","3","9"]
# grid_Ll = file_grid1
# grid_Sl = file_grid3
# grid_Ss = file_grid9
# grid_m_Ll = file_pm_grid1
# grid_m_Sl = file_pm_grid3
# grid_m_Ss = file_pm_grid9

name_splicing = "splicing_D"
size_small = 125
size_large = 250
N_small = 640
N_large = 1280
grid_to_plot = ["5","6","9"]
grid_verification = "2"
grid_Ll = file_grid5
grid_Sl = file_grid6
grid_Ss = file_grid9
grid_m_Ll = file_pm_grid5
grid_m_Sl = file_pm_grid6
grid_m_Ss = file_pm_grid9
use_nyquist = False
power_weighted = False
error_estimator="uncorrelated"
epsilon = 0.05


### Plot options

mu_bins = [[0.0],[0.25],[0.5],[0.75]]

name_out_splice = f"{name_splicing}_z{redshift}_spliced_spectra.pdf"

kwargs = {"x_min_lim" : 0.001,
          "x_max_lim" : 100,
          "y_max_lim" : 0.1,
          "y_min_lim" : 10**-7,
          "linestyle" : ["-","--"],
          "k_multiplication" : True,
          "y_label" : r"$\Delta(k) = \frac{k^{3} P(k)}{2\pi^{2}}$",
          "y_unit" : False,
          "y_label_comparison" : r"$\frac{\Delta(Spliced) - \Delta(i)}{\Delta(Spliced)}$",
          "y_unit_comparison" : False,
          "y_max_lim_comparison" : 0.20,
          "y_min_lim_comparison" : -0.20,
          "error_bar_comparison" : False
          }
kwargs_pm = {
          "y_max_lim_comparison" : 0.05,
          "y_min_lim_comparison" : -0.05,
          "error_bar_comparison" : False
          }

kwargs2 = {"k_multiplication" : True,
           "y_label" : r"$\Delta(k) = \frac{k^{3} P(k)}{2\pi^{2}}$",
           "y_unit" : False,
          }

os.makedirs("splicing_plots",exist_ok=True)
name_out_splice = os.path.join(os.getcwd(),"splicing_plots",name_out_splice)


base_out = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/GPU_grid"
dir_out_splice = os.path.join(base_out,name_splicing)

os.makedirs(dir_out_splice,exist_ok=True)

if __name__ == "__main__":

    # L : 250 Mpc  # 500 Mpc
    # S : 125 Mpc
    # l : 2560
    # s : 1280  # 640
    power_Ll = fitter.read_pfkmu(grid_Ll,
                                 power_weighted=power_weighted,
                                 error_estimator=error_estimator,
                                 epsilon=epsilon)
    power_Sl = fitter.read_pfkmu(grid_Sl,
                                 power_weighted=power_weighted,
                                 error_estimator=error_estimator,
                                 epsilon=epsilon)
    power_Ss = fitter.read_pfkmu(grid_Ss,
                                 power_weighted=power_weighted,
                                 error_estimator=error_estimator,
                                 epsilon=epsilon)



    power_spliced = power_spectra.splice_3D(power_Ll,
                                            power_Sl,
                                            power_Ss,
                                            size_small,
                                            size_large,
                                            N_small,
                                            N_large,
                                            use_nyquist=use_nyquist)

    splice_file_out = os.path.join(dir_out_splice,f"{name_splicing}_spliced_spectra_z{redshift}.txt")
    power_spliced.write_to_gimlet(splice_file_out,power_weighted=False)


    power_m_Ll = fitter.read_pk(grid_m_Ll,
                                power_weighted=power_weighted,
                                error_estimator=error_estimator,
                                epsilon=epsilon)
    power_m_Sl = fitter.read_pk(grid_m_Sl,
                                power_weighted=power_weighted,
                                error_estimator=error_estimator,
                                epsilon=epsilon)
    power_m_Ss = fitter.read_pk(grid_m_Ss,
                                power_weighted=power_weighted,
                                error_estimator=error_estimator,
                                epsilon=epsilon)


    power_m_spliced = power_spectra.splice_1D(power_m_Ll,
                                              power_m_Sl,
                                              power_m_Ss,
                                              size_small,
                                              size_large,
                                              N_small,
                                              N_large,
                                              use_nyquist=use_nyquist)

    splice_file_pm_out = os.path.join(dir_out_splice,f"{name_splicing}_spliced_matter_power_spectra.txt")
    power_m_spliced.write_to_gimlet(splice_file_pm_out,power_weighted=power_weighted)



    for i in range(len(mu_bins)):
        mu = mu_bins[i]
        name_out = f"{name_splicing}_z{redshift}_mu{mu[0]}.pdf"

        name_out = os.path.join(os.getcwd(),"splicing_plots",name_out)

        kwargs["comparison"] = power_spliced
        legend = []
        for j in range(len(grid_to_plot)):
            color = f"C{j}"
            legend.append(eval(f"legend_grid{grid_to_plot[j]}"))
            power_f = fitter.read_pfkmu(eval(f"file_grid{grid_to_plot[j]}"),
                                        power_weighted=power_weighted,
                                        error_estimator=error_estimator,
                                        epsilon=epsilon)


            if(j == 0):power_f.open_subplot(figsize =(8,8))
            power_f.plot_2d_pk(mu,color=[color,color],**kwargs)


        legend.append("Spliced spectra")
        kwargs["comparison"] = None

        power_spliced.plot_2d_pk(mu,
                                color=[f"C{j+1}",f"C{j+1}"],
                                legend=legend,
                                **kwargs)
        power_spliced.save_plot(name_out)


    for i in range(len(mu_bins)):
        mu = mu_bins[i]
        name_out = f"{name_splicing}_z{redshift}_mu{mu[0]}_verification.pdf"

        name_out = os.path.join(os.getcwd(),"splicing_plots",name_out)

        kwargs["comparison"] = power_spliced
        legend = []
        color = f"C{0}"
        legend.append(eval(f"legend_grid{grid_verification}"))
        power_f = fitter.read_pfkmu(eval(f"file_grid{grid_verification}"),
                                    power_weighted=power_weighted,
                                    error_estimator=error_estimator,
                                    epsilon=epsilon)


        power_f.open_subplot(figsize =(8,8))
        power_f.plot_2d_pk(mu,color=["C0","C0"],**kwargs)


        legend.append("Spliced spectra")
        kwargs["comparison"] = None

        power_spliced.plot_2d_pk(mu,
                                color=["C1","C1"],
                                legend=legend,
                                **kwargs)
        power_spliced.save_plot(name_out)



    l = [r"0.0 < $|\mu|$ < 0.25",r"0.25 < $|\mu|$ < 0.5",r"0.5 < $|\mu|$ < 0.75",r"0.75 < $|\mu|$ < 1.0"]

    power_spliced.open_plot()
    power_spliced.plot_2d_pk([0.0,0.25,0.5,0.75],
                             legend=l,
                             **kwargs2)
    power_spliced.save_plot(name_out_splice)


    legend = []

    kwargs_pm["comparison"] = power_m_spliced

    for j in range(len(grid_to_plot)):
        color = f"C{j}"
        legend.append(eval(f"legend_grid{grid_to_plot[j]}"))
        power_m = fitter.read_pk(eval(f"file_pm_grid{grid_to_plot[j]}"),
                                 power_weighted=power_weighted,
                                 error_estimator=error_estimator,
                                 epsilon=epsilon)


        if(j == 0):power_m.open_subplot(figsize =(8,8))
        power_m.plot_1d_pk(color=color,**kwargs_pm)

    legend.append("Spliced spectra")
    kwargs_pm["comparison"] = None

    name_out_pm = f"{name_splicing}_pm_z{redshift}.pdf"

    name_out_pm = os.path.join(os.getcwd(),"splicing_plots",name_out_pm)

    power_m_spliced.plot_1d_pk(color=f"C{j+1}",legend=legend,**kwargs_pm)
    power_m_spliced.save_plot(name_out_pm)
