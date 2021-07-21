from rsd_fitter import fitter
from matplotlib.lines import Line2D
import glob,os





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


grid_michael = "/local/home/cravoux/Documents/Simulations/RSD_Nyx_GPU/Fitter/Michael_grid/all_model_outputs_copy.hdf5"
grid_michael_field_name = 'fiducial/redshift_2.2/rescale_Fbar_fiducial/3d power kmu/default binning'

legend_grid_michael = "4096_80mpc"


grid_to_plot = ["3","7"]
grid_to_plot = ["2","6"]
grid_to_plot = ["3","6","9"]



h = 0.7

power_weighted = False
error_estimator="uncorrelated"
epsilon = 0.05
rebin = None # #16

# legend = None
# l = [r"0.0 < $|\mu|$ < 0.25",r"0.75 < $|\mu|$ < 1.0","2560_250mpc","2560_125mpc","1280_125mpc","4096_80mpc"]
# legend_elements = [Line2D([0], [0], color='k', lw=1, label=l[0]),
#                    Line2D([0], [0], color='k',linestyle ="--", lw=1, label=l[1]),
#                    Line2D([0], [0], color='C0',linestyle ="-", lw=1, label=l[2]),
#                    Line2D([0], [0], color='C1',linestyle ="-", lw=1, label=l[3]),
#                    Line2D([0], [0], color='C2',linestyle ="-", lw=1, label=l[4]),
#                    Line2D([0], [0], color='C3',linestyle ="-", lw=1, label=l[5]),
#                    Line2D([0], [0], color='C4',linestyle ="-", lw=1, label=l[5]),
#                    Line2D([0], [0], color='C4',linestyle ="-", lw=1, label=l[5])
#                    ]
# mu_bins = [[0.0,0.25],[0.75,1.0]]
# name_out = "p3d_x_z2.2_mu0_mu1.pdf"



plot_michael_grid = True
legend = []
legend_elements = None
mu_bins = [[0.0,0.2]]  # [[0.7,0.8]]

name_out = "p3d_x_z2.2_mu0_large_grid.pdf"  # "p3d_x_z2.2_mu1_large_grid.pdf"

kwargs = {"x_min_lim" : 0.001,
          "x_max_lim" : 100,
          "y_max_lim" : 0.1,
          "y_min_lim" : 10**-7,
          "linestyle" : ["-","--"],
          "k_multiplication" : True,
          }



if __name__ == "__main__":

    if(plot_michael_grid):
        legend.append(legend_grid_michael)
        power_f_6 = fitter.read_pfkmu_hdf5(grid_michael,
                                           grid_michael_field_name,
                                           power_weighted=power_weighted,
                                           error_estimator=error_estimator,
                                           epsilon=epsilon)
        power_f_6.open_plot()
        power_f_6.k_array[0] = power_f_6.k_array[0]/h
        power_f_6.power_array = power_f_6.power_array*(h**3)
        power_f_6.plot_2d_pk(mu_bins,color=[f"C{0}",f"C{0}"],**kwargs)


    for i in range(len(grid_to_plot)):
        color = f"C{i}"
        if(plot_michael_grid):
            color = f"C{i+1}"
        legend.append(eval(f"legend_grid{grid_to_plot[i]}"))
        power_f = fitter.read_pfkmu(eval(f"file_grid{grid_to_plot[i]}"),
                                    power_weighted=power_weighted,
                                    error_estimator=error_estimator,
                                    epsilon=epsilon)


        if((i == 0)&(not(plot_michael_grid))):power_f.open_plot()
        if(i == len(grid_to_plot) - 1):
            power_f.plot_2d_pk(mu_bins,color=[color,color],legend=legend,legend_elements=legend_elements,**kwargs)
        else:
            power_f.plot_2d_pk(mu_bins,color=[color,color],**kwargs)







    power_f.save_plot(name_out)
