import matplotlib.pyplot as plt
import numpy as np

z = np.array([2.0,2.2,2.4,2.6,2.8,3.0,3.2])

b = - np.array([91.6e-3,116.3e-3,  0.1456, 0.1784 ,0.2158 ,0.2570 , 0.3046  ])
b_err = [0.7e-3, 1.0e-3 , 0.0013,0.0017  , 0.0022,0.0028, 0.0025      ]
beta = [1.94, 1.89,1.80,1.69  ,1.56 ,  1.44, 1.25   ]
beta_err = [0.05,0.05, 0.05,0.05, 0.05, 0.05,  0.05  ]


log_F_bar = - 0.0025 * (1 + z)**3.7
b_tau = b /log_F_bar
b_tau_err = b_err /log_F_bar


z_arr = np.array([2.2,2.4,2.6,2.8,3.0])
log_F_bar_arr = - 0.0023 * (1 + z_arr)**3.65

b_tau_arr = np.array([0.6249,0.6202,0.6092,0.5962,0.5808])
b_tau_arr_err = np.array([0.0105,0.0100,0.0131,0.0104,0.0098])
b_arr = b_tau_arr * log_F_bar_arr
b_arr_err = b_tau_arr_err * log_F_bar_arr
beta_arr = [1.431,1.409,1.365,1.296,1.204]
beta_arr_err = [0.036,0.052,0.053,0.028,0.032]


z_dr16 = [2.334]
beta_dr16 = [1.669 ]
beta_dr16_err = [0.071]
b_dr16 = [-0.2014]
b_dr16_err = [0.0032]


plt.style.use("/local/home/cravoux/Software/desi_ec/ec_style.mplstyle")



# fig = plt.figure(figsize=(7,6))
# ax = plt.gca()
# ax.grid()
# ax.errorbar(z,beta,b_err,linestyle='None',marker = ".")
#
# label_size = 23
# size = 15
# ax.set_ylabel(r"$\beta_{\alpha}$", fontsize=label_size,color="C0")
# ax.set_xlabel(r"$z$", fontsize=label_size)
# ax.tick_params(axis='y', labelsize=size, labelcolor="C0")
# ax.tick_params(axis='x', labelsize=size)
#
#
# ax2 = ax.twinx()
# ax2.grid()
# ax2.errorbar(z,b,b_err,linestyle='None',marker = ".",color="C3")
# ax2.set_ylabel(r"$b_{\alpha}$", fontsize=label_size,color="C3")
# ax2.tick_params(axis='y', labelsize=size, labelcolor="C3")
#
# plt.tight_layout()
# plt.savefig("b_beta_measurement.pdf",format="pdf")


label_size = 23
size = 16
markersize = 16
legend_size= 17


fig,ax = plt.subplots(1,2,figsize=(13,6))
ax[0].errorbar(z,beta,beta_err,marker = ".",markersize=markersize)
ax[0].errorbar(z_arr,beta_arr,beta_arr_err,marker = ".",markersize=markersize)
ax[0].errorbar(z_dr16,beta_dr16,beta_dr16_err,linestyle='None',marker = ".",markersize=markersize)
ax[0].set_ylabel(r"$\beta_{\alpha}$", fontsize=label_size)
ax[0].set_xlabel(r"$z$", fontsize=label_size)
ax[0].tick_params(axis='y', labelsize=size)
ax[0].tick_params(axis='x', labelsize=size)
ax[0].legend(["This study","AR15","eBOSS DR16"], fontsize=legend_size)


ax[1].errorbar(z,b,b_err,marker = ".",markersize=markersize)
ax[1].errorbar(z_arr,b_arr,b_arr_err,marker = ".",markersize=markersize)
ax[1].errorbar(z_dr16,b_dr16,b_dr16_err,linestyle='None',marker = ".",markersize=markersize)

# ax[1].errorbar(z,b_tau,b_tau_err,marker = ".",markersize=markersize)
# ax[1].errorbar(z_arr,b_tau_arr,b_tau_arr_err,marker = ".",markersize=markersize)

ax[1].set_ylabel(r"$b_{\alpha}$", fontsize=label_size)
ax[1].set_xlabel(r"$z$", fontsize=label_size)
ax[1].tick_params(axis='y', labelsize=size)
ax[1].tick_params(axis='x', labelsize=size)

plt.tight_layout()
plt.savefig("b_beta_measurement_separated.pdf",format="pdf")
