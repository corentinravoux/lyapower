#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from lyapower import fitter
import glob, os
from matplotlib.lines import Line2D


power_in = "./plt00781_fB14meandirection_flux_pkmu.txt"
name = "P3D"
name_out = "./p3d_plot.pdf"

power_comparison_in = "plt00781_fB14x_flux_pkmu.txt"
name_comparison = "Comparison"


power_weighted = False
error_estimator = "uncorrelated"
epsilon = 0.0


mu_bins = [0.0, 0.25, 0.5, 0.75]
legend = [
    r"0.0 < $|\mu|$ < 0.25",
    r"0.25 < $|\mu|$ < 0.5",
    r"0.5 < $|\mu|$ < 0.75",
    r"0.75 < $|\mu|$ < 1.0",
]

legend_elements = [
    Line2D([0], [0], color=f"C{i}", lw=1, label=legend[i]) for i in range(len(legend))
]
legend_elements.append(Line2D([0], [0], color=f"k", ls="-", lw=1, label=name))
legend_elements.append(
    Line2D([0], [0], color=f"k", ls="--", lw=1, label=name_comparison)
)

kwargs = {
    "x_min_lim": 0.01,
    "x_max_lim": 100,
    "y_max_lim": None,
    "y_min_lim": None,
    "k_multiplication": True,
    "y_label": r"$\Delta^2_{\alpha}(k,\mu) = \frac{k^{3} P_{\alpha}(k,\mu)}{2\pi^{2}}$",
    "y_unit": False,
    "style": "~/2_Software/desi_ec/style.mplstyle",
    "fontsize": 18,
    "labelsize_x": 18,
    "labelsize_y": 18,
    "color": [f"C{i}" for i in range(len(mu_bins))],
    "y_label_comparison": r"$\frac{\Delta - \Delta(Comparison)}{\Delta}$",
    "y_unit_comparison": False,
    "y_max_lim_comparison": 0.20,
    "y_min_lim_comparison": -0.20,
    "error_bar_comparison": False,
    "fontsize_comparison": 18,
    "labelsize_x_comparison": 18,
    "labelsize_y_comparison": 18,
    "legend_elements": legend_elements,
}

if __name__ == "__main__":
    power_f = fitter.read_pfkmu(
        power_in,
        power_weighted=power_weighted,
        error_estimator=error_estimator,
        epsilon=epsilon,
    )

    power_f_comparison = fitter.read_pfkmu(
        power_comparison_in,
        power_weighted=power_weighted,
        error_estimator=error_estimator,
        epsilon=epsilon,
    )

    power_f.open_subplot(figsize=(10, 10))
    power_f.plot_2d_pk(mu_bins, **kwargs)

    kwargs.update({"linestyle": ["--" for i in range(len(mu_bins))]})
    kwargs.update({"comparison": power_f})
    power_f_comparison.plot_2d_pk(mu_bins, **kwargs)

    power_f.save_plot(name_out)
