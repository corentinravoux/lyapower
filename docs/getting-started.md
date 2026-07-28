# Getting started

## Installation

```bash
git clone https://github.com/corentinravoux/lyapower.git
cd lyapower
pip install .            # or: pip install -e ".[dev]" for tests
```

Core dependencies (`fitsio`, `numpy`, `scipy`, `matplotlib`, `iminuit`) install
automatically. A few features need extra packages:

- [`CLASS`](api/CLASS.md) uses `classy` (CLASS) and/or `cosmoprimo` for the
  linear power spectrum;
- [`utils_fitter`](api/utils_fitter.md) uses `h5py` to rebin Nyx/HDF5 fields;
- some readers use `fitsio` / `h5py` for gimlet HDF5 outputs.

The `data_example/` directory ships a matter P₃D and a flux P(k, μ) gimlet file
so the readers can be tried immediately.

## Reading a power spectrum

```python
from lyapower import power_spectra

# 3D flux power spectrum P(k, mu) from a gimlet text file
flux = power_spectra.FluxPowerSpectrum.init_3D_from_gimlet(
    "data_example/plt00216z_flux_pkmu.txt", type_file="txt", kmu=True,
)
k, mu = flux.k_array          # stacked [k, mu]
power = flux.power_array

# matter P(k) from a gimlet file
matter = power_spectra.MatterPowerSpectrum.init_from_gimlet(
    "data_example/plt00216rhom_ps3d.txt", specie="matter",
)
```

## Rebinning and cuts

```python
matter.cut_extremum(kmin=0.1, kmax=5.0)     # keep a k range
matter.rebin_arrays(nb_bin=20, operation="mean")
```

## The linear power spectrum

```python
from lyapower.CLASS import CosmoprimoInterface
import numpy as np

cosmo = CosmoprimoInterface(pwd=".", settings={"h": 0.675, "omega_cdm": 0.12})
k = np.logspace(-3, 1, 200)
pk_lin = cosmo.Pl_class_cosmoprimo(k, z=2.4)
```

## Fitting the Arinyo model

The fit combines a linear `P(k)` with the Arinyo non-linear terms (`D0`/`D1`,
optionally with BAO broadening) using iminuit — see [`fitter`](api/fitter.md):

```python
from lyapower import fitter

# the non-linear correction terms are plain functions of (k, mu, params)
d0 = fitter.D0(k=0.5, mu=0.6, k_nl=6.4, a_nl=0.4, k_p=13.0, a_p=1.6,
               k_v0=1.2, a_v0=1.2, k_v1=1.5, a_v1=1.5)
```
