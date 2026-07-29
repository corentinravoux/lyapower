# lyapower

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21653030.svg)](https://doi.org/10.5281/zenodo.21653030)

Fit the **3D flux power spectrum P₃D(k, μ)** of the Lyman-alpha forest measured
from **Nyx** hydrodynamical simulations (post-processed by
[gimlet](https://github.com/DrJRyan/gimlet)), implementing the non-linear models
of **Arinyo-i-Prats et al. 2015** and the **splicing** method for combining
simulation boxes of different resolution.

📖 **Full documentation (API reference + user guide):** **<https://corentinravoux.github.io/lyapower/>**

---

## What it does

```
Nyx simulation ─▶ gimlet (P3D / P1D / P(k))
               ─▶ power_spectra   read / rebin / splice
               ─▶ CLASS           linear matter P(k)
               ─▶ fitter          Arinyo model fit (iminuit)  ─▶ b, beta, non-linear params
```

## Install

```bash
pip install .              # core deps: fitsio, numpy, scipy, matplotlib, iminuit
pip install -e ".[dev]"    # + pytest, to run the test suite
```

Optional extras: `classy` / `cosmoprimo` (linear `P(k)` in [`CLASS`](lyapower/CLASS.py)),
`h5py` (gimlet HDF5 outputs and Nyx-field rebinning). The `data_example/`
directory ships a matter P₃D and a flux P(k, μ) gimlet file.

## Quick start

```python
from lyapower import power_spectra

# 3D flux power spectrum P(k, mu) from a gimlet text file
flux = power_spectra.FluxPowerSpectrum.init_3D_from_gimlet(
    "data_example/plt00216z_flux_pkmu.txt", type_file="txt", kmu=True,
)
k, mu = flux.k_array          # stacked [k, mu], k in h/Mpc, mu = k_par/k
power = flux.power_array

# matter P(k) from a gimlet file
matter = power_spectra.MatterPowerSpectrum.init_from_gimlet(
    "data_example/plt00216rhom_ps3d.txt", specie="matter",
)
```

## Package map

| Module | Role |
|---|---|
| `power_spectra` | `PowerSpectrum` / `MatterPowerSpectrum` / `FluxPowerSpectrum` objects; gimlet / genpk / ascii I/O; rebinning; splicing (`splice_1D`, `splice_3D`) |
| `fitter` | The Arinyo non-linear model (`D0`, `D1`, BAO variants), the flux model `Pf_model`, iminuit cost functions and fit drivers, gimlet readers and plots |
| `CLASS` | Linear matter power spectra and transfer functions via CLASS (`MyClass`) or cosmoprimo (`CosmoprimoInterface`) |
| `utils_fitter` | Logging, pluggable power-spectrum error models, and n-dimensional / on-disk array rebinning |

`launchers/` holds edit-then-run drivers (fit, plot, splice, BAO, compute-mean).

## Conventions

- Wavenumber `k` in **h/Mpc** (`change_k_normalization` converts to/from 1/Mpc);
  `mu = k_parallel / k`.
- Power spectra are `P(k)` (1D) or `P(k, μ)` (3D, with `k_array` stacked as
  `[k, mu]`).
- `b` and `beta` are the linear bias and redshift-space-distortion parameter.

## Documentation

📖 **Online: <https://corentinravoux.github.io/lyapower/>**

Built from the in-code docstrings with
[MkDocs Material](https://squidfunk.github.io/mkdocs-material/) +
[mkdocstrings](https://mkdocstrings.github.io/) and deployed to GitHub Pages on
every push to `master`. Build locally with:

```bash
pip install -r requirements-docs.txt
mkdocs serve
```

## Tests

```bash
pip install -e ".[dev]"
pytest -q
```

See [`test/README.md`](test/README.md) for the layout of the suite.
