# lyapower test suite

Run from the package root:

```bash
pip install -e ".[dev]"   # installs pytest alongside the runtime deps
pytest -q                 # or: pytest test/
```

The suite exercises the pure model / utility functions plus the gimlet readers
against the two files shipped in `data_example/` (a matter P₃D and a flux
P(k, μ)) — no simulation outputs or network access are needed.

| File | Scope | What it covers |
|---|---|---|
| `test_utils_fitter.py` | unit | `return_key`, `latex_float`, `bin_ndarray` (against its docstring example) and the four power-spectrum error estimators (incl. dispatch + `KeyError` guards) |
| `test_power_spectra.py` | unit + e2e | gimlet matter / flux readers on the bundled data; ascii reader; `center_wavenumbers_1d`, `rebin_arrays`, `cut_extremum` (2D flux) and `change_k_normalization` |
| `test_fitter.py` | unit + e2e | the Arinyo `D0` / `D1` non-linear terms (against their closed forms and the k→0 limit), `rebin_matter_power` (interpolation + out-of-range NaN), `mask_data`, and the `read_pfkmu` gimlet reader |
| `test_imports_smoke.py` | smoke | every shipped module (and launcher) byte-compiles; the four package modules import |

## Notes

- The **CLASS / cosmoprimo** linear-power-spectrum paths and the iminuit fit
  drivers need a Boltzmann engine (`classy`) and are documented but not run
  here; the pure model terms they feed (`D0`, `D1`) are covered directly.
- `cut_extremum` indexes `k_array[0, :]`, so it applies to the 2D (flux)
  spectra; the 1D matter spectrum is exercised through the other transforms.
