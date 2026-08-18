# ww_sig1000 turbulence validation

Scratch scripts that reproduce the validation of the ported HR beam-5 spectral
dissipation (`ww_sig1000/turbulence.py`) against the published paper product. They are
**not part of the pipeline** and are kept as a reproducibility record. Each runs
standalone (`python ww_sig1000/validation/<script>.py`) and adds the repo root to
`sys.path` automatically.

They depend on external reference data that is **not in the repo** — Devon Northcott's
Dryad set (`doi:10.5061/dryad.8sf7m0d44`) and the raw TLC 2023 `.ad2cp`. Edit the paths
at the top of each script to point at your copies.

| Script | What it checks |
|--------|----------------|
| `rawprofile2_validate.py` | Single gold-standard profile: runs `process_cast_turbulence` on `RawProfile2_20231005-1132.nc` and compares ε(z) to the nearest `NortekTurbulenceData.nc` profile (the Fig-3 `plot_inds`). Writes `RawProfile2_validation.png`. |
| `turb_window_validate.py` | A window of the raw `.ad2cp`: detects upcasts, runs `process_cast_turbulence` per cast, and reports the bin-wise offset / correlation vs the published product. Args: `START N` (ensemble range). |
| `turb_deployment_validate.py` | The full assembled product (`build_turbulence_streaming` output) vs `NortekTurbulenceData.nc`: bin-wise offset/correlation + a Fig-4c-style ε(time, depth) section. Writes `TLC_turb_deployment.png`. |

Result of record: median offset **−0.00 dex**, rms ~0.08, corr(log ε) **0.994** over the
full TLC 2023 deployment (2116 casts × 67 depths).
