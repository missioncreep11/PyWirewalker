# PyWirewalker

Processing pipelines for the RBR family loggers and Nortek family Doppler Sonars carried onboard the
[Wirewalker](https://www.delmarocean.com/) profiler. Raw instrument files are converted through a staged
L0–L3 archive of self-describing NetCDF products; every derived quantity is reproducible from
the raw data and version-controlled configuration files.

This README summarizes the methods and usage.

---

## Quick start

**Coming from the MATLAB Wirewalker toolbox?** The workflow maps over one-to-one. The constants you
used to edit at the top of a processing script now live in one JSON **config** per instrument, and
each instrument has a single **driver** you run. You point the driver straight at the raw `.rsk`
(CTD) or `.ad2cp` (ADCP) — there is **no `.mat` export and no `merge_signature` / `sort_file` step**;
DOLfYN reads the binary directly. Output is self-describing **NetCDF** — dimensions, units, and full
provenance live in the file attributes — not `.mat` with loose script constants.

| In the MATLAB toolbox you… | In PyWirewalker you… |
|---|---|
| edit constants at the top of the `.m` script | edit a `config_*.json` (kept next to the deployment's data) |
| export to `.mat`, then `merge_signature` / `sort_file` | nothing — the driver reads the raw `.rsk` / `.ad2cp` directly |
| run the processing `.m` script | run `process_wirewalker_rbr.py` (CTD) / `process_ww_sig1000.py` (ADCP) |
| get `.mat` files + whatever the script printed | get NetCDF `L1`–`L3`, with the settings recorded in the attributes |

```bash
# 1. one-time: build the isolated `wirewalker` conda environment (like a self-contained toolbox path)
conda env create -f environment.yml
conda activate wirewalker
python -m pytest ww_rbr/tests ww_sig1000/tests -q            # sanity check

# 2. configure — copy a template and edit it for your deployment (this is your "script header")
cp config_ctd.example.json  config_ctd.json                  # CTD
cp config_adcp.example.json config_adcp.json                 # ADCP
#    edit at minimum: the raw-file path, output_dir, basename, latitude/longitude, mooring
#    CTD  also: atmospheric_pressure_dbar, sampling_hz, grid sizes
#    ADCP also: velocity.motion / attitude / sail   (see "Configuration" below — v3 is recommended)

# 3a. CTD:  raw .rsk   -> L1 (convert) -> L2 (gridded upcasts) -> L3 (regular depth-time grid)
python process_wirewalker_rbr.py --level all                 # or --level 1 / 2 / 3

# 3b. ADCP: raw .ad2cp -> (depth, cast) velocity and/or turbulence products
python process_ww_sig1000.py --product velocity
python process_ww_sig1000.py --product turbulence

# 4. peek at a product (self-describing NetCDF):
ncdump -h  <output_dir>/L2/<basename>_L2_*.nc                 # header: dims, vars, units, provenance
python -c "import xarray as xr; print(xr.open_dataset('<path>.nc'))"
```

**One config per deployment, kept with its data.** Copy the template beside the raw file and pass it
explicitly — `--config /path/to/<deployment>/config_adcp.json` (an absolute path avoids the
working-directory ambiguity a bare filename would have). Any single setting is overridable on the
command line, and paths may use `~` or the `WW_RSK` / `WW_AD2CP` / `WW_OUTPUT_DIR` / `WW_CONFIG` /
`WW_ADCP_CONFIG` environment variables. Data products and figures are **not tracked in git** (see
`.gitignore`); rebuild them from the raw files at any time — reprocessing is never silent, since
every product records the exact settings that made it.

---

## Processing architecture (L0–L3)

Each level is a self-describing NetCDF product derived deterministically from the level below
it. Raw instrument files are never modified.

| Level | Content | Character |
|-------|---------|-----------|
| **L0** | Raw instrument file (RBR `.rsk`, Nortek `.ad2cp`) | Vendor binary; read-only, never edited |
| **L1** | Full-resolution converted time series | Physical units, per-sample cast/profile flags |
| **L2** | Gridded per-cast product, dims `(depth, cast)` | QC applied, derived variables, one profile per column |
| **L3** | Regular `(depth, time)` grid | Uniform time base for spectral/temporal analysis |

Two design principles are common to both processing pipelines:

1. **Strict tiered dependence.** `L(n)` reads only `L(n−1)` — never the raw file — so a level
   inherits only the variables present beneath it. Provenance is explicit and reprocessing from
   source is never silent.
2. **Configuration, not code, carries deployment specifics.** All deployment- and
   machine-specific settings (paths, metadata, calibration constants, grid parameters) live in a
   JSON config — `config_ctd.json` / `config_adcp.json`. No paths or constants are hard-coded;
   resolved values are written into the output NetCDF attributes for provenance. Every
   deployment-specific value below is a configuration item, named in `code font` where it appears.

The **`(depth, cast)`** L2 convention is shared by the CTD and ADCP: each column is one upcast,
each row a depth bin. This lets CTD, velocity, and turbulence products be co-registered by cast and
depth without interpolation. Because a buoyant upcast is **slanted in time** (the vehicle samples
the deep bins minutes before the shallow ones), the CTD L2 stores `time` — and `pressure` — as full
2-D `(depth, cast)` fields like every other variable, rather than one timestamp per profile. That
matters when regridding to L3: each depth cell is placed into the time bin of *its own* sample time,
so a slow/deep profile spanning several L3 time steps is spread across them instead of being
collapsed to a single cast time (a real bias for long casts).

---

## Repository layout

| Path | What |
|------|------|
| `process_wirewalker_rbr.py` | CTD driver (`--level 1/2/3/all`) |
| `process_ww_sig1000.py` | ADCP driver (`--product velocity/turbulence`) |
| `ww_rbr/` | CTD package (`config`, `rsk`, `levels`, `derive`) + tests |
| `ww_sig1000/` | ADCP package (`transforms`, `geometry`, `casts`, `velocity`, `motion`, `l2`, `turbulence`, `turb_product`, `index`) + tests |
| `ww_sig1000/validation/` | turbulence reproducibility scripts (not part of the pipeline) |
| `config_ctd.example.json`, `config_adcp.example.json` | generic config templates — copy and edit |
| `docs/processing_report.tex` | formal methods report |
| `environment.yml` | conda environment specification |

Deployment configs (`config_*.json`), analysis notebooks (`*.ipynb`), data products (`*.nc`),
raw files, and figures are gitignored — kept local, not distributed. The MATLAB reference
toolbox (from [`modscripps/wirewalker`](https://github.com/modscripps/wirewalker)) is likewise
kept local under `WW_Velocity_Processing_SWOT/`.

---

## CTD processing — RBR Concerto

The CTD chain (`ww_rbr/`, driver `process_wirewalker_rbr.py`) converts a raw RBR Concerto
`.rsk` file into an L1–L3 archive. Reference deployment: mooring **NOPP-Aleutians** (RBR
Concerto³ S/N 213752, ~5.3 × 10⁷ scans).

**Cast detection.** Profiles are detected from the **CTD pressure record itself**, not taken
from the instrument file — a time-aware, vectorised port of the historical MATLAB
`get_upcastRBR.m` (`ww_rbr.rsk.detect_casts`; set `cast_detection.method` to `"ruskin"` to fall
back to the `region` / `regionCast` tables instead). It classifies every sample as rising or
sinking from the sign of the local pressure slope, debounces the flag with a majority vote to
remove brief flips at the turnarounds, then splits the flag into casts, dropping runs too small
in pressure span to be real (surface dwell, telemetry stops). Everything is specified in physical
units — profiling **speed** (`min_slope_dbar_per_s`) and **time** (`slope_window_s`,
`debounce_window_s`) — and the sampling rate is read from the record (median sample interval), so
the same config works across instruments logging at different rates. The record's **time
continuity is checked first**: an RBR cannot skip samples onboard, but a real-time telemetered
file can drop data, so any interval longer than `gap_factor` × the median is treated as a drop and
detection runs independently within each continuous segment — no analysis window ever spans a gap.
On the TLC gold-standard deployment (11.8 M scans, 8 Hz) this recovers 2116 upcasts, matching the
Ruskin segmentation to within 3 casts (all 2116 within 0.5 s median of a Ruskin upcast).

Each profile splits into a *down cast* (slow ratcheting descent) and an *up cast* (fast buoyant
ascent). Only **upcasts** are gridded — the CTD sits on top of the vehicle and is in its wake
during descent, so downcasts are contaminated; up/down agreement is never used as a metric.

**L1 — full-resolution conversion.** `build_L1` reads the raw `.rsk` (SQLite) and writes the
full-rate time series in physical units at whatever rate the instrument recorded (`sampling_hz`;
2 Hz for the reference deployment, arbitrary in general). Measured channels are **discovered
from the instrument metadata tables**, not hard-coded, so auxiliary sensors (backscatter,
fluorescence, dissolved oxygen, PAR, …) carry through automatically. Every sample is tagged with
`cast_number`, `profile_number`, `cast_direction`. Conductivity is stored raw at L1.

**L2 — gridded upcasts, thermodynamics, de-spiking.** `build_L2` reads the **L1 NetCDF**,
selects upcasts, and bin-averages onto a vertical grid (`l2_dz_m`, `zmin_m`, `zmax_m`; 0.5 m over
0–500 m by default). Two per-cast corrections precede gridding:

- **Conductivity-cell thermal-mass correction** (Lueck & Picklo 1990), using the config
  `thermal_mass` block (defaulting to the RBR `pyRSKtools` values α = 0.04, β = 0.1 s⁻¹, γ = 1.0).
  Fitting τ freely is ill-posed on these profiles, so the fixed parameters are used and recorded.
- **C–T alignment.** The optimal lag was measured at ≈ 0 s (the Concerto's C and T are already
  aligned), so no shift is applied.

TEOS-10 variables are then derived with `gsw`: sea pressure `p = p_total − atmospheric_pressure_dbar`,
depth, practical/absolute salinity, conservative temperature, σ₀, and sound speed (latitude and
longitude from `latitude`/`longitude`). The L2 product carries these plus the bin-averaged
auxiliary channels and `n_obs`, all dimensioned `(depth, cast)` — including `pressure` and `time`,
which are stored as full 2-D fields because sampling within a cast is irregular in both.

**L3 — regular depth–time grid and stratification.** `build_L3` grids L2 onto a regular depth–time
matrix whose spacing is set in the config (`l3_dz_m`, `l3_dt`; 1 m × 30 min for the reference
deployment). Being a regular grid, `time` and `pressure` here collapse to **1-D** vectors (the
constant-Δt axis and the ≈1-D-in-depth pressure), while the data variables are 2-D `(depth, time)`.
Each L2 depth cell is binned by its **own** 2-D sample time (see the slant note above), so `n_casts`
— the number of distinct upcasts contributing to a time bin — can exceed one even where cadence is
sparse, and a long profile is spread across the bins it truly spans. L3 is the **continuous**
product — that is its purpose — so whole-empty time bins are linearly interpolated across short gaps
(≤ `l3_interp_max_gap_bins`); longer gaps stay NaN, and `n_casts == 0` marks an interpolated bin (so
`where(n_casts > 0)` recovers the observed-only grid). The build reports how sparse the matrix was
**before** interpolation and stores it on the product (`pre_interpolation_matrix_sparsity_percent`,
`pre_interpolation_empty_time_bins_percent`, `matrix_sparsity_percent`, `time_coverage_percent`): on
the TLC gold-standard (1 m × 15 min) the raw matrix is 13% empty (0.9% of time bins had no upcast),
which the single-bin gap-fill takes to 12% empty / 100% time coverage. The squared buoyancy
frequency is

$$N^2 = \frac{g}{\rho_0} \frac{\partial \sigma_0}{\partial z}$$

($z$ positive down, $g$ = `gravity`) from the gridded σ₀ after a NaN-aware boxcar smooth of length
`n2_vertical_smoothing_m` (5 m default).

---

## ADCP processing — Nortek Signature1000

The ADCP chain (`ww_sig1000/`, driver `process_ww_sig1000.py`) processes the raw `.ad2cp` into
two `(depth, cast)` L2 products: **velocity** (motion-corrected ENU currents and vertical shear
from the four slant beams) and **turbulence** (spectral and structure-function dissipation ε from
the fifth, high-resolution beam). This chain is a Python port of a MATLAB toolbox
([`modscripps/wirewalker`](https://github.com/modscripps/wirewalker); Zheng, Lucas, Le Boyer,
Northcott, Griffin). Improvements introduced *during the port* are marked **★**.

### Ingest and streaming

**★ Direct raw ingest.** The original workflow required a Nortek `.mat` export plus bespoke
`sort_file` / `merge_signature` reassembly. The port reads the raw `.ad2cp` directly via
MHKiT/DOLfYN (`from mhkit import dolfyn`), exposing the slant-beam burst, the HR fifth-beam
burst, the IMU/AHRS records, and the instrument geometry — eliminating the manual export step.

**★ Bounded-memory streaming.** Deployment files routinely exceed 10 GB (the reference NOPP file
is ~38 GB / 4 × 10⁷ ensembles). The driver streams the file in ensemble chunks (`chunk`),
detecting casts on a rolling pressure buffer and carrying boundary-spanning casts into the next
read, so memory is bounded by chunk size. Optional bounds (`start_time`/`end_time`) trim
deployment/recovery transit. Instrument geometry is read from the file, never configured.

### Velocity and shear (slant beams)

Currents are formed by the standard beam → XYZ → ENU transformation of the four slant beams,
after a correlation mask (`corr_min`), depth-alignment of the tilt-corrected per-beam cells, and
box-averaging onto a depth grid (`boxsize_m`, `z_max_m`; 1 m default). Casts are selected by
direction and pressure span (`kind`, `min_span_dbar`). Three motion-correction models are provided
(`motion`, `motion_correct`, `attitude`, `sail`):

- **v1** — a port of `WWcorr_beam`: bandpass-integrated IMU acceleration plus the `dp/dt` ascent,
  rotated by the AHRS attitude. *(legacy)*
- **v2** — a buoyant-ascent model (`dp/dt` vertical, depth-gain-weighted horizontal motion,
  low-passed accelerometer tilt + tilt-compensated magnetometer heading), immune to AHRS attitude
  faults, with the **sail correction** always on. *(legacy — frozen for reproducibility)*
- **v3** *(recommended)* — the same buoyant-ascent engine with two independent, per-deployment
  flags: the **attitude source** (`attitude`: `ahrs` = the instrument's AHRS solution;
  `reconstructed` = low-passed accelerometer tilt + magnetometer heading) and the **sail
  correction** (`sail`: on/off). Spike handling (interpolate + exclude), the `dp/dt` vertical, and
  the depth-gain-weighted horizontal are as in v2.

**When to turn the sail correction off.** The sail term removes the horizontal velocity the vehicle
gains travelling *along an inclined wire*; its magnitude scales with `sin(tilt)` and its direction
follows the vehicle heading. That is correct when the **whole mooring leans** (drawn over by
current), but wrong when the ADCP carries a large **fixed mounting tilt** — there the vehicle still
ascends vertically, so the sail term fabricates a spurious per-cast horizontal velocity that spins
with the (rotating) vehicle heading and stripes the section. The **TLC** deployment, whose ADCP is
bolted at a fixed **25° tilt**, is processed with `sail: false` for exactly this reason; a
near-vertical instrument on a wire that leans under current keeps it on.

**Motion-immune shear (ported).** The L2 product also carries a beam-differenced vertical shear
(`shearE`, `shearN`), ported from the toolbox's `beamshear`. Centred cell differences
`(v[c+1]−v[c−1])/(z[c+1]−z[c−1])` are formed **along each raw beam before rotation**; anything
common to a ping's cells (platform translation, attitude-error leakage, the sail term) cancels
exactly, so the shear is immune to the whole motion/attitude error family and needs no motion
correction. Each velocity and shear bin carries a Doppler-noise standard error (`_sem`).

### Turbulence — spectral method (HR beam 5)

ε is estimated from the pulse-coherent vertical fifth beam. Per cast: mask low-correlation
samples (`corr_min`); angular-demean/unwrap the wrapped velocity; subtract the deep **stagnation
profile**; de-spike and detrend; form band-averaged vertical-wavenumber spectra over depth bins
(`dep_res_m`, `max_dep_m`) and fit the Kolmogorov model

$$S(k) = N + A \cdot k^{-5/3}, \qquad \varepsilon = (A/C_K)^{3/2}, \quad C_K = 0.53,$$

with $N$ the white noise floor and $A$ the inertial-subrange amplitude.

**★ Correct reference implementation.** The vendored `WWturb_upward.m` had diverged from the
published method and produced a systematic **+0.28 dex** high bias. The port reproduces
`ProcessSingleProfile.m` — the final paper code from Northcott et al.
([`doi:10.5061/dryad.8sf7m0d44`](https://doi.org/10.5061/dryad.8sf7m0d44)) — restoring
band-averaged spectra, the data-dependent wavenumber grid, the stagnation subtraction (absent in
the vendored code), the cutoff-nearest-1 m rule, and removal of an extraneous mask.

**Validation.** Against the published `NortekTurbulenceData.nc` for the full TLC 2023 deployment
(2116 profiles × 67 depths @ 3 m), the ported spectral ε matches to **median offset −0.00 dex,
RMS ≈ 0.08 dex, corr(log ε) = 0.994**, with exactly the published profile count.

### Turbulence — structure-function method

**★ A second, gap-robust estimator.** The port ports the structure-function branch of
`WWturb_upward.m` (present but unused in the original) into a parallel estimate from the same
preprocessed velocity, accumulating

$$D(r) = N + A \cdot r^{2/3}, \qquad \varepsilon = (A/C_{SF})^{3/2},$$

over valid sample pairs only. $C_{SF} = 1.476$ is cross-calibrated to the validated spectral ε on
the high-scattering TLC deployment (literature value ≈ 2.0), so `epsilon` and `epsilon_sf` share
one absolute scale. Its value is specific to **low-scattering water**: the spectral method
zero-fills correlation-masked samples before the FFT (negligible at TLC, ~8 % masked; significant
at NOPP, ~30–45 % masked, where it inflates the apparent noise floor), whereas the structure
function skips gaps rather than filling them.

### ★ Quality framework and diagnostics

The port adds quantitative quality measures absent from the original toolbox:

- **ε–noise-floor coupling.** Real turbulence varies independently of the instrument noise floor,
  so `r(log ε, log N)` is a contamination diagnostic. On validated (high-scattering) TLC it is
  ≈ −0.3 (the intrinsic fit trade-off) with no diel structure; on low-scattering NOPP it rises to
  +0.6…+0.9 with depth, with a near-surface diel ε cycle locked to `N` (r ≈ 0.95) — a biological
  (diel-vertical-migration) contamination. The structure-function estimator removes the gap-driven
  component in mid-water, extending the trustworthy range.
- **Transfer-function shear noise floor.** A flat SEM-derived floor is the *wrong shape* for
  beam-differenced shear: differencing plus the box-average shape the noise as
  $|H(m)|^2 \propto \sin^2(2\pi m c_z)\mathrm{sinc}^2(mL)$, which rises as $m^2$ at low
  wavenumber (differencing removes the mean). Against this floor the shear is signal-dominated to
  ≈ 6 m vertical scales, not the ≈ 50 m a flat floor implies. Rotary (CW/CCW) spectra of the
  complex shear resolve the up- and down-propagating internal-wave field.

### ★ Software engineering

- **Config-driven**, symmetric with the CTD; resolved settings stamped into product attributes.
- **Unit-tested numerics** (transforms, geometry, cast detection, spectral fit, structure-function
  fit) run alongside the CTD tests.
- **Validated** against the published product (scripts under `ww_sig1000/validation/`).
- **Crash-safe output** — products are written to a temp file and atomically renamed, so a failed
  or interrupted write never truncates an existing product and succeeds even while a reader (e.g.
  an open notebook) holds the old file.

---

## Configuration

Copy a tracked template to the working name and edit it:

```bash
cp config_ctd.example.json  config_ctd.json
cp config_adcp.example.json config_adcp.json
```

The real `config_*.json` are gitignored (kept local). Every processing value is set here; the
CTD and ADCP sections above name each key at its point of use.

**CTD** (`config_ctd.example.json`):

```json
{
  "rsk_file": "/path/to/deployment.rsk",
  "output_dir": "/path/to/output",
  "basename": "MOORING_INSTRUMENT",
  "mooring": "MOORING_NAME",
  "instrument": "RBR Concerto SN000000",
  "latitude": 0.0,
  "longitude": 0.0,
  "atmospheric_pressure_dbar": 10.1325,
  "sampling_hz": 2.0,
  "thermal_mass": { "alpha": 0.04, "beta_per_s": 0.1, "gamma": 1.0 },
  "grid": {
    "l2_dz_m": 0.5, "zmin_m": 0.0, "zmax_m": 500.0,
    "l3_dz_m": 1.0, "l3_dt": "30min", "l3_interp_max_gap_bins": 1
  },
  "cast_detection": {
    "method": "pressure",
    "slope_window_s": 5.0, "debounce_window_s": 7.5,
    "min_slope_dbar_per_s": 0.04, "min_span_dbar": 5.0, "gap_factor": 4.0
  },
  "n2_vertical_smoothing_m": 5.0,
  "gravity": 9.81
}
```

**ADCP** (`config_adcp.example.json`): instrument geometry is *not* configured — it is read from
the `.ad2cp`. Record-trim (`cast.start_time`/`end_time`) and the velocity `motion`/`attitude`
selections are optional and default when absent.

```json
{
  "ad2cp_file": "/path/to/deployment.ad2cp",
  "output_dir": "/path/to/output",
  "basename": "MOORING_INSTRUMENT",
  "mooring": "MOORING_NAME",
  "instrument": "Nortek Signature1000 SN000000",
  "latitude": null,
  "longitude": null,
  "velocity":   { "boxsize_m": 1.0, "z_max_m": null, "motion_correct": true, "motion": "v3", "attitude": "reconstructed", "sail": true },
  "turbulence": { "dep_res_m": 3.0, "max_dep_m": 100.0 },
  "cast":       { "kind": null, "min_span_dbar": 40.0, "corr_min": 50, "chunk": 500000 }
}
```

---

## Summary of ADCP port improvements

| Area | Original MATLAB toolbox | PyWirewalker port |
|------|-------------------------|-------------------|
| Ingest | `.mat` export + `sort_file`/`merge_signature` | Direct `.ad2cp` read via DOLfYN; bounded-memory streaming |
| Turbulence (spectral) | `WWturb_upward.m` (+0.28 dex bias) | `ProcessSingleProfile.m` method; validated to −0.00 dex / corr 0.994 |
| Turbulence (2nd estimator) | Present but unused | Structure-function ε in every product; gap-robust in low scattering |
| Motion correction | Single AHRS-based model | Selectable v1/v2/v3 models; v3 has independent attitude-source (`ahrs`/`reconstructed`) and `sail` flags, robust to AHRS faults |
| Quality control | — | ε–noise-floor coupling; transfer-function shear noise floor; rotary spectra |
| Reproducibility | Script constants | Config-driven, unit-tested, provenance in attributes, atomic writes |

Velocity, beam-differenced shear, and the sail correction were **ported** from the original
toolbox, not introduced during the port.

---

## References

- Lueck, R. G. & Picklo, J. J. (1990). Thermal inertia of conductivity cells. *J. Atmos. Oceanic Technol.*, 7, 756–768.
- Northcott, D. et al. (2026). Wirewalker Signature1000 turbulence dataset. Dryad, [`doi:10.5061/dryad.8sf7m0d44`](https://doi.org/10.5061/dryad.8sf7m0d44).
- McDougall, T. J. & Barker, P. M. (2011). *Getting started with TEOS-10 and the GSW Oceanographic Toolbox.*
- Wiles, P. J. et al. (2006). A novel technique for measuring turbulent dissipation. *Geophys. Res. Lett.*, 33, L21608.
