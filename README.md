# PyWirewalker

Processing pipelines for the RBR family loggers and Nortek family Doppler Sonars carried onboard the
[Wirewalker](https://www.delmarocean.com/) profiler. Raw instrument files are converted through a staged
L0–L3 archive of self-describing NetCDF products; every derived quantity is reproducible from
the raw data and version-controlled configuration files.

This README summarizes the methods and usage.

---

## Getting the code

The project lives on GitHub at [`missioncreep11/PyWirewalker`](https://github.com/missioncreep11/PyWirewalker).
First install [git](https://git-scm.com/) (check with `git --version`):

- **Windows** — install [Git for Windows](https://git-scm.com/download/win); it includes **Git Bash**,
  a Unix-like shell. Run every command in this README from **Git Bash** or the **Anaconda Prompt**
  (installed with Miniconda/Anaconda) — not `cmd.exe` — so the `bash`-style snippets work as written.
- **macOS** — `xcode-select --install` (or `brew install git`).
- **Linux** — your package manager, e.g. `sudo apt install git`.

Then clone the repository (identical on every platform):

```bash
git clone https://github.com/missioncreep11/PyWirewalker.git
cd PyWirewalker
```

That is everything needed to run the pipelines — continue to **Quick start** below to build the
conda environment (Miniconda/Anaconda works the same on Windows, macOS, and Linux). Pull later
updates with `git pull` from inside the folder. The optional [GitHub CLI](https://cli.github.com/)
(`gh`) simplifies sign-in and pull requests on all three: `gh auth login`.

> On Windows, a couple of Unix-only commands used elsewhere in this README have equivalents: in the
> Anaconda Prompt use `copy` instead of `cp`; `ncdump` is optional (the `xarray` one-liner works
> everywhere). In Git Bash, `cp` and the other snippets work unchanged.

### Installing Python

This project runs on Python via **conda**. The simplest route on any system is to install
[Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install) (a minimal Python +
`conda`), which the *Quick start* then uses to build the `wirewalker` environment for you.

- **Windows** — the Miniconda [Windows installer](https://www.anaconda.com/docs/getting-started/miniconda/install#windows-installation) (adds the "Anaconda Prompt").
- **macOS / Linux (Unix)** — the Miniconda [macOS](https://www.anaconda.com/docs/getting-started/miniconda/install#macos-installation) / [Linux](https://www.anaconda.com/docs/getting-started/miniconda/install#linux-installation) installer, or a system Python from [python.org/downloads](https://www.python.org/downloads/) or your package manager.

### Contributing

Work on a branch and merge through a pull request — never commit straight to `main`:

```bash
git checkout -b short-topic-name                    # a branch for your change
# ... edit code ...
python -m pytest ww_rbr/tests ww_sig1000/tests -q   # keep the tests green
git add -p && git commit -m "Short description of the change"
git push -u origin short-topic-name                 # then open a PR on GitHub (or: gh pr create)
```

A few conventions:

- **Only source, tests, and the `*.example.json` templates are versioned.** Deployment configs,
  data products, notebooks, and figures are gitignored (see `.gitignore`) — they stay local with
  each deployment's data. Never `git add -f` them.
- **Add or update a test** for any change to the numerics (`ww_rbr/tests/`, `ww_sig1000/tests/`),
  and run the full suite before pushing.
- Keep new processing knobs **config-driven** (not hard-coded) and document them in the
  *Configuration* tables below.

---

## Quick start

**Coming from the MATLAB Wirewalker toolbox?** The workflow maps over one-to-one. The constants you
used to edit at the top of a processing script now live in one JSON **config** per instrument, and
each instrument has a single **driver** you run. You point the driver straight at the raw `.rsk`
(CTD) or `.ad2cp` (Doppler Sonar) — there is **no `.mat` export and no `merge_signature` / `sort_file` step**;
DOLfYN reads the binary directly. Output is self-describing **NetCDF** — dimensions, units, and full
provenance live in the file attributes — not `.mat` with loose script constants.

| In the MATLAB toolbox you… | In PyWirewalker you… |
|---|---|
| edit constants at the top of the `.m` script | edit a `config_*.json` (kept next to the deployment's data) |
| export to `.mat`, then `merge_signature` / `sort_file` | nothing — the driver reads the raw `.rsk` / `.ad2cp` directly |
| run the processing `.m` script | run `process_wirewalker_rbr.py` (CTD) / `process_ww_sig1000.py` (Doppler Sonar) |
| get `.mat` files + whatever the script printed | get NetCDF `L1`–`L3`, with the settings recorded in the attributes |

```bash
# 1. one-time: build the isolated `wirewalker` conda environment (like a self-contained toolbox path)
conda env create -f environment.yml
conda activate wirewalker
python -m pytest ww_rbr/tests ww_sig1000/tests -q            # sanity check

# 2. configure — copy a template and edit it for your deployment (this is your "script header")
cp config_ctd.example.json  config_ctd.json                  # CTD
cp config_adcp.example.json config_adcp.json                 # Doppler Sonar
#    edit at minimum: the raw-file path, output_dir, basename, latitude/longitude, mooring
#    CTD  also: atmospheric_pressure_dbar, grid sizes (sampling rate is read from the .rsk)
#    Doppler Sonar also: velocity.motion / attitude / sail   (see "Configuration" below — v3 is recommended)

# 3a. CTD:  raw .rsk   -> L1 (convert) -> L2 (gridded upcasts) -> L3 (regular depth-time grid)
python process_wirewalker_rbr.py --level all                 # or --level 1 / 2 / 3

# 3b. Doppler Sonar: raw .ad2cp -> (depth, cast) velocity and/or turbulence products
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

1. **Tiered dependence.** `L(n)` reads only `L(n−1)` — never the raw file — so a level
   inherits only the variables present beneath it. Provenance is explicit and reprocessing from
   source is never silent.
2. **Configuration, not code, carries deployment specifics.** All deployment- and
   machine-specific settings (paths, metadata, calibration constants, grid parameters) live in a
   JSON config — `config_ctd.json` / `config_adcp.json`. No paths or constants are hard-coded;
   resolved values are written into the output NetCDF attributes for provenance. Every
   deployment-specific value below is a configuration item, named in `code font` where it appears.

The **`(depth, cast)`** L2 convention is shared by the CTD and Doppler Sonar: each column is one upcast,
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
| `process_ww_sig1000.py` | Doppler Sonar driver (`--product velocity/turbulence`) |
| `ww_rbr/` | CTD package (`config`, `rsk`, `levels`, `derive`) + tests |
| `ww_sig1000/` | Doppler Sonar package (`transforms`, `geometry`, `casts`, `velocity`, `motion`, `l2`, `turbulence`, `turb_product`, `index`) + tests |
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
`.rsk` file into an L1–L3 archive.

**Cast detection.** Profiles are detected from the **CTD pressure record itself**, not taken
from the instrument file — a time-aware, vectorised port of the historical MATLAB
`get_upcastRBR.m` (`ww_rbr.rsk.detect_casts`; set `cast_detection.method` to `"ruskin"` to fall
back to the `region` / `regionCast` tables instead). It classifies every sample as rising or
sinking from the sign of the local pressure slope, debounces the flag with a majority vote to
remove false turnarounds, then splits the flag into casts, dropping runs too small
in pressure span to be full profiles. Everything is specified in physical
units — profiling **speed** (`min_slope_dbar_per_s`) and **time** (`slope_window_s`,
`debounce_window_s`) — and the sampling rate is read from the record (median sample interval), so
the same config works across instruments logging at different rates. The record's **time
continuity is checked first**, so any interval longer than `gap_factor` × the median is treated as a drop and
detection runs independently within each continuous segment — no analysis window ever spans a gap.

Each profile splits into a *down cast* (slow ratcheting descent) and an *up cast* (fast buoyant
ascent). Only free ascent **upcasts** are used for the L2 and L3 products.

**L1 — full-resolution conversion.** `build_L1` reads the raw `.rsk` (SQLite) and writes the
full-rate time series in physical units at whatever rate the instrument recorded. The sampling
rate is **read from the record** (median sample interval — the RBR logs at a steady rate), so it is
not a config item; L2 records the value used and its source in the attributes. Measured channels
are **discovered from the instrument metadata tables**, not hard-coded, so auxiliary sensors (backscatter,
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
Each L2 depth cell is binned by its sample time (see the slant note above), so `n_casts`
— the number of distinct upcasts contributing to a time bin — can exceed one even where cadence is
sparse, and a long profile is spread across the bins it truly spans. L3 is the **continuous**
product — that is its purpose — so whole-empty time bins are linearly interpolated across short gaps
(≤ `l3_interp_max_gap_bins`); longer gaps stay NaN, and `n_casts == 0` marks an interpolated bin (so
`where(n_casts > 0)` recovers the observed-only grid). The build reports how sparse the matrix was
**before** interpolation and stores it on the product (`pre_interpolation_matrix_sparsity_percent`,
`pre_interpolation_empty_time_bins_percent`, `matrix_sparsity_percent`, `time_coverage_percent`). 

The squared buoyancy frequency is calculated as

$$N^2 = \frac{g}{\rho_0} \frac{\partial \sigma_0}{\partial z}$$

($z$ positive down, $g$ = `gravity`) from the gridded σ₀ after a NaN-aware boxcar smooth of length
`n2_vertical_smoothing_m` (5 m default).

---

## Doppler Sonar processing — Nortek Signature1000

The Doppler Sonar chain (`ww_sig1000/`, driver `process_ww_sig1000.py`) processes the raw `.ad2cp` into
two `(depth, cast)` L2 products: **velocity** (motion-corrected ENU currents and vertical shear
from the four slant beams) and **turbulence** (spectral and structure-function dissipation ε from
the fifth, high-resolution beam). This chain is a Python port of a MATLAB toolbox
([`modscripps/wirewalker`](https://github.com/modscripps/wirewalker); Zheng, Lucas, Le Boyer,
Northcott, Griffin), implementing the Wirewalker Doppler Sonar velocity method of Zheng et al. (2021) and the
fifth-beam turbulence method behind Northcott et al. (see [References](#references)). 

Note that at the time of puclishing, turbulence calculation on the slant beams is not supported, this will be incorporated on a future release.


### Ingest and streaming

**★ Direct raw ingest.** The original workflow required a Nortek `.mat` export plus bespoke
`sort_file` / `merge_signature` reassembly. The port reads the raw `.ad2cp` directly via
MHKiT/DOLfYN (`from mhkit import dolfyn`), exposing the slant-beam burst, the HR fifth-beam
burst, the IMU/AHRS records, and the instrument geometry — eliminating the manual export step.

**★ Bounded-memory streaming.** Since the data files are often large, routinely exceeding 10 GB, the driver streams the file in ensemble chunks (`chunk`),
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
current), but wrong when the Doppler Sonar carries a large **fixed mounting tilt** — there the vehicle is not aligned with the Doppler Sonar, so the sail term fabricates a spurious per-cast horizontal velocity that spins
with the (rotating) vehicle heading and stripes the section. The **TLC** deployment, whose Doppler Sonar is
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

Every processing value lives in a per-deployment JSON file kept **with the data** (gitignored);
only the `*.example.json` templates are tracked. Copy a template to the working name and edit it:

```bash
cp config_ctd.example.json  config_ctd.json
cp config_adcp.example.json config_adcp.json
```

**Resolution and overrides.** The driver looks for the config via `--config`, then `$WW_CONFIG`
(CTD) / `$WW_ADCP_CONFIG` (Doppler Sonar), then a `config_*.json` in the working directory or repo root; an
ambiguous (relative) path must be confirmed. Paths can be overridden without editing the file:
`$WW_RSK` / `$WW_AD2CP` (raw input) and `$WW_OUTPUT_DIR` (output). For the Doppler Sonar, **any CLI flag
overrides the corresponding config value** (`--boxsize`, `--motion`, `--attitude`, `--start-time`,
…). Product filenames encode the grid parameters (e.g. `_L2_grid0.5m`, `_L3_grid1m_15min`), so
changing a grid size writes a new file rather than overwriting the old one, and the resolved
settings are stamped into every product's NetCDF attributes.

Required keys have no default (the run errors without them); everything else defaults as shown. The
detailed behaviour behind the processing choices is in the CTD and Doppler Sonar sections above — the tables
here are the exhaustive key reference.

### CTD options (`config_ctd.json`)

| Block | Key | Default | Description |
|-------|-----|---------|-------------|
| *(top)* | `rsk_file` | *(required)* | Raw RBR `.rsk` (SQLite). `$WW_RSK` overrides. |
| *(top)* | `output_dir` | *(required)* | Directory for the `L1/ L2/ L3/` products. `$WW_OUTPUT_DIR` overrides. |
| *(top)* | `basename` | *(required)* | Prefix for every product filename. |
| *(top)* | `mooring` | `""` | Deployment name (written to product attributes). |
| *(top)* | `instrument` | `""` | Instrument description (attributes). |
| *(top)* | `latitude` | *(required)* | Deployment latitude — TEOS-10/`gsw` needs it for depth and absolute salinity. |
| *(top)* | `longitude` | *(required)* | Deployment longitude (as above). |
| *(top)* | `atmospheric_pressure_dbar` | `10.1325` | Subtracted from total pressure to get sea pressure. |
| *(top)* | `sampling_hz` | *(derived)* | Sampling rate (Hz). Optional — omitted, it is read from the record (median sample interval); set only to override. |
| *(top)* | `n2_vertical_smoothing_m` | `5.0` | Boxcar length (m) applied to σ₀ before differencing for N². |
| *(top)* | `gravity` | `9.81` | g used in N² = (g/ρ₀) ∂σ₀/∂z. |
| `thermal_mass` | `alpha` | `0.04` | Lueck & Picklo cell-thermal-mass amplitude. |
| `thermal_mass` | `beta_per_s` | `0.1` | Lueck & Picklo inverse time constant (s⁻¹). |
| `thermal_mass` | `gamma` | `1.0` | Overall correction scale (1.0 = full). |
| `grid` | `l2_dz_m` | `0.5` | L2 depth-bin size (m). |
| `grid` | `zmin_m` | `0.0` | L2 grid top (m). |
| `grid` | `zmax_m` | `500.0` | L2 grid bottom (m); set to the deployment's max depth. |
| `grid` | `l3_dz_m` | `1.0` | L3 depth-bin size (m); a multiple of `l2_dz_m`. |
| `grid` | `l3_dt` | `"30min"` | L3 time-bin width — any pandas offset (`"10min"`, `"1h"`, …). |
| `grid` | `l3_interp_max_gap_bins` | `1` | Gap-fill whole-empty L3 time bins across runs up to this many bins. |
| `cast_detection` | `method` | `"pressure"` | `"pressure"` = detect casts from CTD pressure (port of `get_upcastRBR`); `"ruskin"` = reuse the `.rsk` region tables. |
| `cast_detection` | `slope_window_s` | `5.0` | Window (s) for the centred pressure slope that classifies rising vs sinking. |
| `cast_detection` | `debounce_window_s` | `7.5` | Majority-vote window (s) removing brief flips at the apex/nadir turnarounds. |
| `cast_detection` | `min_slope_dbar_per_s` | `0.04` | Profiling-speed threshold (dbar s⁻¹) separating a cast from a dwell. |
| `cast_detection` | `min_span_dbar` | `5.0` | Discard detected runs spanning less than this (surface dwell, telemetry stops). |
| `cast_detection` | `gap_factor` | `4.0` | Split the record at any interval > this × the median (a data drop) and detect within each continuous segment. |

*Notes.* Only **upcasts** reach L2/L3. The sampling rate and time continuity are read from the
record, so `sampling_hz` is optional and a telemetered file with drops is handled automatically
(see *Cast detection*). `l3_dt` should be chosen relative to the cast cadence — wide enough that
most time bins catch a cast, but the wider it is the more a long cast's time slant matters (each
depth is binned by its own sample time). `zmax_m` only sets the grid extent; deeper-than-profiled
bins are simply never populated.

### Doppler Sonar options (`config_adcp.json`)

Instrument geometry (cell size, blanking, ambiguity velocity, sample rate) is **read from the
`.ad2cp`**, never configured.

| Block | Key | Default | Description |
|-------|-----|---------|-------------|
| *(top)* | `ad2cp_file` | *(required)* | Raw Nortek `.ad2cp`. `$WW_AD2CP` overrides. |
| *(top)* | `output_dir` | *(required)* | Directory for the products. `$WW_OUTPUT_DIR` overrides. |
| *(top)* | `basename` | `"adcp"` | Prefix for product filenames. |
| *(top)* | `mooring` | `""` | Deployment name (attributes). |
| *(top)* | `instrument` | `""` | Instrument description (attributes). |
| *(top)* | `latitude` | `null` | Metadata only (ENU comes from the heading, not position). |
| *(top)* | `longitude` | `null` | Metadata only. |
| `velocity` | `boxsize_m` | `1.0` | L2 depth-bin size (m). |
| `velocity` | `z_max_m` | `null` | L2 grid max depth (m); `null` → auto from max pressure. Set to the profiling depth to drop always-empty deep bins. |
| `velocity` | `motion_correct` | `true` | Apply IMU platform-motion correction; `false` leaves raw beam→ENU velocities. |
| `velocity` | `motion` | `"v1"` | Motion model: `v1`/`v2` (frozen legacy) or `v3` (flexible). See *Velocity and shear*. |
| `velocity` | `attitude` | `"ahrs"` | Pitch/roll source (v3): `ahrs`, `reconstructed` (accel-derived), or `auto` (reconstruct only on AHRS-faulted casts). |
| `velocity` | `sail` | `true` | v3 along-wire "sail" correction. Turn **off** for a large fixed mount tilt (see *When to turn the sail correction off*). |
| `velocity` | `bin_average` | `"boxcar"` | Depth-bin estimator: `boxcar` mean, or `notch` (ridged constant + wave-band fit above 60 m, suppressing residual surface-wave contamination). |
| `velocity` | `l3_dz_m` | *(= boxsize)* | L3 depth-bin size (m); omit/`null` → `boxsize_m`. |
| `velocity` | `l3_dt` | `"15min"` | L3 time-bin width (pandas offset). |
| `velocity` | `l3_interp_max_gap_bins` | `1` | L3 whole-empty-bin gap-fill run length. |
| `turbulence` | `dep_res_m` | `3.0` | Dissipation depth-bin resolution (m). |
| `turbulence` | `max_dep_m` | `100.0` | Deepest ε bin (m). |
| `cast` | `kind` | `null` | `up`/`down`/`both`; `null` → per-product default (velocity `both`, turbulence `up`). |
| `cast` | `min_span_dbar` | `40.0` | Discard casts whose pressure span is under this. |
| `cast` | `corr_min` | `50` | Beam-correlation threshold (%); cells below it are masked before averaging. |
| `cast` | `chunk` | `500000` | Ensembles per streaming read (memory vs. speed). |
| `cast` | `start_ensemble` | `0` | First ensemble to process (trims deployment transit). |
| `cast` | `end_ensemble` | `null` | Stop before this ensemble; `null` → end of record. |
| `cast` | `start_time` | `null` | ISO time; first ensemble at/after it. **Overrides** `start_ensemble`. |
| `cast` | `end_time` | `null` | ISO time; stop here. **Overrides** `end_ensemble`. |

*Notes.* The velocity and turbulence products share the `cast` block but default `kind`
differently (velocity keeps both directions; turbulence uses upcasts). Only velocity has an **L3**
(regular depth–time grid, upcasts only, gap-filled) — there is no turbulence L3. For a mount with a
large fixed tilt (e.g. the TLC Doppler Sonar at 25°), set `motion: "v3"`, `attitude: "ahrs"`, `sail: false`;
the sail term is correct only for a leaning wire, not a fixed tilt. `bin_average: "notch"` costs a
little time for 7–17 % less near-surface velocity noise and is identical to `boxcar` below 60 m.

### Example templates

```json
// config_ctd.example.json
{
  "rsk_file": "/path/to/deployment.rsk",
  "output_dir": "/path/to/output",
  "basename": "MOORING_INSTRUMENT",
  "mooring": "MOORING_NAME",
  "instrument": "RBR Concerto SN000000",
  "latitude": 0.0, "longitude": 0.0,
  "atmospheric_pressure_dbar": 10.1325,
  "thermal_mass": { "alpha": 0.04, "beta_per_s": 0.1, "gamma": 1.0 },
  "grid": { "l2_dz_m": 0.5, "zmin_m": 0.0, "zmax_m": 500.0,
            "l3_dz_m": 1.0, "l3_dt": "30min", "l3_interp_max_gap_bins": 1 },
  "cast_detection": { "method": "pressure", "slope_window_s": 5.0, "debounce_window_s": 7.5,
                      "min_slope_dbar_per_s": 0.04, "min_span_dbar": 5.0, "gap_factor": 4.0 },
  "n2_vertical_smoothing_m": 5.0, "gravity": 9.81
}
```

```json
// config_adcp.example.json
{
  "ad2cp_file": "/path/to/deployment.ad2cp",
  "output_dir": "/path/to/output",
  "basename": "MOORING_INSTRUMENT",
  "mooring": "MOORING_NAME",
  "instrument": "Nortek Signature1000 SN000000",
  "latitude": null, "longitude": null,
  "velocity":   { "boxsize_m": 1.0, "z_max_m": null, "motion_correct": true,
                  "motion": "v3", "attitude": "reconstructed", "sail": true,
                  "l3_dz_m": 2.0, "l3_dt": "15min", "l3_interp_max_gap_bins": 1 },
  "turbulence": { "dep_res_m": 3.0, "max_dep_m": 100.0 },
  "cast":       { "kind": null, "min_span_dbar": 40.0, "corr_min": 50, "chunk": 500000 }
}
```

---

## Summary of Doppler Sonar port improvements

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

## Known issues and limitations

- **Salinity spiking (CTD) is not resolved.** Sharp thermal gradients still produce salinity
  spikes from imperfect conductivity–temperature response matching. L2 applies the Lueck & Picklo
  thermal-mass correction and the measured (≈ 0 s) C–T lag, but performs **no** dedicated de-spiking
  or adaptive C–T alignment — so spikes persist through the interfaces where they are worst, and
  bin-averaging only partially smooths them. Treat near-interface salinity and σ₀ (and anything
  derived from them, e.g. N²) with caution. Improving this is open work.
- **Turbulence is fifth-beam (HR) only.** ε is estimated solely from the pulse-coherent vertical
  fifth beam (spectral + structure-function). The **four-beam HR mode is not implemented**, so there
  is no slant-beam dissipation estimate or independent cross-check, and no ε is produced when the
  fifth beam is absent, saturated, or below the correlation threshold.

---

## References

*Wirewalker instrument*

- Rainville, L. & Pinkel, R. (2001). Wirewalker: An autonomous wave-powered vertical profiler. *J. Atmos. Oceanic Technol.*, 18, 1048–1051. `doi:10.1175/1520-0426(2001)018<1048:WAAWPV>2.0.CO;2`
- Pinkel, R., Goldin, M. A., Sun, O. M., Aja, A. A., Bui, M. N. & Hughen, T. (2011). The Wirewalker: A vertically profiling instrument carrier powered by ocean waves. *J. Atmos. Oceanic Technol.*, 28, 426–446. [`doi:10.1175/2010JTECHO805.1`](https://doi.org/10.1175/2010JTECHO805.1).
- Lucas, A. J., Pinkel, R. & Alford, M. (2017). Ocean wave energy for long endurance, broad bandwidth ocean monitoring. *Oceanography*, 30. [`doi:10.5670/oceanog.2017.232`](https://doi.org/10.5670/oceanog.2017.232).

*Velocity and turbulence on the Wirewalker*

- Zheng, B., Lucas, A. J., Pinkel, R. & Le Boyer, A. (2021). Fine-scale velocity measurement on the Wirewalker wave-powered profiler. *J. Atmos. Oceanic Technol.*, 38. [`doi:10.1175/JTECH-D-21-0048.1`](https://doi.org/10.1175/JTECH-D-21-0048.1). *(the velocity / motion-correction method ported here.)*
- Le Boyer, A., Alford, M. H., Couto, N., Goldin, M., Lastuka, S., Goheen, S., Nguyen, S., Lucas, A. J. & Hennon, T. D. (2021). Modular, flexible, low-cost microstructure measurements: The Epsilometer. *J. Atmos. Oceanic Technol.*, 38(3), 657–668. [`doi:10.1175/JTECH-D-20-0116.1`](https://doi.org/10.1175/JTECH-D-20-0116.1).
- Wiles, P. J. et al. (2006). A novel technique for measuring turbulent dissipation. *Geophys. Res. Lett.*, 33, L21608. *(structure-function method.)*
- Northcott, D. et al. (2026). Wirewalker Signature1000 turbulence dataset. Dryad, [`doi:10.5061/dryad.8sf7m0d44`](https://doi.org/10.5061/dryad.8sf7m0d44). *(fifth-beam spectral ε; `ProcessSingleProfile.m` reference implementation.)*

*CTD / thermodynamics*

- Lueck, R. G. & Picklo, J. J. (1990). Thermal inertia of conductivity cells. *J. Atmos. Oceanic Technol.*, 7, 756–768.
- McDougall, T. J. & Barker, P. M. (2011). *Getting started with TEOS-10 and the GSW Oceanographic Toolbox.*
