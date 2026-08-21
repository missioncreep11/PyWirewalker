# PyWirewalker — Wirewalker CTD + ADCP processing

Processing pipelines for the two instruments carried on a [Wirewalker](https://www.delmarocean.com/)
wave-powered profiling mooring:

- **RBR Concerto CTD** — an L0→L1→L2→L3 archive chain for temperature/salinity/density
  and any auxiliary channels (`ww_rbr/`, driven by `process_wirewalker_rbr.py`).
- **Nortek Signature1000 ADCP** — motion-corrected ENU currents from the 4 slant beams,
  and HR beam-5 spectral turbulent dissipation ε (`ww_sig1000/`, driven by
  `process_ww_sig1000.py`).

The two halves are independent — different instruments, different raw formats, different
drivers — but share one conda environment (`wirewalker`) and the same `(depth, cast)`
gridded-product convention. Both read the raw instrument files directly (RBR `.rsk`,
Nortek `.ad2cp`); the ADCP turbulence ε reproduces the published paper product to
**−0.00 dex / corr 0.994** (see [Validation](#validation-adcp)).

## Repo layout

| Path | What |
|------|------|
| `process_wirewalker_rbr.py` | CTD driver (`--level 1/2/3/all`), reads `config_ctd.json` |
| `process_ww_sig1000.py` | ADCP driver (`--product velocity/turbulence`), args on the CLI |
| `ww_rbr/` | CTD package (`config`, `rsk`, `levels`, `derive`) + tests |
| `ww_sig1000/` | ADCP package (`transforms`, `geometry`, `casts`, `velocity`, `motion`, `l2`, `turbulence`, `turb_product`, `index`) + tests |
| `ww_sig1000/validation/` | turbulence reproducibility scripts (not part of the pipeline) |
| `WW_Velocity_Processing_SWOT/` | MATLAB reference toolbox — kept local, not in the repo (from [`modscripps/wirewalker`](https://github.com/modscripps/wirewalker)) |
| `config_ctd.json` | CTD deployment/machine settings (no paths hardcoded in code) |
| `config_adcp.json` | ADCP deployment settings (paths, metadata, per-product params) |
| `wirewalker_ctd_processing.ipynb`, `wirewalker_ctd_plots.ipynb` | CTD diagnostics |

## Environment

```bash
conda env create -f environment.yml          # creates the `wirewalker` env
conda activate wirewalker
python -m ipykernel install --user --name wirewalker --display-name "Python (wirewalker)"
python -m pytest ww_rbr/tests ww_sig1000/tests -q   # sanity check
```

Data products and figures are **not tracked in git** (see `.gitignore`); rebuild them
from the raw files with the drivers below.

---

# CTD — RBR Concerto

A configurable L0→L1→L2→L3 chain and diagnostic notebook for Wirewalker-mounted RBR
Concerto CTDs. Reference deployment: mooring **NOPP-Aleutians**, RBR Concerto³ S/N 213752,
2025-07 → 2026-05, 2 Hz continuous (~52.8 M scans, max ~518 dbar).

All deployment- and machine-specific settings live in **`config_ctd.json`** — no paths are
hardcoded. Point it at your `.rsk` and an output directory.

```bash
# 1. configure: edit config_ctd.json -> rsk_file, output_dir, basename, lat/lon, atm pressure
#    (paths may use ~; or override with env vars WW_RSK / WW_OUTPUT_DIR / WW_CONFIG)
# 2. build the products  (L1 from the .rsk, then L2, L3 from the level below)
python process_wirewalker_rbr.py --level all               # or --level 1 / 2 / 3
python process_wirewalker_rbr.py --level all --config /path/to/other.json
python process_wirewalker_rbr.py --level 1 --max-casts 50  # quick test subset
# 3. explore: open wirewalker_ctd_plots.ipynb with the "Python (wirewalker)" kernel
```

Each level reads the product below it: `--level 2` reads the L1 NetCDF, `--level 3` reads
L2; both error if the input is missing. L2 builds in ~13 s, L3 in ~2 s.

## Archive levels

Processing chain is strictly **L0 → L1 → L2**: `build_L1` reads the raw `.rsk`; `build_L2`
reads the **L1 NetCDF** (not the `.rsk`), so L2 can only inherit L1's channels.

| Level | File | Content |
|-------|------|---------|
| **L0** | `…_DeploymentData.rsk` | raw RBR binary (SQLite). The original — not modified. |
| **L1** | `L1/…_L1_converted.nc` | full 2 Hz time series, minimal converted set. Every sample tagged with `cast_number`, `profile_number`, `cast_direction` (0=down, 1=up). |
| **L2** | `L2/…_L2_upcast_grid0.5m.nc` | **upcasts only**, derived from L1, bin-averaged to a 0.5 m depth grid (0–500 m, 1000 bins). dims `(depth, cast)`. |
| **L3** | `L3/…_L3_grid1m_30min.nc` | regular **1 m × 30 min** depth × time grid, derived from L2. dims `(depth, time)`. Empty bins NaN. |
| **L3 (interp)** | `L3/…_L3_grid1m_30min_interp.nc` | companion to L3 with single empty 30-min bins linearly interpolated. `n_casts==0` flags filled bins. |

### Variables

The measured channels are **discovered from the `.rsk`** (the `channels` +
`instrumentChannels` tables), not hardcoded, so any deployment's extra sensors are carried
through automatically. Conductivity, the `temp14` C-T cell thermistor and the measured
pressure get the canonical names below (for the TEOS-10 step); every other measured channel
is passed through under a name slugged from its long name (e.g. `backscatter`, `chlorophyll`,
`rhodamine`, `dissolved_o2_concentration`, `irradiance`, `par`, `temperature_2`).

- **L1**: core `conductivity` (mS/cm), `temperature` (°C, `temp14`), `pressure` (dbar,
  total), `depth` (m) + **all other measured channels** (verbatim) + flags. *Raw*
  conductivity (no thermal-mass correction at L1).
- **L2**: `conductivity` (thermal-mass corrected), `temperature`, `practical_salinity`,
  `absolute_salinity`, `conservative_temperature`, `sigma0`, `sound_speed`, the extra
  channels (bin-averaged raw), and `n_obs`.
- **L3**: same variables as L2 (minus `n_obs`) on a regular `(depth, time)` grid, plus
  `buoyancy_frequency_squared` (N², s⁻²) and `n_casts` = upcasts averaged into each time bin.

## Processing notes (CTD)

- **Profiles/casts** reuse Ruskin's instrument-generated detection (`region` / `regionCast`
  tables): each profile splits into a DOWN cast (slow ratcheting descent) and an UP cast
  (fast buoyant ascent). **Only upcasts go to L2** — the CTD sits on top of the Wirewalker
  and is in the vehicle wake on the descent, so downcasts are contaminated.
- **TEOS-10 conversion** (`gsw`, at L2): sea pressure = `pressure` − `atmospheric_pressure_dbar`;
  depth from `gsw.z_from_p`; practical salinity `SP_from_C`; absolute salinity, conservative
  temperature, σ₀, sound speed.
- **Salinity de-spiking** (L2): conductivity-cell thermal-mass correction (Lueck & Picklo
  1990), per upcast, `pyRSKtools` defaults **α = 0.04, β = 0.1 s⁻¹, γ = 1.0**; C-T alignment
  lag = 0 s (measured ≈0). Parameters recorded in the L2 NetCDF attributes.
- **L3 gridding**: vertical 0.5 m → 1 m (adjacent-pair nan-mean); temporal **30 min** bins.
  Empty bins left **NaN — no temporal interpolation** in the primary; the `_interp`
  companion fills single 30-min gaps only.
- **Buoyancy frequency** (L3): `buoyancy_frequency_squared` = (g/ρ₀)·dσ₀/dz after a **5 m**
  nan-aware boxcar vertical smooth; z positive down. Length in `n2_vertical_smoothing_m`.

### Configuration (`config_ctd.json`)

- **paths** — `rsk_file`, `output_dir` (where `L1/ L2/ L3/` are written), `basename`.
  May use `~`; override with `WW_RSK` / `WW_OUTPUT_DIR`.
- **metadata** — `mooring`, `instrument`, `latitude`, `longitude`, `atmospheric_pressure_dbar`.
- **processing** — `sampling_hz`, `thermal_mass` (α/β/γ), `grid` (L2/L3 bin sizes, gap-fill),
  `n2_vertical_smoothing_m`, `gravity`.

`config_ctd.json` is the default (NOPP-Aleutians reference deployment); `config_astral_ctd.json`
is a second tracked example (ASTRAL_1) — process it with
`process_wirewalker_rbr.py --level all --config config_astral_ctd.json`.

The plots notebook (`wirewalker_ctd_plots.ipynb`, loads L1/L2/L3/L3i) produces cast-flag
checks, time–depth sections, T–S diagrams, deployment-mean profiles, N² sections, an
isopycnal-depth series, and a single-isopycnal depth spectrum. Static figures go to `figs/`.

---

# ADCP — Nortek Signature1000

A Python port of the Wirewalker Signature1000 velocity/turbulence toolbox (original MATLAB
from [`modscripps/wirewalker`](https://github.com/modscripps/wirewalker); Zheng/Lucas/Le
Boyer/Northcott — kept locally under `WW_Velocity_Processing_SWOT/`, not distributed here).
Reads the raw `.ad2cp` **directly** via MHKiT/DOLfYN (`from mhkit import dolfyn`) — the Nortek
`.mat` export step is obsolete. The 5-beam head gives two independent products, both streamed
from the raw file in ensemble chunks (casts crossing a chunk boundary are carried), both
`(depth, cast)` NetCDF.

Deployment settings live in **`config_adcp.json`** (parallel to the CTD's `config_ctd.json`):
paths, metadata, and the per-product processing choices. Any CLI flag overrides the config,
and the resolved values are written into the output NetCDF attributes for provenance.
Instrument *geometry* (cell size, blanking, ambiguity velocity, sample rate) is **not** in
the config — it is read straight from the `.ad2cp` at run time.

```bash
# 1. configure: edit config_adcp.json -> ad2cp_file, output_dir, basename, mooring, params
#    (paths may use ~; override with env vars WW_AD2CP / WW_OUTPUT_DIR / WW_ADCP_CONFIG)

# 2. build a product (config-driven; output name derived from basename + grid)
python process_ww_sig1000.py --product velocity      # motion-corrected ENU currents (slant beams)
python process_ww_sig1000.py --product turbulence    # HR beam-5 spectral dissipation eps
python process_ww_sig1000.py --product turbulence --config /path/to/other.json

# 3. or override any config value on the CLI (e.g. a one-off file)
python process_ww_sig1000.py --product velocity \
    --file ww_sig1000/test_data/S101913A013_ASTRAL_1_U.ad2cp \
    --out  out/ASTRAL_1_U_L2.nc --mooring ASTRAL_1_U --boxsize 1.0
```

Config sections (all optional beyond the paths): **top-level** `ad2cp_file`, `output_dir`,
`basename`, `mooring`, `instrument`, `latitude`/`longitude`; **`velocity`** `boxsize_m`,
`z_max_m` (null → auto), `motion_correct`; **`turbulence`** `dep_res_m`, `max_dep_m`;
**`cast`** (shared) `kind` (null → `both` for velocity, `up` for turbulence), `min_span_dbar`,
`corr_min`, `chunk`. The matching CLI overrides are `--file/--out/--mooring`, `--boxsize`,
`--z-max`, `--no-motion`, `--dep-res`, `--max-dep`, `--cast-kind`, `--min-span-dbar`,
`--corr-min`, `--chunk`. With no config file present, built-in defaults apply so the tool
still runs purely from CLI flags.

> **Runtime:** the *first* dolfyn read of a raw `.ad2cp` builds a `.ad2cp.index` sidecar
> next to it (one-time, can take minutes on a large file — this is normal, not a hang).
> With the index in place, a full ~10 GB deployment streams in ~25–30 min; progress prints
> per chunk. The `.index` files are gitignored.

## Products

- **velocity** (`ww_sig1000/{transforms,geometry,casts,velocity,motion,l2}.py`) — beam →
  XYZ → ENU on the 4 slant beams, IMU motion correction, bin-averaged to a `--boxsize`
  depth grid. Reproduces the paper's Fig-9f internal-tide currents.

- **turbulence** (`ww_sig1000/turbulence.py` + streaming assembler `turb_product.py`) —
  spectral dissipation ε from the pulse-coherent HR beam-5. An **exact port of
  `ProcessSingleProfile.m`** (Devon Northcott's final paper code from Dryad
  `doi:10.5061/dryad.8sf7m0d44`), *not* the vendored `WWturb_upward.m` — the vendored script
  diverged and produced a +0.28 dex bias. Pipeline: angular demean/unwrap of the wrapped HR
  velocity → deep-profile (`pressure>50`) stagnation subtraction → cutoff-nearest-1 m →
  despike/detrend → band-averaged wavenumber spectra → Kolmogorov fit `S(k)=N+A·k^(−5/3)`,
  ε = (A/0.53)^(3/2). Output vars: `epsilon`, `N`, `SNR`, `A`, `corr`, `num_spectra`.

## Ambiguous config paths need confirmation

A relative path resolves against the process working directory, so
`--config config_adcp.json` names a different deployment depending on where the shell
happens to be — and the run then completes against the wrong raw file, reporting
success. Both drivers therefore **warn and require agreement** before using a relative
`--config` (or a relative `$WW_ADCP_CONFIG` / `$WW_CONFIG`). The warning names the
mooring and raw file it is about to load, which is what reveals a wrong one:

```
WARNING: --config 'config_adcp.json' is a relative path, resolved against the
  current working directory (/Users/drew/PyWirewalker).
  It resolves to : /Users/drew/PyWirewalker/config_adcp.json
  which contains : TLC_23  ->  ~/TLC/2023/Wirewalker/Nortek/S101913A008_TLC_23.ad2cp
  Proceed with this config? [y/N]
```

Anything but `y`/`yes` aborts, as do Ctrl-C and Ctrl-D. With no terminal attached
(a script, a queued job) there is no way to agree, so the run stops unless `-y` /
`--yes` is passed. An **absolute** path — `~` counts — is unambiguous and never
prompts, so it stays the right habit for anything scripted.

With no `--config`, a `config_<ctd|adcp>.json` in the working directory or at the repo
root is used silently. If **both** exist, the same confirmation applies, showing which
deployment each one names; the working-directory copy is the one offered.

## Trimming deployment/recovery transit

A record usually starts before the mooring is in the water and ends after it comes
back, and that transit can produce pressure excursions large enough to pass
`min_span_dbar` and enter the product as spurious casts. Bound the ensemble range with
`cast.start_time` / `cast.end_time` in the config (ISO times, resolved against the
dolfyn index) or `cast.start_ensemble` / `cast.end_ensemble`; the CLI equivalents are
`--start-time` / `--end-time` / `--start-ensemble` / `--end-ensemble`, and times win
over ensembles. The resolved range is written to the `ensemble_range` attribute.

`ww_sig1000/index.py` reads the `.ad2cp.index` sidecar directly, which also gives the
**exact** ensemble count — previously the streaming builders probed the reader by
bisection, which cost many reads and rounded *down* to a 5000-ensemble tolerance,
silently dropping up to ~10 min of data at 8 Hz. It doubles as a CLI:

```bash
python -m ww_sig1000.index raw.ad2cp                          # build/inspect the index
python -m ww_sig1000.index raw.ad2cp --at-time 2023-09-19T18:30:00
```

`summarize` reports record counts, ensemble count, start/end times and whether beam-5
(record 0x18) is present — useful for deciding up front whether a deployment can
produce the turbulence product at all.

> Note: `nortek2_lib.get_index` defaults to `eof=2**32` and silently stops indexing at
> 4 GB; `ww_sig1000.index` passes the true file length, as dolfyn's own reader does.
> Index dates use the raw year byte (year − 1900) and a **1-based** month — not the
> zero-based month of the data record.

## Duty-cycled (burst) deployments

Cast detection is burst-aware. `casts.detect_bursts` splits the record into contiguous
sampling blocks (any time step > `gap_s`, default 30 s, starts a new burst), and the
low-pass + turning-point detection runs **within each burst independently** — filtering
across a multi-hour gap smears the pressure discontinuity into the burst edges and
orphans samples there. A continuously sampled record is a single burst, so this is a
no-op for one: detection is bit-identical to the pre-burst-aware code (checked against
ASTRAL_1_U, 72 casts, exact match).

Because a burst rarely starts or ends on a profile turning point, casts at burst
boundaries are **clipped** and cover only part of the water column. Every cast carries a
`profile_complete` flag (1 = bounded by turning points, 0 = truncated by a burst edge or
the record edge), alongside `pressure_min` / `pressure_max`, so gridding can weight or
exclude partial profiles; `n_casts_truncated` is in the global attributes. Truncated
casts are still valid velocity data over the depth range they do cover — they are
flagged, not dropped.

On the NOPP1-California record (40-min bursts every 2h20m, 0–500 dbar profiles) each
burst yields ~2 complete profiles bracketed by 2 truncated ones: **48% of casts are
truncated**, and they cover 47% of the depth bins on average against 95% for complete
casts. The streaming reader carries a boundary cast — and its truncation state — across
chunks, so the product is bit-identical regardless of `--chunk`.

## Platform kinematics and AHRS validation

`ww_sig1000/platform.py` reconstructs what the vehicle was doing during a cast — tilt,
rotation, climb, wave heave — and **checks the instrument's orientation solution against
its own raw sensors**. That check is not optional hygiene: on NOPP_d2 the AHRS emitted
physically impossible attitudes on **15.7% of casts**, including 11 casts reporting tilt
past horizontal, while the accelerometer and gyro stayed healthy in every one of them.

```bash
python -m ww_sig1000.platform --config /abs/path/config_adcp.json --cast-kind both
```

writes one row per cast to `ahrs_scan.csv` (~40 min for a 41 GB file), with
`ahrs_error_deg`, `tilt_accel_deg` vs `tilt_ahrs_deg`, `gyro_dps` vs `ahrs_dps`, and
`ahrs_ok`. It needs only `orientmat`, `accel`, `angrt` and `pressure`, so any Signature
record with the AHRS enabled can be screened without rebuilding its velocity product.

`ahrs_error` is the **per-ping** angle between where the AHRS puts earth-up (third row
of `orientmat`) and where the accelerometer measures it, compared in the instrument
frame. Per-ping matters: rotating acceleration into earth and measuring how far the
*mean* lands from vertical collapses when the false attitude is rotating, because the
horizontal residuals average away — that method reports 3.2° for a 43° error.

Why it matters for velocity: an attitude error leaks the platform's own ~0.45 m/s ascent
into the horizontal as `0.45·sin(δθ)`. Doppler beam noise dominates *per ping* by ~100×,
but it averages down as `1/√n_obs` while a correlated attitude error does not — so after
binning, attitude is the term that survives. `classify()` therefore reports
`ahrs_fault` before any attitude-derived label, since none of them mean anything on a
bad solution.

## Validation (ADCP)

Turbulence ε is validated against the published `NortekTurbulenceData.nc` (Dryad DA_final).
Full TLC 2023 deployment (2116 casts × 67 depths @ `dep-res=3 m`): **exactly 2116 casts**,
median offset **−0.00 dex**, rms ~0.08, corr(log ε) **0.994**; the gold-standard single
profile (`RawProfile2_20231005-1132.nc`) matches to −0.00 dex. Reproduce with the scripts in
[`ww_sig1000/validation/`](ww_sig1000/validation/) (they need the external Dryad data — see
that folder's README). Beware `Turbulence.mat` in the Nortek folder: it is a coarser earlier
product (2263×25); the correct reference is the Dryad `NortekTurbulenceData.nc`.
