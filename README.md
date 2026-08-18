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
| `process_wirewalker_rbr.py` | CTD driver (`--level 1/2/3/all`), reads `config.json` |
| `process_ww_sig1000.py` | ADCP driver (`--product velocity/turbulence`), args on the CLI |
| `ww_rbr/` | CTD package (`config`, `rsk`, `levels`, `derive`) + tests |
| `ww_sig1000/` | ADCP package (`transforms`, `geometry`, `casts`, `velocity`, `motion`, `l2`, `turbulence`, `turb_product`) + tests |
| `ww_sig1000/validation/` | turbulence reproducibility scripts (not part of the pipeline) |
| `WW_Velocity_Processing_SWOT/` | MATLAB reference toolbox — kept local, not in the repo (from [`modscripps/wirewalker`](https://github.com/modscripps/wirewalker)) |
| `config.json` | CTD deployment/machine settings (no paths hardcoded in code) |
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

All deployment- and machine-specific settings live in **`config.json`** — no paths are
hardcoded. Point it at your `.rsk` and an output directory.

```bash
# 1. configure: edit config.json -> rsk_file, output_dir, basename, lat/lon, atm pressure
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

### Configuration (`config.json`)

- **paths** — `rsk_file`, `output_dir` (where `L1/ L2/ L3/` are written), `basename`.
  May use `~`; override with `WW_RSK` / `WW_OUTPUT_DIR`.
- **metadata** — `mooring`, `instrument`, `latitude`, `longitude`, `atmospheric_pressure_dbar`.
- **processing** — `sampling_hz`, `thermal_mass` (α/β/γ), `grid` (L2/L3 bin sizes, gap-fill),
  `n2_vertical_smoothing_m`, `gravity`.

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

Deployment settings live in **`config_adcp.json`** (parallel to the CTD's `config.json`):
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

## Validation (ADCP)

Turbulence ε is validated against the published `NortekTurbulenceData.nc` (Dryad DA_final).
Full TLC 2023 deployment (2116 casts × 67 depths @ `dep-res=3 m`): **exactly 2116 casts**,
median offset **−0.00 dex**, rms ~0.08, corr(log ε) **0.994**; the gold-standard single
profile (`RawProfile2_20231005-1132.nc`) matches to −0.00 dex. Reproduce with the scripts in
[`ww_sig1000/validation/`](ww_sig1000/validation/) (they need the external Dryad data — see
that folder's README). Beware `Turbulence.mat` in the Nortek folder: it is a coarser earlier
product (2263×25); the correct reference is the Dryad `NortekTurbulenceData.nc`.
