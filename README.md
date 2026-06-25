# PyWirewalker — Wirewalker RBR Concerto CTD processing

A configurable L0→L1→L2→L3 processing chain and diagnostic notebook for Wirewalker-mounted
RBR Concerto CTDs. Reference deployment: mooring **NOPP-Aleutians**, RBR Concerto³ S/N 213752,
2025-07 → 2026-05, 2 Hz continuous (~52.8 M scans, max ~518 dbar).

All deployment- and machine-specific settings live in **`config.json`** — no paths are
hardcoded. Point it at your `.rsk` and an output directory; the data and figures are not
tracked in git.

## Quick start

```bash
# 1. environment
conda env create -f environment.yml          # creates the `wirewalker` env
conda activate wirewalker
python -m ipykernel install --user --name wirewalker --display-name "Python (wirewalker)"

# 2. configure: edit config.json -> rsk_file, output_dir, basename, lat/lon, atm pressure
#    (paths may use ~; or override with env vars WW_RSK / WW_OUTPUT_DIR)

# 3. build the products  (L1 from the .rsk, then L2, L3 from the level below)
python process_wirewalker.py --level all      # or --level 1 / 2 / 3
python process_wirewalker.py --level all --config /path/to/other.json   # another deployment

# 4. explore: open wirewalker_ctd_plots.ipynb with the "Python (wirewalker)" kernel
```

The script and notebook both read `config.json` (found beside the script / by walking up
from the notebook's directory). Outputs (`L1/ L2/ L3/`, `*.nc`) go to `output_dir`; figures
go to `figs/` in the repo. Both are gitignored.

## Archive levels

Processing chain is strictly **L0 → L1 → L2**: `build_L1` reads the raw `.rsk`; `build_L2`
reads the **L1 NetCDF** (not the `.rsk`), so L2 can only inherit L1's channels.

| Level | File | Content |
|-------|------|---------|
| **L0** | `../202605_NOPP-Aleutians_RBR-Concerto-Serial213752_DeploymentData.rsk` | raw RBR binary (SQLite). The original — not modified. |
| **L1** | `L1/NOPP-Aleutians_RBR-Concerto-213752_L1_converted.nc` | full 2 Hz time series, minimal converted set. Every sample tagged with `cast_number`, `profile_number`, `cast_direction` (0=down, 1=up). |
| **L2** | `L2/NOPP-Aleutians_RBR-Concerto-213752_L2_upcast_grid0.5m.nc` | **upcasts only**, derived from L1, bin-averaged to a 0.5 m depth grid (0–500 m, 1000 bins). dims `(depth, cast)`. |
| **L3** | `L3/NOPP-Aleutians_RBR-Concerto-213752_L3_grid1m_30min.nc` | regular **1 m × 30 min** depth × time grid, derived from L2. dims `(depth, time)`, 500 × ~14,634. Empty bins NaN. |
| **L3 (interp)** | `L3/..._L3_grid1m_30min_interp.nc` | companion to L3 with single empty 30-min bins linearly interpolated (coverage 92→99%). Same grid; `n_casts==0` flags filled bins. |

### Variables

- **L1**: `conductivity` (mS/cm), `temperature` (°C, `temp14` C-T cell thermistor),
  `pressure` (dbar, total), `depth` (m) + flags `cast_number`, `profile_number`,
  `cast_direction`. *Raw* conductivity (no thermal-mass correction at L1). The extra
  thermistors (`temp22`, `temp10`) and all derived salinities/density are **not** kept.
- **L2**: `conductivity` (thermal-mass corrected), `temperature`, `practical_salinity`,
  `absolute_salinity`, `conservative_temperature`, `sigma0`, `sound_speed`, `n_obs`.
  No extra thermistors and no sea-pressure (absent by construction, since L2 comes from L1).
- **L3**: same variables as L2 (minus `n_obs`) on a regular `(depth, time)` grid, plus
  `buoyancy_frequency_squared` (N², s⁻²) and `n_casts` = number of upcasts averaged into
  each 30-min time bin.

## Processing notes

- **Profiles/casts** reuse Ruskin's instrument-generated detection (`region` / `regionCast`
  tables): 14,173 profiles, each split into a DOWN cast (slow ratcheting descent) and an
  UP cast (fast buoyant ascent). **Only upcasts go to L2** — the CTD sits on top of the
  Wirewalker and is in the vehicle wake on the descent, so downcasts are contaminated.
- **TEOS-10 conversion** (`gsw`, computed at L2): sea pressure = `pressure` − 10.1325 dbar
  (the `ATMOSPHERE` value stored in the `.rsk`); depth from `gsw.z_from_p`; practical
  salinity `SP_from_C`; absolute salinity, conservative temperature, σ₀, sound speed.
  Lat/lon = 49.5, −159.
- **Salinity de-spiking** (L2): conductivity-cell thermal-mass correction (Lueck & Picklo
  1990), applied per upcast using RBR `pyRSKtools` defaults **α = 0.04, β = 0.1 s⁻¹
  (τ = 10 s), γ = 1.0**. **C-T alignment lag = 0 s** — the optimal lag was measured at
  ≈0 across the deployment (the Concerto's C and T are already aligned). Parameters are
  recorded in the L2 NetCDF attributes.
- **L3 gridding** (from L2): vertical 0.5 m → **1 m** (adjacent-pair nan-mean, centres
  0.5…499.5); temporal **30 min** bins (upcasts bin-averaged by cast time). Native upcast
  cadence is ~31 min, so ≈1 upcast/bin → **Nyquist = 1 cph**. ~88% of bins hold one
  upcast, ~4% two or more, ~8% empty. Empty bins are left **NaN — no temporal
  interpolation** (preserves independence). Parameters in the L3 NetCDF attributes.
- **Buoyancy frequency** (L3): `buoyancy_frequency_squared` = (g/ρ₀)·dσ₀/dz, computed in
  `build_L3` from the gridded σ₀ field after a **5 m** nan-aware boxcar vertical smooth
  (`N2_SMOOTH_M`), z positive down so N²>0 is stable; ρ₀ = 1000+σ₀ (local). The 5 m length
  was chosen from a dσ₀/dz smoothing comparison (1→5 m removes most gradient noise:
  5.7%→3.2% spurious unstable points, while preserving the pycnocline peak). Present in
  both L3 and L3i; smoothing length recorded in the `n2_vertical_smoothing_m` attribute.
- **L3 interp companion**: a second file where empty bins are linearly interpolated in
  time across gaps ≤ `L3_INTERP_MAXGAP` bins (default **1** = single 30-min bins, the
  ~31-min cadence phase slip). Fills 1,017 isolated bins → 99% coverage; real gaps (2–9
  bins) stay NaN. Filled bins are flagged by `n_casts == 0`. Use the NaN version as the
  archival primary; the interp version for plotting / spectra that dislike gaps.

## Rebuild

```bash
python3 process_wirewalker.py --level all          # L1 -> L2 -> L3
python3 process_wirewalker.py --level 1             # just L1 (from the .rsk)
python3 process_wirewalker.py --level 2             # just L2 (requires L1)
python3 process_wirewalker.py --level 3             # just L3 (requires L2)
python3 process_wirewalker.py --level 1 --max-casts 50   # quick test subset
```

Each level reads the product below it: `--level 2` reads the L1 NetCDF, `--level 3` reads
L2; both error if the input is missing. L2 builds in ~13 s, L3 in ~2 s (no SQLite re-query).
Use `--config <file>` (or `$WW_CONFIG`) to process a different deployment.

### Environment

Create the `wirewalker` conda env from `environment.yml` (see **Quick start** above). In
VS Code / Jupyter select the **Python (wirewalker)** kernel.

### Configuration (`config.json`)

All deployment- and machine-specific settings live in `config.json`:

- **paths** — `rsk_file`, `output_dir` (where `L1/ L2/ L3/` are written), `basename`
  (output filename prefix). May use `~`; override with `WW_RSK` / `WW_OUTPUT_DIR` env vars.
- **metadata** — `mooring`, `instrument`, `latitude`, `longitude`,
  `atmospheric_pressure_dbar`.
- **processing** — `sampling_hz`, `thermal_mass` (α/β/γ), `grid` (L2/L3 bin sizes, gap-fill),
  `n2_vertical_smoothing_m`, `gravity`.

## Plots

Open `wirewalker_ctd_plots.ipynb` (loads L1 / L2 / L3 / L3i): L1 vehicle depth-vs-time
cast-flag check (two-week panels + full-res zoom), L2 time–depth sections, individual
profiles, T–S diagram (all data + seasonal), single-profile up-vs-down from the L1 flags
(salinity on the fly with `gsw`), deployment-mean profile, L3 deployment sections
(primary + interpolated + N²), an isopycnal-depth time series, an interactive single-isopycnal
viewer (Plotly, drag-to-zoom in VS Code), and a single-isopycnal depth spectrum
(spike removal → 30-day Welch PSD, 1/30 d–1 cph, with tidal/inertial reference lines).

All static figures are written to `figs/`. The base folder holds only the notebook, the
`process_wirewalker.py` script, the `L1/ L2/ L3/` data dirs, `figs/`, and this README.

### Derived product

`L3/NOPP-Aleutians_RBR-Concerto-213752_isopycnal_depths.nc` — depth of target isopycnals
vs time, dims `(isopycnal, time)`. Targets = deployment-mean σ₀ (from L3i) sampled every
5 m; depth at each time step by inverting σ₀(z). Coords `target_sigma0`, `reference_depth`
(mean-profile depth of each isopycnal), `time`. Written by notebook section 7.
