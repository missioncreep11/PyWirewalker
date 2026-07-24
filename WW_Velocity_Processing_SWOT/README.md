# WW_Velocity_Processing_SWOT — MATLAB usage guide

MATLAB toolbox for processing **Nortek Signature1000 ADCP** velocity and
turbulence data from a Wirewalker, developed by
[`modscripps/wirewalker`](https://github.com/modscripps/wirewalker/tree/master/WW_Velocity_Processing_SWOT)
(Zheng / Lucas / Le Boyer / Northcott / Griffin).

This guide walks through running the MATLAB scripts **as written** — i.e. from
raw `.ad2cp` files, converted to `.mat`, to a gridded velocity (and optional
turbulence) product.

---

## 1. Prerequisites

- **MATLAB** R2017a or newer (uses `smoothdata`, `fillmissing`).
- **Signal Processing Toolbox** (`butter`, `filtfilt`).
- **Statistics and Machine Learning Toolbox** (`nanmean`, `nanstd`, etc.).
- **cbrewer** colormap toolbox — *only* needed for the final quicklook plot
  (`plot_result_adcp.m`). If you don't have it, skip that step or swap in a
  built-in colormap.
- **Nortek Signature Deployment** software (Windows; "MIDAS" / "Signature
  Deployment") to convert the raw binary `.ad2cp` into MATLAB `.mat`.

---

## 2. Convert `.ad2cp` → `.mat`

The scripts do **not** read `.ad2cp` directly. First convert with Nortek's
desktop software:

1. Open **Signature Deployment** and load your deployment / `.ad2cp` file(s).
2. Export to MATLAB. This produces one or more `.mat` files, each containing two
   structs:
   - `Data` — all the burst fields (`Data.Burst_Time`, `Data.Burst_Pressure`,
     `Data.Burst_VelBeam1..4` or `Data.Burst_Velocity_Beam`,
     `Data.Burst_CorBeam1..4`, `Data.Burst_AmpBeam1..4`,
     `Data.Burst_Pitch/Roll/Heading`, `Data.Burst_Accelerometer*`, and — if HR
     mode is enabled — `Data.IBurstHR_*`).
   - `Config` — deployment settings (blanking distance, cell size, sampling
     rate, number of cells).
3. Put all the exported `.mat` files for one deployment into a single folder.
   This folder is your **`aqdpath`**.

> A long deployment is usually exported as many sequential `.mat` files. The
> toolbox sorts and stitches them automatically (steps below).

### Find your acquisition settings

Open one `.mat` in MATLAB and inspect `Config` (and `Data`) to read off the
values you'll need in step 4:

- **blanking distance** (`blockdis`, m)
- **cell size** (`cellsize`, m)
- **sampling rate** (`saprate`, Hz)
- number of cells (used automatically via `Data.Burst_NCells` /
  `Burst_NumberofCells`)

---

## 3. Set the toolbox path

In MATLAB, add this folder to the path (or `cd` into it):

```matlab
addpath('/path/to/PyWirewalker/WW_Velocity_Processing_SWOT')
```

`SetupPath.m` will create the output directory tree for you (next step).

---

## 4. Edit the driver: `process_WW_ADCP_main.m`

Everything is driven by `process_WW_ADCP_main.m`. Open it and edit the header
block for your deployment:

```matlab
MainPath   = '/path/to/proc/';             % root for ALL outputs
Wirewalker = 'WW';                         % vehicle id  (becomes a subfolder)
Deployment = 'D1';                         % deployment  (becomes a subfolder)

WWmeta.aqdpath     = '/path/to/raw/';      % folder of .mat from step 2
WWmeta.root_script = '/path/to/WW_Velocity_Processing_SWOT';
WWmeta.name_aqd    = ['SN23_' Wirewalker '_' Deployment];  % output base name (freeform)
```

`SetupPath(WWmeta,MainPath,Wirewalker,Deployment)` then creates, under
`MainPath`:

| Folder | Contents |
|---|---|
| `Combined/<WW>/<Dep>/` | merged raw groups |
| `Profile/<WW>/<Dep>/` | per-cast structs (up/down) |
| `ReOrdered/<WW>/<Dep>/` | profiles after cross-file stitching |
| `Grid/<WW>/<Dep>/` | **final gridded products** |
| `Fig/<WW>/<Dep>/` | cast-detection QC figures |

### Set the processing `variables`

Fill in the `variables` struct from your `Config` (see step 2) and choices:

```matlab
variables.NUM_combining_files = 1;    % # raw .mat files merged per group (e.g. 20)
variables.blockdis = 0.1;             % blanking distance (m)  -- from Config
variables.cellsize = 0.5;             % cell size (m)          -- from Config
variables.saprate  = 8;               % sampling rate (Hz)     -- from Config
variables.boxsize  = 1;               % vertical bin for averaging (typically ≤ 1 m)
variables.z_max    = 510;             % max profile depth (m); set slightly > actual max depth
variables.k        = 0;               % 1 = also process/save downcasts
variables.thhold   = 2;               % min samples to count as a profile
variables.direction = 'up';           % 'up' or 'down' facing ADCP
variables.sail_corr = 1;              % 1 = correct WW horizontal motion
variables.z_unit   = [0,0,1];         % unit vector of Nortek z-axis vs WW z-axis
                                      % (e.g. tilted head: [-1/sqrt(2)*sind(22.5), -1/sqrt(2)*sind(22.5), cosd(22.5)])

variables.HRturb = 0;                 % 1 = also compute HR-mode turbulence (ε)
if variables.HRturb == 1
    variables.HRbeams    = [1,2,3,4]; % beam(s) running HR mode
    variables.HRblockdis = 0.1;       % HR blanking distance (m)
    variables.HRcellsize = 0.06;      % HR cell size (m)
    variables.HRboxsize  = 50*variables.HRcellsize;   % ε depth resolution
end
```

> `NUM_combining_files`: bigger = fewer, larger intermediate files. 20 is a
> typical default; if combined files get too large to save, lower it.
>
> The driver also filters out macOS AppleDouble sidecar files (names starting
> with `._`) from the raw `.mat` list — harmless on other platforms.

---

## 5. Run the pipeline

Run `process_WW_ADCP_main.m` **section by section** (use MATLAB's cell mode,
Ctrl/Cmd+Enter) so you can watch each stage:

1. **Set path / variables** — runs `SetupPath`, lists the raw `.mat` files into
   `WWmeta.dd0`.
2. **Sort files** (`sort_file`) — orders raw `.mat` by burst time. *Slow* (loads
   every file once).
3. **Merge + split into profiles** — loops over groups of `NUM_combining_files`:
   - `merge_signature` concatenates the group into `Combined/`.
   - `create_profiles` → `get_aqd_2G` low-passes the pressure to separate
     up/down casts (and splits on >30 s time gaps from duty cycling), saving
     per-cast structs to `Profile/`. **A QC figure is written to `Fig/`** for
     every group — open a few and confirm the red (up) / blue (down) separation
     looks right.
4. **Combine cut-off profiles** (`combine_cutoff`) — copies `Profile/` →
   `ReOrdered/` and stitches casts that span two groups.
5. **Velocity analysis** (`WWvel_upward`) — the core: correlation + sidelobe
   masking, IMU motion correction, beam→ENU rotation, optional sail correction,
   box-averaging onto a `0 : -boxsize : -z_max` grid. Results are packed into an
   `ADCP` struct and saved to `Grid/<name>_<n>.mat`. Long deployments are
   processed in chunks of `splitfiles=25` groups (`_1`, `_2`, …).
6. **Turbulence** (`WWturb_upward`, only if `HRturb==1`) — dissipation ε from
   structure-function and Kolmogorov spectral fits; saved to
   `Grid/<name>_<n>_HR_Turbulence.mat`.
7. **Quicklook** (`plot_result_adcp`) — pcolor of E/N velocity, shear, and
   backscatter (needs cbrewer). **Commented out in the driver by default** —
   uncomment the last line of `process_WW_ADCP_main.m` to run it.

---

## 6. Outputs

All final products are written to `Grid/<WW>/<Dep>/`. The `Combined/`,
`Profile/`, and `ReOrdered/` `.mat` files are bulky intermediates (often >1 GB
total) and can be deleted once the Grid products look good; `Fig/` holds the
cast-detection QC PNGs.

**What you get depends on `HRturb`:**

| `variables.HRturb` | Files produced in `Grid/` |
|---|---|
| `0` (velocity only) | `<name>_<n>.mat` (`ADCP` struct) |
| `1` (velocity + turbulence) | `<name>_<n>.mat` (`ADCP`) **and** `<name>_<n>_HR_Turbulence.mat` (`turb`) |

`<n>` is the chunk index (`splitfiles` groups per chunk). A short deployment is a
single chunk, `_1`.

---

### 6a. Velocity product — `<name>_<n>.mat` → struct `ADCP`

**Always produced** (HR mode on or off). Most fields are a 2-D grid of
**depth × cast**, where the depth axis is `0 : -boxsize : -z_max` (so `size 1` =
`z_max/boxsize + 1` rows) and each column is one upcast (only upcasts unless
`k==1`).

| Field | Dims | Units | Meaning |
|---|---|---|---|
| `time` | depth × cast | MATLAB datenum | box-averaged sample time per bin |
| `dz` | depth × cast | m | depth grid (positive down); identical in every column |
| `velE`, `velN`, `velU` | depth × cast | m/s | ENU velocity (East, North, Up), motion-corrected |
| `velE_corr`, `velN_corr` | depth × cast | m/s | horizontal velocity with "sail" (wire drift) correction — **only if `sail_corr==1`** |
| `shearE`, `shearN` | depth × cast | s⁻¹ | E/N velocity shear (∂u/∂z, ∂v/∂z) |
| `velE_var`, `velN_var`, `velU_var` | depth × cast | m/s | within-bin std of each velocity component |
| `amp` | depth × cast × 4 | dB (rel.) | range-normalized backscatter amplitude, per beam |
| `amp_var` | depth × cast × 4 | dB | within-bin std of amplitude, per beam |
| `N` | depth × cast | count | number of valid samples averaged into each bin |
| `surf_vel` | 2 × cast | m/s | surface-echo velocity [E; N] per cast |
| `Nav` | struct | — | per-bin binned navigation (see below) |
| `Notes` | char | — | free-text notes string (empty by default) |

`ADCP.Nav` is a struct; each field is a **depth × cast** grid of the bin-averaged
sensor value: `Burst_WaterTemperature`, `Burst_Temperature`, `Burst_Heading`,
`Burst_Pitch`, `Burst_Roll`.

> Bins with fewer than the minimum valid samples are left `NaN` (velocity
> requires >10 good samples per bin).

#### Sail (wire-drift) correction — `velE_corr` / `velN_corr`

*Contributed by **Caeli Griffin**.* Enabled with `variables.sail_corr = 1`,
implemented in `WWvel_upward.m`.

As the Wirewalker climbs the wire it does not rise straight up: the package
leans ("sails") and translates horizontally, so the ADCP frame itself is moving
through the water. That platform motion aliases directly into the measured
horizontal water velocity. The sail correction estimates and removes it:

1. **Rise rate** — differentiate pressure to get the vertical speed
   `dp/dt`, then smooth it (Gaussian, ~24 samples) to suppress wave jitter.
2. **Package tilt** — rotate the instrument's body z-axis (`variables.z_unit`)
   into ENU using the measured pitch/roll/heading (`XYZ2ENU`), giving the
   horizontal (`b_h`) and vertical (`b_z`) components of the tilt unit vector.
3. **Horizontal drift speed** — geometry gives the horizontal translation rate
   from the climb: `v_h = (dp/dt) · b_h / b_z`, decomposed into East/North
   components `v_x`, `v_y` along the lean direction.
4. **Apply** — add the drift back to the ENU velocity:
   `velE_corr = velE + v_x`, `velN_corr = velN + v_y`.

The uncorrected `velE`/`velN` are retained alongside, so you can compare. The
correction grows with rise rate and tilt, so it matters most on fast, leaning
casts. `variables.z_unit` must describe how the Nortek z-axis is mounted
relative to the Wirewalker z-axis (e.g. a 22.5°-tilted head:
`[-1/sqrt(2)*sind(22.5), -1/sqrt(2)*sind(22.5), cosd(22.5)]`).

---

### 6b. Turbulence product — `<name>_<n>_HR_Turbulence.mat` → struct `turb`

**Only produced when `HRturb==1`** (requires `IBurstHR_*` HR-mode fields in the
export). Dissipation ε is estimated two independent ways — a **wavenumber
spectral** Kolmogorov fit and a **second-order structure function** fit — for
each HR beam in `variables.HRbeams`.

Coordinate vectors:

| Field | Dims | Meaning |
|---|---|---|
| `z` | depth × 1 | depth bin centers (m) for the **spectral** estimates |
| `z_struct` | depth × 1 | depth grid (m) for the **structure-function** estimates |
| `time` | 1 × prof | profile time (MATLAB datenum), de-duplicated |
| `k` | nk × 1 | along-beam wavenumber (cpm) for `spec` |
| `r` | nr × 1 | separation distance (m) for `struct_fun` |
| `beam_number` | 1×1×nbeam | which HR beam each slice corresponds to |

Spectral-method estimates (dims **depth × prof × beam**):

| Field | Units | Meaning |
|---|---|---|
| `ep` | W kg⁻¹ | **dissipation rate ε** (spectral; `(A/0.53)^{3/2}`) |
| `N` | — | spectral noise floor |
| `A` | — | Kolmogorov (k⁻⁵ᐟ³) amplitude coefficient |
| `SNR` | — | signal-to-noise (`A/N`) |
| `corr` | counts | mean beam correlation in the bin (QC) |
| `slope`, `N_slope` | — | fitted log-log spectral slope and intercept (QC) |
| `spec` | (depth × prof × nk × beam) | the averaged wavenumber spectra |
| `spec_num` | same as `spec` | # of spectra averaged per estimate |

Structure-function estimates (dims **depth × prof × beam**, on the separate
decimated depth grid `turb.z_struct`):

| Field | Units | Meaning |
|---|---|---|
| `ep_struct` | W kg⁻¹ | **dissipation rate ε** (structure-function; `(A/2)^{3/2}`) |
| `A_struct` | — | structure-function amplitude (∝ ε²ᐟ³) |
| `N_struct` | — | structure-function noise/offset term |
| `struct_fun` | (depth × prof × nr × beam) | the second-order structure functions D(r) |

> Two ε estimates are provided so you can cross-check; they should broadly agree.
> Plot either on a **log x-axis** (see `plot_uv_tke_profile.m`). The spectral
> `ep` uses `turb.z`; `ep_struct` uses `turb.z_struct`.

---

Merge per-chunk `Grid` files with `Combine_Grid_Files.m` (edit the folder and
file-range at the top) to get one continuous `ADCP` record.

---

## 7. Tips & gotchas

- **Run cell-by-cell the first time.** The cast-detection figures in `Fig/` are
  your main sanity check; if up/down are mislabeled, revisit `thhold` and the
  pressure record.
- **`direction`** must match the instrument orientation (`'up'` vs `'down'`); it
  controls how per-bin depth is computed.
- **HR turbulence** requires `IBurstHR_*` fields in the export — only present if
  HR mode was enabled on the listed `HRbeams`. If absent, set `HRturb = 0`.
- **Beam geometry is hard-coded** for the Signature1000 (beam angle 25°,
  φ = 65°, azimuths [0, −90, 180, 90]°) in `Beam2ENU`, `WWvel_upward`,
  `WWcorr_beam`, `WWturb_upward`. Change there if you use a different head.
- **Correlation threshold** for masking bad samples is 50 (hard-coded in
  `WWvel_upward`/`WWturb_upward`).
- **SWOT duty-cycle emulation**: to test 1/3 duty cycling on a continuous
  record, swap `merge_signature` for `merge_signature_SWOT_emulator` in the
  driver loop.
- **Fixed/moored ADCP** (not on a Wirewalker): use `ProcessFixedADCP.m` instead
  of this driver.
- See the bundled manuals for more detail: `Manual for WW_ADCP.docx`,
  `VelocityDataReadMe.docx`, `NortekTurbulenceDataReadMe.docx`.

---

## 8. Authorship & credits

- **Core velocity processing** (file sorting, merging, cast splitting,
  beam→ENU, IMU motion correction, box-averaging) — **Bofu Zheng**, **Arnaud
  Le Boyer**, and **Drew Lucas**. Method described in **Zheng et al. (2022)**
  [[1]](#references).
- **Turbulence processing** (`WWturb_upward.m`; HR-mode dissipation ε via
  structure-function and Kolmogorov spectral fits) — **Devon Northcott**.
  Method described in **Northcott et al. (2026)** [[2]](#references).
- **Sail (wire-drift) correction** (`sail_corr`; horizontal-motion removal in
  `WWvel_upward.m`, see §6a) — **Caeli Griffin**.

Please credit the appropriate authors and cite the references below when using or
adapting these components.

---

## References

1. Zheng, B., A. J. Lucas, R. Pinkel, and A. Le Boyer, 2022: Fine-Scale Velocity
   Measurement on the Wirewalker Wave-Powered Profiler. *Journal of Atmospheric
   and Oceanic Technology*, **39**(2), 133–147,
   <https://doi.org/10.1175/JTECH-D-21-0048.1>.

2. Northcott, D., A. Le Boyer, J. MacKinnon, M. H. Alford, and A. J. Lucas, 2026:
   Observations of the Internal Wave to Turbulence Cascade. *Journal of Physical
   Oceanography*, **56**(4), 839–853,
   <https://doi.org/10.1175/JPO-D-25-0114.1>.
