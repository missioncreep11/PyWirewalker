# Python port plan — WW_Velocity_Processing_SWOT

This directory holds the **original MATLAB** Wirewalker ADCP velocity/turbulence
toolbox, pulled verbatim from
[`modscripps/wirewalker/WW_Velocity_Processing_SWOT`](https://github.com/modscripps/wirewalker/tree/master/WW_Velocity_Processing_SWOT)
(Zheng / Lucas / Le Boyer / Northcott). It is kept as the reference
implementation for a future Python port that follows the conventions already
established by [`process_wirewalker_rbr.py`](../process_wirewalker_rbr.py) (config-driven,
numpy/xarray/gsw/scipy, level-based NetCDF archive).

The MATLAB `.m` files and `.docx` manuals are reference only — nothing here is
wired into the Python pipeline yet.

## What the pipeline does

Processes **Nortek Signature1000** ADCP data collected on an upward- (or
downward-) looking Wirewalker. The MATLAB reference reads `.mat` files exported
by Nortek's *Signature Deployment* software (which converts the raw `.ad2cp`).
**The Python port skips that step and reads the raw `.ad2cp` directly** (see
"Ingestion" below). Output is a depth × time grid of ENU velocity, shear,
backscatter, and (HR mode) turbulent dissipation.

## Ingestion — read `.ad2cp` directly with DOLfYN

The port does **not** depend on the Nortek Signature Deployment MATLAB export.
Instead it ingests the raw binary `.ad2cp` files directly using
[**DOLfYN**](https://dolfyn.readthedocs.io/) (`pip install dolfyn`), the
industry-standard open-source library for reading and processing binary Nortek
and TRDI ADCP/ADV files. DOLfYN natively parses `.ad2cp`, handles the burst /
IBurstHR record layouts, and returns an **`xarray.Dataset`** — so it slots
straight into the PyWirewalker stack (numpy/xarray/gsw/scipy) and the rest of the
pipeline operates on xarray from the first stage.

What this changes vs. the MATLAB reference:

- **Replaces `merge_signature` / `sort_file`** as the loader: `dolfyn.read()`
  (or `read_example` / globbed multi-file read) gives a time-ordered Dataset, so
  the manual `.mat` sort/concat-by-burst-time logic is unnecessary. Keep the
  duty-cycle emulator idea (`merge_signature_SWOT_emulator`) as a post-load mask.
- **Field-name fallbacks largely go away**: DOLfYN normalizes variable names
  (`vel`, `amp`, `corr`, `pressure`, `heading`/`pitch`/`roll`, `accel`, beam
  ranges, time) regardless of firmware, so the `Burst_VelBeam1` vs
  `Burst_Velocity_Beam` and `Burst_Time` vs `Burst_MatlabTimeStamp` branching
  collapses into the DOLfYN schema.
- **Time** comes back as `datetime64` — no MATLAB-datenum conversion needed.
- DOLfYN also provides vetted **beam↔inst↔earth (ENU) rotations**
  (`dolfyn.rotate2`) and a config-aware orientation matrix. Decide per function
  whether to use DOLfYN's rotations or port the toolbox's own
  (`Beam2ENU`/`XYZ2ENU`/`GetUnitVectors`) — the toolbox versions bake in the
  Wirewalker-specific surface-echo / sidelobe / motion-correction handling, so we
  likely keep those but cross-check against DOLfYN on a sample.
- Raw `.ad2cp` retains **all** records (every beam, IBurstHR, full attitude/IMU),
  so nothing is lost relative to the `.mat` export.

Stages, in `process_WW_ADCP_main.m` order:

1. **Ingest** `.ad2cp` with DOLfYN → one time-ordered `xarray.Dataset`
   (replaces `sort_file` + `merge_signature`).
2. *(optional)* **Duty-cycle mask** to emulate SWOT 1/3 on-time
   (`merge_signature_SWOT_emulator`), now applied as a post-load time mask.
3. **Split into profiles** (`create_profiles` → `get_aqd_2G`): low-pass the
   pressure, find turning points to separate up/down casts, split on large time
   gaps (duty cycling), save each cast as a struct.
4. **Stitch boundary casts** that straddle two file-groups (`combine_cutoff`).
5. **Velocity** (`WWvel_upward`): the heavy lifting — beam geometry, per-ping
   depth coordinates, correlation + sidelobe masking, amplitude normalization,
   surface-echo velocity, beamwise shear, IMU motion correction
   (`WWcorr_beam`), beam→ENU rotation (`Beam2ENU`), optional "sail"
   (horizontal-drift) correction, then box-average onto a uniform depth grid.
6. **Turbulence** (`WWturb_upward`, HR mode): dissipation ε from the
   second-order structure function and from Kolmogorov spectral fits
   (`FitKolmogorov`).
7. **Quicklook** plots (`plot_result_adcp`) and grid concatenation
   (`Combine_Grid_Files`).

## File-by-file port map

| MATLAB file | Role | Python target | Notes |
|---|---|---|---|
| `process_WW_ADCP_main.m` | driver script + `variables` struct | `process_ww_adcp.py` CLI + `config.json` keys | mirror `process_wirewalker_rbr.py` arg/level style; move all `variables.*` into config |
| `SetupPath.m` | make Combined/Profile/ReOrdered/Grid/Fig dirs | `pathlib`-based helper | derive from `output_dir`, levels as subdirs |
| `sort_file.m` | sort raw `.mat` by time | **obsolete** — DOLfYN | `dolfyn.read()` returns a time-ordered Dataset |
| `merge_signature.m` | concat N raw files | **obsolete** — DOLfYN | read `.ad2cp` directly; no `.mat` merge step |
| `merge_signature_SWOT_emulator.m` | duty-cycle emulator | optional flag | 30-min on/off, keep 1 of 3 windows; post-load time mask |
| `get_aqd_2G.m` | up/down cast detection | `split_casts()` | Butterworth low-pass (`scipy.signal.butter`/`filtfilt`), turning points, 30 s time-gap split |
| `create_profiles.m` | save per-cast structs | profile builder | replace cell-array-of-structs with list of `xr.Dataset` |
| `combine_cutoff.m` | stitch cross-file casts | profile builder | pressure jump <10 dbar & time gap <20 s |
| `Beam2XYZ.m` / `Beam2ENU.m` / `XYZ2ENU.m` | coordinate transforms | `transforms.py` | vectorize the per-ping loop; θ=25°, heading−90° for ENU |
| `GetUnitVectors.m` | beam unit vectors vs attitude | `transforms.py` | Klymak rotation; φ=65°, azi=[0,−90,180,90]° |
| `WWcorr.m` / `WWcorr_beam.m` | IMU motion correction | `motion_correction.py` | de-tilt accel, bandpass 0.1–1.2 Hz, integrate to velocity, combine with dp/dt; `_beam` returns correction in beam coords |
| `WWvel_upward.m` | main velocity processing | `velocity.py` | the bulk of the work — see below |
| `WWvel_downward_2.m` | downward-looking variant | `velocity.py` (mode) | fold in as `direction='down'` path |
| `WWturb_upward.m` | HR-mode dissipation | `turbulence.py` | structure-function + spectral ε; needs phase-unwrap / ambiguity handling |
| `FitKolmogorov.m` | k^(−5/3) least-squares fit | `turbulence.py` | small linear-algebra helper |
| `plot_result_adcp.m` | quicklook pcolor | `plot_adcp.py` | matplotlib `pcolormesh`; drop cbrewer dep |
| `Combine_Grid_Files.m` | concat grid outputs | merge util | xarray `concat` over time |
| `ProcessFixedADCP.m` (+`.asv`) | fixed/moored ADCP variant | `process_fixed_adcp.py` | separate driver, lower priority |
| `Plot_vel_downward_info.m` | downward quicklook | `plot_adcp.py` | low priority |
| `*.docx` | manuals / readmes | keep as docs | `Manual for WW_ADCP.docx`, velocity & turbulence ReadMes |

## `WWvel_upward` — the core, broken down

Per cast, per ping group:
- Gather beam velocity / correlation / amplitude (handle both
  `Burst_VelBeam{1..4}` and `Burst_Velocity_Beam` 3-D layouts).
- Mask samples with any beam correlation < 50.
- Beam unit vectors from attitude (`GetUnitVectors`); per-bin depth
  `z = -pressure + range·bZ`, range corrected by `/cos(25°)`.
- Surface echo: bin nearest z=0 → surface beam velocity.
- Amplitude normalization: add transmission loss `10·log10((2r)²)+2·0.37·r`.
- Sidelobe mask near the surface bounce.
- Beamwise shear via centered differences along range.
- Interpolate beam velocity & shear onto a nominal even-depth grid (`interp1`).
- Motion-correct in beam coordinates (`WWcorr_beam`), then rotate to ENU
  (`Beam2ENU`).
- Optional **sail correction**: horizontal velocity from rise rate `dp/dt` and
  the platform z-axis tilt (`z_unit`).
- Box-average E/N/U velocity, shear, amplitude (+variances), nav vars, time onto
  `z = 0:-boxsize:-z_max`; require >10 valid samples per bin.

Output (currently a 17-element cell array `out{}`) should become a single
`xr.Dataset` (dims `depth`, `cast`) with variables `velE, velN, velU, shearE,
shearN, amp, amp_var, *_var, surf_vel, Nav.*, N` and, when `sail_corr`,
`velE_corr/velN_corr`.

## Config keys to add (mirroring `variables` struct)

```
adcp:
  ad2cp_glob (path/glob to raw .ad2cp files),  # ingested with DOLfYN
  num_combining_files, blockdis_m, cellsize_m, sample_rate_hz,
  boxsize_m, z_max_m, process_downcast, profile_threshold,
  direction (up|down), sail_corr, z_unit [3],
  corr_min (50), beam_angle_deg (25), beam_phi_deg (65),
  beam_azi_deg [0,-90,180,90],
  motion_bandpass_hz [0.1, 1.2],
  hr_turb: { enabled, beams, blockdis_m, cellsize_m, boxsize_m }
```

## Suggested port order

0. **Ingestion spike**: `dolfyn.read()` on a sample `.ad2cp`; map DOLfYN
   variable/coord names → what the rest of the pipeline needs; confirm IBurstHR
   (turbulence) records load.
1. `transforms.py` (`GetUnitVectors`, `Beam2XYZ/ENU`, `XYZ2ENU`) + unit tests
   against MATLAB outputs on a saved sample; cross-check vs `dolfyn.rotate2`.
2. Loaders on the DOLfYN Dataset: cast split / stitch (sort + merge are now DOLfYN).
3. `motion_correction.py` (`WWcorr_beam`).
4. `velocity.py` (`WWvel_upward`) — validate a gridded profile against MATLAB.
5. `plot_adcp.py` quicklook.
6. `turbulence.py` (`FitKolmogorov`, `WWturb_upward`).
7. Fixed-ADCP and downward variants last.

## Gotchas

- The MATLAB reference uses **MATLAB datenums** (days; `·86400` → seconds) and
  per-firmware field-name fallbacks. With DOLfYN ingestion both largely go away
  (DOLfYN gives `datetime64` and normalized variable names) — but watch any time
  arithmetic ported verbatim from `.m` (e.g. `·86400`, `dp/dt`) and convert to
  seconds explicitly.
- `.ad2cp` are large binary files; DOLfYN can read incrementally — consider
  chunking / a `.nc` cache of the raw ingest so re-runs don't reparse.
- Cross-check DOLfYN's coordinate-system convention (beam/inst/earth,
  up vs down config, declination) against the toolbox's assumptions before
  trusting either rotation path.
- MATLAB is 1-indexed and column-major; watch `squeeze`/`permute` when moving to
  numpy.
- `filtfilt` edge handling here uses manual flip-pad-flip thirds — replicate or
  use `scipy.signal.filtfilt`'s padding and compare.
- **Two depth grids in the turbulence output**: spectral ε (`ep`) is on
  `turb.z` (bin centers); structure-function ε (`ep_struct`) is on a *different*
  decimated grid. The stock upstream `WWturb_upward.m` never saved that second
  axis — our local copy adds `turb.z_struct` (see below). The Python port must
  carry **both** depth coordinates explicitly (e.g. two coords or two Datasets),
  not assume `ep` and `ep_struct` share a depth axis.

## Local modifications to the MATLAB reference

These diverge from the verbatim `modscripps/wirewalker` copy and should be
carried into the Python port:

- **`WWturb_upward.m`**: added `turb.z_struct = transpose(z_grid_dec);` so the
  structure-function estimates (`ep_struct`, `A_struct`, `N_struct`,
  `struct_fun`) have an explicit depth axis. Upstream only saved `turb.z` (the
  spectral grid), leaving `ep_struct` effectively un-georeferenced in depth.
  → Port: emit a `z_struct` coordinate alongside `z`.
- **Turbulence output filename**: dropped the upstream `_Test2` suffix; the file
  is now `<name>_<n>_HR_Turbulence.mat`. → Port: use a clean, stable name.
- **`HRbeams`**: this instrument runs HR mode on **beam 5 only**
  (`IBurstHR_*Beam5`); the upstream/Caeli drivers variously default to `[5]` or
  `[1,2,3,4]`. → Port: detect available `IBurstHR_*Beam<n>` fields rather than
  hardcoding the beam list.
