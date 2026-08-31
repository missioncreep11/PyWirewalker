"""Assemble per-cast gridded velocity profiles into an L2 (depth x cast) product.

`build_l2` grids an in-memory dolfyn Dataset (beam coords); `build_l2_streaming`
does the same over a raw `.ad2cp` too large to hold in memory, by reading the
file in chunks and carrying the boundary cast across chunks. Both stack per-cast
profiles (`velocity.process_cast`) into an `xarray.Dataset` with dims
(depth, cast), mirroring the CTD L2 conventions.
"""
from __future__ import annotations

import os
import time as _time

import numpy as np
import xarray as xr

from .attitude import DEFAULT_CUTOFF_HZ, reconstruct
from .casts import detect_casts
from .geometry import look_direction
from .platform import AHRS_BAD_DEG, ahrs_error
from .velocity import NOTCH_MAX_DEPTH_M, process_cast, output_grid

ATTITUDE_MODES = ("ahrs", "reconstructed", "auto")
MOTION_VERSIONS = ("v1", "v2", "v3")
BIN_AVERAGE_MODES = ("boxcar", "notch")
_TILT_SOURCE_FLAG = {"ahrs": 0, "lp_accel": 1, "ahrs_fallback": 2}
_HEADING_SOURCE_FLAG = {"ahrs": 0, "mag": 1}


def _cast_attitude(dsc, mode):
    """Choose the attitude source for one cast.

    Returns ``(dsc, source_flag, ahrs_error_deg)`` where the flag is
    0 = AHRS pitch/roll used as-is, 1 = replaced by the accel reconstruction,
    2 = replacement wanted but the reconstruction was not usable, AHRS kept.

    ``auto`` substitutes only on detector-flagged casts (`ahrs_error` >=
    `platform.AHRS_BAD_DEG`); ``reconstructed`` substitutes everywhere it is valid.
    The error is computed in every mode so the product always carries the QC.
    """
    err = np.nan
    if "accel" in dsc and "orientmat" in dsc:
        err = ahrs_error(dsc["accel"].values, dsc["orientmat"].values)
    if mode == "ahrs" or "accel" not in dsc:
        return dsc, 0, err
    if mode == "auto" and not (np.isfinite(err) and err >= AHRS_BAD_DEG):
        return dsc, 0, err
    rec = reconstruct(dsc)
    if not rec.usable:
        return dsc, 2, err
    return dsc.assign(pitch=("time", rec.pitch_deg),
                      roll=("time", rec.roll_deg)), 1, err


def _motion_desc(motion, attitude, sail):
    """Human-readable motion_correction provenance. v1/v2 are frozen legacy models;
    v3 records the user-selected attitude source and sail state."""
    if motion == "v1":
        return "v1 (legacy): WWcorr_beam (bandpass-integrated IMU + dp/dt), AHRS attitude"
    if motion == "v2":
        return ("v2 (legacy): buoyant-ascent; LP-accel tilt + tilt-compensated mag heading; "
                "sail ON; depth-gain-weighted bandpass IMU horizontal; dp/dt vertical; spike-interp")
    att = "AHRS" if attitude == "ahrs" else "LP-accel tilt + tilt-compensated mag heading"
    return (f"v3: buoyant-ascent; {att}; sail {'ON' if sail else 'OFF'}; "
            "depth-gain-weighted bandpass IMU horizontal; dp/dt vertical; spike-interp")


def _assemble(results, *, boxsize, z_max, look, cast_kind, corr_min, min_bin_samples,
              motion_correct, mooring, source, attitude="ahrs", motion="v1",
              bin_average="boxcar", sail=True):
    """Stack a list of (cast_result, Cast, att_src, ahrs_err) into an L2 xarray.Dataset."""
    zc = output_grid(boxsize, z_max)
    nz, ncast = zc.size, len(results)
    G = {k: np.full((nz, ncast), np.nan, np.float32)
         for k in ("velE", "velN", "velU", "amp", "shearE", "shearN",
                   "velE_sem", "velN_sem", "velU_sem", "shearE_sem", "shearN_sem")}
    nobs = np.zeros((nz, ncast), np.int32)
    # per-(depth, cast) sample time as float ms since epoch (NaN where empty). Stored
    # with a CF units attr rather than as datetime64 so xarray writes it as a plain
    # float column: a 2-D datetime64 coord with NaT triggers a fill-value that overflows
    # on decode. (Mirrors the CTD L2, which stores its 2-D time the same way.)
    gtime = np.full((nz, ncast), np.nan)
    ctime = np.empty(ncast, "datetime64[ns]")         # per-cast mid-time (for ordering)
    cpmax = np.zeros(ncast)
    cpmin = np.zeros(ncast)
    cdir = np.zeros(ncast, np.int8)
    ccomplete = np.zeros(ncast, np.int8)
    csrc = np.zeros(ncast, np.int8)
    chead = np.zeros(ncast, np.int8)
    cerr = np.full(ncast, np.nan)
    for j, (g, cast, src, err) in enumerate(results):
        for k in G:
            G[k][:, j] = g[k]
        nobs[:, j] = g["n_obs"]
        _tt = np.asarray(g["time"]).astype("datetime64[ms]")     # (nz,) per-bin time
        _tc = _tt.view("int64").astype(float)
        _tc[np.isnat(_tt)] = np.nan
        gtime[:, j] = _tc
        ctime[j] = np.datetime64(g["cast_time"], "ns")
        cpmax[j] = g["pressure_max"]
        cpmin[j] = g["pressure_min"]
        cdir[j] = 1 if cast.direction == "up" else 0
        ccomplete[j] = 0 if cast.truncated else 1
        # motion v2 chooses the tilt inside process_cast; its report wins over
        # the _cast_attitude decision when present
        csrc[j] = _TILT_SOURCE_FLAG[g["tilt_source"]] if "tilt_source" in g else src
        chead[j] = _HEADING_SOURCE_FLAG.get(g.get("heading_source"), 0)
        cerr[j] = err

    order = np.argsort(ctime)                      # ensure time order
    for k in G:
        G[k] = G[k][:, order]
    gtime = gtime[:, order]
    nobs, ctime, cdir = nobs[:, order], ctime[order], cdir[order]
    cpmax, cpmin, ccomplete = cpmax[order], cpmin[order], ccomplete[order]
    csrc, chead, cerr = csrc[order], chead[order], cerr[order]
    n_trunc = int((ccomplete == 0).sum())

    return xr.Dataset(
        {"velE": (("depth", "cast"), G["velE"], {"units": "m s-1", "long_name": "eastward velocity"}),
         "velN": (("depth", "cast"), G["velN"], {"units": "m s-1", "long_name": "northward velocity"}),
         "velU": (("depth", "cast"), G["velU"], {"units": "m s-1", "long_name": "upward velocity"}),
         "amp": (("depth", "cast"), G["amp"], {"units": "dB", "long_name": "beam-mean backscatter amplitude"}),
         **{f"shear{c}": (("depth", "cast"), G[f"shear{c}"],
                          {"units": "s-1",
                           "long_name": f"vertical shear of vel{c} (beam-differenced)",
                           "comment": "centred cell differences along each beam "
                                      "before rotation (WWvel_upward beamshear), "
                                      "so per-ping common-mode errors - platform "
                                      "motion, attitude leakage, the sail term - "
                                      "cancel exactly; no motion correction is "
                                      "applied or needed. Gaussian-smoothed ~1 s "
                                      "along pings; two cells nearest the "
                                      "transducer masked (stagnation)"})
            for c in ("E", "N")},
         **{f"shear{c}_sem": (("depth", "cast"), G[f"shear{c}_sem"],
                              {"units": "s-1",
                               "long_name": f"standard error of shear{c} "
                                            f"(Doppler noise only)"})
            for c in ("E", "N")},
         **{f"vel{c}_sem": (("depth", "cast"), G[f"vel{c}_sem"],
                            {"units": "m s-1",
                             "long_name": f"standard error of vel{c} (Doppler noise only)",
                             "comment": "from the measured beam-correlation noise "
                                        "relation (ww_sig1000.velocity.beam_sigma); "
                                        "excludes correlated errors (attitude "
                                        "systematics, residual surface-wave "
                                        "contamination), so a floor, not a total "
                                        "error bar"})
            for c in ("E", "N", "U")},
         "n_obs": (("depth", "cast"), nobs, {"long_name": "samples averaged per bin"})},
        coords={"depth": ("depth", zc.astype(np.float32), {"units": "m", "positive": "down"}),
                "cast": ("cast", np.arange(ncast, dtype=np.int32)),
                "time": (("depth", "cast"), gtime,
                         {"units": "milliseconds since 1970-01-01 00:00:00 UTC",
                          "standard_name": "time", "calendar": "standard",
                          "long_name": "mean sample time in each (depth, cast) bin",
                          "comment": "2-D because a buoyant upcast is slanted in time "
                                     "(deep sampled before shallow)"}),
                "pressure_max": ("cast", cpmax.astype(np.float32), {"units": "dbar"}),
                "pressure_min": ("cast", cpmin.astype(np.float32), {"units": "dbar"}),
                "cast_direction": ("cast", cdir,
                                   {"flag_values": np.array([0, 1], np.int8),
                                    "flag_meanings": "down up"}),
                "profile_complete": ("cast", ccomplete,
                                     {"long_name": "profile completeness flag",
                                      "flag_values": np.array([0, 1], np.int8),
                                      "flag_meanings": "truncated complete",
                                      "comment": "0 = clipped by a duty-cycle burst boundary "
                                                 "or the record edge, so the cast covers only "
                                                 "part of the profile; 1 = bounded by pressure "
                                                 "turning points"}),
                "attitude_source": ("cast", csrc,
                                    {"long_name": "pitch/roll source for this cast",
                                     "flag_values": np.array([0, 1, 2], np.int8),
                                     "flag_meanings": "ahrs accel_reconstructed ahrs_fallback",
                                     "comment": "0 = AHRS pitch/roll; 1 = replaced by the "
                                                "low-passed accelerometer tilt (heading "
                                                "remains AHRS); 2 = replacement requested "
                                                "but the reconstruction was not usable "
                                                "(vehicle too tilted or |a| != g), AHRS "
                                                "kept - treat this cast's velocity as "
                                                "suspect"}),
                "heading_source": ("cast", chead,
                                   {"long_name": "heading source for this cast",
                                    "flag_values": np.array([0, 1], np.int8),
                                    "flag_meanings": "ahrs magnetometer",
                                    "comment": "1 = tilt-compensated magnetometer "
                                               "compass (Stage 2); 0 = AHRS heading "
                                               "(v1, fallback, or mag field unusable). "
                                               "The AHRS has heading-only fault modes "
                                               "invisible to ahrs_error_deg"}),
                "ahrs_error_deg": ("cast", cerr.astype(np.float32),
                                   {"units": "degrees",
                                    "long_name": "median angle between AHRS up and measured gravity",
                                    "comment": f"per-cast AHRS attitude-fault statistic "
                                               f"(ww_sig1000.platform.ahrs_error); healthy "
                                               f"casts sit at 4-7 deg, >= {AHRS_BAD_DEG:g} "
                                               f"deg means the orientation solution is not "
                                               f"trustworthy"})},
        attrs={"title": f"{mooring} Wirewalker ADCP L2 (gridded velocity)",
               "mooring": mooring, "source_file": source,
               "instrument_look": look, "cast_kind": cast_kind,
               "n_casts_truncated": n_trunc,
               "attitude_mode": attitude,
               "bin_average": ("boxcar mean" if bin_average == "boxcar"
                               else f"notch above {NOTCH_MAX_DEPTH_M:g} m (ridged "
                                    f"constant + wave-band fit per bin dwell), "
                                    f"boxcar mean below"),
               "n_casts_attitude_reconstructed": int((csrc == 1).sum()),
               "n_casts_attitude_fallback": int((csrc == 2).sum()),
               "n_casts_heading_mag": int((chead == 1).sum()),
               "attitude_reconstruction": (f"accel lowpass {DEFAULT_CUTOFF_HZ:g} Hz "
                                           f"(ww_sig1000.attitude)"),
               "grid_boxsize_m": boxsize, "grid_z_max_m": z_max,
               "corr_min": corr_min, "min_bin_samples": min_bin_samples,
               "motion_version": motion,
               "motion_correction": (
                   "none" if not motion_correct
                   else _motion_desc(motion, attitude, sail)),
               "processing": "ww_sig1000 port of WW_Velocity_Processing_SWOT",
               "date_created": _time.strftime("%Y-%m-%dT%H:%M:%S")},
    )


def _select_casts(press, t_s, fs, kinds, thhold_s, min_span_dbar):
    return [c for c in detect_casts(press, t_s, thhold=int(thhold_s * fs))
            if c.direction in kinds and np.ptp(press[c.start:c.stop + 1]) >= min_span_dbar]


def build_l2(ds, *, boxsize=1.0, z_max=None, cast_kind="both", min_span_dbar=40.0,
             corr_min=50, min_bin_samples=10, thhold_s=30.0, motion_correct=True,
             attitude="ahrs", motion="v1", bin_average="boxcar", sail=True,
             mooring="", source=""):
    """Grid substantial casts of an in-memory Dataset into an L2 Dataset."""
    if attitude not in ATTITUDE_MODES:
        raise ValueError(f"attitude must be one of {ATTITUDE_MODES}, got {attitude!r}")
    if motion not in MOTION_VERSIONS:
        raise ValueError(f"motion must be one of {MOTION_VERSIONS}, got {motion!r}")
    if bin_average not in BIN_AVERAGE_MODES:
        raise ValueError(f"bin_average must be one of {BIN_AVERAGE_MODES}, "
                         f"got {bin_average!r}")
    fs = float(ds.attrs["fs"])
    t_s = ds["time"].values.astype("datetime64[ns]").astype("int64") / 1e9
    press = ds["pressure"].values
    look = look_direction(ds["pitch"].values, ds["roll"].values)
    if z_max is None:
        z_max = float(np.ceil(np.nanmax(press) / boxsize) * boxsize)
    kinds = ("up", "down") if cast_kind == "both" else (cast_kind,)

    casts = _select_casts(press, t_s, fs, kinds, thhold_s, min_span_dbar)
    if not casts:
        raise RuntimeError(f"no {cast_kind}-casts with span >= {min_span_dbar} dbar found")
    results = []
    for c in casts:
        dsc, src, err = _cast_attitude(ds.isel(time=slice(c.start, c.stop + 1)), attitude)
        g = process_cast(dsc, corr_min=corr_min, boxsize=boxsize, z_max=z_max,
                         direction=look, min_bin_samples=min_bin_samples,
                         motion_correct=motion_correct, motion=motion,
                         bin_average=bin_average, sail=sail, attitude=attitude)
        results.append((g, c, src, err))
    return _assemble(results, boxsize=boxsize, z_max=z_max, look=look, cast_kind=cast_kind,
                     corr_min=corr_min, min_bin_samples=min_bin_samples,
                     motion_correct=motion_correct, mooring=mooring, source=source,
                     attitude=attitude, motion=motion, bin_average=bin_average, sail=sail)


def _count_ensembles(fn, reader):
    """Total burst ensembles in `fn`.

    Read exactly from the dolfyn index, which dolfyn has to build anyway before it
    can read a single ensemble. Falls back to an exponential probe + bisection if
    the index cannot be read — that fallback rounds *down* to a 5000-ensemble
    tolerance, so it can drop ~10 min of data at 8 Hz; the index path does not.
    """
    try:
        from .index import count_ensembles
        return count_ensembles(fn)
    except Exception as e:                          # pragma: no cover - fallback path
        print(f"  [warn] could not read the dolfyn index ({e}); "
              f"falling back to bisection probe", flush=True)

    def ok(start):
        try:
            return reader(fn, nens=[start, start + 4]).sizes.get("time", 0) > 0
        except Exception:
            return False
    hi = 100_000
    while ok(hi):
        hi *= 2
    lo = hi // 2
    while hi - lo > 5000:
        mid = (lo + hi) // 2
        if ok(mid):
            lo = mid
        else:
            hi = mid
    return lo


def build_l2_streaming(fn, reader, *, chunk=500_000, total=None, ens_start=0, boxsize=1.0,
                       z_max=None, cast_kind="both", min_span_dbar=40.0, corr_min=50,
                       min_bin_samples=10, thhold_s=30.0, gap_s=30.0, motion_correct=True,
                       attitude="ahrs", motion="v1", bin_average="boxcar", sail=True,
                       mooring="", source="", progress=True):
    """Grid a raw `.ad2cp` too large for memory. `reader` is dolfyn.read.

    Reads the file in `chunk`-ensemble windows, detects/grids casts on a rolling
    buffer, and carries the last (possibly boundary-spanning) cast into the next
    chunk. z_max, if not given, is set from the first chunk's max pressure + range.

    `ens_start` / `total` bound the ensemble range, so the deployment/recovery
    transit — when the vehicle is on deck or being lowered and is not profiling —
    can be trimmed off rather than entering the product as spurious casts.

    `attitude` selects the pitch/roll source per cast (see `_cast_attitude`):
    "ahrs" uses the instrument's fused solution as-is, "reconstructed" replaces it
    with the low-passed accelerometer tilt wherever that is valid, and "auto"
    replaces it only on casts the AHRS-fault detector flags. The choice and the
    per-cast outcome are recorded in the product.
    """
    if attitude not in ATTITUDE_MODES:
        raise ValueError(f"attitude must be one of {ATTITUDE_MODES}, got {attitude!r}")
    if motion not in MOTION_VERSIONS:
        raise ValueError(f"motion must be one of {MOTION_VERSIONS}, got {motion!r}")
    if bin_average not in BIN_AVERAGE_MODES:
        raise ValueError(f"bin_average must be one of {BIN_AVERAGE_MODES}, "
                         f"got {bin_average!r}")
    if total is None:
        total = _count_ensembles(fn, reader)
    kinds = ("up", "down") if cast_kind == "both" else (cast_kind,)
    results, look, buf = [], None, None
    carry_trunc = None      # start-truncation of the cast carried in from the last chunk
    start, t0 = int(ens_start), _time.time()

    while start < total:
        stop = min(start + chunk, total)
        ds_chunk = reader(fn, nens=[start, stop])
        # drop HR beam-5 data (own time_b5/range_b5 axes) — not used for velocity L2,
        # and it breaks concatenation along 'time'.
        ds_chunk = ds_chunk.drop_dims([d for d in ("time_b5", "range_b5")
                                       if d in ds_chunk.dims])
        buf = ds_chunk if buf is None else xr.concat([buf, ds_chunk], dim="time",
                                                     data_vars="minimal", coords="minimal")

        fs = float(buf.attrs["fs"])
        press = buf["pressure"].values
        t_s = buf["time"].values.astype("datetime64[ns]").astype("int64") / 1e9
        if look is None:
            look = look_direction(buf["pitch"].values, buf["roll"].values)
        if z_max is None:
            z_max = float(np.ceil((np.nanmax(press) + buf.sizes["range"] * buf.attrs["cell_size"])
                                  / boxsize) * boxsize)

        casts = detect_casts(press, t_s, thhold=int(thhold_s * fs), gap_s=gap_s,
                             first_is_continuation=carry_trunc is not None)
        if carry_trunc and casts:
            casts[0].truncated = True      # it was already clipped in the previous chunk
        last = start + chunk < total       # more data coming -> hold back the final cast
        to_do = casts[:-1] if (last and casts) else casts
        for c in to_do:
            if c.direction in kinds and np.ptp(press[c.start:c.stop + 1]) >= min_span_dbar:
                dsc, src, err = _cast_attitude(buf.isel(time=slice(c.start, c.stop + 1)),
                                               attitude)
                g = process_cast(dsc, corr_min=corr_min,
                                 boxsize=boxsize, z_max=z_max, direction=look,
                                 min_bin_samples=min_bin_samples,
                                 motion_correct=motion_correct, motion=motion,
                                 bin_average=bin_average, sail=sail, attitude=attitude)
                results.append((g, c, src, err))
        # carry the tail (from the start of the held-back cast) into the next chunk,
        # remembering whether a real gap clipped its start (unknowable next chunk)
        if last and casts:
            held = casts[-1]
            carry_trunc = bool(held.start > 0
                               and t_s[held.start] - t_s[held.start - 1] > gap_s) \
                or (held.start == 0 and bool(carry_trunc))
            buf = buf.isel(time=slice(int(held.start), None))
        elif last:
            # no cast at all in this window (vehicle not profiling): keep the tail so a
            # cast straddling the boundary survives, bounded so an idle stretch can't
            # grow the buffer without limit
            carry_trunc = None
            buf = buf.isel(time=slice(-min(buf.sizes["time"], chunk), None))
        else:
            carry_trunc = None
            buf = None
        if progress:
            print(f"  [{start:>10,},{stop:>10,}) casts_total={len(results):>4d} "
                  f"p={press.min():.0f}..{press.max():.0f} {_time.time()-t0:.0f}s", flush=True)
        start = stop

    ds = _assemble(results, boxsize=boxsize, z_max=z_max, look=look, cast_kind=cast_kind,
                   corr_min=corr_min, min_bin_samples=min_bin_samples,
                   motion_correct=motion_correct, mooring=mooring, source=source,
                   attitude=attitude, motion=motion, bin_average=bin_average, sail=sail)
    ds.attrs["ensemble_range"] = f"{int(ens_start)}:{int(total)}"
    return ds


def _atomic_to_netcdf(ds, path, enc):
    """Write to a temp file then atomically replace `path`. A rebuild then never
    truncates the existing product on failure, and the swap succeeds even if a
    reader (e.g. a Jupyter kernel) still holds the old file open (POSIX rename
    keeps the reader on the old inode)."""
    tmp = f"{path}.tmp{os.getpid()}"
    try:
        ds.to_netcdf(tmp, encoding=enc)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise


def save_l2(ds_l2, path):
    enc = {v: {"zlib": True, "complevel": 4} for v in ds_l2.data_vars}
    _atomic_to_netcdf(ds_l2, str(path), enc)
    return path
