"""L1/L2/L3 product builders for the Wirewalker RBR Concerto CTD.

Archive levels
--------------
L1 : full-resolution (native Hz) converted time series (measured channels +
     depth + cast flags). One NetCDF with an unlimited time dimension.
L2 : UPCASTS ONLY, bin-averaged onto a `grid_dz` m depth grid; thermal-mass
     correction + TEOS-10 quantities computed here. dims (depth, cast).
L3 : regular, gap-filled depth x time grid (l3_dz m x l3_dt); n_casts==0 marks
     interpolated time bins.

Casts are detected from the CTD pressure record (rsk.detect_casts, a time-aware
port of the MATLAB get_upcastRBR), splitting each profile into a DOWN and an UP
cast. Set cast_detection.method="ruskin" to reuse the region/regionCast tables.

Each builder takes a `Config` (see ww_rbr.config) — no module-level state.
"""
from __future__ import annotations

import sqlite3
import time as _time
import warnings

import gsw
import netCDF4
import numpy as np
import xarray as xr

from .rsk import (discover_channels, load_casts, casts_from_pressure,
                  read_cast_data, VAR_META, global_attrs)
from .derive import correct_thermal_mass, convert, buoyancy_n2


# --------------------------------------------------------------------------- #
# L1 : full-resolution converted time series
# --------------------------------------------------------------------------- #
def build_L1(cfg, max_casts: int | None = None):
    """Stream casts and append to a NetCDF with an unlimited time dimension.

    Low memory: one cast in RAM at a time. float32 data vars + zlib compression.
    """
    cfg.l1_path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{cfg.rsk_path}?mode=ro", uri=True)
    chans = discover_channels(con)
    if cfg.cast_method == "ruskin":
        casts = load_casts(con)
        cast_src = "Ruskin region tables"
    else:
        casts = casts_from_pressure(con, chans, cfg)
        cast_src = (f"pressure (slope {cfg.cast_slope_window_s:g}s, debounce "
                    f"{cfg.cast_debounce_window_s:g}s, min span {cfg.cast_min_span_dbar:g} dbar)")
    if max_casts:
        casts = casts[:max_casts]

    extra = [ch for ch in chans if ch["var"] not in ("conductivity", "temperature", "pressure")]
    nup = sum(c.direction == 1 for c in casts)
    print(f"[L1] {len(casts)} casts ({nup} up) detected from {cast_src} -> {cfg.l1_path.name}")
    print(f"[L1] {len(chans)} measured channels from {cfg.rsk_path.name}:")
    for ch in chans:
        tag = " (core CTD)" if ch["var"] in ("conductivity", "temperature", "pressure") else ""
        print(f"       {ch['col']} -> {ch['var']:<28s} [{ch['units']}]  {ch['long_name']}{tag}")

    base_vars = ["conductivity", "temperature", "pressure", "depth"]

    cfg.l1_path.unlink(missing_ok=True)   # avoid EACCES if a stale reader still holds the file
    nc = netCDF4.Dataset(cfg.l1_path, "w", format="NETCDF4")
    for k, v in global_attrs(cfg, "L1").items():
        nc.setncattr(k, v)
    nc.createDimension("time", None)  # unlimited

    vtime = nc.createVariable("time", "f8", ("time",), zlib=True, complevel=4)
    vtime.units = "milliseconds since 1970-01-01 00:00:00 UTC"
    vtime.standard_name = "time"
    vtime.calendar = "standard"

    ncv = {}
    for name in base_vars:
        comment, std, units = VAR_META[name]
        var = nc.createVariable(name, "f4", ("time",), zlib=True, complevel=4,
                                chunksizes=(65536,))
        var.units = units
        var.comment = comment
        if std:
            var.standard_name = std
        ncv[name] = var
    for ch in extra:
        var = nc.createVariable(ch["var"], "f4", ("time",), zlib=True, complevel=4,
                                chunksizes=(65536,))
        var.units = ch["units"]
        var.comment = ch["long_name"]
        var.rsk_channel = ch["short"]
        ncv[ch["var"]] = var

    v_castn = nc.createVariable("cast_number", "i4", ("time",), zlib=True, complevel=4)
    v_castn.comment = "sequential cast index over the deployment (down + up)"
    v_profn = nc.createVariable("profile_number", "i4", ("time",), zlib=True, complevel=4)
    v_profn.comment = "Ruskin instrument-generated profile index"
    v_dir = nc.createVariable("cast_direction", "i1", ("time",), zlib=True, complevel=4)
    v_dir.comment = "cast direction flag"
    v_dir.flag_values = np.array([0, 1], dtype=np.int8)
    v_dir.flag_meanings = "down up"

    pos = 0
    t0 = _time.time()
    for c in casts:
        d = read_cast_data(con, c, chans)
        if not d:
            continue
        n = d["tstamp"].size
        depth = -gsw.z_from_p(d["pressure"] - cfg.atm_dbar, cfg.lat)   # m, positive down
        sl = slice(pos, pos + n)
        vtime[sl] = d["tstamp"]
        ncv["conductivity"][sl] = d["conductivity"]
        ncv["temperature"][sl] = d["temperature"]
        ncv["pressure"][sl] = d["pressure"]
        ncv["depth"][sl] = depth
        for ch in extra:
            ncv[ch["var"]][sl] = d[ch["var"]]
        v_castn[sl] = c.cast_number
        v_profn[sl] = c.profile_number
        v_dir[sl] = c.direction
        pos += n
        if c.cast_number % 1000 == 0:
            print(f"  cast {c.cast_number:>6d}/{len(casts)}  samples={pos:>10d}  {_time.time()-t0:5.1f}s")

    nc.setncattr("n_samples", pos)
    nc.setncattr("n_casts", len(casts))
    nc.setncattr("cast_detection", cast_src)
    nc.close()
    con.close()
    print(f"[L1] done: {pos} samples in {_time.time()-t0:.1f}s -> {cfg.l1_path}")


# --------------------------------------------------------------------------- #
# L2 : upcasts gridded to grid_dz m
# --------------------------------------------------------------------------- #
def build_L2(cfg, max_casts: int | None = None):
    """Bin-average every UPCAST from L1 onto the depth grid -> NetCDF (depth, cast).

    Derived from L1 (not the raw .rsk): inherits L1's channels (core C/T/P plus any
    extra measured sensors, bin-averaged raw). Thermal-mass correction is applied
    per upcast before deriving salinity; TEOS-10 quantities are computed here.
    """
    if not cfg.l1_path.exists():
        raise FileNotFoundError(
            f"L2 is built from L1, but {cfg.l1_path} is missing. Run --level 1 first.")
    cfg.l2_path.parent.mkdir(parents=True, exist_ok=True)

    l1 = xr.open_dataset(cfg.l1_path)
    cond_all = l1["conductivity"].values
    temp_all = l1["temperature"].values
    pres_all = l1["pressure"].values
    cn_all = l1["cast_number"].values
    cd_all = l1["cast_direction"].values
    pn_all = l1["profile_number"].values
    tms_all = l1["time"].values.astype("datetime64[ms]").astype(np.int64)  # ms since epoch
    _core = {"conductivity", "temperature", "pressure", "depth",
             "cast_number", "profile_number", "cast_direction"}
    extra_vars = [v for v in l1.data_vars if v not in _core]
    extra_all = {v: l1[v].values for v in extra_vars}
    extra_meta = {v: (l1[v].attrs.get("comment", ""), None, l1[v].attrs.get("units", ""))
                  for v in extra_vars}
    l1.close()

    # Sampling rate for the thermal-mass correction: the RBR logs at a steady rate, so
    # take it from the record (median sample interval) rather than requiring it in the
    # config. cfg.fs, if set, is an explicit override. (Inter-cast gaps are rare, so the
    # median is the within-cast rate.)
    if cfg.fs:
        fs = float(cfg.fs)
    else:
        fs = 1000.0 / float(np.median(np.diff(tms_all)))
    print(f"[L2] sampling rate {fs:.4g} Hz ({'config override' if cfg.fs else 'from record'})")

    # Each cast is a contiguous block of equal cast_number (L1 is time-ordered).
    bnds = np.flatnonzero(np.diff(cn_all) != 0) + 1
    starts = np.concatenate(([0], bnds))
    ends = np.concatenate((bnds, [cn_all.size]))
    casts = [(s, e) for s, e in zip(starts, ends) if cd_all[s] == 1]  # upcasts only
    if max_casts:
        casts = casts[:max_casts]
    print(f"[L2] {len(casts)} upcasts from {cfg.l1_path.name} -> {cfg.l2_path.name}")

    edges = np.arange(cfg.grid_zmin, cfg.grid_zmax + cfg.grid_dz, cfg.grid_dz)
    zc = 0.5 * (edges[:-1] + edges[1:])     # bin centres
    nz = zc.size
    ncast = len(casts)

    gridvars = ["conductivity", "temperature", "pressure", "practical_salinity",
                "absolute_salinity", "conservative_temperature", "sigma0",
                "sound_speed"] + extra_vars
    grids = {k: np.full((nz, ncast), np.nan, np.float32) for k in gridvars}
    nobs = np.zeros((nz, ncast), np.int32)
    # time is 2-D like every other variable: a buoyant upcast is slanted in time
    # (deep bins sampled before shallow ones), so each (depth, cast) bin carries its
    # own mean sample time. L3, being a regular grid, later collapses time to 1-D.
    gtime = np.full((nz, ncast), np.nan)    # ms, bin-mean sample time (depth, cast)
    pcastn = np.zeros(ncast, np.int32)
    pprofn = np.zeros(ncast, np.int32)

    t0 = _time.time()
    for j, (s, e) in enumerate(casts):
        temp = temp_all[s:e]
        pres = pres_all[s:e]
        cond = correct_thermal_mass(cond_all[s:e], temp,
                                    cfg.tm_alpha, cfg.tm_beta, cfg.tm_gamma, fs)
        der = convert(cond, temp, pres, cfg)
        depth = der["depth"]
        col = {
            "conductivity": cond, "temperature": temp, "pressure": pres,
            "practical_salinity": der["practical_salinity"],
            "absolute_salinity": der["absolute_salinity"],
            "conservative_temperature": der["conservative_temperature"],
            "sigma0": der["sigma0"], "sound_speed": der["sound_speed"],
        }
        for v in extra_vars:
            col[v] = extra_all[v][s:e]
        ib = np.digitize(depth, edges) - 1            # bin index per sample
        valid = (ib >= 0) & (ib < nz) & np.isfinite(depth)
        ibv = ib[valid]
        nobs[:, j] = np.bincount(ibv, minlength=nz)
        with np.errstate(invalid="ignore"):
            for k, v in col.items():
                vv = v[valid]
                m = np.isfinite(vv)
                ssum = np.bincount(ibv[m], weights=vv[m], minlength=nz)
                cnt = np.bincount(ibv[m], minlength=nz).astype(float)
                cnt[cnt == 0] = np.nan
                grids[k][:, j] = (ssum / cnt).astype(np.float32)
            # mean sample time in each depth bin (same bins as the data above)
            tsum = np.bincount(ibv, weights=tms_all[s:e].astype(float)[valid], minlength=nz)
            tcnt = nobs[:, j].astype(float)
            tcnt[tcnt == 0] = np.nan
            gtime[:, j] = tsum / tcnt
        pcastn[j] = cn_all[s]
        pprofn[j] = pn_all[s]
        if j % 1000 == 0:
            print(f"  upcast {j:>6d}/{ncast}  {_time.time()-t0:5.1f}s")

    data_vars = {}
    for k in gridvars:
        comment, std, units = VAR_META.get(k, extra_meta.get(k, ("", None, "")))
        attrs = {"units": units, "comment": comment}
        if std:
            attrs["standard_name"] = std
        data_vars[k] = (("depth", "cast"), grids[k], attrs)
    data_vars["n_obs"] = (("depth", "cast"), nobs,
                          {"comment": "number of samples averaged in each bin"})

    ds = xr.Dataset(
        data_vars,
        coords={
            "depth": ("depth", zc.astype(np.float32),
                      {"units": "m", "standard_name": "depth", "positive": "down"}),
            "cast": ("cast", np.arange(ncast, dtype=np.int32)),
            "time": (("depth", "cast"), gtime,
                     {"units": "milliseconds since 1970-01-01 00:00:00 UTC",
                      "standard_name": "time", "calendar": "standard",
                      "comment": "mean sample time in each (depth, cast) bin; 2-D because a "
                                 "buoyant upcast is slanted in time (deep sampled before shallow)"}),
            "cast_number": ("cast", pcastn, {"comment": "sequential cast index in L1"}),
            "profile_number": ("cast", pprofn, {"comment": "Ruskin profile index"}),
        },
        attrs={**global_attrs(cfg, "L2"),
               "derived_from": cfg.l1_path.name,
               "cast_direction": "up (buoyant ascent)",
               "grid_dz_m": cfg.grid_dz, "grid_zmin_m": cfg.grid_zmin, "grid_zmax_m": cfg.grid_zmax,
               "sampling_hz": round(fs, 4),
               "sampling_hz_source": "config override" if cfg.fs else "derived from record",
               "thermal_mass_correction": "Lueck & Picklo (1990), RBR pyRSKtools correctTM",
               "thermal_mass_alpha": cfg.tm_alpha, "thermal_mass_beta_per_s": cfg.tm_beta,
               "thermal_mass_gamma": cfg.tm_gamma,
               "ct_alignment_lag_s": 0.0,
               "conductivity_note": "thermal-mass-corrected; salinity derived from it"},
    )
    enc = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars}
    cfg.l2_path.unlink(missing_ok=True)   # avoid EACCES if a stale reader still holds the file
    ds.to_netcdf(cfg.l2_path, encoding=enc)
    print(f"[L2] done: {ncast} upcasts x {nz} levels in {_time.time()-t0:.1f}s -> {cfg.l2_path}")


# --------------------------------------------------------------------------- #
# L3 : regular depth x time grid
# --------------------------------------------------------------------------- #
def _fill_empty_col_gaps(grids, empty, max_gap):
    """Linearly interpolate empty time bins across short gaps, per depth.

    Only whole-empty time bins (`empty[t]` True) are candidates, and only runs of
    length <= `max_gap` bracketed by non-empty bins are filled (no extrapolation at
    the ends, no bridging of longer gaps). A deep bin that is NaN in a bracketing
    column stays NaN. Returns a new dict of arrays.
    """
    out = {k: v.copy() for k, v in grids.items()}
    empty = np.asarray(empty, bool)
    nt = empty.size
    i = 0
    while i < nt:
        if empty[i]:
            j = i
            while j < nt and empty[j]:
                j += 1
            a, b = i - 1, j           # filled bins bracketing the run [i, j-1]
            run = j - i
            if run <= max_gap and a >= 0 and b < nt:
                for t in range(i, j):
                    w = (t - a) / (b - a)
                    for k in out:
                        out[k][:, t] = grids[k][:, a] * (1 - w) + grids[k][:, b] * w
            i = j
        else:
            i += 1
    return out


def build_L3(cfg):
    """Regular, gap-filled depth x time grid from L2 (l3_dz m x l3_dt).

    Vertical: L2's grid_dz bins are averaged in consecutive pairs to l3_dz.
    Temporal: upcasts bin-averaged into fixed l3_dt bins by cast time. L3 is the
    continuous product, so whole-empty time bins are linearly interpolated across
    short gaps (<= l3_interp_max_gap_bins); longer gaps stay NaN. The pre-
    interpolation sparsity is reported and stored on the product; n_casts==0 marks
    interpolated bins (mask with n_casts>0 for observed-only).
    """
    import pandas as pd
    if not cfg.l2_path.exists():
        raise FileNotFoundError(
            f"L3 is built from L2, but {cfg.l2_path} is missing. Run --level 2 first.")
    cfg.l3_path.parent.mkdir(parents=True, exist_ok=True)

    l2 = xr.open_dataset(cfg.l2_path)
    # In L2 every field is 2-D (depth, cast); in the regular L3 grid pressure collapses
    # to a 1-D vertical vector and time to the 1-D regular axis, so both are handled
    # separately from the gridded data variables.
    gridvars = [v for v in l2.data_vars if v not in ("n_obs", "pressure")]
    l2meta = {v: (l2[v].attrs.get("comment", ""),
                  l2[v].attrs.get("standard_name"),
                  l2[v].attrs.get("units", "")) for v in gridvars}
    ncast = l2.sizes["cast"]
    nz2 = l2.sizes["depth"]
    assert nz2 % 2 == 0, "expected an even number of grid_dz bins"
    nz3 = nz2 // 2
    z3 = l2["depth"].values.reshape(nz3, 2).mean(axis=1).astype(np.float32)  # l3_dz centres

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "Mean of empty slice")
        warnings.filterwarnings("ignore", "invalid value encountered")
        V = {k: np.nanmean(l2[k].values.reshape(nz3, 2, ncast), axis=1)     # (nz3, ncast)
             for k in gridvars}
        # 1-D pressure vector for the regular grid: mean bin pressure over all casts
        pres_1d = np.nanmean(l2["pressure"].values.reshape(nz3, 2, ncast),
                             axis=(1, 2)).astype(np.float32)
        # per-(depth, cast) sample time as ms since epoch (float, NaT -> NaN), paired to
        # L3 depths. Work in ms (not ns) so magnitudes stay well within float64 precision.
        tv = l2["time"].values.astype("datetime64[ms]")                     # (nz2, ncast)
        T2 = tv.view("int64").astype("float64")
        T2[np.isnat(tv)] = np.nan
        T3 = np.nanmean(T2.reshape(nz3, 2, ncast), axis=1)                  # (nz3, ncast) ms
    l2.close()

    # Regular time axis spanning the actual sample times.
    t0 = pd.Timestamp(int(np.nanmin(T3)), unit="ms").floor(cfg.l3_dt)
    t1 = pd.Timestamp(int(np.nanmax(T3)), unit="ms").ceil(cfg.l3_dt)
    edges = pd.date_range(t0, t1, freq=cfg.l3_dt)
    centers = (edges[:-1] + (edges[1] - edges[0]) / 2)
    ntime = len(edges) - 1
    ee = edges.values.astype("datetime64[ms]").astype("int64")             # ms since epoch

    # Bin each (depth, cast) cell into the time bin of *its own* sample time. A buoyant
    # upcast is slanted in time, so a slow/deep profile that spans several l3_dt bins is
    # spread across them correctly instead of being dropped whole into one cast-time bin.
    sums = {k: np.zeros((nz3, ntime)) for k in gridvars}
    cnts = {k: np.zeros((nz3, ntime)) for k in gridvars}
    ncast_bin = np.zeros(ntime, np.int32)          # distinct upcasts contributing per bin
    rows = np.arange(nz3)
    t0w = _time.time()
    for j in range(ncast):
        tj = T3[:, j]
        finT = np.isfinite(tj)
        bins = np.full(nz3, -1, np.int64)
        bins[finT] = np.searchsorted(ee, tj[finT].astype("int64"), side="right") - 1
        ok = finT & (bins >= 0) & (bins < ntime)
        if ok.any():
            ncast_bin[np.unique(bins[ok])] += 1
        for k in gridvars:
            col = V[k][:, j]
            m = ok & np.isfinite(col)
            np.add.at(sums[k], (rows[m], bins[m]), col[m])
            np.add.at(cnts[k], (rows[m], bins[m]), 1)
    grids = {k: np.where(cnts[k] > 0, sums[k] / np.maximum(cnts[k], 1), np.nan).astype(np.float32)
             for k in gridvars}

    # Sparsity of the raw (depth, time) matrix, before any temporal interpolation.
    # Two views: cell-level (fraction of the whole grid that is empty, using a core
    # variable as reference) and column-level (time bins with no observation at all).
    ref = "practical_salinity" if "practical_salinity" in gridvars else gridvars[0]
    col_empty = cnts[ref].sum(axis=0) == 0
    sparsity_before = float(np.mean(~np.isfinite(grids[ref])) * 100.0)
    empty_cols_before = float(np.mean(col_empty) * 100.0)

    # L3's purpose is a *continuous* product: linearly interpolate whole-empty time
    # bins across short gaps (n_casts==0 marks them); longer gaps stay NaN.
    gridsi = _fill_empty_col_gaps(grids, col_empty, cfg.l3_interp_maxgap)
    sparsity_after = float(np.mean(~np.isfinite(gridsi[ref])) * 100.0)
    coverage_after = float(np.mean(np.isfinite(gridsi[ref]).any(axis=0)) * 100.0)

    print(f"[L3] {ntime} time bins ({cfg.l3_dt}) x {nz3} depths ({cfg.l3_dz:g} m); "
          f"built in {_time.time()-t0w:.1f}s")
    print(f"[L3] before interpolation: matrix {sparsity_before:.1f}% empty "
          f"({empty_cols_before:.1f}% of time bins had no upcast)")
    print(f"[L3] gap-filled across <= {cfg.l3_interp_maxgap} empty bin(s): "
          f"matrix {sparsity_after:.1f}% empty; time coverage {coverage_after:.1f}%")

    tcen_ms = centers.values.astype("datetime64[ms]").astype(np.int64)

    def _write(path, gr, extra_attrs):
        data_vars = {}
        for k in gridvars:
            comment, std, units = VAR_META.get(k, l2meta.get(k, ("", None, "")))
            attrs = {"units": units, "comment": comment}
            if std:
                attrs["standard_name"] = std
            data_vars[k] = (("depth", "time"), gr[k], attrs)
        n2 = buoyancy_n2(gr["sigma0"], z3, cfg.n2_smooth_m, cfg.gravity)
        data_vars["buoyancy_frequency_squared"] = (
            ("depth", "time"), n2.astype(np.float32),
            {"units": "s-2", "long_name": "square of buoyancy (Brunt-Vaisala) frequency",
             "standard_name": "square_of_brunt_vaisala_frequency_in_sea_water",
             "comment": f"(g/rho0) d(sigma0)/dz; sigma0 vertically smoothed with a "
                        f"{cfg.n2_smooth_m:g} m boxcar before differencing; z positive down"})
        data_vars["n_casts"] = (("time",), ncast_bin,
                                {"comment": "number of distinct upcasts contributing to each time "
                                            "bin (0 marks a bin filled by interpolation, if any)"})
        ds = xr.Dataset(
            data_vars,
            coords={
                "depth": ("depth", z3, {"units": "m", "standard_name": "depth",
                                        "positive": "down", "comment": f"{cfg.l3_dz:g} m bin centres"}),
                "pressure": ("depth", pres_1d,
                             {"units": "dbar", "standard_name": "sea_water_pressure",
                              "comment": "mean in-bin pressure over all casts (regular grid: "
                                         "pressure is ~1-D in depth)"}),
                "time": ("time", tcen_ms,
                         {"units": "milliseconds since 1970-01-01 00:00:00 UTC",
                          "standard_name": "time", "calendar": "standard",
                          "comment": f"{cfg.l3_dt} bin centres"}),
            },
            attrs={**global_attrs(cfg, "L3"), "derived_from": cfg.l2_path.name,
                   "cast_direction": "up (buoyant ascent)",
                   "grid_dz_m": cfg.l3_dz, "grid_dt": cfg.l3_dt, "nyquist_frequency_cph": 1.0,
                   "n2_vertical_smoothing_m": cfg.n2_smooth_m,
                   **extra_attrs},
        )
        enc = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars}
        path.unlink(missing_ok=True)   # avoid EACCES if a stale reader still holds the file
        ds.to_netcdf(path, encoding=enc)

    _write(cfg.l3_path, gridsi, {
        "gap_handling": (f"empty time bins linearly interpolated across gaps <= "
                         f"{cfg.l3_interp_maxgap} bin(s); longer gaps left as NaN. "
                         f"n_casts==0 marks interpolated bins (mask with n_casts>0 for "
                         f"observed-only)."),
        "interpolation_max_gap_bins": cfg.l3_interp_maxgap,
        "pre_interpolation_matrix_sparsity_percent": round(sparsity_before, 2),
        "pre_interpolation_empty_time_bins_percent": round(empty_cols_before, 2),
        "matrix_sparsity_percent": round(sparsity_after, 2),
        "time_coverage_percent": round(coverage_after, 2),
    })
    print(f"[L3] done -> {cfg.l3_path}")
