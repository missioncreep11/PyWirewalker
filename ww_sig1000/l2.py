"""Assemble per-cast gridded velocity profiles into an L2 (depth x cast) product.

`build_l2` grids an in-memory dolfyn Dataset (beam coords); `build_l2_streaming`
does the same over a raw `.ad2cp` too large to hold in memory, by reading the
file in chunks and carrying the boundary cast across chunks. Both stack per-cast
profiles (`velocity.process_cast`) into an `xarray.Dataset` with dims
(depth, cast), mirroring the CTD L2 conventions.
"""
from __future__ import annotations

import time as _time

import numpy as np
import xarray as xr

from .casts import detect_casts
from .geometry import look_direction
from .velocity import process_cast, output_grid


def _assemble(results, *, boxsize, z_max, look, cast_kind, corr_min, min_bin_samples,
              motion_correct, mooring, source):
    """Stack a list of (cast_result, Cast) into an L2 xarray.Dataset."""
    zc = output_grid(boxsize, z_max)
    nz, ncast = zc.size, len(results)
    G = {k: np.full((nz, ncast), np.nan, np.float32) for k in ("velE", "velN", "velU", "amp")}
    nobs = np.zeros((nz, ncast), np.int32)
    ctime = np.empty(ncast, "datetime64[ns]")
    cpmax = np.zeros(ncast)
    cpmin = np.zeros(ncast)
    cdir = np.zeros(ncast, np.int8)
    ccomplete = np.zeros(ncast, np.int8)
    for j, (g, cast) in enumerate(results):
        for k in G:
            G[k][:, j] = g[k]
        nobs[:, j] = g["n_obs"]
        ctime[j] = np.datetime64(g["time"], "ns")
        cpmax[j] = g["pressure_max"]
        cpmin[j] = g["pressure_min"]
        cdir[j] = 1 if cast.direction == "up" else 0
        ccomplete[j] = 0 if cast.truncated else 1

    order = np.argsort(ctime)                      # ensure time order
    for k in G:
        G[k] = G[k][:, order]
    nobs, ctime, cdir = nobs[:, order], ctime[order], cdir[order]
    cpmax, cpmin, ccomplete = cpmax[order], cpmin[order], ccomplete[order]
    n_trunc = int((ccomplete == 0).sum())

    return xr.Dataset(
        {"velE": (("depth", "cast"), G["velE"], {"units": "m s-1", "long_name": "eastward velocity"}),
         "velN": (("depth", "cast"), G["velN"], {"units": "m s-1", "long_name": "northward velocity"}),
         "velU": (("depth", "cast"), G["velU"], {"units": "m s-1", "long_name": "upward velocity"}),
         "amp": (("depth", "cast"), G["amp"], {"units": "dB", "long_name": "beam-mean backscatter amplitude"}),
         "n_obs": (("depth", "cast"), nobs, {"long_name": "samples averaged per bin"})},
        coords={"depth": ("depth", zc.astype(np.float32), {"units": "m", "positive": "down"}),
                "cast": ("cast", np.arange(ncast, dtype=np.int32)),
                "time": ("cast", ctime, {"long_name": "cast mid-time"}),
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
                                                 "turning points"})},
        attrs={"title": f"{mooring} Wirewalker ADCP L2 (gridded velocity)",
               "mooring": mooring, "source_file": source,
               "instrument_look": look, "cast_kind": cast_kind,
               "n_casts_truncated": n_trunc,
               "grid_boxsize_m": boxsize, "grid_z_max_m": z_max,
               "corr_min": corr_min, "min_bin_samples": min_bin_samples,
               "motion_correction": ("WWcorr_beam (bandpass-integrated IMU + dp/dt)"
                                     if motion_correct else "none"),
               "processing": "ww_sig1000 port of WW_Velocity_Processing_SWOT",
               "date_created": _time.strftime("%Y-%m-%dT%H:%M:%S")},
    )


def _select_casts(press, t_s, fs, kinds, thhold_s, min_span_dbar):
    return [c for c in detect_casts(press, t_s, thhold=int(thhold_s * fs))
            if c.direction in kinds and np.ptp(press[c.start:c.stop + 1]) >= min_span_dbar]


def build_l2(ds, *, boxsize=1.0, z_max=None, cast_kind="both", min_span_dbar=40.0,
             corr_min=50, min_bin_samples=10, thhold_s=30.0, motion_correct=True,
             mooring="", source=""):
    """Grid substantial casts of an in-memory Dataset into an L2 Dataset."""
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
    results = [(process_cast(ds.isel(time=slice(c.start, c.stop + 1)), corr_min=corr_min,
                             boxsize=boxsize, z_max=z_max, direction=look,
                             min_bin_samples=min_bin_samples, motion_correct=motion_correct),
                c) for c in casts]
    return _assemble(results, boxsize=boxsize, z_max=z_max, look=look, cast_kind=cast_kind,
                     corr_min=corr_min, min_bin_samples=min_bin_samples,
                     motion_correct=motion_correct, mooring=mooring, source=source)


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
                       mooring="", source="", progress=True):
    """Grid a raw `.ad2cp` too large for memory. `reader` is dolfyn.read.

    Reads the file in `chunk`-ensemble windows, detects/grids casts on a rolling
    buffer, and carries the last (possibly boundary-spanning) cast into the next
    chunk. z_max, if not given, is set from the first chunk's max pressure + range.

    `ens_start` / `total` bound the ensemble range, so the deployment/recovery
    transit — when the vehicle is on deck or being lowered and is not profiling —
    can be trimmed off rather than entering the product as spurious casts.
    """
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
                g = process_cast(buf.isel(time=slice(c.start, c.stop + 1)), corr_min=corr_min,
                                 boxsize=boxsize, z_max=z_max, direction=look,
                                 min_bin_samples=min_bin_samples, motion_correct=motion_correct)
                results.append((g, c))
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
                   motion_correct=motion_correct, mooring=mooring, source=source)
    ds.attrs["ensemble_range"] = f"{int(ens_start)}:{int(total)}"
    return ds


def save_l2(ds_l2, path):
    enc = {v: {"zlib": True, "complevel": 4} for v in ds_l2.data_vars}
    ds_l2.to_netcdf(path, encoding=enc)
    return path
