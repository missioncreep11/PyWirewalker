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
    """Stack a list of (cast_result, direction) into an L2 xarray.Dataset."""
    zc = output_grid(boxsize, z_max)
    nz, ncast = zc.size, len(results)
    G = {k: np.full((nz, ncast), np.nan, np.float32) for k in ("velE", "velN", "velU", "amp")}
    nobs = np.zeros((nz, ncast), np.int32)
    ctime = np.empty(ncast, "datetime64[ns]")
    cpmax = np.zeros(ncast)
    cdir = np.zeros(ncast, np.int8)
    for j, (g, direction) in enumerate(results):
        for k in G:
            G[k][:, j] = g[k]
        nobs[:, j] = g["n_obs"]
        ctime[j] = np.datetime64(g["time"], "ns")
        cpmax[j] = g["pressure_max"]
        cdir[j] = 1 if direction == "up" else 0

    order = np.argsort(ctime)                      # ensure time order
    for k in G:
        G[k] = G[k][:, order]
    nobs, ctime, cpmax, cdir = nobs[:, order], ctime[order], cpmax[order], cdir[order]

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
                "cast_direction": ("cast", cdir,
                                   {"flag_values": np.array([0, 1], np.int8),
                                    "flag_meanings": "down up"})},
        attrs={"title": f"{mooring} Wirewalker ADCP L2 (gridded velocity)",
               "mooring": mooring, "source_file": source,
               "instrument_look": look, "cast_kind": cast_kind,
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
                c.direction) for c in casts]
    return _assemble(results, boxsize=boxsize, z_max=z_max, look=look, cast_kind=cast_kind,
                     corr_min=corr_min, min_bin_samples=min_bin_samples,
                     motion_correct=motion_correct, mooring=mooring, source=source)


def _count_ensembles(fn, reader):
    """Total burst ensembles in `fn`, via exponential probe + bisection (tiny reads)."""
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


def build_l2_streaming(fn, reader, *, chunk=500_000, total=None, boxsize=1.0, z_max=None,
                       cast_kind="both", min_span_dbar=40.0, corr_min=50, min_bin_samples=10,
                       thhold_s=30.0, motion_correct=True, mooring="", source="", progress=True):
    """Grid a raw `.ad2cp` too large for memory. `reader` is dolfyn.read.

    Reads the file in `chunk`-ensemble windows, detects/grids casts on a rolling
    buffer, and carries the last (possibly boundary-spanning) cast into the next
    chunk. z_max, if not given, is set from the first chunk's max pressure + range.
    """
    if total is None:
        total = _count_ensembles(fn, reader)
    kinds = ("up", "down") if cast_kind == "both" else (cast_kind,)
    results, look, buf = [], None, None
    start, t0 = 0, _time.time()

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

        casts = detect_casts(press, t_s, thhold=int(thhold_s * fs))
        last = start + chunk < total       # more data coming -> hold back the final cast
        to_do = casts[:-1] if (last and casts) else casts
        for c in to_do:
            if c.direction in kinds and np.ptp(press[c.start:c.stop + 1]) >= min_span_dbar:
                g = process_cast(buf.isel(time=slice(c.start, c.stop + 1)), corr_min=corr_min,
                                 boxsize=boxsize, z_max=z_max, direction=look,
                                 min_bin_samples=min_bin_samples, motion_correct=motion_correct)
                results.append((g, c.direction))
        # carry the tail (from the start of the held-back cast) into the next chunk
        if last and casts:
            buf = buf.isel(time=slice(int(casts[-1].start), None))
        else:
            buf = None
        if progress:
            print(f"  [{start:>10,},{stop:>10,}) casts_total={len(results):>4d} "
                  f"p={press.min():.0f}..{press.max():.0f} {_time.time()-t0:.0f}s", flush=True)
        start = stop

    return _assemble(results, boxsize=boxsize, z_max=z_max, look=look, cast_kind=cast_kind,
                     corr_min=corr_min, min_bin_samples=min_bin_samples,
                     motion_correct=motion_correct, mooring=mooring, source=source)


def save_l2(ds_l2, path):
    enc = {v: {"zlib": True, "complevel": 4} for v in ds_l2.data_vars}
    ds_l2.to_netcdf(path, encoding=enc)
    return path
