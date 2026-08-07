"""Assemble per-cast gridded velocity profiles into an L2 (depth x cast) product.

`build_l2` takes an in-memory dolfyn Dataset (beam coords), detects casts, grids
each selected cast (`velocity.process_cast`), and stacks them into an
`xarray.Dataset` with dims (depth, cast). Mirrors the CTD L2 conventions.
"""
from __future__ import annotations

import time as _time

import numpy as np
import xarray as xr

from .casts import detect_casts
from .geometry import look_direction
from .velocity import process_cast, output_grid


def build_l2(ds, *, boxsize=1.0, z_max=None, cast_kind="both", min_span_dbar=40.0,
             corr_min=50, min_bin_samples=10, thhold_s=30.0, motion_correct=True,
             mooring="", source=""):
    """Grid substantial casts into an L2 (depth x cast) Dataset.

    Parameters
    ----------
    ds : dolfyn Dataset in beam coords (full record or a window).
    cast_kind : 'both' (default), 'up', or 'down'. Each cast is tagged by direction.
    min_span_dbar : only keep casts spanning at least this pressure range.
    """
    fs = float(ds.attrs["fs"])
    t_s = ds["time"].values.astype("datetime64[ns]").astype("int64") / 1e9
    press = ds["pressure"].values
    look = look_direction(ds["pitch"].values, ds["roll"].values)
    if z_max is None:
        z_max = float(np.ceil(np.nanmax(press) / boxsize) * boxsize)
    kinds = ("up", "down") if cast_kind == "both" else (cast_kind,)

    casts = [c for c in detect_casts(press, t_s, thhold=int(thhold_s * fs))
             if c.direction in kinds and np.ptp(press[c.start:c.stop + 1]) >= min_span_dbar]
    if not casts:
        raise RuntimeError(f"no {cast_kind}-casts with span >= {min_span_dbar} dbar found")

    zc = output_grid(boxsize, z_max)
    nz, ncast = zc.size, len(casts)
    G = {k: np.full((nz, ncast), np.nan, np.float32) for k in ("velE", "velN", "velU", "amp")}
    nobs = np.zeros((nz, ncast), np.int32)
    ctime = np.empty(ncast, "datetime64[ns]")
    cpmax = np.zeros(ncast)
    cdir = np.zeros(ncast, np.int8)          # 0 = down, 1 = up

    for j, c in enumerate(casts):
        g = process_cast(ds.isel(time=slice(c.start, c.stop + 1)), corr_min=corr_min,
                         boxsize=boxsize, z_max=z_max, direction=look,
                         min_bin_samples=min_bin_samples, motion_correct=motion_correct)
        for k in G:
            G[k][:, j] = g[k]
        nobs[:, j] = g["n_obs"]
        ctime[j] = np.datetime64(g["time"], "ns")
        cpmax[j] = g["pressure_max"]
        cdir[j] = 1 if c.direction == "up" else 0

    out = xr.Dataset(
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
               "processing": "ww_adcp port of WW_Velocity_Processing_SWOT",
               "date_created": _time.strftime("%Y-%m-%dT%H:%M:%S")},
    )
    return out


def save_l2(ds_l2, path):
    enc = {v: {"zlib": True, "complevel": 4} for v in ds_l2.data_vars}
    ds_l2.to_netcdf(path, encoding=enc)
    return path
