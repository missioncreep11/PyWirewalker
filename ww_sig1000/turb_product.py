"""Assemble per-cast HR beam-5 dissipation into a (depth x cast) turbulence product.

`build_turbulence_streaming` streams a raw `.ad2cp` too large to hold in memory,
detecting upcasts on a rolling buffer and running `turbulence.process_cast_turbulence`
(the ProcessSingleProfile.m spectral method) on each. Unlike the velocity L2 it keeps
the HR beam-5 fields (vel_b5/corr_b5 on time_b5/range_b5, 1:1 aligned to the burst
pings) and drops the 4-beam data. Output mirrors NortekTurbulenceData.nc: dims
(depth, cast) with epsilon, N (noise floor), SNR, A, corr and num_spectra.
"""
from __future__ import annotations

import time as _time

import numpy as np
import xarray as xr

from .casts import detect_casts
from .l2 import _count_ensembles
from .turbulence import process_cast_turbulence, hr_bins

_VARS = ("eps", "N", "SNR", "A", "corr", "num")


def _assemble(results, *, cellsize, dep_res, max_dep, corr_min, mooring, source):
    """Stack a list of per-cast turbulence dicts into a (depth, cast) Dataset."""
    zc, _ = hr_bins(cellsize, dep_res, max_dep)
    nz, ncast = zc.size, len(results)
    G = {k: np.full((nz, ncast), np.nan, np.float32) for k in _VARS}
    ctime = np.empty(ncast, "datetime64[ns]")
    for j, r in enumerate(results):
        for k in _VARS:
            G[k][:, j] = r[k]
        ctime[j] = r["time"]

    order = np.argsort(ctime)
    for k in G:
        G[k] = G[k][:, order]
    ctime = ctime[order]

    return xr.Dataset(
        {"epsilon": (("depth", "cast"), G["eps"],
                     {"units": "W kg-1", "long_name": "TKE dissipation rate (spectral)"}),
         "N": (("depth", "cast"), G["N"],
               {"units": "(m s-1)2 rad-1 m", "long_name": "spectral noise floor"}),
         "SNR": (("depth", "cast"), G["SNR"], {"long_name": "A/N signal-to-noise ratio"}),
         "A": (("depth", "cast"), G["A"],
               {"long_name": "inertial-subrange amplitude (Ck eps^2/3)"}),
         "corr": (("depth", "cast"), G["corr"],
                  {"units": "percent", "long_name": "beam-5 correlation, bin mean"}),
         "num_spectra": (("depth", "cast"), G["num"].astype(np.float32),
                         {"long_name": "spectra averaged per bin"})},
        coords={"depth": ("depth", zc.astype(np.float32), {"units": "m", "positive": "down"}),
                "cast": ("cast", np.arange(ncast, dtype=np.int32)),
                "time": ("cast", ctime, {"long_name": "cast mid-time"})},
        attrs={"title": f"{mooring} Wirewalker ADCP HR turbulence (spectral eps)",
               "mooring": mooring, "source_file": source,
               "method": "wavenumber-spectrum, S(k)=N+A k^-5/3, eps=(A/0.53)^3/2",
               "dep_res_m": dep_res, "max_dep_m": max_dep, "corr_min": corr_min,
               "processing": "ww_sig1000 port of ProcessSingleProfile.m (Northcott et al. 2026)",
               "date_created": _time.strftime("%Y-%m-%dT%H:%M:%S")},
    )


def build_turbulence_streaming(fn, reader, *, chunk=500_000, total=None,
                               dep_res=3.0, max_dep=100.0, cast_kind="up",
                               min_span_dbar=40.0, corr_min=50, thhold_s=30.0,
                               mooring="", source="", progress=True):
    """Stream a raw `.ad2cp`, producing a (depth, cast) turbulence Dataset.

    `reader` is dolfyn.read. Casts are detected on buffered pressure; the last
    (possibly boundary-spanning) cast is carried into the next chunk. Beam-5
    geometry (cell_size_b5, blank_dist_b5, ambig_vel_b5) is read from the file.
    """
    if total is None:
        total = _count_ensembles(fn, reader)
    kinds = ("up", "down") if cast_kind == "both" else (cast_kind,)
    results, buf, geom = [], None, None
    start, t0 = 0, _time.time()

    while start < total:
        stop = min(start + chunk, total)
        ds = reader(fn, nens=[start, stop])
        if ds.sizes["time"] != ds.sizes.get("time_b5", ds.sizes["time"]):
            raise RuntimeError("time and time_b5 lengths differ; beam-5 not 1:1 aligned")
        if geom is None:
            a = ds.attrs
            geom = dict(fs=float(a["fs"]), cs=float(a["cell_size_b5"]),
                        bd=float(a["blank_dist_b5"]), va=float(a["ambig_vel_b5"]))
        cur = {"press": ds["pressure"].values, "pitch": ds["pitch"].values,
               "roll": ds["roll"].values,
               "t": ds["time"].values.astype("datetime64[ns]"),
               "vel": ds["vel_b5"].values, "corr": ds["corr_b5"].values}   # (cells, pings)
        if buf is None:
            buf = cur
        else:
            buf = {k: np.concatenate([buf[k], cur[k]], axis=(1 if cur[k].ndim == 2 else 0))
                   for k in cur}

        press = buf["press"]
        t_s = buf["t"].astype("int64") / 1e9
        casts = detect_casts(press, t_s, thhold=int(thhold_s * geom["fs"]))
        last = stop < total                       # hold back final cast if more coming
        to_do = casts[:-1] if (last and casts) else casts
        for c in to_do:
            if c.direction in kinds and np.ptp(press[c.start:c.stop + 1]) >= min_span_dbar:
                sl = slice(c.start, c.stop + 1)
                r = process_cast_turbulence(
                    buf["vel"][:, sl].T, buf["corr"][:, sl].T, press[sl],
                    buf["pitch"][sl], buf["roll"][sl], geom["va"],
                    cellsize=geom["cs"], blockdis=geom["bd"], dep_res=dep_res,
                    max_dep=max_dep, corr_min=corr_min,
                    time=buf["t"][sl][(c.stop - c.start) // 2])
                if r is not None:
                    results.append(r)
        if last and casts:
            s0 = int(casts[-1].start)
            buf = {k: (buf[k][:, s0:] if buf[k].ndim == 2 else buf[k][s0:]) for k in buf}
        else:
            buf = None
        if progress:
            print(f"  [{start:>10,},{stop:>10,}) casts_total={len(results):>4d} "
                  f"p={press.min():.0f}..{press.max():.0f} {_time.time()-t0:.0f}s", flush=True)
        start = stop

    if not results:
        raise RuntimeError(f"no {cast_kind}-casts with span >= {min_span_dbar} dbar found")
    return _assemble(results, cellsize=geom["cs"], dep_res=dep_res, max_dep=max_dep,
                     corr_min=corr_min, mooring=mooring, source=source)


def save_turbulence(ds, path):
    enc = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars}
    ds.to_netcdf(path, encoding=enc)
    return path
