"""Regular depth-time L3 grid for the ADCP **velocity** product (no turbulence L3).

Derived from the (depth, cast) velocity L2: upcasts only, each L2 depth cell binned
onto a regular ``l3_dz_m`` x ``l3_dt`` grid **by its own 2-D sample time** — a buoyant
upcast is slanted in time, so a slow/deep profile is spread across the time bins it
truly spans rather than dropped whole into one cast-time bin. Whole-empty time columns
are then linearly interpolated across short gaps, so L3 is a continuous product.

Fields carried: velE, velN, velU, shearE, shearN, amp (+ n_casts). SEMs are dropped.
Mirrors the CTD `ww_rbr.levels.build_L3`; the two packages stay independent.
"""
from __future__ import annotations

import time as _time

import numpy as np
import pandas as pd
import xarray as xr

from .l2 import _atomic_to_netcdf

L3_VARS = ("velE", "velN", "velU", "shearE", "shearN", "amp")


def _fill_empty_col_gaps(grids, empty, max_gap):
    """Linearly interpolate whole-empty time bins across short gaps, per depth.

    Only bins with ``empty[t]`` True are candidates, and only runs of length
    <= ``max_gap`` bracketed by non-empty bins are filled (no extrapolation at the
    ends, no bridging of longer gaps). A depth that is NaN in a bracketing column
    stays NaN. Returns a new dict of arrays.
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
            a, b = i - 1, j                       # filled bins bracketing the run [i, j-1]
            if (j - i) <= max_gap and a >= 0 and b < nt:
                for t in range(i, j):
                    w = (t - a) / (b - a)
                    for k in out:
                        out[k][:, t] = grids[k][:, a] * (1 - w) + grids[k][:, b] * w
            i = j
        else:
            i += 1
    return out


def build_velocity_l3(l2_path, *, l3_dz_m, l3_dt, interp_max_gap_bins=1):
    """Build the regular depth-time velocity L3 Dataset from a velocity L2 file."""
    l2 = xr.open_dataset(l2_path)
    up = np.flatnonzero(l2["cast_direction"].values == 1)     # upcasts only
    if up.size == 0:
        l2.close()
        raise RuntimeError("no upcasts in the velocity L2; nothing to grid")
    sub = l2.isel(cast=up)
    ncast = sub.sizes["cast"]
    z2 = sub["depth"].values.astype(float)                    # (nz2,) L2 depth centres
    V = {k: sub[k].values for k in L3_VARS}                   # each (nz2, ncast)
    units = {k: l2[k].attrs.get("units", "") for k in L3_VARS}
    lname = {k: l2[k].attrs.get("long_name", k) for k in L3_VARS}
    # per-bin sample time -> ms since epoch (float, NaT -> NaN). A current L2 stores a
    # 2-D (depth, cast) time (slant-correct); a legacy L2 has a 1-D cast mid-time, which
    # we broadcast across depth (whole cast into one time bin, no slant correction).
    tvals = sub["time"].values
    if tvals.ndim == 1:
        tvals = np.broadcast_to(tvals[None, :], (z2.size, ncast))
        slant = False
    else:
        slant = True
    tv = np.ascontiguousarray(tvals).astype("datetime64[ms]")  # (nz2, ncast)
    T2 = tv.view("int64").astype("float64")
    T2[np.isnat(tv)] = np.nan
    boxsize = float(l2.attrs.get("grid_boxsize_m", np.median(np.diff(z2)) if z2.size > 1 else l3_dz_m))
    z_max = float(l2.attrs.get("grid_z_max_m", z2.max() + boxsize / 2))
    l2.close()

    t0w = _time.time()

    # regular L3 depth grid; each L2 depth maps to one L3 depth bin
    dedges = np.arange(0.0, z_max + l3_dz_m, l3_dz_m)
    z3 = 0.5 * (dedges[:-1] + dedges[1:])
    nz3 = z3.size
    dib = np.digitize(z2, dedges) - 1                         # (nz2,) L3 depth bin per L2 depth
    zok = (dib >= 0) & (dib < nz3)

    # regular L3 time axis spanning the actual sample times
    t0 = pd.Timestamp(int(np.nanmin(T2)), unit="ms").floor(l3_dt)
    t1 = pd.Timestamp(int(np.nanmax(T2)), unit="ms").ceil(l3_dt)
    edges = pd.date_range(t0, t1, freq=l3_dt)
    centers = edges[:-1] + (edges[1] - edges[0]) / 2
    ntime = len(edges) - 1
    ee = edges.values.astype("datetime64[ms]").astype("int64")

    # bin each (depth, cast) cell by its own sample time (slant-correct)
    sums = {k: np.zeros((nz3, ntime)) for k in L3_VARS}
    cnts = {k: np.zeros((nz3, ntime)) for k in L3_VARS}
    ncast_bin = np.zeros(ntime, np.int32)                     # distinct upcasts per time bin
    for j in range(ncast):
        tj = T2[:, j]
        finT = np.isfinite(tj)
        tib = np.full(z2.size, -1, np.int64)
        tib[finT] = np.searchsorted(ee, tj[finT].astype("int64"), side="right") - 1
        ok = zok & finT & (tib >= 0) & (tib < ntime)
        if ok.any():
            ncast_bin[np.unique(tib[ok])] += 1
        for k in L3_VARS:
            col = V[k][:, j]
            m = ok & np.isfinite(col)
            np.add.at(sums[k], (dib[m], tib[m]), col[m])
            np.add.at(cnts[k], (dib[m], tib[m]), 1)
    grids = {k: np.where(cnts[k] > 0, sums[k] / np.maximum(cnts[k], 1), np.nan).astype(np.float32)
             for k in L3_VARS}

    # pre-interpolation sparsity (velE as reference), then gap-fill whole-empty columns
    ref = "velE"
    col_empty = cnts[ref].sum(axis=0) == 0
    sparsity_before = float(np.mean(~np.isfinite(grids[ref])) * 100.0)
    empty_cols_before = float(np.mean(col_empty) * 100.0)
    gridsi = _fill_empty_col_gaps(grids, col_empty, interp_max_gap_bins)
    sparsity_after = float(np.mean(~np.isfinite(gridsi[ref])) * 100.0)
    coverage_after = float(np.mean(np.isfinite(gridsi[ref]).any(axis=0)) * 100.0)

    print(f"[L3] {ntime} time bins ({l3_dt}) x {nz3} depths ({l3_dz_m:g} m) from "
          f"{ncast} upcasts; built in {_time.time()-t0w:.1f}s")
    print(f"[L3] before interpolation: matrix {sparsity_before:.1f}% empty "
          f"({empty_cols_before:.1f}% of time bins had no upcast)")
    print(f"[L3] gap-filled across <= {interp_max_gap_bins} empty bin(s): "
          f"matrix {sparsity_after:.1f}% empty; time coverage {coverage_after:.1f}%")

    tcen_ms = centers.values.astype("datetime64[ms]").astype(np.int64)
    data_vars = {k: (("depth", "time"), gridsi[k],
                     {"units": units[k], "long_name": lname[k]}) for k in L3_VARS}
    data_vars["n_casts"] = (("time",), ncast_bin,
                            {"long_name": "distinct upcasts contributing to each time bin",
                             "comment": "0 marks a bin filled by interpolation, if any"})
    return xr.Dataset(
        data_vars,
        coords={
            "depth": ("depth", z3.astype(np.float32),
                      {"units": "m", "positive": "down", "comment": f"{l3_dz_m:g} m bin centres"}),
            "time": ("time", tcen_ms,
                     {"units": "milliseconds since 1970-01-01 00:00:00 UTC",
                      "standard_name": "time", "calendar": "standard",
                      "comment": f"{l3_dt} bin centres"}),
        },
        attrs={
            "title": "Wirewalker ADCP L3 (regular depth-time velocity grid)",
            "derived_from": str(l2_path).split("/")[-1],
            "cast_direction": "up (buoyant ascent)",
            "grid_dz_m": l3_dz_m, "grid_dt": l3_dt,
            "gap_handling": (f"whole-empty time bins linearly interpolated across gaps <= "
                             f"{interp_max_gap_bins} bin(s); longer gaps left as NaN. "
                             f"n_casts==0 marks interpolated bins (mask with n_casts>0 for "
                             f"observed-only)."),
            "interpolation_max_gap_bins": interp_max_gap_bins,
            "pre_interpolation_matrix_sparsity_percent": round(sparsity_before, 2),
            "pre_interpolation_empty_time_bins_percent": round(empty_cols_before, 2),
            "matrix_sparsity_percent": round(sparsity_after, 2),
            "time_coverage_percent": round(coverage_after, 2),
            "time_binning": ("each L2 depth cell binned by its own 2-D sample time "
                             "(a buoyant upcast is slanted in time)") if slant else
                            ("legacy L2 with 1-D cast time: whole cast binned into one "
                             "time bin (no within-cast slant correction)"),
            "processing": "ww_sig1000.l3.build_velocity_l3",
            "date_created": _time.strftime("%Y-%m-%dT%H:%M:%S"),
        },
    )


def save_l3(ds, path):
    enc = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars}
    _atomic_to_netcdf(ds, str(path), enc)
    return path
