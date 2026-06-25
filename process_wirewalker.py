#!/usr/bin/env python3
"""
Wirewalker RBR Concerto CTD processing workflow
===============================================

Archive levels
--------------
L0 : raw RBR binary  ->  the .rsk file (SQLite). Nothing to build; it is the original.
L1 : NetCDF of the full-resolution (2 Hz) converted time series.
     C/T/P converted to practical & absolute salinity, conservative temperature,
     potential density (sigma0), sound speed and depth via TEOS-10 (gsw).
     EVERY sample is tagged with:
         profile_number   - Ruskin instrument-generated profile index
         cast_number      - sequential index over all casts (down + up)
         cast_direction   - 0 = down, 1 = up   (also stored as a string attr)
L2 : NetCDF of UPCASTS ONLY, each cast bin-averaged onto a 0.5 m depth grid.
     dims (depth, cast); time / direction / cast_number carried as coords.
L3 : regular depth x time grid -- TODO (see build_L3 stub).

The Ruskin software already detected profiles and split each into a DOWN and an
UP cast (region / regionCast / regionProfile tables). We reuse that detection
rather than re-detecting, so cast boundaries match the manufacturer tool.

Usage
-----
    python3 process_wirewalker.py --level 1
    python3 process_wirewalker.py --level 2
    python3 process_wirewalker.py --level 1 --max-casts 50      # quick test
    python3 process_wirewalker.py --level all

Author: NOPP-Aleutians processing
"""

from __future__ import annotations

import argparse
import sqlite3
import time as _time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import gsw
import netCDF4
import xarray as xr

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
HERE = Path(__file__).resolve().parent
DATA_DIR = HERE.parent  # /Users/drew/NOPP/data

RSK_PATH = DATA_DIR / "202605_NOPP-Aleutians_RBR-Concerto-Serial213752_DeploymentData.rsk"

L1_PATH = HERE / "L1" / "NOPP-Aleutians_RBR-Concerto-213752_L1_converted.nc"
L2_PATH = HERE / "L2" / "NOPP-Aleutians_RBR-Concerto-213752_L2_upcast_grid0.5m.nc"
L3_PATH = HERE / "L3" / "NOPP-Aleutians_RBR-Concerto-213752_L3_grid1m_30min.nc"
L3I_PATH = HERE / "L3" / "NOPP-Aleutians_RBR-Concerto-213752_L3_grid1m_30min_interp.nc"

# Deployment metadata
LAT = 49.5            # deg N
LON = -159.0          # deg E (i.e. 159 W)
ATM_DBAR = 10.1325    # atmospheric pressure offset (from rsk parameterKeys ATMOSPHERE)
INSTRUMENT = "RBRconcerto3 S/N 213752"
MOORING = "NOPP-Aleutians"

# L2 vertical grid (0.5 m bins, 0-500 m); max pressure in file is ~518 dbar
GRID_DZ = 0.5
GRID_ZMIN = 0.0
GRID_ZMAX = 500.0

# L3 regular depth x time grid: 1 m vertical, 30 min temporal.
# Native upcast cadence ~31 min, so a 30-min grid is ~1 upcast/bin (Nyquist 1 cph).
# Empty time bins (the few gaps) are left as NaN -- no temporal interpolation.
L3_DZ = 1.0
L3_DT = "30min"
# Companion gap-filled L3: linearly interpolate empty time bins (whole NaN columns)
# across gaps no longer than this many bins. Default 1 = fill only isolated single
# missing 30-min bins (the ~31-min-cadence phase slip); never bridges real gaps.
L3_INTERP_MAXGAP = 1
# Buoyancy frequency N^2 = (g/rho0) d(sigma0)/dz: sigma0 is vertically smoothed with a
# boxcar of this length before differentiating (chosen from the dsigma0/dz comparison).
N2_SMOOTH_M = 5.0
G_GRAV = 9.81

# Sampling rate (continuous, 2 Hz) and conductivity-cell thermal-mass correction.
# RBR pyRSKtools `correctTM` defaults (Lueck & Picklo 1990): alpha=0.04, beta=0.1 s^-1
# (relaxation time tau = 1/beta = 10 s), gamma = dC/dT = 1.0 (conductivity in mS/cm).
# C-T alignment is NOT applied: optimal C-T lag was found to be ~0 s for this logger.
FS = 2.0
TM_ALPHA = 0.04
TM_BETA = 0.1
TM_GAMMA = 1.0

# Mapping of stored `data` columns -> physical channels (from `channels` table)
#   channel01 = cond19  Conductivity   mS/cm
#   channel02 = temp14  Temperature    degC   (C-T cell thermistor -> used for salinity)
#   channel03 = pres24  Pressure       dbar   (total, incl. atmosphere)
#   channel09 = temp22  Temperature    degC   (extra thermistor)
#   channel10 = temp10  Temperature    degC   (extra thermistor)
DATA_COLS = ["tstamp", "channel01", "channel02", "channel03", "channel09", "channel10"]


# --------------------------------------------------------------------------- #
# Cast / profile bookkeeping
# --------------------------------------------------------------------------- #
@dataclass
class Cast:
    cast_number: int       # sequential 0..N-1 over all casts
    profile_number: int    # Ruskin profile index (regionProfileID)
    direction: int         # 0 down, 1 up
    t1: int                # start tstamp (ms)
    t2: int                # end tstamp (ms)


def load_casts(con: sqlite3.Connection) -> list[Cast]:
    """Read cast intervals from the Ruskin region tables, ordered in time."""
    q = """
        SELECT r.tstamp1, r.tstamp2, rc.type, rc.regionProfileID
        FROM regionCast rc
        JOIN region r ON r.regionID = rc.regionID
        ORDER BY r.tstamp1 ASC
    """
    rows = con.execute(q).fetchall()
    casts: list[Cast] = []
    for i, (t1, t2, ctype, prof) in enumerate(rows):
        direction = 1 if ctype.upper() == "UP" else 0
        casts.append(Cast(cast_number=i, profile_number=int(prof),
                           direction=direction, t1=int(t1), t2=int(t2)))
    return casts


def read_cast_data(con: sqlite3.Connection, cast: Cast) -> dict[str, np.ndarray]:
    """Pull the raw samples for one cast as float64 arrays."""
    q = (f"SELECT {', '.join(DATA_COLS)} FROM data "
         f"WHERE tstamp >= ? AND tstamp < ? ORDER BY tstamp ASC")
    rows = con.execute(q, (cast.t1, cast.t2)).fetchall()
    if not rows:
        return {}
    arr = np.asarray(rows, dtype=np.float64)
    return {
        "tstamp": arr[:, 0].astype(np.int64),
        "cond": arr[:, 1],     # mS/cm
        "temp": arr[:, 2],     # degC  (temp14)
        "pres": arr[:, 3],     # dbar  (total)
        "temp22": arr[:, 4],
        "temp10": arr[:, 5],
    }


# --------------------------------------------------------------------------- #
# TEOS-10 conversion
# --------------------------------------------------------------------------- #
def correct_thermal_mass(cond, temp, alpha=TM_ALPHA, beta=TM_BETA, gamma=TM_GAMMA, fs=FS):
    """Conductivity cell thermal-mass correction (Lueck & Picklo 1990), RBR form.

    Matches RBR pyRSKtools `RSK.correctTM`. `cond` in mS/cm, `temp` in degC, evenly
    sampled at `fs` Hz. Returns corrected conductivity (mS/cm). The recursive term is
    a first-order IIR filter, evaluated here with scipy.signal.lfilter per cast.
    """
    from scipy.signal import lfilter
    a = (4 * fs / 2) * (alpha / beta) / (1 + 4 * fs / 2 / beta)
    b = 1 - 2 * a / alpha
    dT = np.diff(temp, prepend=temp[0])
    corr = lfilter([1.0], [1.0, b], gamma * a * dT)
    return cond + corr


def convert(cond, temp, pres):
    """Raw conductivity/temperature/total-pressure -> derived CTD quantities.

    Returns a dict of float arrays. Uses gsw / TEOS-10 with the configured
    deployment lat/lon. `cond` mS/cm, `temp` ITS-90 degC, `pres` dbar (total).
    """
    sea_p = pres - ATM_DBAR                         # sea pressure, dbar
    depth = -gsw.z_from_p(sea_p, LAT)               # m, positive down
    SP = gsw.SP_from_C(cond, temp, sea_p)           # practical salinity
    SA = gsw.SA_from_SP(SP, sea_p, LON, LAT)        # absolute salinity g/kg
    CT = gsw.CT_from_t(SA, temp, sea_p)             # conservative temperature
    sigma0 = gsw.sigma0(SA, CT)                     # potential density anomaly kg/m3
    svel = gsw.sound_speed(SA, CT, sea_p)           # sound speed m/s
    return {
        "sea_pressure": sea_p,
        "depth": depth,
        "practical_salinity": SP,
        "absolute_salinity": SA,
        "conservative_temperature": CT,
        "sigma0": sigma0,
        "sound_speed": svel,
    }


# --------------------------------------------------------------------------- #
# Common NetCDF attributes
# --------------------------------------------------------------------------- #
def global_attrs(level: str) -> dict:
    return {
        "title": f"{MOORING} Wirewalker RBR Concerto CTD - {level}",
        "instrument": INSTRUMENT,
        "mooring": MOORING,
        "source_file": RSK_PATH.name,
        "processing_level": level,
        "geospatial_lat": LAT,
        "geospatial_lon": LON,
        "atmospheric_pressure_dbar": ATM_DBAR,
        "TEOS10_note": "Salinity uses C-T cell thermistor (temp14). gsw " + gsw.__version__,
        "date_created": _time.strftime("%Y-%m-%dT%H:%M:%S"),
        "Conventions": "CF-1.8",
    }


VAR_META = {
    "conductivity":             ("S m-1?mS/cm", "sea_water_electrical_conductivity", "mS/cm"),
    "temperature":              ("temp14", "sea_water_temperature", "degC"),
    "temperature_2":            ("temp22", "sea_water_temperature", "degC"),
    "temperature_3":            ("temp10", "sea_water_temperature", "degC"),
    "pressure":                 ("total pressure", "sea_water_pressure", "dbar"),
    "sea_pressure":             ("pressure - atmosphere", "sea_water_pressure_due_to_sea_water", "dbar"),
    "depth":                    ("depth, positive down", "depth", "m"),
    "practical_salinity":       ("PSS-78", "sea_water_practical_salinity", "1"),
    "absolute_salinity":        ("TEOS-10 SA", "sea_water_absolute_salinity", "g kg-1"),
    "conservative_temperature": ("TEOS-10 CT", "sea_water_conservative_temperature", "degC"),
    "sigma0":                   ("potential density anomaly ref 0 dbar", "sea_water_sigma_theta", "kg m-3"),
    "sound_speed":              ("TEOS-10 sound speed", "speed_of_sound_in_sea_water", "m s-1"),
}


# --------------------------------------------------------------------------- #
# L1 : full-resolution converted time series
# --------------------------------------------------------------------------- #
def build_L1(max_casts: int | None = None):
    """Stream casts and append to a NetCDF with an unlimited time dimension.

    Low memory: one cast in RAM at a time. float32 data vars + zlib compression.
    """
    L1_PATH.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{RSK_PATH}?mode=ro", uri=True)
    casts = load_casts(con)
    if max_casts:
        casts = casts[:max_casts]
    print(f"[L1] {len(casts)} casts to process -> {L1_PATH.name}")

    # L1 is a minimal converted product: measured channels + depth + cast flags.
    # All TEOS-10 derived quantities (salinity/CT/sigma0/sound speed) are produced
    # at L2 instead.
    f32vars = ["conductivity", "temperature", "pressure", "depth"]

    nc = netCDF4.Dataset(L1_PATH, "w", format="NETCDF4")
    for k, v in global_attrs("L1").items():
        nc.setncattr(k, v)
    nc.createDimension("time", None)  # unlimited

    vtime = nc.createVariable("time", "f8", ("time",), zlib=True, complevel=4)
    vtime.units = "milliseconds since 1970-01-01 00:00:00 UTC"
    vtime.standard_name = "time"
    vtime.calendar = "standard"

    ncv = {}
    for name in f32vars:
        comment, std, units = VAR_META[name]
        var = nc.createVariable(name, "f4", ("time",), zlib=True, complevel=4,
                                chunksizes=(65536,))
        var.units = units
        var.comment = comment
        if std:
            var.standard_name = std
        ncv[name] = var

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
        d = read_cast_data(con, c)
        if not d:
            continue
        n = d["tstamp"].size
        depth = -gsw.z_from_p(d["pres"] - ATM_DBAR, LAT)   # m, positive down
        sl = slice(pos, pos + n)
        vtime[sl] = d["tstamp"]
        ncv["conductivity"][sl] = d["cond"]
        ncv["temperature"][sl] = d["temp"]
        ncv["pressure"][sl] = d["pres"]
        ncv["depth"][sl] = depth
        v_castn[sl] = c.cast_number
        v_profn[sl] = c.profile_number
        v_dir[sl] = c.direction
        pos += n
        if c.cast_number % 1000 == 0:
            el = _time.time() - t0
            print(f"  cast {c.cast_number:>6d}/{len(casts)}  samples={pos:>10d}  {el:5.1f}s")

    nc.setncattr("n_samples", pos)
    nc.setncattr("n_casts", len(casts))
    nc.close()
    con.close()
    print(f"[L1] done: {pos} samples in {_time.time()-t0:.1f}s -> {L1_PATH}")


# --------------------------------------------------------------------------- #
# L2 : upcasts gridded to 0.5 m
# --------------------------------------------------------------------------- #
def build_L2(max_casts: int | None = None):
    """Bin-average every UPCAST from L1 onto the 0.5 m depth grid -> NetCDF (depth, cast).

    L2 is derived from the L1 product (NOT the raw .rsk), so it inherits only L1's
    channels: conductivity, temperature, pressure and depth, plus the cast flags.
    The extra thermistors and sea pressure are therefore absent by construction.
    The conductivity thermal-mass correction is applied per upcast before deriving
    salinity; the TEOS-10 quantities are computed here.
    """
    if not L1_PATH.exists():
        raise FileNotFoundError(
            f"L2 is built from L1, but {L1_PATH} is missing. Run --level 1 first.")
    L2_PATH.parent.mkdir(parents=True, exist_ok=True)

    l1 = xr.open_dataset(L1_PATH)
    cond_all = l1["conductivity"].values
    temp_all = l1["temperature"].values
    pres_all = l1["pressure"].values
    cn_all = l1["cast_number"].values
    cd_all = l1["cast_direction"].values
    pn_all = l1["profile_number"].values
    tms_all = l1["time"].values.astype("datetime64[ms]").astype(np.int64)  # ms since epoch
    l1.close()

    # Each cast is a contiguous block of equal cast_number (L1 is time-ordered).
    bnds = np.flatnonzero(np.diff(cn_all) != 0) + 1
    starts = np.concatenate(([0], bnds))
    ends = np.concatenate((bnds, [cn_all.size]))
    casts = [(s, e) for s, e in zip(starts, ends) if cd_all[s] == 1]  # upcasts only
    if max_casts:
        casts = casts[:max_casts]
    print(f"[L2] {len(casts)} upcasts from {L1_PATH.name} -> {L2_PATH.name}")

    edges = np.arange(GRID_ZMIN, GRID_ZMAX + GRID_DZ, GRID_DZ)
    zc = 0.5 * (edges[:-1] + edges[1:])     # bin centres
    nz = zc.size
    ncast = len(casts)

    gridvars = ["conductivity", "temperature", "practical_salinity",
                "absolute_salinity", "conservative_temperature", "sigma0",
                "sound_speed"]
    grids = {k: np.full((nz, ncast), np.nan, np.float32) for k in gridvars}
    nobs = np.zeros((nz, ncast), np.int32)
    ptime = np.full(ncast, np.nan)          # ms, mean sample time
    pcastn = np.zeros(ncast, np.int32)
    pprofn = np.zeros(ncast, np.int32)

    t0 = _time.time()
    for j, (s, e) in enumerate(casts):
        temp = temp_all[s:e]
        pres = pres_all[s:e]
        cond = correct_thermal_mass(cond_all[s:e], temp)   # RBR cell thermal-mass
        der = convert(cond, temp, pres)
        depth = der["depth"]
        col = {
            "conductivity": cond, "temperature": temp,
            "practical_salinity": der["practical_salinity"],
            "absolute_salinity": der["absolute_salinity"],
            "conservative_temperature": der["conservative_temperature"],
            "sigma0": der["sigma0"], "sound_speed": der["sound_speed"],
        }
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
        ptime[j] = float(np.mean(tms_all[s:e]))
        pcastn[j] = cn_all[s]
        pprofn[j] = pn_all[s]
        if j % 1000 == 0:
            print(f"  upcast {j:>6d}/{ncast}  {_time.time()-t0:5.1f}s")

    # Assemble xarray dataset
    data_vars = {}
    for k in gridvars:
        comment, std, units = VAR_META[k]
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
            "time": ("cast", ptime,
                     {"units": "milliseconds since 1970-01-01 00:00:00 UTC",
                      "standard_name": "time", "calendar": "standard"}),
            "cast_number": ("cast", pcastn,
                            {"comment": "sequential cast index in L1"}),
            "profile_number": ("cast", pprofn,
                               {"comment": "Ruskin profile index"}),
        },
        attrs={**global_attrs("L2"),
               "derived_from": L1_PATH.name,
               "cast_direction": "up (buoyant ascent)",
               "grid_dz_m": GRID_DZ, "grid_zmin_m": GRID_ZMIN, "grid_zmax_m": GRID_ZMAX,
               "thermal_mass_correction": "Lueck & Picklo (1990), RBR pyRSKtools correctTM",
               "thermal_mass_alpha": TM_ALPHA, "thermal_mass_beta_per_s": TM_BETA,
               "thermal_mass_gamma": TM_GAMMA,
               "ct_alignment_lag_s": 0.0,
               "conductivity_note": "thermal-mass-corrected; salinity derived from it"},
    )
    enc = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars}
    ds.to_netcdf(L2_PATH, encoding=enc)
    print(f"[L2] done: {ncast} upcasts x {nz} levels in {_time.time()-t0:.1f}s -> {L2_PATH}")


# --------------------------------------------------------------------------- #
# L3 : regular depth x time grid  (TODO)
# --------------------------------------------------------------------------- #
def _buoyancy_n2(sigma0, z, Lm=N2_SMOOTH_M):
    """Buoyancy frequency squared from potential density: N^2 = (g/rho0) d(sigma0)/dz.

    `sigma0` (nz, nt) on a uniform depth grid `z` (m, positive down). The density field
    is first smoothed vertically with a nan-aware boxcar of length `Lm` m, then
    differentiated in depth; rho0 = 1000 + smoothed sigma0 (local). Returns N^2 in s^-2
    (positive = stable). Because z is positive down, N^2 = (g/rho0) d(sigma0)/d(depth).
    """
    from scipy.ndimage import uniform_filter1d
    dz = float(np.median(np.diff(z)))
    win = max(1, int(round(Lm / dz)))
    if win > 1:
        fin = np.isfinite(sigma0).astype(float)
        c = np.where(np.isfinite(sigma0), sigma0, 0.0)
        num = uniform_filter1d(c, win, axis=0, mode="nearest")
        den = uniform_filter1d(fin, win, axis=0, mode="nearest")
        sm = np.where(den > 0, num / den, np.nan)
        sm[den < 0.5] = np.nan
    else:
        sm = sigma0
    dpdz = np.gradient(sm, dz, axis=0)              # d(sigma0)/d(depth), kg m^-4
    return G_GRAV * dpdz / (1000.0 + sm)            # N^2, s^-2


def _fill_empty_col_gaps(grids, ncast_bin, max_gap):
    """Linearly interpolate empty time bins across short gaps, per depth.

    Only whole-empty time bins (ncast_bin == 0) are candidates, and only runs of
    length <= `max_gap` that are bracketed by non-empty bins are filled (no
    extrapolation at the record ends, no bridging of longer/real gaps). A deep bin
    that is NaN in a bracketing column stays NaN. Returns a new dict of arrays.
    """
    out = {k: v.copy() for k, v in grids.items()}
    empty = ncast_bin == 0
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


def build_L3():
    """Regular depth x time grid from L2: 1 m vertical, 30 min temporal.

    Vertical: L2's 0.5 m bins are averaged in consecutive pairs to 1 m
    (centres 0.5, 1.5, ... 499.5). Temporal: upcasts are bin-averaged into fixed
    30 min bins by their cast time; bins with no upcast stay NaN (no interpolation).
    """
    import pandas as pd
    if not L2_PATH.exists():
        raise FileNotFoundError(
            f"L3 is built from L2, but {L2_PATH} is missing. Run --level 2 first.")
    L3_PATH.parent.mkdir(parents=True, exist_ok=True)

    l2 = xr.open_dataset(L2_PATH)
    gridvars = ["conductivity", "temperature", "practical_salinity",
                "absolute_salinity", "conservative_temperature", "sigma0",
                "sound_speed"]
    ncast = l2.sizes["cast"]
    nz2 = l2.sizes["depth"]
    assert nz2 % 2 == 0, "expected an even number of 0.5 m bins"
    nz3 = nz2 // 2
    z3 = l2["depth"].values.reshape(nz3, 2).mean(axis=1).astype(np.float32)  # 1 m centres

    # vertical 0.5 m -> 1 m (nan-aware average of each adjacent pair)
    import warnings
    V = {}
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "Mean of empty slice")
        for k in gridvars:
            a = l2[k].values.reshape(nz3, 2, ncast)
            V[k] = np.nanmean(a, axis=1)            # (nz3, ncast)

    # 30 min time bins spanning the deployment
    ctime = l2["time"].values.astype("datetime64[ns]")
    t0 = pd.Timestamp(ctime.min()).floor(L3_DT)
    t1 = pd.Timestamp(ctime.max()).ceil(L3_DT)
    edges = pd.date_range(t0, t1, freq=L3_DT)
    centers = (edges[:-1] + (edges[1] - edges[0]) / 2)
    ntime = len(edges) - 1
    ee = edges.values.astype("datetime64[ns]").astype("int64")
    ce = ctime.astype("int64")
    tbin = np.searchsorted(ee, ce, side="right") - 1   # bin index per cast
    tbin = np.clip(tbin, 0, ntime - 1)

    sums = {k: np.zeros((nz3, ntime)) for k in gridvars}
    cnts = {k: np.zeros((nz3, ntime)) for k in gridvars}
    ncast_bin = np.zeros(ntime, np.int32)
    t0w = _time.time()
    for j in range(ncast):
        b = tbin[j]
        ncast_bin[b] += 1
        for k in gridvars:
            col = V[k][:, j]
            fin = np.isfinite(col)
            sums[k][fin, b] += col[fin]
            cnts[k][fin, b] += 1
    grids = {k: np.where(cnts[k] > 0, sums[k] / np.maximum(cnts[k], 1), np.nan).astype(np.float32)
             for k in gridvars}
    l2.close()

    filled = np.mean(ncast_bin > 0) * 100
    print(f"[L3] {ntime} time bins ({L3_DT}) x {nz3} depths (1 m); "
          f"{filled:.1f}% of bins have >=1 upcast; built in {_time.time()-t0w:.1f}s")

    tcen_ms = centers.values.astype("datetime64[ms]").astype(np.int64)

    def _write(path, gr, extra_attrs):
        data_vars = {}
        for k in gridvars:
            comment, std, units = VAR_META[k]
            attrs = {"units": units, "comment": comment}
            if std:
                attrs["standard_name"] = std
            data_vars[k] = (("depth", "time"), gr[k], attrs)
        n2 = _buoyancy_n2(gr["sigma0"], z3, N2_SMOOTH_M)
        data_vars["buoyancy_frequency_squared"] = (
            ("depth", "time"), n2.astype(np.float32),
            {"units": "s-2", "long_name": "square of buoyancy (Brunt-Vaisala) frequency",
             "standard_name": "square_of_brunt_vaisala_frequency_in_sea_water",
             "comment": f"(g/rho0) d(sigma0)/dz; sigma0 vertically smoothed with a "
                        f"{N2_SMOOTH_M:g} m boxcar before differencing; z positive down"})
        data_vars["n_casts"] = (("time",), ncast_bin,
                                {"comment": "number of upcasts averaged into each time bin "
                                            "(0 marks a bin filled by interpolation, if any)"})
        ds = xr.Dataset(
            data_vars,
            coords={
                "depth": ("depth", z3, {"units": "m", "standard_name": "depth",
                                        "positive": "down", "comment": "1 m bin centres"}),
                "time": ("time", tcen_ms,
                         {"units": "milliseconds since 1970-01-01 00:00:00 UTC",
                          "standard_name": "time", "calendar": "standard",
                          "comment": f"{L3_DT} bin centres"}),
            },
            attrs={**global_attrs("L3"), "derived_from": L2_PATH.name,
                   "cast_direction": "up (buoyant ascent)",
                   "grid_dz_m": L3_DZ, "grid_dt": L3_DT, "nyquist_frequency_cph": 1.0,
                   "n2_vertical_smoothing_m": N2_SMOOTH_M,
                   **extra_attrs},
        )
        enc = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars}
        ds.to_netcdf(path, encoding=enc)

    # primary: no temporal interpolation
    _write(L3_PATH, grids,
           {"gap_handling": "empty 30-min bins left as NaN (no temporal interpolation)"})
    print(f"[L3] done -> {L3_PATH}")

    # companion: short-gap linear interpolation of empty bins
    gridsi = _fill_empty_col_gaps(grids, ncast_bin, L3_INTERP_MAXGAP)
    filledi = np.mean([np.isfinite(gridsi["practical_salinity"]).any(0)]) * 100
    _write(L3I_PATH, gridsi,
           {"gap_handling": f"empty bins linearly interpolated across gaps <= "
                            f"{L3_INTERP_MAXGAP} bin(s) ({L3_INTERP_MAXGAP*30} min); "
                            f"longer gaps left as NaN. n_casts==0 marks interpolated bins.",
            "interpolation_max_gap_bins": L3_INTERP_MAXGAP})
    print(f"[L3] done -> {L3I_PATH}  (time coverage now {filledi:.1f}%)")


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--level", choices=["1", "2", "3", "all"], required=True)
    ap.add_argument("--max-casts", type=int, default=None,
                    help="process only the first N casts (for testing)")
    args = ap.parse_args()

    if args.level in ("1", "all"):
        build_L1(args.max_casts)
    if args.level in ("2", "all"):
        build_L2(args.max_casts)
    if args.level in ("3", "all"):
        build_L3()


if __name__ == "__main__":
    main()
