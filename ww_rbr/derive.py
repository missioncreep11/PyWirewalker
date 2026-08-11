"""Derived-quantity physics for the CTD pipeline.

Pure functions (no config globals): conductivity thermal-mass correction,
TEOS-10 conversion, and buoyancy frequency from potential density.
"""
from __future__ import annotations

import gsw
import numpy as np


def correct_thermal_mass(cond, temp, alpha, beta, gamma, fs):
    """Conductivity cell thermal-mass correction (Lueck & Picklo 1990), RBR form.

    Matches RBR pyRSKtools `RSK.correctTM`. `cond` in mS/cm, `temp` in degC, evenly
    sampled at `fs` Hz. Returns corrected conductivity (mS/cm). The recursive term is
    a first-order IIR filter, evaluated with scipy.signal.lfilter per cast.
    """
    from scipy.signal import lfilter
    a = (4 * fs / 2) * (alpha / beta) / (1 + 4 * fs / 2 / beta)
    b = 1 - 2 * a / alpha
    dT = np.diff(temp, prepend=temp[0])
    corr = lfilter([1.0], [1.0, b], gamma * a * dT)
    return cond + corr


def convert(cond, temp, pres, cfg):
    """Raw conductivity/temperature/total-pressure -> derived CTD quantities (TEOS-10).

    Uses gsw with the deployment lat/lon and atmospheric pressure from `cfg`.
    `cond` mS/cm, `temp` ITS-90 degC, `pres` dbar (total). Returns a dict of arrays.
    """
    sea_p = pres - cfg.atm_dbar                     # sea pressure, dbar
    depth = -gsw.z_from_p(sea_p, cfg.lat)           # m, positive down
    SP = gsw.SP_from_C(cond, temp, sea_p)           # practical salinity
    SA = gsw.SA_from_SP(SP, sea_p, cfg.lon, cfg.lat)  # absolute salinity g/kg
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


def buoyancy_n2(sigma0, z, Lm, gravity):
    """Buoyancy frequency squared: N^2 = (g/rho0) d(sigma0)/dz.

    `sigma0` (nz, nt) on a uniform depth grid `z` (m, positive down). The density
    field is first smoothed vertically with a nan-aware boxcar of length `Lm` m,
    then differentiated in depth; rho0 = 1000 + smoothed sigma0. Returns N^2 (s^-2,
    positive = stable). Because z is positive down, N^2 = (g/rho0) d(sigma0)/d(depth).
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
    return gravity * dpdz / (1000.0 + sm)           # N^2, s^-2
