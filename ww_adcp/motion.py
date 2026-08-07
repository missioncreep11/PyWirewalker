"""Wirewalker platform-motion correction for a downward-looking Signature ADCP.

Port of ``WWcorr_beam.m`` (Bofu Zheng). Estimates the platform's translational
velocity from the IMU (bandpass-integrated dynamic acceleration) plus the
pressure-derived vertical rate, and returns the per-ping, per-beam correction to
add to the measured beam velocities so they become water-relative-to-earth.
"""
from __future__ import annotations

import numpy as np
from scipy.signal import butter, filtfilt, detrend

from .transforms import xyz2enu, get_unit_vectors
from .geometry import BEAM_PHI_DEG, BEAM_AZI_DEG

G = 9.81  # m s-2


def _bandpass_reflect(x, fs, lo=0.1, hi=1.2, order=1):
    """Butterworth band-pass with flip-pad-flip edge handling, detrend + demean,
    matching WWcorr_beam's acceleration filtering."""
    b, a = butter(order, [lo / (fs / 2), hi / (fs / 2)], "bandpass")
    xf = np.concatenate([x[::-1], x, x[::-1]])
    y = filtfilt(b, a, xf)
    y = detrend(y, type="constant")
    n = x.size
    y = y[n:2 * n]
    return y - np.nanmean(y)


def _lowpass_reflect(x, fs, fc=0.3, order=1):
    b, a = butter(order, fc / (fs / 2), "low")
    xf = np.concatenate([x[::-1], x, x[::-1]])
    y = filtfilt(b, a, xf)
    n = x.size
    return y[n:2 * n]


def beam_motion_correction(time_s, pressure, accel_xyz, pitch, roll, heading, fs,
                           beam_angle=25.0):
    """Per-ping platform-velocity correction projected onto the four beams.

    Parameters
    ----------
    time_s : (nping,) seconds.
    pressure : (nping,) dbar (~m).
    accel_xyz : (3, nping) instrument-frame acceleration, m s-2 (dolfyn `accel`).
    pitch, roll, heading : (nping,) degrees.
    fs : sample rate (Hz).

    Returns
    -------
    corr_beam : (4, nping) velocity to ADD to each measured beam velocity.
    platform_enu : (3, nping) estimated platform velocity [E, N, Up] (diagnostic).
    """
    accel_xyz = np.asarray(accel_xyz, float)
    # acceleration -> ENU, remove static gravity on the up-axis
    aENU = xyz2enu(accel_xyz[:, None, :], heading, pitch, roll)[:, 0, :]  # (3, nping)
    aENU[2] = aENU[2] - G
    # spike mask: zero any ping whose |a| exceeds 0.3 g on any axis
    bad = (np.abs(aENU) > 0.3 * G).any(axis=0)
    aENU[:, bad] = 0.0

    acu = _bandpass_reflect(aENU[0], fs)
    acv = _bandpass_reflect(aENU[1], fs)
    acw = _bandpass_reflect(aENU[2], fs)

    # pressure-derived vertical rate (dp/dt), low-passed
    dt = 1.0 / fs
    dpw = _lowpass_reflect(np.gradient(np.asarray(pressure, float), dt), fs, fc=0.3)

    # integrate bandpassed acceleration -> platform translational velocity (m/s)
    def _integrate(ac):
        w = np.zeros_like(ac)
        w[1:] = np.cumsum(ac[:-1] * dt)
        return w
    WWu, WWv, WWw = _integrate(acu), _integrate(acv), _integrate(acw)

    # platform velocity in ENU: vertical from pressure + IMU high-freq residual
    platform_enu = np.vstack([WWu, WWv, -dpw + WWw])          # (3, nping)

    # rotate platform velocity ENU -> instrument XYZ
    vel_xyz = xyz2enu(platform_enu[:, None, :], heading, pitch, roll, reverse=True)[:, 0, :]

    # project onto the (level) beam unit vectors
    phi = np.deg2rad(BEAM_PHI_DEG)
    azi = np.deg2rad(np.asarray(BEAM_AZI_DEG, float))
    corr_beam = np.empty((4, time_s.size))
    for b in range(4):
        bX, bY, bZ = (float(v) for v in get_unit_vectors(phi, azi[b], 0.0, 0.0))
        corr_beam[b] = vel_xyz[0] * bX + vel_xyz[1] * bY + vel_xyz[2] * bZ
    return corr_beam, platform_enu
