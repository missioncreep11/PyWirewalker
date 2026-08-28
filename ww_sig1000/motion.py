"""Wirewalker platform-motion correction for a downward-looking Signature ADCP.

Two generations:

``beam_motion_correction`` (v1) is a port of ``WWcorr_beam.m`` (Bofu Zheng):
platform velocity from bandpass-integrated IMU acceleration (all three axes,
rotated with the per-ping AHRS attitude) plus the pressure-derived vertical rate.

``beam_motion_correction_v2`` is optimized for the Wirewalker's buoyant-ascent
upcast, following the NOPP_d2 diagnosis (21 raw upcasts, depth-banded wave-band
variance accounting):

- **Vertical from pressure alone.** On the upcast the vehicle is mechanically
  decoupled from the heaving surface buoy, so there is no wave-band vertical
  platform motion for the IMU to measure; the v1 IMU-vertical term only injected
  noise (wave-band velU rms 0.156 m/s with it vs 0.084 without, in the top 50 m).
- **Rotation with low-passed accelerometer tilt**, never the AHRS fusion. Real
  wave-band tilt at depth is < 0.1 deg (gyro), while the AHRS carries ~0.3 deg of
  in-band tilt noise - and fails outright on 15.7% of casts. Using the accel tilt
  makes the correction immune to AHRS faults by construction: on the 2024-05-02
  fault cast it cut wave-band velE rms from 0.13-0.14 to 0.036-0.040 m/s.
- **Horizontal IMU correction weighted by a depth gain.** Wire-transmitted lateral shaking
  is real near the surface (in-band accel 0.22 m/s2 at 0-50 m) and gone by 150 m
  (flat 0.01 m/s2 sensor floor below), so the bandpass-integrated horizontal
  term tapers to zero across `H_GAIN_FULL_M`..`H_GAIN_ZERO_M`.
- **Spike pings are interpolated over and flagged** (`ping_ok`) instead of the
  v1 zeroing, which turned each spike into a step in the integrated velocity.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import butter, filtfilt, detrend

from .transforms import xyz2enu, get_unit_vectors
from .geometry import BEAM_PHI_DEG, BEAM_AZI_DEG

G = 9.81  # m s-2

# v2 horizontal-correction depth gain (dbar ~ m), from the measured decay of
# wire-transmitted wave-band acceleration on NOPP_d2.
H_GAIN_FULL_M = 100.0    # full weight above this depth
H_GAIN_ZERO_M = 180.0    # zero weight below this depth
SPIKE_TOL = 0.3 * G      # |a| further than this from g = not a usable ping


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


# --------------------------------------------------------------------------- #
# v2: buoyant-ascent model
# --------------------------------------------------------------------------- #
@dataclass
class MotionV2:
    """v2 correction plus the attitude it used and per-ping validity."""
    corr_beam: np.ndarray       # (4, n) velocity to ADD to measured beam velocities
    platform_enu: np.ndarray    # (3, n) estimated platform velocity [E, N, Up]
    pitch_deg: np.ndarray       # LP-accel tilt actually used for every rotation -
    roll_deg: np.ndarray        # the caller must use the same for beam2enu/cell_depths
    heading_deg: np.ndarray     # heading used for every rotation (mag when usable)
    heading_source: str         # "mag" | "ahrs"
    ping_ok: np.ndarray         # (n,) False where |a| was a spike (interpolated over)
    h_gain: np.ndarray          # (n,) depth gain applied to the horizontal term
    usable: bool                # tilt guard + gravity check passed


def _depth_gain(pressure, full_m=H_GAIN_FULL_M, zero_m=H_GAIN_ZERO_M):
    """Cosine taper 1 -> 0 across the depth band where wire-transmitted wave-band
    motion dies out."""
    p = np.asarray(pressure, float)
    x = np.clip((p - full_m) / (zero_m - full_m), 0.0, 1.0)
    return 0.5 * (1.0 + np.cos(np.pi * x))


def beam_motion_correction_v2(time_s, pressure, accel_xyz, heading, fs, *, mag=None,
                              beam_angle=25.0, tilt_cutoff_hz=None, sail=True,
                              attitude_source="lp_accel", pitch_ahrs=None, roll_ahrs=None,
                              z_unit=(0.0, 0.0, 1.0),
                              h_full_m=H_GAIN_FULL_M, h_zero_m=H_GAIN_ZERO_M) -> MotionV2:
    """Per-ping platform-velocity correction for the buoyant-ascent upcast.

    Note the signature: the AHRS pitch/roll are *not* inputs. Tilt comes from the
    low-passed accelerometer (``attitude.gravity_direction``), so an AHRS attitude
    fault cannot enter the correction. When ``mag`` (the raw magnetometer, (3, n))
    is given and the field is sane, heading likewise comes from the
    tilt-compensated compass (Stage 2) instead of the passed-in AHRS heading — the
    AHRS has heading-only fault modes that the tilt detector cannot see, while the
    mag heading tracks the gyro at r ~ 0.98 straight through them. The heading
    actually used is returned in ``heading_deg``/``heading_source`` and the caller
    must use it for beam2enu.

    Returns a `MotionV2`; when ``usable`` is False (vehicle too tilted for the
    low-pass to be legal, or the accelerometer is not measuring gravity) the caller
    should fall back to the v1 path rather than trust these fields.
    """
    from .attitude import (DEFAULT_CUTOFF_HZ, LOWPASS_SMEAR_MAX_DEG,
                           gravity_direction, lowpass_smear_deg, mag_field_ok,
                           mag_heading, pitch_roll_from_up)
    from .platform import ACCEL_TOL

    cutoff = DEFAULT_CUTOFF_HZ if tilt_cutoff_hz is None else tilt_cutoff_hz
    acc = np.asarray(accel_xyz, float).copy()
    n = acc.shape[1]
    dt = 1.0 / fs

    # spikes: interpolate over them (v1 zeroed them, turning each into a velocity
    # step after integration) and flag the pings for the caller to exclude
    amag = np.linalg.norm(acc, axis=0)
    spike = np.abs(amag - G) > SPIKE_TOL
    ping_ok = ~spike
    if spike.any() and ping_ok.any():
        idx = np.arange(n)
        for i in range(3):
            acc[i, spike] = np.interp(idx[spike], idx[ping_ok], acc[i, ping_ok])

    heading = np.asarray(heading, float)
    if attitude_source == "ahrs":
        # trust the instrument's AHRS attitude as-is (handles an arbitrary fixed tilt
        # the LP-accel/mag path is not tuned for); heading is the AHRS heading.
        pitch = np.asarray(pitch_ahrs, float)
        roll = np.asarray(roll_ahrs, float)
        heading_source = "ahrs"
        usable = True
    else:
        # tilt from the low-passed gravity direction; guard as in attitude.reconstruct
        u = gravity_direction(acc, fs, cutoff)
        pitch, roll = pitch_roll_from_up(u)
        usable = bool(lowpass_smear_deg(acc, fs, cutoff, u_lp=u) < LOWPASS_SMEAR_MAX_DEG
                      and abs(np.median(amag[ping_ok]) - G) < ACCEL_TOL) if ping_ok.any() else False
        # heading from the tilt-compensated compass (Stage 2), if the field is sane
        heading_source = "ahrs"
        if mag is not None:
            h_mag, mt = mag_heading(mag, pitch, roll)
            if mag_field_ok(mag, mt):
                heading = h_mag
                heading_source = "mag"

    # horizontal: bandpass-integrated earth-frame acceleration, depth-gain weighted.
    # (no IMU vertical term at all - on the upcast dp/dt is the vertical motion)
    aENU = xyz2enu(acc[:, None, :], heading, pitch, roll)[:, 0, :]
    acu = _bandpass_reflect(aENU[0], fs)
    acv = _bandpass_reflect(aENU[1], fs)
    WWu = np.zeros(n)
    WWv = np.zeros(n)
    WWu[1:] = np.cumsum(acu[:-1] * dt)
    WWv[1:] = np.cumsum(acv[:-1] * dt)
    gain = _depth_gain(pressure, h_full_m, h_zero_m)

    dpw = _lowpass_reflect(np.gradient(np.asarray(pressure, float), dt), fs, fc=0.3)

    # "sail" term (MATLAB WWvel_upward sail_corr): the vehicle travels along the
    # wire, i.e. along its own axis, so on an angled wire the ascent has a real
    # horizontal component: v_platform = s * zhat with s set by dp/dt through the
    # axis's vertical component (s * b_z = -dp/dt). At the 8-13 deg lean NOPP_d2
    # carried after June this is 0.06-0.11 m/s along the lean azimuth - first
    # order against the currents. The MATLAB rotates with the AHRS attitude; here
    # zhat comes from the LP tilt + compass heading, so a fault cannot inject
    # 0.46*tan(43 deg) of fiction.
    if sail:
        zed = np.tile(np.asarray(z_unit, float).reshape(3, 1, 1), (1, 1, n))
        zhat = xyz2enu(zed, heading, pitch, roll)[:, 0, :]        # mount axis in ENU
        b_z = np.clip(zhat[2], 0.5, 1.0)                          # guard: tilt < 60 deg
        s = -dpw / b_z
        sail_e, sail_n = s * zhat[0], s * zhat[1]
    else:
        sail_e = sail_n = 0.0
    platform_enu = np.vstack([gain * WWu + sail_e, gain * WWv + sail_n, -dpw])

    vel_xyz = xyz2enu(platform_enu[:, None, :], heading, pitch, roll, reverse=True)[:, 0, :]
    phi = np.deg2rad(90.0 - beam_angle)     # geometry phi is elevation from horizontal
    azi = np.deg2rad(np.asarray(BEAM_AZI_DEG, float))
    corr_beam = np.empty((4, n))
    for b in range(4):
        bX, bY, bZ = (float(v) for v in get_unit_vectors(phi, azi[b], 0.0, 0.0))
        corr_beam[b] = vel_xyz[0] * bX + vel_xyz[1] * bY + vel_xyz[2] * bZ
    return MotionV2(corr_beam=corr_beam, platform_enu=platform_enu,
                    pitch_deg=pitch, roll_deg=roll,
                    heading_deg=heading, heading_source=heading_source,
                    ping_ok=ping_ok, h_gain=gain, usable=usable)
