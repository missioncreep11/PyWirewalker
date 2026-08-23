"""Reconstruct instrument tilt from the raw accelerometer.

For casts where the AHRS orientation solution has failed (see `platform.cast_qc`),
the raw sensors are still good — on NOPP_d2 the accelerometer measured clean gravity
in 100% of the 941 faulted casts. This module turns that measurement back into the
`pitch`/`roll` the velocity pipeline consumes.

Why low-pass, and when it is legal
----------------------------------
The accelerometer measures gravity *plus* the vehicle's own acceleration. Wave-band
motion (0.05-0.5 Hz) is the dominant contaminant, so low-passing well below it isolates
the gravity direction and sharply reduces tilt noise.

That is only valid while the filter is not smearing *real* attitude motion. On a
leaning vehicle, spin precesses gravity around the instrument frame at the spin rate —
but observed Wirewalker spin (4-30 mHz) sits at or below the cutoff, so the precession
passes through the filter and even a 13-deg steady lean low-passes cleanly. What must
be caught is motion of the gravity direction *above* the cutoff. `reconstruct`
therefore measures the error directly: `lowpass_smear_deg`, the median angle between
the raw and low-passed gravity directions, gates `lowpass_valid` (with
`highfreq_fraction` kept as a diagnostic), so the assumption is checked per cast
rather than assumed.

Frame
-----
`pitch`/`roll` follow the convention of `transforms._tilt_heading_matrix`, whose row 2
(earth-up in instrument coordinates) is ``[sin p, sin r cos p, cos p cos r]``. Inverting
gives ``pitch = asin(u_x)``, ``roll = atan2(u_y, u_z)``; the round trip is exact to
float precision, and is pinned by a test.

The pipeline's own `pitch`/`roll` sit ~2 deg from the raw accelerometer (closer than
`orientmat` does, at ~4 deg). Reprocessing healthy casts with accel-derived tilt shifts
mean velU by 0.00025 m/s — below what 3 casts can resolve — so no frame correction is
applied. That is a measured decision, not an assumption; revisit it with a larger
sample if the frame question is ever settled properly.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .platform import ACCEL_TOL, G, _filt

# Validity guard: the low-pass is legal when the error it introduces is small,
# measured directly as the median angle between the raw and low-passed gravity
# directions. On NOPP_d2 every real upcast sits at 0.2-1.8 deg (including casts
# leaning a genuine 13 deg with spin at ~20 mHz, which the old tilt-based guard
# wrongly rejected - 29% of the record after the June current increase leaned
# >10 deg); a genuinely fast-coning vehicle would show tens of degrees. The spin
# rates observed (4-30 mHz) sit at or below the cutoff, so a steady lean's
# precession *passes* the filter rather than being smeared by it.
LOWPASS_SMEAR_MAX_DEG = 5.0
DEFAULT_CUTOFF_HZ = 0.03          # ~33 s; well below the 0.05-0.5 Hz wave band


def pitch_roll_from_up(u) -> tuple[np.ndarray, np.ndarray]:
    """Pipeline-convention pitch/roll (deg) from an instrument-frame up-vector (3, n)."""
    u = np.asarray(u, float)
    return (np.rad2deg(np.arcsin(np.clip(u[0], -1.0, 1.0))),
            np.rad2deg(np.arctan2(u[1], u[2])))


def up_from_pitch_roll(pitch_deg, roll_deg) -> np.ndarray:
    """Instrument-frame up-vector (3, n) from pitch/roll — inverse of the above."""
    p, r = np.deg2rad(np.asarray(pitch_deg, float)), np.deg2rad(np.asarray(roll_deg, float))
    return np.stack([np.sin(p), np.sin(r) * np.cos(p), np.cos(p) * np.cos(r)])


def gravity_direction(accel, fs, cutoff_hz=DEFAULT_CUTOFF_HZ) -> np.ndarray:
    """Low-passed unit gravity direction in the instrument frame, shape (3, n)."""
    a = np.asarray(accel, float)
    lp = np.stack([_filt(a[i], fs, cutoff_hz, "low") for i in range(3)])
    return lp / np.maximum(np.linalg.norm(lp, axis=0), 1e-9)


def highfreq_fraction(accel, fs, cutoff_hz=DEFAULT_CUTOFF_HZ) -> float:
    """Fraction of gravity-direction variance sitting above `cutoff_hz`.

    Large values mean the low-pass is discarding real attitude change, not just noise —
    the signal that `reconstruct`'s validity guard is built on.
    """
    a = np.asarray(accel, float)
    raw = a / np.maximum(np.linalg.norm(a, axis=0), 1e-9)
    lp = gravity_direction(a, fs, cutoff_hz)
    v_hi = np.var(raw - lp, axis=1).sum()
    v_all = np.var(raw, axis=1).sum()
    return float(v_hi / v_all) if v_all > 0 else np.nan


@dataclass
class TiltReconstruction:
    """Accel-derived tilt for one cast, with the checks that say whether to trust it."""
    pitch_deg: np.ndarray
    roll_deg: np.ndarray
    tilt_deg: np.ndarray
    up: np.ndarray                  # (3, n) instrument-frame gravity direction
    accel_mag: np.ndarray
    cutoff_hz: float
    highfreq_fraction: float
    lowpass_smear_deg: float        # median angle(raw, low-passed) gravity direction
    lowpass_valid: bool             # the low-pass introduces little error, so it is safe
    accel_is_gravity: bool          # |a| ~ g, so the accelerometer is measuring gravity
    n_bad_mag: int

    @property
    def usable(self) -> bool:
        return self.lowpass_valid and self.accel_is_gravity

    def summary(self) -> dict:
        return {"tilt_median_deg": float(np.median(self.tilt_deg)),
                "pitch_median_deg": float(np.median(self.pitch_deg)),
                "roll_median_deg": float(np.median(self.roll_deg)),
                "accel_mag_median": float(np.median(self.accel_mag)),
                "highfreq_fraction": self.highfreq_fraction,
                "lowpass_smear_deg": self.lowpass_smear_deg,
                "lowpass_valid": self.lowpass_valid,
                "accel_is_gravity": self.accel_is_gravity,
                "usable": self.usable,
                "cutoff_hz": self.cutoff_hz}


# --------------------------------------------------------------------------- #
# Stage 2: heading from the tilt-compensated magnetometer
# --------------------------------------------------------------------------- #
def mag_heading(mag, pitch_deg, roll_deg, hard_iron=(0.0, 0.0, 0.0)):
    """Tilt-compensated magnetic heading (deg), in the pipeline's convention.

    Inverts ``transforms._tilt_heading_matrix``: the earth's field has zero *east*
    component in magnetic coordinates, so after tilt-correcting the measured field
    with P(pitch, roll), ``heading = atan2(-m_x, m_y) + 90``. With ``DECL = 0`` in
    the deployment configuration the AHRS heading is also magnetic, so the two are
    directly comparable (measured agreement on healthy NOPP_d2 casts: <= 1.6 deg
    with no calibration).

    Why trust it: on the 2024-05-02 AHRS fault the mag-derived heading rate tracks
    the gyro at r = 0.985 while the AHRS claims -47.7 deg/s; and on two nominally
    healthy casts (2024-04-29, 2024-05-06) the AHRS heading drifted at -4 deg/s
    against both gyro and magnetometer - a *heading-only* AHRS fault mode that the
    tilt-based `ahrs_error` detector cannot see.

    Returns ``(heading_deg, m_tilt)`` where ``m_tilt`` is the (3, n) tilt-corrected
    field - its horizontal magnitude is the per-ping quality check.
    """
    from .transforms import _tilt_heading_matrix
    m = np.asarray(mag, float) - np.asarray(hard_iron, float)[:, None]
    n = m.shape[1]
    P = _tilt_heading_matrix(np.full(n, 90.0), pitch_deg, roll_deg)   # H(0) = I
    mt = np.einsum("nij,jn->in", P, m)
    head = (np.rad2deg(np.arctan2(-mt[0], mt[1])) + 90.0) % 360.0
    return head, mt


def fit_hard_iron(mt_x, mt_y) -> tuple[float, float, float]:
    """Kasa circle fit of the tilt-compensated horizontal field over spinning data.

    Returns (cx, cy, radius): the hard-iron offset in the tilt-corrected horizontal
    plane and the field radius. On NOPP_d2 the offset is ~2.7% of the radius
    (<= 1.6 deg of heading), so calibration is optional there - but the Wirewalker's
    spin makes the fit geometry ideal wherever it is needed.
    """
    x = np.asarray(mt_x, float)
    y = np.asarray(mt_y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    A = np.c_[2 * x, 2 * y, np.ones_like(x)]
    c = np.linalg.lstsq(A, x ** 2 + y ** 2, rcond=None)[0]
    return float(c[0]), float(c[1]), float(np.sqrt(c[2] + c[0] ** 2 + c[1] ** 2))


def mag_field_ok(mag, mt, *, min_horiz_frac=0.2, max_cv=0.15) -> bool:
    """Is the magnetometer field sane enough to steer by?

    The tilt-corrected horizontal magnitude must be a stable, substantial fraction
    of the total field (dip ~58 deg in California leaves ~53% horizontal; a mooring
    with severe iron or a failed sensor does not).
    """
    h = np.hypot(mt[0], mt[1])
    hm = float(np.median(h))
    total = float(np.median(np.linalg.norm(np.asarray(mag, float), axis=0)))
    return bool(hm > min_horiz_frac * total and np.std(h) / max(hm, 1e-9) < max_cv)


def lowpass_smear_deg(accel, fs, cutoff_hz=DEFAULT_CUTOFF_HZ, u_lp=None) -> float:
    """Median angle (deg) between the raw and low-passed gravity directions.

    The error the low-pass actually introduces, real motion and noise combined —
    the quantity that decides whether filtering is legal. Real upcasts sit at
    0.2-1.8 deg; fast large-angle coning would sit at tens of degrees.
    """
    a = np.asarray(accel, float)
    u_raw = a / np.maximum(np.linalg.norm(a, axis=0), 1e-9)
    if u_lp is None:
        u_lp = gravity_direction(a, fs, cutoff_hz)
    cos = np.clip((u_raw * u_lp).sum(axis=0), -1.0, 1.0)
    return float(np.median(np.rad2deg(np.arccos(cos))))


def reconstruct(ds, sl=None, *, cutoff_hz=DEFAULT_CUTOFF_HZ) -> TiltReconstruction:
    """Reconstruct tilt for one cast from `ds['accel']`."""
    sl = slice(None) if sl is None else sl
    fs = float(ds.attrs["fs"])
    acc = ds["accel"].values[:, sl]
    mag = np.linalg.norm(acc, axis=0)
    u = gravity_direction(acc, fs, cutoff_hz)
    pitch, roll = pitch_roll_from_up(u)
    tilt = np.rad2deg(np.arccos(np.clip(u[2], -1.0, 1.0)))
    smear = lowpass_smear_deg(acc, fs, cutoff_hz, u_lp=u)
    return TiltReconstruction(
        pitch_deg=pitch, roll_deg=roll, tilt_deg=tilt, up=u, accel_mag=mag,
        cutoff_hz=cutoff_hz,
        highfreq_fraction=highfreq_fraction(acc, fs, cutoff_hz),
        lowpass_smear_deg=smear,
        lowpass_valid=bool(smear < LOWPASS_SMEAR_MAX_DEG),
        accel_is_gravity=bool(abs(np.median(mag) - G) < ACCEL_TOL),
        n_bad_mag=int((np.abs(mag - G) > ACCEL_TOL).sum()))


def apply_to(ds, sl=None, *, cutoff_hz=DEFAULT_CUTOFF_HZ):
    """A copy of `ds` (optionally one cast) with pitch/roll replaced by the
    reconstruction. Heading is left untouched — that is Stage 2."""
    sub = ds if sl is None else ds.isel(time=sl)
    rec = reconstruct(sub, cutoff_hz=cutoff_hz)
    out = sub.copy(deep=True)
    out["pitch"] = ("time", rec.pitch_deg)
    out["roll"] = ("time", rec.roll_deg)
    out.attrs["attitude_source"] = f"accel-reconstructed (lowpass {cutoff_hz} Hz)"
    return out, rec


# --------------------------------------------------------------------------- #
# validation harness
# --------------------------------------------------------------------------- #
def _ping_noise(x):
    """Per-ping noise from the first difference: sd(diff)/sqrt(2)."""
    return float(np.std(np.diff(np.asarray(x, float))) / np.sqrt(2))


def noise_vs_cutoff(ds, sl=None, cutoffs=(0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.3, 1.0)):
    """Tilt noise as a function of low-pass cutoff, for choosing one on evidence.

    Returns a list of dicts with the per-ping pitch/roll noise and the fraction of
    gravity-direction variance removed. The useful cutoff is where noise stops
    improving while `highfreq_fraction` is still small.
    """
    sl = slice(None) if sl is None else sl
    fs = float(ds.attrs["fs"])
    acc = ds["accel"].values[:, sl]
    rows = []
    for fc in cutoffs:
        u = gravity_direction(acc, fs, fc)
        p, r = pitch_roll_from_up(u)
        rows.append({"cutoff_hz": fc,
                     "pitch_noise_deg": _ping_noise(p),
                     "roll_noise_deg": _ping_noise(r),
                     "tilt_noise_deg": float(np.hypot(_ping_noise(p), _ping_noise(r))),
                     "highfreq_fraction": highfreq_fraction(acc, fs, fc)})
    # unfiltered reference
    u0 = acc / np.maximum(np.linalg.norm(acc, axis=0), 1e-9)
    p0, r0 = pitch_roll_from_up(u0)
    rows.append({"cutoff_hz": np.inf, "pitch_noise_deg": _ping_noise(p0),
                 "roll_noise_deg": _ping_noise(r0),
                 "tilt_noise_deg": float(np.hypot(_ping_noise(p0), _ping_noise(r0))),
                 "highfreq_fraction": 0.0})
    return rows


def compare_with_pipeline(ds, sl=None, *, cutoff_hz=DEFAULT_CUTOFF_HZ) -> dict:
    """Reconstruction vs the pipeline's own pitch/roll, for a cast.

    On a *healthy* cast the pipeline attitude is trustworthy, so this measures the
    reconstruction's accuracy. On a faulted one the same numbers measure the size of
    the fault instead — which is why callers must screen with `platform.cast_qc` first
    and know which case they are in.
    """
    sl = slice(None) if sl is None else sl
    rec = reconstruct(ds, sl, cutoff_hz=cutoff_hz)
    pitch, roll = ds["pitch"].values[sl], ds["roll"].values[sl]
    u_pipe = up_from_pitch_roll(pitch, roll)
    ang = np.rad2deg(np.arccos(np.clip((u_pipe * rec.up).sum(axis=0), -1.0, 1.0)))
    return {
        "pitch_bias_deg": float(np.median(rec.pitch_deg - pitch)),
        "roll_bias_deg": float(np.median(rec.roll_deg - roll)),
        "pitch_rms_deg": float(np.sqrt(np.mean((rec.pitch_deg - pitch) ** 2))),
        "roll_rms_deg": float(np.sqrt(np.mean((rec.roll_deg - roll) ** 2))),
        "up_angle_median_deg": float(np.median(ang)),
        "noise_pipeline_deg": float(np.hypot(_ping_noise(pitch), _ping_noise(roll))),
        "noise_reconstructed_deg": float(np.hypot(_ping_noise(rec.pitch_deg),
                                                  _ping_noise(rec.roll_deg))),
        "usable": rec.usable,
        "highfreq_fraction": rec.highfreq_fraction,
    }
