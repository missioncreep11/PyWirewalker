"""Platform kinematics and AHRS validation for a Wirewalker Signature1000.

Reconstructs what the vehicle was doing during a cast — tilt, rotation, climb — and,
critically, **checks the instrument's own orientation solution against its raw
sensors**.

Why the AHRS is the thing under test
------------------------------------
On the NOPP_d2 record a cast at 2024-05-02 09:51 produced velocity ~6x rougher than
its neighbours. The AHRS reported a steady 43 deg tilt rotating at 48 deg/s. That was
false. Three independent lines say the vehicle was upright and normal:

- the **accelerometer** measured |a| = 9.76 m/s2 with a tight spread (gravity, with no
  kinematic acceleration to corrupt it) lying **2.2 deg** off the instrument z-axis;
- the **gyro** (``angrt``) reported ~3 deg/s of rotation, not 48;
- the **pressure** record showed a full 500 dbar span at a normal 0.46 m/s, whereas a
  genuinely leaning wire would shorten the vertical travel by cos(tilt).

So the raw sensors are the reference and ``orientmat`` is validated against them. The
decisive statistic is `ahrs_error`: the per-ping angle between where the AHRS puts
earth-up and where the accelerometer measures it. Healthy casts sit at 4-7 deg, the
floor set by accelerometer noise; the bad cast sits at 43.8 deg.

This matters because a wrong attitude corrupts velocity twice over: ``beam2enu``
rotates beam velocities with AHRS heading/pitch/roll, and ``geometry.cell_depths``
maps cell depths with the same pitch/roll.

Sign conventions
----------------
Earth frame is (East, North, Up). ``lean_azimuth`` is a compass bearing — degrees
clockwise from north — so it advances *opposite* to ``omega_up``, which is positive
counterclockwise by the right-hand rule: ``d(lean_azimuth)/dt == -omega_up``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.signal import butter, detrend, filtfilt, welch

G = 9.81                       # m s-2

# thresholds (degrees, deg/s)
TILT_QUIET_DEG = 8.0           # below this the platform is effectively upright
SPIN_QUIET_DPS = 6.0           # below this there is no sustained rotation
CIRCULARITY_CONING = 0.5       # (pitch, roll) covariance eigenvalue ratio
# Calibrated on NOPP_d2: healthy casts sit at a 4-7 deg median error (p90 <= 8.6),
# which is the floor set by accelerometer noise and real vehicle acceleration
# (~0.7 m/s2 rms / g ~ 4 deg). The 2024-05-02 09:51 fault sits at 43.8 deg. 15 deg is
# comfortably clear of the healthy population without approaching the fault.
AHRS_BAD_DEG = 15.0            # attitude error above this = solution not trustworthy
ACCEL_TOL = 0.6                # |a| may differ from g by this much and still be gravity


# --------------------------------------------------------------------------- #
# attitude: the AHRS solution, and the raw-sensor reference
# --------------------------------------------------------------------------- #
def attitude(orientmat) -> tuple[np.ndarray, np.ndarray]:
    """Tilt magnitude and lean azimuth from the AHRS orientation matrix.

    `orientmat` is (3, 3, nping), dims (earth, inst, time), so column 2 is the
    instrument z-axis in earth (E, N, U) coordinates.
    """
    om = np.asarray(orientmat, float)
    z_earth = om[:, 2, :]
    tilt = np.rad2deg(np.arccos(np.clip(z_earth[2], -1.0, 1.0)))
    azi = np.rad2deg(np.arctan2(z_earth[0], z_earth[1])) % 360.0
    return tilt, azi


def attitude_from_accel(accel) -> tuple[np.ndarray, np.ndarray]:
    """Tilt magnitude and |a|, from the accelerometer alone.

    With the vehicle in steady profiling the accelerometer measures gravity, so the
    angle between the measured vector and the instrument z-axis *is* the tilt. This
    needs no fusion, no magnetometer and no gyro integration, which is what makes it a
    valid reference for the AHRS. `accel_mag` near g is the check that no kinematic
    acceleration is contaminating it.
    """
    a = np.asarray(accel, float)
    mag = np.linalg.norm(a, axis=0)
    tilt = np.rad2deg(np.arccos(np.clip(a[2] / np.maximum(mag, 1e-9), -1.0, 1.0)))
    return tilt, mag


def ahrs_error(accel, orientmat, per_ping=False) -> float | np.ndarray:
    """Angle (deg) between the AHRS's idea of "up" and measured gravity.

    Both directions are compared *in the instrument frame*: the AHRS places earth-up
    along ``R^T z``, i.e. the third row of the orientation matrix, while the
    accelerometer measures it along ``a/|a|``. The angle between them is the attitude
    error. Returns the per-ping median by default; `per_ping=True` gives the series.

    Comparing them per ping matters. The obvious alternative — rotate acceleration
    into earth coordinates and measure how far the *mean* lands from vertical — fails
    exactly when the false attitude is rotating, because the horizontal residuals then
    average to nearly zero. On the NOPP_d2 fault that method reported 3.2 deg for a
    43 deg error; this one reports the 43 deg.
    """
    a = np.asarray(accel, float)
    a = a / np.maximum(np.linalg.norm(a, axis=0), 1e-9)      # (3, n) unit
    up_inst = np.asarray(orientmat, float)[2, :, :]          # earth-up in inst coords
    cos = np.clip((a * up_inst).sum(axis=0), -1.0, 1.0)
    err = np.rad2deg(np.arccos(cos))
    return err if per_ping else float(np.median(err))


def angular_velocity(orientmat, fs) -> np.ndarray:
    """Body-frame angular velocity (nping, 3) in rad/s, from Omega = R^T dR/dt.

    This is the rotation rate *the AHRS solution implies*. Compare it against the gyro
    (``angrt``) rather than trusting it: on a failed solution the two diverge, and the
    gyro is the honest one.
    """
    R = np.moveaxis(np.asarray(orientmat, float), -1, 0)
    dR = np.gradient(R, 1.0 / fs, axis=0)
    S = np.einsum("nji,njk->nik", R, dR)
    S = 0.5 * (S - np.transpose(S, (0, 2, 1)))
    return np.stack([S[:, 2, 1], S[:, 0, 2], S[:, 1, 0]], axis=1)


def to_earth(vec_body, orientmat) -> np.ndarray:
    """Rotate a (nping, 3) body-frame vector into earth (E, N, U)."""
    R = np.moveaxis(np.asarray(orientmat, float), -1, 0)
    return np.einsum("nij,nj->ni", R, np.asarray(vec_body, float))


def circularity(pitch, roll) -> float:
    """Eigenvalue ratio of the (pitch, roll) covariance, in [0, 1].

    ~1 = the attitude traces a circle (a steady lean on a rotating platform);
    ~0 = a straight line (planar rocking). Medians of pitch and roll cannot tell
    these apart — both sit near zero.
    """
    p, r = np.asarray(pitch, float), np.asarray(roll, float)
    m = np.isfinite(p) & np.isfinite(r)
    if m.sum() < 3:
        return np.nan
    ev = np.linalg.eigvalsh(np.cov(p[m], r[m]))
    return float(ev[0] / ev[1]) if ev[1] > 0 else np.nan


# --------------------------------------------------------------------------- #
# vertical motion
# --------------------------------------------------------------------------- #
def _filt(x, fs, cutoff, btype, order=2):
    nyq = fs / 2.0
    wn = np.clip(np.atleast_1d(cutoff) / nyq, 1e-6, 0.999)
    b, a = butter(order, wn if wn.size > 1 else wn[0], btype)
    pad = np.concatenate([x[::-1], x, x[::-1]])
    return filtfilt(b, a, pad)[x.size:2 * x.size]


def vertical(pressure, fs, wave_band=(0.05, 0.5), climb_cutoff=0.02):
    """Split pressure into mean climb and wave-band heave.

    Returns (depth_m, ascent_rate_ms, heave_m); ascent is positive upward and
    low-passed below `climb_cutoff` Hz so it is the profiling rate, not orbital motion.
    """
    p = np.asarray(pressure, float)
    depth = p.copy()
    ascent = -_filt(np.gradient(depth, 1.0 / fs), fs, climb_cutoff, "low")
    heave = _filt(detrend(depth, type="linear"), fs, np.array(wave_band), "bandpass")
    return depth, ascent, heave


def spectrum(x, fs, nperseg=None):
    """Welch PSD of a kinematic series. Returns (freq_hz, psd)."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < 16:
        return np.array([]), np.array([])
    nperseg = nperseg or int(min(x.size, max(256, 2 ** int(np.log2(x.size / 8)))))
    return welch(detrend(x, type="constant"), fs=fs, nperseg=nperseg)


# --------------------------------------------------------------------------- #
# per-cast QC — the detector
# --------------------------------------------------------------------------- #
def cast_qc(ds, sl=None) -> dict:
    """One row of AHRS health for a cast. Cheap enough to run over a whole record.

    Keys: `ahrs_error_deg` (the headline number), `tilt_accel_deg` (truth),
    `tilt_ahrs_deg` (reported), `accel_mag`, `gyro_dps` vs `ahrs_dps` (rotation rate
    from the gyro against the rate the AHRS solution implies), and `ahrs_ok`.
    """
    sl = slice(None) if sl is None else sl
    fs = float(ds.attrs["fs"])
    om = ds["orientmat"].values[:, :, sl]
    acc = ds["accel"].values[:, sl]

    tilt_ahrs, _ = attitude(om)
    tilt_acc, mag = attitude_from_accel(acc)
    err = ahrs_error(acc, om)
    ahrs_dps = np.rad2deg(np.median(np.abs(to_earth(angular_velocity(om, fs), om)[:, 2])))
    gyro_dps = (np.rad2deg(np.median(np.abs(ds["angrt"].values[:, sl])))
                if "angrt" in ds else np.nan)
    accel_mag = float(np.median(mag))
    return {
        "ahrs_error_deg": err,
        "tilt_accel_deg": float(np.median(tilt_acc)),
        "tilt_ahrs_deg": float(np.median(tilt_ahrs)),
        "tilt_disagree_deg": float(np.median(tilt_ahrs) - np.median(tilt_acc)),
        "accel_mag": accel_mag,
        "accel_is_gravity": bool(abs(accel_mag - G) < ACCEL_TOL),
        "ahrs_dps": float(ahrs_dps),
        "gyro_dps": float(gyro_dps),
        "ahrs_ok": bool(err < AHRS_BAD_DEG),
        "n_pings": int(om.shape[-1]),
    }


# --------------------------------------------------------------------------- #
# assembled kinematics
# --------------------------------------------------------------------------- #
@dataclass
class Kinematics:
    """Per-cast platform kinematics. Angles in degrees, rates in deg/s."""
    time: np.ndarray
    fs: float
    tilt_deg: np.ndarray                # from the AHRS - only meaningful if ahrs_ok
    tilt_accel_deg: np.ndarray          # from the accelerometer - the reference
    ahrs_error_deg: float
    lean_azimuth_deg: np.ndarray
    omega_body_dps: np.ndarray
    omega_up_dps: np.ndarray
    spin_rate_dps: np.ndarray
    pitch_deg: np.ndarray
    roll_deg: np.ndarray
    depth_m: np.ndarray
    ascent_rate_ms: np.ndarray
    heave_m: np.ndarray
    accel_mag: np.ndarray
    angrt_raw: np.ndarray = field(default=None, repr=False)

    @property
    def ahrs_ok(self) -> bool:
        return self.ahrs_error_deg < AHRS_BAD_DEG

    @property
    def regime(self) -> str:
        return classify(self)

    def summary(self) -> dict:
        omega_mag = np.linalg.norm(self.omega_body_dps, axis=1)
        spin = float(np.median(self.spin_rate_dps))
        return {
            "regime": self.regime,
            "ahrs_ok": self.ahrs_ok,
            "ahrs_error_deg": self.ahrs_error_deg,
            "tilt_accel_median_deg": float(np.median(self.tilt_accel_deg)),
            "tilt_ahrs_median_deg": float(np.median(self.tilt_deg)),
            "spin_median_dps": spin,
            "period_s": float(abs(360.0 / spin)) if abs(spin) > 1e-6 else np.inf,
            "omega_mag_median_dps": float(np.median(omega_mag)),
            "circularity": circularity(self.pitch_deg, self.roll_deg),
            "ascent_median_ms": float(np.median(self.ascent_rate_ms)),
            "heave_rms_m": float(np.sqrt(np.nanmean(self.heave_m ** 2))),
            "accel_mag_median": float(np.median(self.accel_mag)),
        }


def classify(k: Kinematics) -> str:
    """Label the platform's behaviour over a cast.

    ``ahrs_fault`` is checked **first**: if the orientation solution disagrees with
    measured gravity, none of the attitude-derived labels mean anything, and the
    velocity built from that attitude is suspect.

    - ``ahrs_fault`` : orientation solution inconsistent with the accelerometer.
    - ``upright``    : small tilt, no sustained rotation — the normal case.
    - ``spinning``   : sustained rotation while near upright; harmless for geometry.
    - ``coning``     : sustained lean rotating about the vertical.
    - ``rocking``    : tilt oscillating in a plane, no sustained rotation.
    """
    if not k.ahrs_ok:
        return "ahrs_fault"
    tilt = float(np.median(k.tilt_accel_deg))
    spin = abs(float(np.median(k.spin_rate_dps)))
    circ = circularity(k.pitch_deg, k.roll_deg)
    if tilt < TILT_QUIET_DEG:
        return "spinning" if spin >= SPIN_QUIET_DPS else "upright"
    if spin >= SPIN_QUIET_DPS and (np.isnan(circ) or circ >= CIRCULARITY_CONING):
        return "coning"
    return "rocking"


def scan_casts(fn, reader, *, chunk=500_000, ens_start=0, total=None, cast_kind="up",
               min_span_dbar=40.0, thhold_s=30.0, gap_s=30.0, progress=True) -> list[dict]:
    """Stream a raw `.ad2cp` and return one AHRS-QC row per cast.

    Same chunking and cast detection as the L2 builder, but it keeps only the AHRS and
    pressure fields, so a whole deployment can be screened for attitude faults without
    reproducing the velocity product.
    """
    import time as _time

    from .casts import detect_casts
    from .index import count_ensembles

    if total is None:
        total = count_ensembles(fn)
    kinds = ("up", "down") if cast_kind == "both" else (cast_kind,)
    heavy = ("vel", "amp", "corr", "vel_b5", "amp_b5", "corr_b5", "echo")
    rows, buf = [], None
    start, t0 = int(ens_start), _time.time()

    while start < total:
        stop = min(start + chunk, total)
        ds = reader(fn, nens=[start, stop])
        ds = ds.drop_dims([d for d in ("time_b5", "range_b5") if d in ds.dims])
        ds = ds.drop_vars([v for v in heavy if v in ds], errors="ignore")
        buf = ds if buf is None else xr.concat([buf, ds], dim="time",
                                               data_vars="minimal", coords="minimal")
        press = buf["pressure"].values
        t_s = buf["time"].values.astype("datetime64[ns]").astype("int64") / 1e9
        fs = float(buf.attrs["fs"])

        casts = detect_casts(press, t_s, thhold=int(thhold_s * fs), gap_s=gap_s)
        more = stop < total
        for c in (casts[:-1] if (more and casts) else casts):
            if c.direction in kinds and np.ptp(press[c.start:c.stop + 1]) >= min_span_dbar:
                sl = slice(c.start, c.stop + 1)
                row = cast_qc(buf, sl)
                row["time"] = buf["time"].values[sl][(c.stop - c.start) // 2]
                row["direction"] = c.direction
                row["truncated"] = bool(c.truncated)
                row["pressure_min"] = float(np.nanmin(press[sl]))
                row["pressure_max"] = float(np.nanmax(press[sl]))
                rows.append(row)
        buf = buf.isel(time=slice(int(casts[-1].start), None)) if (more and casts) else None
        if progress:
            bad = sum(not r["ahrs_ok"] for r in rows)
            print(f"  [{start:>11,},{stop:>11,}) casts={len(rows):>5,} "
                  f"ahrs_fault={bad:>4,} {_time.time() - t0:.0f}s", flush=True)
        start = stop
    return rows


def kinematics(ds, sl=None) -> Kinematics:
    """Build `Kinematics` from a dolfyn Dataset (optionally one cast's slice)."""
    sl = slice(None) if sl is None else sl
    fs = float(ds.attrs["fs"])
    om = ds["orientmat"].values[:, :, sl]
    acc = ds["accel"].values[:, sl]
    tilt, azi = attitude(om)
    tilt_acc, mag = attitude_from_accel(acc)
    w_body = angular_velocity(om, fs)
    w_up = to_earth(w_body, om)[:, 2]
    depth, ascent, heave = vertical(ds["pressure"].values[sl], fs)
    return Kinematics(
        time=ds["time"].values[sl], fs=fs,
        tilt_deg=tilt, tilt_accel_deg=tilt_acc,
        ahrs_error_deg=ahrs_error(acc, om),
        lean_azimuth_deg=azi,
        omega_body_dps=np.rad2deg(w_body),
        omega_up_dps=np.rad2deg(w_up),
        spin_rate_dps=-np.rad2deg(w_up),
        pitch_deg=ds["pitch"].values[sl], roll_deg=ds["roll"].values[sl],
        depth_m=depth, ascent_rate_ms=ascent, heave_m=heave, accel_mag=mag,
        angrt_raw=ds["angrt"].values[:, sl] if "angrt" in ds else None,
    )


def _main():
    """CLI: screen a deployment for AHRS attitude faults.

        python -m ww_sig1000.platform --config /abs/path/config_adcp.json
    """
    import argparse
    import pandas as pd
    from mhkit import dolfyn
    from .config import AmbiguousConfigError, load_adcp_config

    ap = argparse.ArgumentParser(description=_main.__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None, help="path to config_adcp.json (absolute)")
    ap.add_argument("--file", default=None, help="raw .ad2cp (overrides the config)")
    ap.add_argument("--out", default=None, help="CSV to write (default: <output_dir>/ahrs_scan.csv)")
    ap.add_argument("--cast-kind", default="up", choices=["both", "up", "down"])
    ap.add_argument("--chunk", type=int, default=None)
    ap.add_argument("-y", "--yes", action="store_true")
    args = ap.parse_args()

    try:
        cfg = load_adcp_config(args.config, assume_yes=args.yes)
    except AmbiguousConfigError as e:
        ap.error(str(e))
    fn = args.file or str(cfg.ad2cp_path)
    ens_start, ens_stop = cfg.resolve_trim()
    out = Path(args.out) if args.out else cfg.output_dir / "ahrs_scan.csv"
    out.parent.mkdir(parents=True, exist_ok=True)

    print(f"[scan] {fn}")
    print(f"[scan] ensembles {ens_start:,}:{ens_stop if ens_stop else 'end'}  "
          f"cast_kind={args.cast_kind}", flush=True)
    rows = scan_casts(fn, dolfyn.read, chunk=args.chunk or cfg.chunk,
                      ens_start=ens_start, total=ens_stop, cast_kind=args.cast_kind,
                      min_span_dbar=cfg.min_span_dbar)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)

    n = len(df)
    bad = int((~df["ahrs_ok"]).sum())
    print(f"\n[scan] {n:,} casts -> {out}")
    print(f"[scan] ahrs_fault: {bad:,} ({bad / max(n, 1) * 100:.2f}%)")
    print(f"[scan] ahrs_error deg: median {df['ahrs_error_deg'].median():.2f}  "
          f"p90 {df['ahrs_error_deg'].quantile(.9):.2f}  max {df['ahrs_error_deg'].max():.2f}")
    if bad:
        print("\n[scan] worst casts:")
        w = df.nlargest(min(10, bad), "ahrs_error_deg")
        for _, r in w.iterrows():
            print(f"    {str(r['time'])[:19]}  err {r['ahrs_error_deg']:6.1f} deg  "
                  f"tilt accel {r['tilt_accel_deg']:5.1f} / ahrs {r['tilt_ahrs_deg']:5.1f}  "
                  f"gyro {r['gyro_dps']:6.2f} vs ahrs {r['ahrs_dps']:6.1f} deg/s")


if __name__ == "__main__":
    _main()
