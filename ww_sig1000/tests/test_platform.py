"""Platform kinematics and AHRS validation, against synthetic motions.

Geometry is built from exact ZYZ Euler rotations ``R = Rz(alpha) Ry(tilt) Rz(beta)``
mapping instrument -> earth (E, N, U), so every quantity has a closed-form expected
value and no data files are needed.

The AHRS-fault cases give the accelerometer one attitude and ``orientmat`` another —
exactly the failure seen on NOPP_d2 at 2024-05-02 09:51, where the accelerometer
measured a 2.2 deg tilt while the orientation solution claimed 43 deg.
"""
import numpy as np
import pytest
import xarray as xr

from ww_sig1000.platform import (
    G, Kinematics, ahrs_error, angular_velocity, attitude, attitude_from_accel,
    cast_qc, circularity, classify, kinematics, spectrum, to_earth, vertical,
)

FS = 8.0


# --------------------------------------------------------------------------- #
# synthetic attitude
# --------------------------------------------------------------------------- #
def _rz(a):
    n, c, s = a.size, np.cos(a), np.sin(a)
    R = np.zeros((n, 3, 3))
    R[:, 0, 0], R[:, 0, 1] = c, -s
    R[:, 1, 0], R[:, 1, 1] = s, c
    R[:, 2, 2] = 1.0
    return R


def _ry(t):
    n, c, s = t.size, np.cos(t), np.sin(t)
    R = np.zeros((n, 3, 3))
    R[:, 0, 0], R[:, 0, 2] = c, s
    R[:, 1, 1] = 1.0
    R[:, 2, 0], R[:, 2, 2] = -s, c
    return R


def _om(alpha_deg=0.0, tilt_deg=0.0, beta_deg=0.0, n=1):
    """(3,3,n) inst->earth from ZYZ Euler angles, in degrees."""
    a, t, b = (np.deg2rad(np.broadcast_to(np.atleast_1d(np.asarray(x, float)), (n,)).copy())
               for x in (alpha_deg, tilt_deg, beta_deg))
    return np.moveaxis(_rz(a) @ _ry(t) @ _rz(b), 0, -1)


def _gravity_in_inst(om):
    """Accelerometer reading for a body at rest with attitude `om`: R^T (0,0,g)."""
    R = np.moveaxis(om, -1, 0)
    return np.einsum("nji,j->in", R, np.array([0.0, 0.0, G]))


def _ds(om, accel=None, pressure=None, angrt=None):
    n = om.shape[-1]
    accel = _gravity_in_inst(om) if accel is None else accel
    pitch = np.rad2deg(np.arcsin(np.clip(-om[2, 0], -1, 1)))
    roll = np.rad2deg(np.arctan2(om[2, 1], om[2, 2]))
    data = {"orientmat": (("earth", "inst", "time"), om),
            "pitch": ("time", pitch), "roll": ("time", roll),
            "pressure": ("time", np.linspace(500, 0, n) if pressure is None else pressure),
            "accel": (("dirIMU", "time"), accel)}
    if angrt is not None:
        data["angrt"] = (("dirIMU", "time"), angrt)
    return xr.Dataset(data, coords={"time": np.arange(n)}, attrs={"fs": FS})


def _kin(tilt_accel, spin_dps, pitch, roll, ahrs_err=0.0):
    """A Kinematics carrying just the fields `classify` reads."""
    n = np.size(pitch)
    return Kinematics(
        time=np.arange(n), fs=FS,
        tilt_deg=np.full(n, tilt_accel), tilt_accel_deg=np.full(n, tilt_accel),
        ahrs_error_deg=ahrs_err, lean_azimuth_deg=np.zeros(n),
        omega_body_dps=np.zeros((n, 3)), omega_up_dps=np.full(n, -spin_dps),
        spin_rate_dps=np.full(n, spin_dps),
        pitch_deg=np.asarray(pitch, float), roll_deg=np.asarray(roll, float),
        depth_m=np.zeros(n), ascent_rate_ms=np.zeros(n), heave_m=np.zeros(n),
        accel_mag=np.full(n, G))


# --------------------------------------------------------------------------- #
# attitude basics
# --------------------------------------------------------------------------- #
def test_attitude_recovers_prescribed_tilt():
    tilt, _ = attitude(_om(tilt_deg=[0.0, 10.0, 43.0, 60.0], n=4))
    assert np.allclose(tilt, [0.0, 10.0, 43.0, 60.0], atol=1e-6)


def test_lean_azimuth_is_a_compass_bearing():
    """alpha rotates the lean counterclockwise, so the compass bearing is 90 - alpha."""
    _, azi = attitude(_om(alpha_deg=[0.0, 90.0, 180.0], tilt_deg=30.0, n=3))
    assert np.allclose(azi, [90.0, 0.0, 270.0], atol=1e-6)


def test_attitude_from_accel_matches_the_orientation_it_was_built_from():
    om = _om(alpha_deg=[0.0, 30.0, 200.0], tilt_deg=[0.0, 5.0, 43.0], n=3)
    tilt, mag = attitude_from_accel(_gravity_in_inst(om))
    assert np.allclose(tilt, [0.0, 5.0, 43.0], atol=1e-6)
    assert np.allclose(mag, G, atol=1e-9)


def test_angular_velocity_matches_a_known_rotation_rate():
    """R = Rz(alpha) Ry(tilt) with constant d(alpha)/dt gives omega = alpha_dot * z."""
    n, rate = 800, 48.0                       # deg/s about the vertical
    om = _om(alpha_deg=rate * np.arange(n) / FS, tilt_deg=43.0, n=n)
    w_up = np.rad2deg(to_earth(angular_velocity(om, FS), om)[:, 2])
    interior = slice(5, -5)
    assert np.median(w_up[interior]) == pytest.approx(rate, rel=0.02)
    mag = np.rad2deg(np.linalg.norm(angular_velocity(om, FS), axis=1))
    assert np.median(np.abs(w_up[interior])) / np.median(mag[interior]) > 0.98


def test_lean_azimuth_advances_opposite_to_omega_up():
    n, rate = 800, 36.0
    om = _om(alpha_deg=rate * np.arange(n) / FS, tilt_deg=30.0, n=n)
    _, azi = attitude(om)
    w_up = to_earth(angular_velocity(om, FS), om)[:, 2]
    interior = slice(5, -5)
    azi_rate = np.gradient(np.unwrap(np.deg2rad(azi)), 1 / FS)
    assert np.median(azi_rate[interior]) == pytest.approx(-np.median(w_up[interior]), rel=0.02)


def test_circularity_separates_circular_from_planar():
    n = 800
    ph = np.linspace(0, 8 * np.pi, n)
    assert circularity(30 * np.cos(ph), 30 * np.sin(ph)) > 0.9      # circle
    assert circularity(30 * np.cos(ph), 30 * np.cos(ph)) < 0.05     # line


# --------------------------------------------------------------------------- #
# the detector
# --------------------------------------------------------------------------- #
def test_ahrs_error_is_zero_for_a_consistent_solution():
    om = _om(alpha_deg=10.0, tilt_deg=43.0, n=200)
    assert ahrs_error(_gravity_in_inst(om), om) == pytest.approx(0.0, abs=1e-6)


@pytest.mark.parametrize("true_tilt, reported_tilt", [(0.0, 43.0), (2.2, 42.8), (5.0, 30.0)])
def test_ahrs_error_recovers_the_mismatch(true_tilt, reported_tilt):
    n = 400
    truth = _om(tilt_deg=true_tilt, n=n)
    claimed = _om(tilt_deg=reported_tilt, n=n)
    assert ahrs_error(_gravity_in_inst(truth), claimed) == pytest.approx(
        abs(reported_tilt - true_tilt), abs=0.5)


def test_ahrs_error_does_not_false_positive_on_wave_accelerations():
    """A healthy cast must stay under the threshold despite real vehicle motion.

    0.7 m/s2 rms is the worst non-gravity acceleration measured across healthy NOPP_d2
    casts; it puts a ~4 deg floor on the per-ping angle, well under AHRS_BAD_DEG.
    """
    from ww_sig1000.platform import AHRS_BAD_DEG
    n = 4000
    om = _om(tilt_deg=3.0, n=n)
    rng = np.random.default_rng(0)
    acc = _gravity_in_inst(om) + 0.7 * rng.standard_normal((3, n))
    err = ahrs_error(acc, om)
    assert err < AHRS_BAD_DEG
    assert err == pytest.approx(4.0, abs=3.0)      # the expected noise floor


def test_cast_qc_passes_a_healthy_cast():
    om = _om(tilt_deg=2.2, n=600)
    q = cast_qc(_ds(om))
    assert q["ahrs_ok"] and q["ahrs_error_deg"] < 1.0
    assert q["accel_is_gravity"]
    assert q["tilt_accel_deg"] == pytest.approx(2.2, abs=0.1)


def test_cast_qc_flags_the_nopp_d2_failure():
    """Accelerometer says 2.2 deg upright; orientmat claims 43 deg rotating at 48 deg/s."""
    n = 600
    acc = _gravity_in_inst(_om(tilt_deg=2.2, n=n))
    bad = _om(alpha_deg=48.0 * np.arange(n) / FS, tilt_deg=43.0, n=n)
    q = cast_qc(_ds(bad, accel=acc))
    assert not q["ahrs_ok"]
    assert q["ahrs_error_deg"] == pytest.approx(43.0, abs=3.0)
    assert q["tilt_accel_deg"] == pytest.approx(2.2, abs=0.5)
    assert q["tilt_ahrs_deg"] == pytest.approx(43.0, abs=0.5)
    assert q["accel_is_gravity"], "the accelerometer itself is healthy"


def test_cast_qc_reports_the_gyro_ahrs_rate_disagreement():
    """The NOPP_d2 tell: the gyro reports ~3 deg/s while the solution implies ~48."""
    n = 600
    acc = _gravity_in_inst(_om(tilt_deg=2.2, n=n))
    bad = _om(alpha_deg=48.0 * np.arange(n) / FS, tilt_deg=43.0, n=n)
    q = cast_qc(_ds(bad, accel=acc, angrt=np.full((3, n), np.deg2rad(3.0))))
    assert q["gyro_dps"] == pytest.approx(3.0, abs=0.5)
    assert q["ahrs_dps"] > 40
    assert not q["ahrs_ok"]


# --------------------------------------------------------------------------- #
# classification
# --------------------------------------------------------------------------- #
def test_classify_checks_ahrs_health_before_anything_else():
    n = 600
    ph = np.linspace(0, 8 * np.pi, n)
    k = _kin(2.2, 48.0, 30 * np.cos(ph), 30 * np.sin(ph), ahrs_err=43.0)
    assert k.regime == "ahrs_fault" and not k.ahrs_ok


def test_classify_labels_the_healthy_regimes():
    n = 600
    ph = np.linspace(0, 8 * np.pi, n)
    circle = (30 * np.cos(ph), 30 * np.sin(ph))
    line = (30 * np.cos(ph), 30 * np.cos(ph))
    assert _kin(3.0, 1.0, *line).regime == "upright"
    assert _kin(3.0, 40.0, *line).regime == "spinning"
    assert _kin(43.0, 48.0, *circle).regime == "coning"
    assert _kin(30.0, 1.0, *line).regime == "rocking"


def test_kinematics_end_to_end_on_a_healthy_record():
    om = _om(tilt_deg=3.0, n=800)
    k = kinematics(_ds(om))
    assert k.ahrs_ok and k.regime == "upright"
    s = k.summary()
    assert s["tilt_accel_median_deg"] == pytest.approx(3.0, abs=0.1)
    assert s["accel_mag_median"] == pytest.approx(G, abs=0.01)


# --------------------------------------------------------------------------- #
# vertical
# --------------------------------------------------------------------------- #
def test_vertical_separates_climb_from_wave_heave():
    n = 4000
    t = np.arange(n) / FS
    climb, wave = 0.5, 0.4 * np.sin(2 * np.pi * 0.12 * t)
    _, ascent, heave = vertical(500 - climb * t + wave, FS)
    interior = slice(200, -200)
    assert np.median(ascent[interior]) == pytest.approx(climb, rel=0.05)
    assert np.sqrt(np.mean(heave[interior] ** 2)) == pytest.approx(0.4 / np.sqrt(2), rel=0.2)


def test_spectrum_finds_the_heave_frequency():
    t = np.arange(4000) / FS
    f, p = spectrum(np.sin(2 * np.pi * 0.12 * t), FS)
    assert f[np.argmax(p)] == pytest.approx(0.12, abs=0.02)
