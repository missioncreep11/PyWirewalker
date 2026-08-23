"""The v2 (buoyant-ascent) motion correction.

Synthetic records with known platform motion, so every claim of the v2 design is
checked against constructed truth: vertical from pressure only, wave-band vertical
acceleration ignored, horizontal depth-gated, spikes interpolated and flagged, and
- structurally - immunity to AHRS attitude faults.
"""
import numpy as np
import pytest
import xarray as xr

from ww_sig1000.attitude import up_from_pitch_roll
from ww_sig1000.motion import (G, H_GAIN_FULL_M, H_GAIN_ZERO_M, _depth_gain,
                               beam_motion_correction_v2)
from ww_sig1000.velocity import process_cast

FS = 8.0
COS25 = np.cos(np.deg2rad(25.0))
INTERIOR = slice(800, -800)          # clear of filter edge effects


def _inputs(n=4000, *, ascent=0.0, p0=300.0, acc_extra=None, heading=0.0):
    """(time_s, pressure, accel, heading) for an upright vehicle."""
    t = np.arange(n) / FS
    press = p0 - ascent * t
    acc = np.zeros((3, n))
    acc[2] = G
    if acc_extra is not None:
        acc = acc + acc_extra
    head = np.full(n, float(heading))
    return t, press, acc, head


# --------------------------------------------------------------------------- #
# the function
# --------------------------------------------------------------------------- #
def test_vertical_correction_comes_from_pressure():
    """Steady 0.4 m/s ascent, quiet IMU: each beam's correction is w*cos(25)."""
    t, press, acc, head = _inputs(ascent=0.4)
    v2 = beam_motion_correction_v2(t, press, acc, head, FS)
    assert v2.usable
    assert np.allclose(v2.platform_enu[2, INTERIOR], 0.4, atol=0.01)
    assert np.allclose(v2.corr_beam[:, INTERIOR], 0.4 * COS25, atol=0.01)


def test_wave_band_vertical_acceleration_is_ignored():
    """v1's core defect: on the upcast, in-band a_z is noise. v2 must not integrate it."""
    n = 4000
    t = np.arange(n) / FS
    az = 0.5 * np.sin(2 * np.pi * 0.2 * t)          # would integrate to 0.4 m/s in v1
    extra = np.zeros((3, n))
    extra[2] = az
    t, press, acc, head = _inputs(n, acc_extra=extra)
    v2 = beam_motion_correction_v2(t, press, acc, head, FS)
    assert np.abs(v2.corr_beam[:, INTERIOR]).max() < 0.01


def test_horizontal_correction_is_depth_gated():
    n = 4000
    t = np.arange(n) / FS
    ax = np.zeros((3, n))
    ax[0] = 0.3 * np.sin(2 * np.pi * 0.25 * t)      # in-band lateral shaking
    shallow = beam_motion_correction_v2(*_inputs(n, p0=10.0, acc_extra=ax), FS)
    deep = beam_motion_correction_v2(*_inputs(n, p0=300.0, acc_extra=ax), FS)
    amp = 0.3 / (2 * np.pi * 0.25)                  # integrated velocity amplitude
    h_sh = np.hypot(shallow.platform_enu[0, INTERIOR], shallow.platform_enu[1, INTERIOR])
    assert np.sqrt(np.mean(h_sh ** 2)) == pytest.approx(amp / np.sqrt(2), rel=0.15)
    h_dp = np.hypot(deep.platform_enu[0], deep.platform_enu[1])
    assert h_dp.max() < 1e-6, "gate must be exactly 0 at depth"


def test_depth_gain_shape():
    assert _depth_gain(np.array([0.0, H_GAIN_FULL_M])).tolist() == [1.0, 1.0]
    assert _depth_gain(np.array([H_GAIN_ZERO_M, 400.0])).tolist() == [0.0, 0.0]
    mid = _depth_gain(np.array([(H_GAIN_FULL_M + H_GAIN_ZERO_M) / 2]))
    assert mid[0] == pytest.approx(0.5)


def test_spike_pings_are_interpolated_and_flagged():
    t, press, acc, head = _inputs(ascent=0.4)
    k = 2000
    acc_spiked = acc.copy()
    acc_spiked[2, k:k + 3] = G + 5.0                # a hit on the wire
    clean = beam_motion_correction_v2(t, press, acc, head, FS)
    spiked = beam_motion_correction_v2(t, press, acc_spiked, head, FS)
    assert not spiked.ping_ok[k] and spiked.ping_ok[k - 1]
    # interpolation keeps the correction continuous - no v1-style step
    assert np.abs(spiked.corr_beam - clean.corr_beam).max() < 5e-3


def _coning_accel(n, amp_deg=15.0, f_hz=0.2):
    """Gravity sweeping a cone above the cutoff - the unfilterable case."""
    t = np.arange(n) / FS
    return up_from_pitch_roll(amp_deg * np.sin(2 * np.pi * f_hz * t),
                              amp_deg * np.cos(2 * np.pi * f_hz * t)) * G


def test_guard_passes_a_steady_lean_but_flags_fast_coning():
    n = 4000
    t, press, _, head = _inputs(n)
    leaning = up_from_pitch_roll(np.full(n, 13.0), np.zeros(n)) * G
    assert beam_motion_correction_v2(t, press, leaning, head, FS).usable
    assert not beam_motion_correction_v2(t, press, _coning_accel(n), head, FS).usable
    assert beam_motion_correction_v2(*_inputs(n), FS).usable


# --------------------------------------------------------------------------- #
# process_cast integration
# --------------------------------------------------------------------------- #
def _cast_ds(*, ahrs_pitch=0.0, n=3000, nc=8, seed=0):
    """Synthetic upcast; the accelerometer says upright, the AHRS says `ahrs_pitch`."""
    rng = np.random.default_rng(seed)
    zeros = np.zeros(n)
    acc = up_from_pitch_roll(zeros, zeros) * G
    step = np.timedelta64(int(1e9 / FS), "ns")
    return xr.Dataset(
        {"vel": (("beam", "range", "time"), 0.05 * rng.standard_normal((4, nc, n))),
         "corr": (("beam", "range", "time"), np.full((4, nc, n), 100, np.int16)),
         "amp": (("beam", "range", "time"), np.zeros((4, nc, n))),
         "accel": (("dirIMU", "time"), acc),
         "pitch": ("time", np.full(n, float(ahrs_pitch))),
         "roll": ("time", zeros.copy()), "heading": ("time", zeros.copy()),
         "pressure": ("time", np.linspace(60.0, 0.5, n))},
        coords={"time": np.datetime64("2024-05-02T09:00:00") + np.arange(n) * step},
        attrs={"fs": FS, "cell_size": 1.0, "blank_dist": 0.5, "beam_angle": 25})


def test_v2_cast_is_immune_to_an_ahrs_fault():
    """Same measurements, AHRS claiming upright vs 43 deg: identical v2 output."""
    good = process_cast(_cast_ds(ahrs_pitch=0.0), motion="v2", z_max=60.0)
    fault = process_cast(_cast_ds(ahrs_pitch=43.0), motion="v2", z_max=60.0)
    for k in ("velE", "velN", "velU"):
        np.testing.assert_array_equal(good[k], fault[k])
    assert good["tilt_source"] == fault["tilt_source"] == "lp_accel"


def test_v2_falls_back_to_v1_when_the_guard_fails():
    ds = _cast_ds()
    ds["accel"] = (("dirIMU", "time"), _coning_accel(ds.sizes["time"]))
    out = process_cast(ds, motion="v2", z_max=60.0)
    ref = process_cast(ds, motion="v1", z_max=60.0)
    assert out["tilt_source"] == "ahrs_fallback"
    for k in ("velE", "velN", "velU"):
        np.testing.assert_array_equal(out[k], ref[k])


def test_v1_output_is_untouched_by_the_new_code_path():
    out = process_cast(_cast_ds(), motion="v1", z_max=60.0)
    assert "tilt_source" not in out


# --------------------------------------------------------------------------- #
# L2 product bookkeeping
# --------------------------------------------------------------------------- #
def test_l2_records_motion_version_and_tilt_source():
    from ww_sig1000.l2 import build_l2
    l2 = build_l2(_cast_ds(ahrs_pitch=43.0), cast_kind="up", min_span_dbar=40.0,
                  motion="v2")
    assert l2.attrs["motion_version"] == "v2"
    assert "v2" in l2.attrs["motion_correction"]
    assert int(l2["attitude_source"].values[0]) == 1
    assert l2.attrs["n_casts_attitude_reconstructed"] == 1
    # no magnetometer in this dataset -> heading stays AHRS
    assert int(l2["heading_source"].values[0]) == 0
    assert l2.attrs["n_casts_heading_mag"] == 0

    with pytest.raises(ValueError, match="motion"):
        build_l2(_cast_ds(), motion="v3")


# --------------------------------------------------------------------------- #
# Stage 2: heading from the magnetometer
# --------------------------------------------------------------------------- #
B_FIELD = np.array([0.0, 100.0, -160.0])      # magnetic ENU: north + down, dip ~58 deg


def _mag_from_attitude(heading, pitch, roll, B=B_FIELD):
    """Forward model: the field the magnetometer would measure, (3, n).

    The earth-fixed field seen in instrument coordinates is R^T B with
    R = _tilt_heading_matrix(heading, pitch, roll) - the exact transform the
    velocity pipeline uses, so recovery through `mag_heading` pins the convention.
    """
    from ww_sig1000.transforms import _tilt_heading_matrix
    R = _tilt_heading_matrix(np.asarray(heading, float), np.asarray(pitch, float),
                             np.asarray(roll, float))
    return np.einsum("nji,j->in", R, B)


def test_mag_heading_round_trips_through_the_pipeline_convention():
    from ww_sig1000.attitude import mag_heading
    n = 100
    for h in (0.0, 45.0, 137.0, 250.0, 359.0):
        for p, r in ((0.0, 0.0), (8.0, -5.0), (-12.0, 7.0)):
            m = _mag_from_attitude(np.full(n, h), np.full(n, p), np.full(n, r))
            hm, _ = mag_heading(m, np.full(n, p), np.full(n, r))
            assert np.allclose((hm - h + 180) % 360 - 180, 0.0, atol=1e-9), (h, p, r)


def test_mag_heading_needs_the_tilt_compensation():
    """At 12 deg of tilt the uncompensated heading is degrees wrong."""
    from ww_sig1000.attitude import mag_heading
    n = 100
    m = _mag_from_attitude(np.full(n, 70.0), np.full(n, 12.0), np.full(n, 6.0))
    good, _ = mag_heading(m, np.full(n, 12.0), np.full(n, 6.0))
    bad, _ = mag_heading(m, np.zeros(n), np.zeros(n))
    assert np.allclose(good, 70.0, atol=1e-9)
    assert np.abs((bad - 70.0 + 180) % 360 - 180).max() > 2.0


def test_fit_hard_iron_recovers_a_known_offset():
    from ww_sig1000.attitude import fit_hard_iron, mag_heading
    n = 720
    heading = np.arange(n) * 0.5 % 360.0        # a spinning vehicle: full circles
    m = _mag_from_attitude(heading, np.zeros(n), np.zeros(n))
    m_iron = m + np.array([7.0, -4.0, 3.0])[:, None]
    _, mt = mag_heading(m_iron, np.zeros(n), np.zeros(n))
    cx, cy, r = fit_hard_iron(mt[0], mt[1])
    assert cx == pytest.approx(7.0, abs=0.01)
    assert cy == pytest.approx(-4.0, abs=0.01)
    assert r == pytest.approx(100.0, abs=0.1)
    # corrected heading is exact again
    h2, _ = mag_heading(m_iron, np.zeros(n), np.zeros(n), hard_iron=(cx, cy, 0.0))
    assert np.abs((h2 - heading + 180) % 360 - 180).max() < 1e-6


def test_v2_uses_the_mag_heading_and_ignores_a_faulted_ahrs_heading():
    """Same measurements, AHRS heading garbage vs truth: identical v2 output when
    the magnetometer is present - the Stage-2 analogue of the tilt immunity test."""
    true_head = 137.0
    ds_good = _cast_ds()
    ds_bad = _cast_ds()
    n = ds_good.sizes["time"]
    m = _mag_from_attitude(np.full(n, true_head), np.zeros(n), np.zeros(n))
    for ds in (ds_good, ds_bad):
        ds["mag"] = (("dirIMU", "time"), m)
    ds_good["heading"] = ("time", np.full(n, true_head))
    ds_bad["heading"] = ("time", (np.arange(n) * 6.0) % 360.0)   # spinning garbage
    good = process_cast(ds_good, motion="v2", z_max=60.0)
    bad = process_cast(ds_bad, motion="v2", z_max=60.0)
    for k in ("velE", "velN", "velU"):
        np.testing.assert_array_equal(good[k], bad[k])
    assert good["heading_source"] == bad["heading_source"] == "mag"


def test_v2_keeps_ahrs_heading_when_the_field_is_junk():
    ds = _cast_ds()
    n = ds.sizes["time"]
    rng = np.random.default_rng(1)
    ds["mag"] = (("dirIMU", "time"), 5.0 * rng.standard_normal((3, n)))  # no stable field
    out = process_cast(ds, motion="v2", z_max=60.0)
    assert out["heading_source"] == "ahrs"
