"""Accel-derived tilt reconstruction.

Synthetic accelerometer records are built from a prescribed attitude, so the truth is
known exactly. The convention round-trip is pinned against `transforms`, which is what
actually defines it for the velocity pipeline.
"""
import numpy as np
import pytest
import xarray as xr

from ww_sig1000.attitude import (
    DEFAULT_CUTOFF_HZ, TILT_LOWPASS_MAX_DEG, apply_to, compare_with_pipeline,
    gravity_direction, highfreq_fraction, noise_vs_cutoff, pitch_roll_from_up,
    reconstruct, up_from_pitch_roll,
)
from ww_sig1000.platform import G
from ww_sig1000.transforms import _tilt_heading_matrix

FS = 8.0


def _ds(pitch, roll, *, n=4000, noise=0.0, seed=0, mag=G, heading=None):
    """Dataset whose accelerometer reports gravity for the given pitch/roll."""
    pitch = np.broadcast_to(np.atleast_1d(np.asarray(pitch, float)), (n,)).copy()
    roll = np.broadcast_to(np.atleast_1d(np.asarray(roll, float)), (n,)).copy()
    acc = up_from_pitch_roll(pitch, roll) * mag
    if noise:
        rng = np.random.default_rng(seed)
        t = np.arange(n) / FS
        acc = acc + noise * np.sin(2 * np.pi * 0.12 * t) * rng.standard_normal((3, n))
    head = np.zeros(n) if heading is None else np.broadcast_to(
        np.atleast_1d(np.asarray(heading, float)), (n,)).copy()
    return xr.Dataset(
        {"accel": (("dirIMU", "time"), acc), "pitch": ("time", pitch),
         "roll": ("time", roll), "heading": ("time", head),
         "pressure": ("time", np.linspace(500, 0, n))},
        coords={"time": np.arange(n)}, attrs={"fs": FS})


# --------------------------------------------------------------------------- #
# the convention
# --------------------------------------------------------------------------- #
def test_pitch_roll_convention_matches_the_transform_the_pipeline_uses():
    """Row 2 of _tilt_heading_matrix is earth-up in instrument coords."""
    pitch = np.array([0.0, 3.0, -5.0, 12.0])
    roll = np.array([0.0, -2.0, 7.0, -9.0])
    R = _tilt_heading_matrix(np.zeros(4), pitch, roll)
    assert np.allclose(R[:, 2, :].T, up_from_pitch_roll(pitch, roll), atol=1e-12)


def test_pitch_roll_round_trip_is_exact():
    pitch = np.array([0.0, 3.0, -5.0, 12.0, -30.0])
    roll = np.array([0.0, -2.0, 7.0, -9.0, 25.0])
    p2, r2 = pitch_roll_from_up(up_from_pitch_roll(pitch, roll))
    assert np.allclose(p2, pitch, atol=1e-10)
    assert np.allclose(r2, roll, atol=1e-10)


def test_heading_does_not_affect_the_up_vector():
    """Tilt must be recoverable from the accelerometer without knowing heading."""
    for hd in (0.0, 90.0, 217.0):
        R = _tilt_heading_matrix(np.full(3, hd), np.full(3, 4.0), np.full(3, -3.0))
        assert np.allclose(R[:, 2, :].T, up_from_pitch_roll(np.full(3, 4.0),
                                                            np.full(3, -3.0)), atol=1e-12)


# --------------------------------------------------------------------------- #
# reconstruction
# --------------------------------------------------------------------------- #
def test_reconstructs_a_clean_attitude_exactly():
    rec = reconstruct(_ds(3.0, -2.0))
    interior = slice(400, -400)                       # filter edges
    assert np.allclose(rec.pitch_deg[interior], 3.0, atol=0.05)
    assert np.allclose(rec.roll_deg[interior], -2.0, atol=0.05)
    assert rec.usable


def test_lowpass_reduces_tilt_noise_substantially():
    ds = _ds(3.0, -2.0, noise=0.7)                    # realistic wave acceleration
    raw = reconstruct(ds, cutoff_hz=3.9)              # effectively unfiltered
    lp = reconstruct(ds, cutoff_hz=DEFAULT_CUTOFF_HZ)
    n_raw = np.std(np.diff(raw.pitch_deg)) / np.sqrt(2)
    n_lp = np.std(np.diff(lp.pitch_deg)) / np.sqrt(2)
    assert n_lp < n_raw / 10, f"expected a large noise reduction, got {n_raw:.4f} -> {n_lp:.4f}"


def test_reconstruction_is_unbiased_under_symmetric_noise():
    rec = reconstruct(_ds(4.0, -3.0, noise=0.7))
    interior = slice(400, -400)
    assert np.median(rec.pitch_deg[interior]) == pytest.approx(4.0, abs=0.3)
    assert np.median(rec.roll_deg[interior]) == pytest.approx(-3.0, abs=0.3)


def test_flags_when_the_accelerometer_is_not_measuring_gravity():
    rec = reconstruct(_ds(3.0, 0.0, mag=G + 2.0))     # sustained extra acceleration
    assert not rec.accel_is_gravity
    assert not rec.usable


def test_flags_when_the_vehicle_is_too_tilted_for_lowpassing():
    """A genuinely tilted vehicle sweeps gravity through the instrument frame."""
    rec = reconstruct(_ds(TILT_LOWPASS_MAX_DEG + 20.0, 0.0))
    assert not rec.lowpass_valid
    assert not rec.usable
    assert reconstruct(_ds(3.0, -2.0)).lowpass_valid


def test_highfreq_fraction_detects_real_attitude_change_above_the_cutoff():
    n = 4000
    t = np.arange(n) / FS
    steady = _ds(3.0, -2.0)
    swinging = _ds(3.0 + 5 * np.sin(2 * np.pi * 0.2 * t), -2.0, n=n)   # 5 s swing
    f_steady = highfreq_fraction(steady["accel"].values, FS)
    f_swing = highfreq_fraction(swinging["accel"].values, FS)
    assert f_swing > f_steady * 5


# --------------------------------------------------------------------------- #
# harness
# --------------------------------------------------------------------------- #
def test_compare_with_pipeline_reports_zero_bias_on_a_consistent_cast():
    cmp = compare_with_pipeline(_ds(3.0, -2.0))
    assert abs(cmp["pitch_bias_deg"]) < 0.05
    assert abs(cmp["roll_bias_deg"]) < 0.05
    assert cmp["up_angle_median_deg"] < 0.1


def test_compare_with_pipeline_measures_a_known_disagreement():
    """Pipeline claims one attitude, accelerometer reports another."""
    ds = _ds(3.0, -2.0)
    ds["pitch"] = ("time", np.full(ds.sizes["time"], 13.0))    # 10 deg wrong
    cmp = compare_with_pipeline(ds)
    assert cmp["pitch_bias_deg"] == pytest.approx(-10.0, abs=0.2)


def test_noise_vs_cutoff_is_monotonic_and_includes_an_unfiltered_row():
    rows = noise_vs_cutoff(_ds(3.0, -2.0, noise=0.7))
    assert np.isinf(rows[-1]["cutoff_hz"])
    finite = [r for r in rows if np.isfinite(r["cutoff_hz"])]
    noise = [r["tilt_noise_deg"] for r in finite]
    assert noise == sorted(noise), "noise should rise with cutoff"
    assert finite[0]["highfreq_fraction"] >= finite[-1]["highfreq_fraction"]


def test_apply_to_replaces_pitch_roll_and_leaves_heading_alone():
    ds = _ds(3.0, -2.0, heading=137.0)
    out, rec = apply_to(ds)
    interior = slice(400, -400)
    assert np.allclose(out["pitch"].values[interior], rec.pitch_deg[interior])
    assert np.array_equal(out["heading"].values, ds["heading"].values)
    assert "accel-reconstructed" in out.attrs["attitude_source"]
    assert np.array_equal(ds["pitch"].values, _ds(3.0, -2.0)["pitch"].values), "input untouched"
