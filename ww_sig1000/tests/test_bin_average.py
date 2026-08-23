"""Per-bin SEM and the depth-gated notch bin-average.

Synthetic casts with uniform correlation make the SEM analytically exact; a
prescribed in-band sinusoid makes the notch's suppression measurable against
the boxcar with known truth.
"""
import numpy as np
import pytest
import xarray as xr

from ww_sig1000.attitude import up_from_pitch_roll
from ww_sig1000.l2 import build_l2
from ww_sig1000.platform import G
from ww_sig1000.velocity import (NOTCH_MAX_DEPTH_M, _SIGMA_CORR_MID, _SIGMA_CORR_VAL,
                                 beam_sigma, process_cast)

FS = 8.0
COS25 = np.cos(np.deg2rad(25.0))


def _ds(*, n=3000, nc=8, p0=60.0, corr=100, vel_fn=None, seed=0):
    """Synthetic upcast; vel_fn(t_s) gives the beam velocity (same on all beams)."""
    rng = np.random.default_rng(seed)
    t_s = np.arange(n) / FS
    beam = vel_fn(t_s) if vel_fn else 0.05 * rng.standard_normal(n)
    zeros = np.zeros(n)
    step = np.timedelta64(int(1e9 / FS), "ns")
    return xr.Dataset(
        {"vel": (("beam", "range", "time"), np.broadcast_to(beam, (4, nc, n)).copy()),
         "corr": (("beam", "range", "time"), np.full((4, nc, n), corr, np.int16)),
         "amp": (("beam", "range", "time"), np.zeros((4, nc, n))),
         "accel": (("dirIMU", "time"), up_from_pitch_roll(zeros, zeros) * G),
         "pitch": ("time", zeros.copy()), "roll": ("time", zeros.copy()),
         "heading": ("time", zeros.copy()),
         "pressure": ("time", np.linspace(p0, 0.5, n))},
        coords={"time": np.datetime64("2024-05-02T09:00:00") + np.arange(n) * step},
        attrs={"fs": FS, "cell_size": 1.0, "blank_dist": 0.5, "beam_angle": 25})


# --------------------------------------------------------------------------- #
# the noise relation and the SEM
# --------------------------------------------------------------------------- #
def test_beam_sigma_is_monotone_decreasing():
    s = beam_sigma(np.arange(30, 101))
    assert np.all(np.diff(s) <= 0)
    assert beam_sigma(97.5) == pytest.approx(0.0430)
    assert beam_sigma(150) == _SIGMA_CORR_VAL[-1], "clamped beyond the table"
    assert _SIGMA_CORR_MID.size == _SIGMA_CORR_VAL.size


def test_sem_matches_the_analytic_value_at_uniform_corr():
    """corr = 100 everywhere: sem = geom_factor * sigma(100) / sqrt(n_obs)."""
    out = process_cast(_ds(corr=100), z_max=60.0, motion_correct=False)
    sig = beam_sigma(100.0)
    th = np.deg2rad(25.0)
    n = out["n_obs"].astype(float)
    ok = n >= 10
    expect_E = sig / (np.sqrt(2) * np.sin(th)) / np.sqrt(n[ok])
    expect_U = sig / (2 * np.cos(th)) / np.sqrt(n[ok])
    assert np.allclose(out["velE_sem"][ok], expect_E, rtol=1e-5)
    assert np.allclose(out["velN_sem"][ok], expect_E, rtol=1e-5)
    assert np.allclose(out["velU_sem"][ok], expect_U, rtol=1e-5)
    assert np.all(np.isnan(out["velE_sem"][~ok]))


def test_sem_grows_when_correlation_drops():
    hi = process_cast(_ds(corr=95), z_max=60.0, motion_correct=False)
    lo = process_cast(_ds(corr=55), z_max=60.0, motion_correct=False)
    ok = (hi["n_obs"] >= 10) & (lo["n_obs"] >= 10)
    assert np.all(lo["velE_sem"][ok] > 1.5 * hi["velE_sem"][ok])


# --------------------------------------------------------------------------- #
# the notch
# --------------------------------------------------------------------------- #
def _gappy_corr(n, nc, seed=5, block=16, frac=0.5):
    """corr field with contiguous low-corr dropouts, like real dwells (which keep
    only 17-32% of near-surface samples)."""
    rng = np.random.default_rng(seed)
    corr = np.full((4, nc, n), 100, np.int16)
    nblk = n // block
    bad = rng.random((nc, nblk)) < frac
    for c in range(nc):
        for b in np.flatnonzero(bad[c]):
            corr[:, c, b * block:(b + 1) * block] = 30
    return corr


def _wave_ds(*, gappy=True, **kw):
    """Constant beam velocity + an in-band 0.2 Hz contamination."""
    ds = _ds(vel_fn=lambda t: 0.10 + 0.30 * np.cos(2 * np.pi * 0.2 * t), **kw)
    if gappy:
        ds["corr"] = (("beam", "range", "time"),
                      _gappy_corr(ds.sizes["time"], ds.sizes["range"]))
    return ds


def test_notch_suppresses_in_band_contamination_on_gappy_dwells():
    """velU truth is 0.10/cos(25). On a gappy dwell (the real case) the boxcar's
    wave cancellation breaks and the notch, which models the wave explicitly on
    the actual sample support, recovers it."""
    truth = 0.10 / COS25
    bx = process_cast(_wave_ds(), z_max=60.0, motion_correct=False)
    nt = process_cast(_wave_ds(), z_max=60.0, motion_correct=False,
                      bin_average="notch")
    ok = bx["n_obs"] >= 20
    e_bx = np.sqrt(np.mean((bx["velU"][ok] - truth) ** 2))
    e_nt = np.sqrt(np.mean((nt["velU"][ok] - truth) ** 2))
    # measured: ~30% rms reduction on this construction (block gaps aliasing a
    # 5 s wave are adversarial; the noisy real-geometry study showed 7-17%)
    assert e_nt < 0.8 * e_bx, (e_nt, e_bx)


def test_notch_is_near_boxcar_on_uniform_dwells():
    """With complete uniform sampling the harmonic columns are ~orthogonal to
    the constant, so the notch constant is essentially the mean (exact only for
    a periodic window; short edge dwells leave a few-mm/s difference). The
    estimator's gain lives entirely on gappy dwells - the real case."""
    bx = process_cast(_wave_ds(gappy=False), z_max=60.0, motion_correct=False)
    nt = process_cast(_wave_ds(gappy=False), z_max=60.0, motion_correct=False,
                      bin_average="notch")
    ok = bx["n_obs"] >= 20
    assert np.allclose(nt["velU"][ok], bx["velU"][ok], atol=5e-3)


def test_notch_is_identical_to_boxcar_below_the_gate():
    ds = _wave_ds(p0=120.0)                       # bins straddle the 60 m gate
    bx = process_cast(ds, z_max=120.0, motion_correct=False)
    nt = process_cast(ds, z_max=120.0, motion_correct=False, bin_average="notch")
    z = bx["z"]
    deep = z > NOTCH_MAX_DEPTH_M
    np.testing.assert_array_equal(nt["velU"][deep], bx["velU"][deep])
    shallow_ok = (~deep) & (bx["n_obs"] >= 20)
    assert not np.allclose(nt["velU"][shallow_ok], bx["velU"][shallow_ok]), \
        "gated bins should actually be refit"


def test_notch_leaves_sem_and_n_obs_untouched():
    bx = process_cast(_wave_ds(), z_max=60.0, motion_correct=False)
    nt = process_cast(_wave_ds(), z_max=60.0, motion_correct=False,
                      bin_average="notch")
    np.testing.assert_array_equal(nt["n_obs"], bx["n_obs"])
    np.testing.assert_array_equal(nt["velE_sem"], bx["velE_sem"])


# --------------------------------------------------------------------------- #
# product plumbing
# --------------------------------------------------------------------------- #
def test_l2_carries_sem_variables_and_bin_average_attr():
    l2 = build_l2(_ds(), cast_kind="up", min_span_dbar=40.0, motion_correct=False)
    for v in ("velE_sem", "velN_sem", "velU_sem"):
        assert v in l2
        assert l2[v].dims == ("depth", "cast")
    assert l2.attrs["bin_average"] == "boxcar mean"

    l2n = build_l2(_ds(), cast_kind="up", min_span_dbar=40.0, motion_correct=False,
                   bin_average="notch")
    assert "notch" in l2n.attrs["bin_average"]
    with pytest.raises(ValueError, match="bin_average"):
        build_l2(_ds(), bin_average="median")
