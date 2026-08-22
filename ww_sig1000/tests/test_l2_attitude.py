"""Attitude-source selection in the L2 builder (--attitude ahrs|reconstructed|auto).

Synthetic datasets carry a known true attitude (via the accelerometer) and a
separately prescribed AHRS claim (via `orientmat` + the pitch/roll scalars), so a
fault is constructed exactly and the selection logic can be checked cast by cast.
"""
import numpy as np
import pytest
import xarray as xr

from ww_sig1000.attitude import up_from_pitch_roll
from ww_sig1000.l2 import ATTITUDE_MODES, _cast_attitude, build_l2
from ww_sig1000.platform import AHRS_BAD_DEG, G
from ww_sig1000.transforms import _tilt_heading_matrix

FS = 4.0


def _orientmat(pitch, roll, heading=None):
    """Dolfyn-shaped (3, 3, n) orientation matrix from pitch/roll/heading (deg)."""
    n = np.size(pitch)
    head = np.zeros(n) if heading is None else np.asarray(heading, float)
    R = _tilt_heading_matrix(head, np.asarray(pitch, float), np.asarray(roll, float))
    return np.moveaxis(R, 0, -1)


def _ds(true_pitch=0.0, *, ahrs_pitch=None, n=2000, mag=G, with_orientmat=True):
    """Cast-like dataset: the accelerometer reports `true_pitch`, the AHRS
    (orientmat + scalars) claims `ahrs_pitch` (defaults to the truth = healthy)."""
    ahrs_pitch = true_pitch if ahrs_pitch is None else ahrs_pitch
    tp = np.full(n, float(true_pitch))
    ap = np.full(n, float(ahrs_pitch))
    zeros = np.zeros(n)
    data = {"accel": (("dirIMU", "time"), up_from_pitch_roll(tp, zeros) * mag),
            "pitch": ("time", ap.copy()), "roll": ("time", zeros.copy()),
            "heading": ("time", zeros.copy()),
            "pressure": ("time", np.linspace(60.0, 0.5, n))}
    if with_orientmat:
        data["orientmat"] = (("earth", "inst", "time"), _orientmat(ap, zeros))
    return xr.Dataset(data, coords={"time": np.arange(n)}, attrs={"fs": FS})


# --------------------------------------------------------------------------- #
# _cast_attitude: the per-cast selection
# --------------------------------------------------------------------------- #
def test_ahrs_mode_keeps_the_pipeline_attitude_but_still_records_the_error():
    ds = _ds(0.0, ahrs_pitch=43.0)                 # badly faulted
    out, src, err = _cast_attitude(ds, "ahrs")
    assert src == 0
    assert err == pytest.approx(43.0, abs=0.5)     # QC recorded even when unused
    assert np.array_equal(out["pitch"].values, ds["pitch"].values)


def test_auto_keeps_ahrs_on_a_healthy_cast():
    out, src, err = _cast_attitude(_ds(3.0), "auto")
    assert src == 0
    assert err < AHRS_BAD_DEG
    assert np.allclose(out["pitch"].values, 3.0)


def test_auto_reconstructs_a_faulted_cast():
    ds = _ds(2.0, ahrs_pitch=43.0)
    out, src, err = _cast_attitude(ds, "auto")
    assert src == 1
    assert err == pytest.approx(41.0, abs=0.5)
    interior = slice(400, -400)                    # filter edges
    assert np.allclose(out["pitch"].values[interior], 2.0, atol=0.1), \
        "pitch should follow the accelerometer, not the faulted AHRS"
    assert np.allclose(out["roll"].values[interior], 0.0, atol=0.1)
    assert np.array_equal(ds["pitch"].values, np.full(ds.sizes["time"], 43.0)), \
        "input untouched"


def test_auto_falls_back_when_the_reconstruction_is_unusable():
    ds = _ds(2.0, ahrs_pitch=43.0, mag=G + 2.0)    # |a| != g: accel not gravity
    out, src, err = _cast_attitude(ds, "auto")
    assert src == 2
    assert err >= AHRS_BAD_DEG
    assert np.array_equal(out["pitch"].values, ds["pitch"].values), "AHRS kept"


def test_reconstructed_mode_substitutes_even_on_healthy_casts():
    _, src, err = _cast_attitude(_ds(3.0), "reconstructed")
    assert src == 1
    assert err < AHRS_BAD_DEG


def test_auto_without_orientmat_cannot_detect_and_keeps_ahrs():
    ds = _ds(3.0, with_orientmat=False)
    _, src, err = _cast_attitude(ds, "auto")
    assert src == 0
    assert np.isnan(err)


def test_build_l2_rejects_an_unknown_mode():
    with pytest.raises(ValueError, match="attitude"):
        build_l2(_ds(0.0), attitude="kalman")
    assert set(ATTITUDE_MODES) == {"ahrs", "reconstructed", "auto"}


# --------------------------------------------------------------------------- #
# end to end: the choice and outcome land in the product
# --------------------------------------------------------------------------- #
def _full_ds(*, ahrs_pitch=0.0, n=3000, nc=8):
    """A single synthetic upcast with the full variable set build_l2 needs."""
    ds = _ds(0.0, ahrs_pitch=ahrs_pitch, n=n)
    ds["vel"] = (("beam", "range", "time"), np.zeros((4, nc, n)))
    ds["corr"] = (("beam", "range", "time"), np.full((4, nc, n), 100, np.int16))
    ds["amp"] = (("beam", "range", "time"), np.zeros((4, nc, n)))
    step = np.timedelta64(int(1e9 / FS), "ns")
    ds = ds.assign_coords(time=np.datetime64("2024-05-02T09:00:00") + np.arange(n) * step)
    ds.attrs.update(cell_size=1.0, blank_dist=0.5, beam_angle=25)
    return ds


def test_l2_product_records_attitude_source_and_error_per_cast():
    l2 = build_l2(_full_ds(ahrs_pitch=43.0), cast_kind="up", min_span_dbar=40.0,
                  motion_correct=False, attitude="auto")
    assert l2.attrs["attitude_mode"] == "auto"
    assert l2.attrs["n_casts_attitude_reconstructed"] == 1
    assert l2.attrs["n_casts_attitude_fallback"] == 0
    assert int(l2["attitude_source"].values[0]) == 1
    assert float(l2["ahrs_error_deg"].values[0]) == pytest.approx(43.0, abs=0.5)


def test_l2_product_default_mode_is_ahrs_and_flags_stay_zero():
    l2 = build_l2(_full_ds(ahrs_pitch=43.0), cast_kind="up", min_span_dbar=40.0,
                  motion_correct=False)
    assert l2.attrs["attitude_mode"] == "ahrs"
    assert l2.attrs["n_casts_attitude_reconstructed"] == 0
    assert int(l2["attitude_source"].values[0]) == 0
    # the QC channel is present regardless of mode
    assert float(l2["ahrs_error_deg"].values[0]) == pytest.approx(43.0, abs=0.5)
