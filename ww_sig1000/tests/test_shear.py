"""Beam-differenced shear (port of WWvel_upward's beamshear).

The forward model builds beam velocities from a prescribed water-velocity
profile through the same beam geometry the estimator inverts, so recovery is
exact-by-construction; the immunity test then proves the estimator's defining
property - per-ping common-mode contamination cancels identically.
"""
import numpy as np
import pytest
import xarray as xr

from ww_sig1000.attitude import up_from_pitch_roll
from ww_sig1000.geometry import cell_depths
from ww_sig1000.platform import G
from ww_sig1000.velocity import process_cast

FS = 8.0
NC = 12
S25, C25 = np.sin(np.deg2rad(25.0)), np.cos(np.deg2rad(25.0))


def _ds(u_of_z, *, n=3000, p0=60.0, ping_offset=None, seed=0):
    """Beam velocities for an eastward flow u(z) seen by an upright vehicle.

    Forward model built from the pipeline's own transform: at heading 0 the
    H(heading-90) rotation maps East to instrument -y, so for water velocity
    (uE, 0, 0): (vx, vy, vz) = (0, -uE, 0), and the beam projections consistent
    with ``beam2xyz_matrix`` are b1 = vx s + vz c, b2 = -vy s + vz c,
    b3 = -vx s + vz c, b4 = vy s + vz c.

    ``ping_offset`` (4, n) adds a per-ping constant to every cell of each beam -
    exactly what platform motion or attitude leakage injects.
    """
    press = np.linspace(p0, 0.5, n)
    zeros = np.zeros(n)
    z, _, _ = cell_depths(press, zeros, zeros, NC, 1.0, 0.5)   # (n, NC, 4), +up
    vel = np.zeros((4, NC, n))
    for b in range(4):
        uE = u_of_z(z[:, :, b]).T                              # (NC, n) at beam cells
        vy = -uE                                               # East -> inst -y
        vel[b] = {0: 0.0 * vy, 1: -vy * S25, 2: 0.0 * vy, 3: vy * S25}[b]
    if ping_offset is not None:
        vel = vel + ping_offset[:, None, :]
    step = np.timedelta64(int(1e9 / FS), "ns")
    return xr.Dataset(
        {"vel": (("beam", "range", "time"), vel),
         "corr": (("beam", "range", "time"), np.full((4, NC, n), 100, np.int16)),
         "amp": (("beam", "range", "time"), np.zeros((4, NC, n))),
         "accel": (("dirIMU", "time"), up_from_pitch_roll(zeros, zeros) * G),
         "pitch": ("time", zeros.copy()), "roll": ("time", zeros.copy()),
         "heading": ("time", zeros.copy()), "pressure": ("time", press)},
        coords={"time": np.datetime64("2024-05-02T09:00:00") + np.arange(n) * step},
        attrs={"fs": FS, "cell_size": 1.0, "blank_dist": 0.5, "beam_angle": 25})


def test_recovers_a_linear_shear_exactly():
    S = 0.02                                   # s-1, eastward
    out = process_cast(_ds(lambda z: S * z), z_max=60.0, motion_correct=False)
    ok = np.isfinite(out["shearE"]) & (out["n_obs"] >= 10)
    assert ok.sum() > 20
    assert np.allclose(out["shearE"][ok], S, atol=1e-3)
    assert np.allclose(out["shearN"][ok], 0.0, atol=1e-3)


def test_shear_is_immune_to_per_ping_common_mode():
    """Add violent fake 'platform motion' (a per-ping constant on every cell):
    the velocities are wrecked, the shear must not move."""
    S = 0.02
    rng = np.random.default_rng(2)
    n = 3000
    contam = 0.5 * rng.standard_normal((1, n)) * np.ones((4, 1))   # same on all beams
    clean = process_cast(_ds(lambda z: S * z), z_max=60.0, motion_correct=False)
    dirty = process_cast(_ds(lambda z: S * z, ping_offset=contam), z_max=60.0,
                         motion_correct=False)
    ok = np.isfinite(clean["shearE"]) & np.isfinite(dirty["shearE"])
    np.testing.assert_allclose(dirty["shearE"][ok], clean["shearE"][ok], atol=1e-10)
    # ...while the velocity product is visibly damaged by the same contamination
    okv = np.isfinite(clean["velU"]) & np.isfinite(dirty["velU"])
    assert np.abs(dirty["velU"][okv] - clean["velU"][okv]).max() > 0.01


def test_curved_profile_and_sem():
    """Quadratic u(z): recovered shear tracks du/dz; SEMs are finite and positive."""
    out = process_cast(_ds(lambda z: 0.001 * z ** 2), z_max=60.0, motion_correct=False)
    zc = out["z"]
    ok = np.isfinite(out["shearE"]) & (out["n_obs"] >= 10) & (zc > 5) & (zc < 55)
    expect = 0.001 * 2 * (-zc[ok])             # z positive up; depth grid positive down
    assert np.corrcoef(out["shearE"][ok], expect)[0, 1] > 0.999
    assert np.all(out["shearE_sem"][ok] > 0)
    assert np.all(np.isfinite(out["shearN_sem"][ok]))


def test_l2_carries_shear_variables():
    from ww_sig1000.l2 import build_l2
    l2 = build_l2(_ds(lambda z: 0.02 * z), cast_kind="up", min_span_dbar=40.0,
                  motion_correct=False)
    for v in ("shearE", "shearN", "shearE_sem", "shearN_sem"):
        assert v in l2
        assert l2[v].dims == ("depth", "cast")
    assert "beam-differenced" in l2["shearE"].attrs["long_name"]
