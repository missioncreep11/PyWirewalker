"""Data-free unit tests for ww_sig1000.transforms (beam<->XYZ<->ENU).

A real-data cross-check against dolfyn.rotate2 lives in dev scripts; these tests
pin the algebra so the ported transforms can't silently drift.
Run: pytest ww_sig1000/tests/
"""
import numpy as np
from ww_sig1000.transforms import beam2xyz_matrix, beam2enu, xyz2enu, BEAM_ANGLE_DEG


def test_beam2xyz_matrix_values():
    T = beam2xyz_matrix(25.0)
    a = 1.0 / np.sin(np.deg2rad(25.0)) / 2.0   # 1.1831...
    c = 1.0 / np.cos(np.deg2rad(25.0)) / 4.0   # 0.2759...
    expected = np.array([[a, 0, -a, 0], [0, -a, 0, a], [c, c, c, c]])
    np.testing.assert_allclose(T, expected, rtol=1e-12)
    # matches the beam2inst coefficients DOLfYN reads from the file (1/(2 sin25))
    np.testing.assert_allclose(a, 1.1831, atol=1e-4)


def test_identity_attitude_gives_xyz():
    """heading=90 (so heading-90=0), pitch=roll=0 -> ENU rotation is identity,
    so beam2enu == beam2xyz."""
    rng = np.random.default_rng(0)
    beam = rng.standard_normal((4, 5, 3))          # (beam, cell, ping)
    head = np.full(3, 90.0); pr = np.zeros(3)
    enu = beam2enu(beam, head, pr, pr)
    T = beam2xyz_matrix()
    xyz = np.einsum("ik,kcn->icn", T, beam)
    np.testing.assert_allclose(enu, xyz, atol=1e-12)


def test_heading_rotates_horizontal():
    """With zero tilt, XYZ->ENU is a pure heading rotation of the horizontal plane.
    A pure +X instrument velocity should rotate into the E/N plane by (heading-90)."""
    xyz = np.zeros((3, 1, 4))
    xyz[0, 0, :] = 1.0                              # unit X
    heading = np.array([90.0, 180.0, 0.0, 45.0])
    enu = xyz2enu(xyz, heading, np.zeros(4), np.zeros(4))
    th = np.deg2rad(heading - 90.0)
    # R = H(th); E = cos*X, N = -sin*X  (per the toolbox H matrix)
    np.testing.assert_allclose(enu[0, 0, :], np.cos(th), atol=1e-12)
    np.testing.assert_allclose(enu[1, 0, :], -np.sin(th), atol=1e-12)
    np.testing.assert_allclose(enu[2, 0, :], 0.0, atol=1e-12)


def test_vertical_is_beam_average_at_no_tilt():
    """At zero tilt, Up = (b1+b2+b3+b4)/(4 cos th)."""
    beam = np.array([0.2, -0.1, 0.05, 0.3]).reshape(4, 1, 1)
    enu = beam2enu(beam, np.array([90.0]), np.array([0.0]), np.array([0.0]))
    c = 1.0 / np.cos(np.deg2rad(BEAM_ANGLE_DEG)) / 4.0
    np.testing.assert_allclose(enu[2, 0, 0], c * beam.sum(), atol=1e-12)
