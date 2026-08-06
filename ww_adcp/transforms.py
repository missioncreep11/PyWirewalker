"""Beam <-> instrument (XYZ) <-> earth (ENU) coordinate transforms.

Python port of the MATLAB toolbox functions
``Beam2XYZ.m`` / ``XYZ2ENU.m`` / ``Beam2ENU.m`` (Bofu Zheng, 2018) for the
Nortek Signature1000 (4-beam Janus, beam angle 25 deg). Vectorized over pings.

Convention notes (matching the MATLAB reference):
- Beam order is 1..4 with beams 1/3 on the instrument X axis and 2/4 on Y.
- ``beam -> XYZ``: X = (b1-b3)/(2 sin th), Y = (b4-b2)/(2 sin th),
  Z = (b1+b2+b3+b4)/(4 cos th)  -- the 4-beam averaged vertical.
- ``XYZ -> ENU``: heading is offset by -90 deg (Nortek ENU convention), then the
  tilt (pitch/roll) matrix is applied; R = H(heading-90) @ P(pitch, roll).
- Angles in **degrees**; velocities in the same units in and out.
"""
from __future__ import annotations

import numpy as np

BEAM_ANGLE_DEG = 25.0


def beam2xyz_matrix(theta_deg: float = BEAM_ANGLE_DEG) -> np.ndarray:
    """3x4 beam->instrument(XYZ) matrix for a 4-beam Signature (MATLAB ``Beam2XYZ``)."""
    t = np.deg2rad(theta_deg)
    a = 1.0 / np.sin(t) / 2.0
    c = 1.0 / np.cos(t) / 4.0
    return np.array([[a, 0.0, -a, 0.0],
                     [0.0, -a, 0.0, a],
                     [c, c, c, c]])


def _tilt_heading_matrix(heading, pitch, roll) -> np.ndarray:
    """Per-ping 3x3 ENU rotation R = H(heading-90) @ P(pitch,roll) (MATLAB ``XYZ2ENU``).

    heading/pitch/roll: array-like (nping,) in degrees. Returns (nping, 3, 3).
    """
    hh = np.deg2rad(np.asarray(heading, float) - 90.0)
    pp = np.deg2rad(np.asarray(pitch, float))
    rr = np.deg2rad(np.asarray(roll, float))
    n = hh.size
    ch, sh = np.cos(hh), np.sin(hh)
    cp, sp = np.cos(pp), np.sin(pp)
    cr, sr = np.cos(rr), np.sin(rr)

    H = np.zeros((n, 3, 3))
    H[:, 0, 0] = ch; H[:, 0, 1] = sh
    H[:, 1, 0] = -sh; H[:, 1, 1] = ch
    H[:, 2, 2] = 1.0

    P = np.zeros((n, 3, 3))
    P[:, 0, 0] = cp; P[:, 0, 1] = -sp * sr; P[:, 0, 2] = -cr * sp
    P[:, 1, 1] = cr; P[:, 1, 2] = -sr
    P[:, 2, 0] = sp; P[:, 2, 1] = sr * cp; P[:, 2, 2] = cp * cr

    return H @ P


def beam2enu(beam, heading, pitch, roll, theta_deg: float = BEAM_ANGLE_DEG) -> np.ndarray:
    """Beam velocities -> ENU (MATLAB ``Beam2ENU``).

    Parameters
    ----------
    beam : (4, ncell, nping) array of the four beam velocities.
    heading, pitch, roll : (nping,) attitude in degrees.

    Returns
    -------
    (3, ncell, nping) ENU velocity: [East, North, Up].
    """
    beam = np.asarray(beam, float)
    T = beam2xyz_matrix(theta_deg)                 # (3, 4)
    R = _tilt_heading_matrix(heading, pitch, roll)  # (nping, 3, 3)
    RT = np.einsum("nij,jk->nik", R, T)            # (nping, 3, 4)
    return np.einsum("nik,kcn->icn", RT, beam)     # (3, ncell, nping)


def xyz2enu(xyz, heading, pitch, roll) -> np.ndarray:
    """Instrument XYZ velocity -> ENU (MATLAB ``XYZ2ENU``).

    xyz : (3, ncell, nping). heading/pitch/roll: (nping,) degrees.
    Returns (3, ncell, nping).
    """
    xyz = np.asarray(xyz, float)
    R = _tilt_heading_matrix(heading, pitch, roll)  # (nping, 3, 3)
    return np.einsum("nij,jcn->icn", R, xyz)
