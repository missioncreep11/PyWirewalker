"""Per-ping cell depth geometry for the upward-looking Wirewalker Signature ADCP.

Ported from the depth-coordinate block of ``WWvel_upward.m``: from pressure,
attitude and the cell layout it computes the vertical coordinate of every cell
along every beam, accounting for instrument tilt.

Conventions (matching the MATLAB reference):
- beam elevation ``phi = 65 deg`` from horizontal (= 90 - 25), azimuths
  ``[0, -90, 180, 90] deg``.
- along-beam cell ranges ``(cellsize*(0..N-1) + blockdis + cellsize) / cos(25 deg)``
  (the /cos(25) converts the instrument's vertically-configured cell spacing to
  slant range, which is then projected back onto the vertical by ``bZ``).
- vertical coordinate ``z = -pressure + range * bZ`` (positive up; z = -pressure
  at the transducer). Pings where any beam is not pointing sufficiently upward
  (``bZ <= 0.1``) are set to NaN.
"""
from __future__ import annotations

import numpy as np

from .transforms import get_unit_vectors, BEAM_ANGLE_DEG

BEAM_PHI_DEG = 65.0
BEAM_AZI_DEG = (0.0, -90.0, 180.0, 90.0)


def beam_ranges(n_cells: int, cellsize: float, blockdis: float,
                beam_angle_deg: float = BEAM_ANGLE_DEG) -> np.ndarray:
    """Along-beam range to each cell centre (m), length ``n_cells`` (MATLAB ``ranges``)."""
    n = np.arange(n_cells)
    return (cellsize * n + blockdis + cellsize) / np.cos(np.deg2rad(beam_angle_deg))


def cell_depths(pressure, pitch, roll, n_cells, cellsize, blockdis,
                phi_deg: float = BEAM_PHI_DEG, azi_deg=BEAM_AZI_DEG,
                beam_angle_deg: float = BEAM_ANGLE_DEG, bz_min: float = 0.1):
    """Vertical coordinate of every cell along every beam.

    Parameters
    ----------
    pressure, pitch, roll : (nping,) arrays; pressure in dbar (~m), attitude in degrees.
    n_cells, cellsize, blockdis : cell layout (from the instrument config).

    Returns
    -------
    z : (nping, n_cells, 4) vertical coordinate (m, positive up); NaN where a ping's
        beams are not all pointing up (``bZ <= bz_min``).
    ranges : (n_cells,) along-beam cell ranges.
    bZ : (nping, 4) vertical component of each beam's unit vector.
    """
    pressure = np.asarray(pressure, float)
    pitch = np.asarray(pitch, float)
    roll = np.asarray(roll, float)
    ranges = beam_ranges(n_cells, cellsize, blockdis, beam_angle_deg)  # (n_cells,)
    phi = np.deg2rad(phi_deg)
    azi = np.deg2rad(np.asarray(azi_deg, float))
    pitch_rad = np.deg2rad(pitch)
    roll_rad = np.deg2rad(roll)

    nping = pressure.size
    bZ = np.empty((nping, 4))
    for b in range(4):
        # WWvel_upward calls GetUnitVectors(phi, azi, roll, pitch): roll into the
        # pitch slot and pitch into the roll slot. Reproduce that swap here.
        _, _, bZ[:, b] = get_unit_vectors(phi, azi[b], roll_rad, pitch_rad)

    # Keep pings whose beams all point the same way (all up OR all down) and are
    # sufficiently off-horizontal. z = -pressure + range*bZ then places cells above
    # (bZ>0, upward-looking) or below (bZ<0, downward-looking) the instrument.
    good = np.all(bZ > bz_min, axis=1) | np.all(bZ < -bz_min, axis=1)
    z = -pressure[:, None, None] + ranges[None, :, None] * bZ[:, None, :]  # (nping,n_cells,4)
    z[~good] = np.nan
    return z, ranges, bZ


def look_direction(pitch, roll, phi_deg: float = BEAM_PHI_DEG, azi_deg=BEAM_AZI_DEG) -> str:
    """Infer whether the instrument looks 'up' or 'down' from the mean beam vertical
    component over the given attitude (sign of mean bZ)."""
    phi = np.deg2rad(phi_deg)
    azi = np.deg2rad(np.asarray(azi_deg, float))
    pr, rr = np.deg2rad(np.asarray(roll, float)), np.deg2rad(np.asarray(pitch, float))
    mean_bz = np.nanmean([get_unit_vectors(phi, azi[b], pr, rr)[2] for b in range(4)])
    return "up" if mean_bz > 0 else "down"
