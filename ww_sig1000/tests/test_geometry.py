"""Data-free unit tests for ww_sig1000.geometry and get_unit_vectors."""
import numpy as np
from ww_sig1000.transforms import get_unit_vectors
from ww_sig1000.geometry import beam_ranges, cell_depths, BEAM_PHI_DEG


def test_unit_vectors_zero_tilt():
    phi = np.deg2rad(65.0)
    for azi_deg in (0.0, -90.0, 180.0, 90.0):
        azi = np.deg2rad(azi_deg)
        bX, bY, bZ = get_unit_vectors(phi, azi, np.zeros(1), np.zeros(1))
        np.testing.assert_allclose(bZ, np.sin(phi), atol=1e-12)              # up component
        np.testing.assert_allclose(bX, np.cos(azi) * np.cos(phi), atol=1e-12)
        np.testing.assert_allclose(bY, np.sin(azi) * np.cos(phi), atol=1e-12)


def test_beam_ranges_layout():
    # cellsize=0.5, blockdis=0.5 -> numerator = 1.0,1.5,...,22.0 (matches dolfyn range),
    # then divided by cos(25 deg).
    r = beam_ranges(43, 0.5, 0.5, 25.0)
    num = r * np.cos(np.deg2rad(25.0))
    np.testing.assert_allclose(num, 0.5 + 0.5 * np.arange(1, 44), atol=1e-9)


def test_cell_depths_zero_tilt():
    # At zero tilt, range*bZ = (num/cos25)*sin(65) = num  (sin65 == cos25),
    # so z = -pressure + (blockdis + cellsize*(n+1)).
    press = np.array([10.0, 20.0])
    z, ranges, bZ = cell_depths(press, np.zeros(2), np.zeros(2),
                                n_cells=5, cellsize=0.5, blockdis=0.5)
    np.testing.assert_allclose(bZ, np.sin(np.deg2rad(BEAM_PHI_DEG)), atol=1e-12)
    num = 0.5 + 0.5 * np.arange(1, 6)                     # (5,)
    for b in range(4):
        np.testing.assert_allclose(z[0, :, b], -10.0 + num, atol=1e-9)
        np.testing.assert_allclose(z[1, :, b], -20.0 + num, atol=1e-9)


def test_cell_depths_downward_beams_are_kept():
    # A downward-looking instrument (phi negative) has bZ<0 but is valid: cells sit
    # BELOW the transducer, so z = -pressure + range*bZ goes deeper (more negative).
    z, ranges, bZ = cell_depths(np.array([5.0]), np.zeros(1), np.zeros(1),
                                n_cells=3, cellsize=0.5, blockdis=0.5, phi_deg=-65.0)
    assert np.all(np.isfinite(z))
    assert np.all(bZ < 0)
    assert np.all(np.diff(z[0, :, 0]) < 0)          # deeper with range


def test_cell_depths_horizontal_beams_are_nan():
    # Near-horizontal beams (phi ~ 0 -> |bZ| < bz_min) can't reference depth -> NaN.
    z, ranges, bZ = cell_depths(np.array([5.0]), np.zeros(1), np.zeros(1),
                                n_cells=3, cellsize=0.5, blockdis=0.5, phi_deg=0.0)
    assert np.all(np.isnan(z))
