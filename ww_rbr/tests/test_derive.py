"""Data-free unit tests for ww_rbr.derive and ww_rbr.config."""
from pathlib import Path

import numpy as np

from ww_rbr.derive import correct_thermal_mass, convert, buoyancy_n2
from ww_rbr.config import Config


def test_thermal_mass_constant_temp_is_identity():
    # No temperature change -> no thermal-mass correction -> conductivity unchanged.
    cond = np.array([3.1, 3.2, 3.15, 3.0, 3.05])
    temp = np.full(5, 12.0)
    out = correct_thermal_mass(cond, temp, alpha=0.04, beta=0.1, gamma=1.0, fs=6.0)
    np.testing.assert_allclose(out, cond, atol=1e-12)


def test_thermal_mass_responds_to_temp_step():
    # A temperature step produces a nonzero, finite correction.
    cond = np.full(20, 3.0)
    temp = np.concatenate([np.full(10, 10.0), np.full(10, 12.0)])
    out = correct_thermal_mass(cond, temp, alpha=0.04, beta=0.1, gamma=1.0, fs=6.0)
    assert np.all(np.isfinite(out))
    assert not np.allclose(out, cond)          # correction actually applied
    assert np.any(np.abs(out - cond)[10:] > 0)  # kicks in at the step


def test_buoyancy_n2_linear_stable():
    dz = 1.0
    z = np.arange(1, 21) * dz - dz / 2         # positive down
    a, b = 25.0, 0.02                          # sigma0 increases with depth (stable)
    sigma0 = (a + b * z)[:, None]
    n2 = buoyancy_n2(sigma0, z, Lm=dz, gravity=9.81)   # Lm=dz -> no smoothing
    expected = 9.81 * b / (1000.0 + sigma0[:, 0])
    np.testing.assert_allclose(n2[:, 0], expected, rtol=1e-6)
    assert np.all(n2 > 0)


def test_buoyancy_n2_unstable_is_negative():
    dz = 1.0
    z = np.arange(1, 11) * dz - dz / 2
    sigma0 = (25.0 - 0.02 * z)[:, None]        # density decreasing downward -> unstable
    n2 = buoyancy_n2(sigma0, z, Lm=dz, gravity=9.81)
    assert np.all(n2 < 0)


def test_convert_sane_seawater():
    cfg = Config(rsk_path=Path("x.rsk"), output_dir=Path("/out"), basename="D",
                 lat=0.0, lon=0.0, atm_dbar=10.1325)
    cond = np.array([50.0, 48.0])              # mS/cm
    temp = np.array([25.0, 20.0])              # degC
    pres = np.array([10.1325 + 5.0, 10.1325 + 50.0])  # total dbar -> sea_p 5, 50
    d = convert(cond, temp, pres, cfg)
    assert set(d) >= {"practical_salinity", "absolute_salinity",
                      "conservative_temperature", "sigma0", "sound_speed", "depth"}
    assert np.all((d["practical_salinity"] > 20) & (d["practical_salinity"] < 40))
    assert np.all(d["depth"] > 0)              # positive down
    np.testing.assert_allclose(d["sea_pressure"], [5.0, 50.0], atol=1e-6)


def test_config_derived_paths():
    cfg = Config(rsk_path=Path("x.rsk"), output_dir=Path("/out"), basename="DEP",
                 grid_dz=0.5, l3_dz=1.0, l3_dt="30min")
    assert cfg.l1_path == Path("/out/L1/DEP_L1_converted.nc")
    assert cfg.l2_path == Path("/out/L2/DEP_L2_upcast_grid0.5m.nc")
    assert cfg.l3_path == Path("/out/L3/DEP_L3_grid1m_30min.nc")
