"""Deployment configuration for the RBR Concerto CTD pipeline.

A `Config` dataclass replaces the old module-level globals: it holds every
deployment- and machine-specific setting plus the derived product paths, and is
passed explicitly to `build_L1/L2/L3`. `load_config` reads the JSON config
(`config.json` beside the repo by default; override with `--config` / `$WW_CONFIG`)
and applies the `$WW_RSK` / `$WW_OUTPUT_DIR` path overrides.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path

_HERE = Path(__file__).resolve().parent


@dataclass
class Config:
    """All settings for one CTD deployment, plus derived product paths."""
    rsk_path: Path
    output_dir: Path
    basename: str
    mooring: str = ""
    instrument: str = ""
    lat: float = 0.0
    lon: float = 0.0
    atm_dbar: float = 10.1325
    fs: float = 2.0
    tm_alpha: float = 0.04
    tm_beta: float = 0.1
    tm_gamma: float = 1.0
    grid_dz: float = 0.5
    grid_zmin: float = 0.0
    grid_zmax: float = 500.0
    l3_dz: float = 1.0
    l3_dt: str = "30min"
    l3_interp_maxgap: int = 1
    n2_smooth_m: float = 5.0
    gravity: float = 9.81

    # --- derived product paths (grid sizes encoded so names track config) ---
    @property
    def l1_path(self) -> Path:
        return self.output_dir / "L1" / f"{self.basename}_L1_converted.nc"

    @property
    def l2_path(self) -> Path:
        return self.output_dir / "L2" / f"{self.basename}_L2_upcast_grid{self.grid_dz:g}m.nc"

    @property
    def _l3name(self) -> str:
        return f"{self.basename}_L3_grid{self.l3_dz:g}m_{self.l3_dt}"

    @property
    def l3_path(self) -> Path:
        return self.output_dir / "L3" / f"{self._l3name}.nc"

    @property
    def l3i_path(self) -> Path:
        return self.output_dir / "L3" / f"{self._l3name}_interp.nc"


def load_config(path=None) -> Config:
    """Load the JSON config into a `Config` (explicit `path`, else `$WW_CONFIG`,
    else `config.json` at the repo root). `$WW_RSK` / `$WW_OUTPUT_DIR` override paths."""
    p = Path(path or os.environ.get("WW_CONFIG") or _HERE.parent / "config.json").expanduser()
    if not p.exists():
        raise FileNotFoundError(
            f"config not found: {p}\nEdit config.json with your deployment paths and metadata.")
    with open(p) as f:
        cfg = json.load(f)
    gr = cfg.get("grid", {})
    tm = cfg.get("thermal_mass", {})
    return Config(
        rsk_path=Path(os.environ.get("WW_RSK", cfg["rsk_file"])).expanduser(),
        output_dir=Path(os.environ.get("WW_OUTPUT_DIR", cfg["output_dir"])).expanduser(),
        basename=cfg["basename"],
        mooring=cfg.get("mooring", ""),
        instrument=cfg.get("instrument", ""),
        lat=cfg["latitude"], lon=cfg["longitude"],
        atm_dbar=cfg["atmospheric_pressure_dbar"],
        fs=cfg.get("sampling_hz", 2.0),
        tm_alpha=tm.get("alpha", 0.04),
        tm_beta=tm.get("beta_per_s", 0.1),
        tm_gamma=tm.get("gamma", 1.0),
        grid_dz=gr.get("l2_dz_m", 0.5),
        grid_zmin=gr.get("zmin_m", 0.0),
        grid_zmax=gr.get("zmax_m", 500.0),
        l3_dz=gr.get("l3_dz_m", 1.0),
        l3_dt=gr.get("l3_dt", "30min"),
        l3_interp_maxgap=gr.get("l3_interp_max_gap_bins", 1),
        n2_smooth_m=cfg.get("n2_vertical_smoothing_m", 5.0),
        gravity=cfg.get("gravity", 9.81),
    )
