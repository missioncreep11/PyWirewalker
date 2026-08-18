"""Deployment configuration for the Nortek Signature1000 ADCP pipeline.

Mirrors the CTD side (`ww_rbr/config.py`): an `AdcpConfig` dataclass holds every
deployment- and machine-specific setting plus derived product paths, and
`load_adcp_config` reads a JSON config (`config_adcp.json` at the repo root by
default; override with `--config` / `$WW_ADCP_CONFIG`), applying the `$WW_AD2CP` /
`$WW_OUTPUT_DIR` path overrides.

Unlike the CTD, the instrument *geometry* (cell size, blanking, ambiguity velocity,
sample rate) is read straight from the `.ad2cp` at run time, so it is deliberately
absent here — the config carries only paths, metadata, and processing choices. The
driver overlays any CLI flags on top of the loaded config, and the resolved values
are written into the output NetCDF attributes for provenance.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

_HERE = Path(__file__).resolve().parent


@dataclass
class AdcpConfig:
    """All settings for one Signature1000 deployment, plus derived product paths."""
    ad2cp_path: Path = Path("")
    output_dir: Path = Path(".")
    basename: str = "adcp"
    mooring: str = ""
    instrument: str = ""
    lat: Optional[float] = None
    lon: Optional[float] = None
    # --- cast detection / QC (shared by both products) ---
    cast_kind: Optional[str] = None            # None -> per-product default (vel: both, turb: up)
    min_span_dbar: float = 40.0
    corr_min: int = 50
    chunk: int = 500_000
    # --- velocity product ---
    boxsize_m: float = 1.0
    z_max_m: Optional[float] = None            # None -> auto from data
    motion_correct: bool = True
    # --- turbulence product ---
    dep_res_m: float = 3.0
    max_dep_m: float = 100.0
    # --- provenance (set by the loader; not from the JSON) ---
    config_path: Optional[Path] = field(default=None)

    @property
    def velocity_path(self) -> Path:
        return self.output_dir / f"{self.basename}_L2_grid{self.boxsize_m:g}m.nc"

    @property
    def turbulence_path(self) -> Path:
        return self.output_dir / f"{self.basename}_turb_dep{self.dep_res_m:g}m.nc"

    def out_path(self, product: str) -> Path:
        return self.turbulence_path if product == "turbulence" else self.velocity_path


def _resolve_config(path=None) -> Optional[Path]:
    """Locate the config file: explicit `path`, else `$WW_ADCP_CONFIG`, else the
    repo-root `config_adcp.json` if it exists (None if none is found)."""
    if path:
        p = Path(path).expanduser()
        if not p.exists():
            raise FileNotFoundError(f"ADCP config not found: {p}")
        return p
    env = os.environ.get("WW_ADCP_CONFIG")
    if env:
        p = Path(env).expanduser()
        if not p.exists():
            raise FileNotFoundError(f"$WW_ADCP_CONFIG points at a missing file: {p}")
        return p
    default = _HERE.parent / "config_adcp.json"
    return default if default.exists() else None


def load_adcp_config(path=None) -> AdcpConfig:
    """Load the ADCP JSON config into an `AdcpConfig`.

    If no config file is found (and none is required via `path`/`$WW_ADCP_CONFIG`),
    returns an `AdcpConfig` of built-in defaults so the driver can run purely from
    CLI flags. `$WW_AD2CP` / `$WW_OUTPUT_DIR` override the file / output paths.
    """
    p = _resolve_config(path)
    if p is None:
        return AdcpConfig()
    with open(p) as f:
        cfg = json.load(f)
    vel = cfg.get("velocity", {})
    turb = cfg.get("turbulence", {})
    cast = cfg.get("cast", {})
    ad2cp = os.environ.get("WW_AD2CP", cfg.get("ad2cp_file", ""))
    outdir = os.environ.get("WW_OUTPUT_DIR", cfg.get("output_dir", "."))
    return AdcpConfig(
        ad2cp_path=Path(ad2cp).expanduser(),
        output_dir=Path(outdir).expanduser(),
        basename=cfg.get("basename", "adcp"),
        mooring=cfg.get("mooring", ""),
        instrument=cfg.get("instrument", ""),
        lat=cfg.get("latitude"),
        lon=cfg.get("longitude"),
        cast_kind=cast.get("kind"),
        min_span_dbar=cast.get("min_span_dbar", 40.0),
        corr_min=cast.get("corr_min", 50),
        chunk=cast.get("chunk", 500_000),
        boxsize_m=vel.get("boxsize_m", 1.0),
        z_max_m=vel.get("z_max_m"),
        motion_correct=vel.get("motion_correct", True),
        dep_res_m=turb.get("dep_res_m", 3.0),
        max_dep_m=turb.get("max_dep_m", 100.0),
        config_path=p,
    )
