"""Deployment configuration for the RBR Concerto CTD pipeline.

A `Config` dataclass replaces the old module-level globals: it holds every
deployment- and machine-specific setting plus the derived product paths, and is
passed explicitly to `build_L1/L2/L3`. `load_config` reads the JSON config
(`config_ctd.json` beside the repo by default; override with `--config` / `$WW_CONFIG`)
and applies the `$WW_RSK` / `$WW_OUTPUT_DIR` path overrides.
"""
from __future__ import annotations

import json
import os
import sys
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


CONFIG_NAME = "config_ctd.json"


class AmbiguousConfigError(RuntimeError):
    """An ambiguous config path was not confirmed, so nothing was loaded."""


def _describe(p: Path) -> str:
    """One-line 'which deployment is this?' summary, for a confirmation prompt."""
    try:
        with open(p) as f:
            cfg = json.load(f)
        return f"{cfg.get('mooring') or '(no mooring)'}  ->  {cfg.get('rsk_file') or '(no raw file)'}"
    except Exception:
        return "(could not be read as JSON)"


def _require_agreement(warning: str, assume_yes: bool) -> None:
    """Warn, then continue only on explicit agreement.

    A relative path is resolved against the *process* working directory, so the same
    `config_ctd.json` names a different deployment depending on where the shell
    happens to be — and the run then succeeds against the wrong raw file with no
    error. (Mirrors `ww_sig1000.config`; the two packages stay independent.)
    """
    print(warning, file=sys.stderr, flush=True)
    if assume_yes:
        print("  Proceeding (--yes).", file=sys.stderr, flush=True)
        return
    if not sys.stdin.isatty():
        raise AmbiguousConfigError(
            "Refusing to guess: no terminal is attached, so this cannot be confirmed.\n"
            "Re-run with an absolute --config, or pass --yes to accept the path above.")
    try:
        answer = input("  Proceed with this config? [y/N] ")
    except (EOFError, KeyboardInterrupt):      # stdin closed, or Ctrl-C/Ctrl-D
        raise AmbiguousConfigError("\nAborted at the config-path confirmation.") from None
    if answer.strip().lower() not in ("y", "yes"):
        raise AmbiguousConfigError("Aborted at the config-path confirmation.")


def _resolve_config(path=None, assume_yes: bool = False) -> Path:
    """Explicit `path`, else `$WW_CONFIG`, else a `config_ctd.json` in the working
    directory or at the repo root. Ambiguous cases need explicit agreement."""
    for raw, source in ((path, "--config"),
                        (None if path else os.environ.get("WW_CONFIG"), "$WW_CONFIG")):
        if not raw:
            continue
        p = Path(raw).expanduser()
        if not p.is_absolute():
            p = (Path.cwd() / p).resolve()
            _require_agreement(
                f"WARNING: {source} {str(raw)!r} is a relative path, resolved against the\n"
                f"  current working directory ({Path.cwd()}).\n"
                f"  If the cwd is not what you expect this loads another deployment.\n"
                f"  It resolves to : {p}\n"
                f"  which contains : {_describe(p) if p.exists() else '(missing)'}",
                assume_yes)
        return p

    here = (Path.cwd() / CONFIG_NAME).resolve()
    repo = (_HERE.parent / CONFIG_NAME).resolve()
    found = [p for p in dict.fromkeys((here, repo)) if p.exists()]
    if len(found) > 1:
        _require_agreement(
            f"WARNING: two different {CONFIG_NAME} files are in play:\n"
            f"  [chosen] {found[0]}\n           {_describe(found[0])}\n"
            f"  [other ] {found[1]}\n           {_describe(found[1])}\n"
            f"  The working-directory one is about to be used.",
            assume_yes)
    return found[0] if found else repo


def load_config(path=None, assume_yes: bool = False) -> Config:
    """Load the JSON config into a `Config` (explicit `path`, else `$WW_CONFIG`,
    else `config_ctd.json` in the cwd or at the repo root). `$WW_RSK` /
    `$WW_OUTPUT_DIR` override paths. Ambiguous paths need explicit agreement."""
    p = _resolve_config(path, assume_yes=assume_yes)
    if not p.exists():
        raise FileNotFoundError(
            f"config not found: {p}\nCopy config_ctd.example.json to config_ctd.json and edit "
            "it with your deployment paths and metadata.")
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
