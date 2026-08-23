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
import sys
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
    # --- record trim: drop deployment/recovery transit, when the vehicle is not
    # profiling. Either ensembles or ISO times; times win and are resolved against
    # the dolfyn index by `resolve_trim`. None -> the natural end of the record. ---
    start_ensemble: int = 0
    end_ensemble: Optional[int] = None
    start_time: Optional[str] = None
    end_time: Optional[str] = None
    # --- velocity product ---
    boxsize_m: float = 1.0
    z_max_m: Optional[float] = None            # None -> auto from data
    motion_correct: bool = True
    attitude: str = "ahrs"                     # ahrs | reconstructed | auto (l2.ATTITUDE_MODES)
    motion: str = "v1"                         # v1 | v2 (l2.MOTION_VERSIONS)
    bin_average: str = "boxcar"                # boxcar | notch (l2.BIN_AVERAGE_MODES)
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

    def resolve_trim(self) -> tuple[int, Optional[int]]:
        """The (start, stop) ensemble range to process.

        `start_time` / `end_time` are resolved against the dolfyn index and take
        precedence over the ensemble forms. A stop of None means "to the end of the
        record"; the builders fill that in from the index.
        """
        start, stop = int(self.start_ensemble or 0), self.end_ensemble
        if self.start_time or self.end_time:
            from .index import ensemble_at_time, read_index
            idx = read_index(self.ad2cp_path)      # read once, reuse for both bounds
            if self.start_time:
                start = ensemble_at_time(self.ad2cp_path, self.start_time, index=idx)
            if self.end_time:
                stop = ensemble_at_time(self.ad2cp_path, self.end_time, index=idx)
        return start, (int(stop) if stop is not None else None)


CONFIG_NAME = "config_adcp.json"


class AmbiguousConfigError(RuntimeError):
    """An ambiguous config path was not confirmed, so nothing was loaded."""


def _describe(p: Path) -> str:
    """One-line 'which deployment is this?' summary, for a confirmation prompt."""
    try:
        with open(p) as f:
            cfg = json.load(f)
        return f"{cfg.get('mooring') or '(no mooring)'}  ->  {cfg.get('ad2cp_file') or '(no raw file)'}"
    except Exception:
        return "(could not be read as JSON)"


def _require_agreement(warning: str, assume_yes: bool) -> None:
    """Warn, then continue only on explicit agreement.

    A relative config path resolves against the *process* working directory, so the
    same `config_adcp.json` names a different deployment depending on where the shell
    happens to be — and the run then succeeds against the wrong raw file with no
    error. Rather than forbid it, surface which deployment is about to be loaded and
    make the user say yes. Without a terminal to ask (a script, a queued job) there
    is no way to agree, so `--yes` is required instead of silently proceeding.
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


def _resolve_config(path=None, assume_yes: bool = False) -> Optional[Path]:
    """Locate the config file.

    Explicit `path`, else `$WW_ADCP_CONFIG`, else a `config_adcp.json` in the working
    directory or at the repo root (None if neither exists). A relative explicit path,
    or two different candidates at the default locations, is resolved only after the
    user agrees to the specific file (see `_require_agreement`).
    """
    for raw, source in ((path, "--config"),
                        (None if path else os.environ.get("WW_ADCP_CONFIG"), "$WW_ADCP_CONFIG")):
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
        if not p.exists():
            raise FileNotFoundError(f"ADCP config not found via {source}: {p}")
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
    return found[0] if found else None


def load_adcp_config(path=None, assume_yes: bool = False) -> AdcpConfig:
    """Load the ADCP JSON config into an `AdcpConfig`.

    If no config file is found (and none is required via `path`/`$WW_ADCP_CONFIG`),
    returns an `AdcpConfig` of built-in defaults so the driver can run purely from
    CLI flags. `$WW_AD2CP` / `$WW_OUTPUT_DIR` override the file / output paths.
    """
    p = _resolve_config(path, assume_yes=assume_yes)
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
        start_ensemble=cast.get("start_ensemble", 0) or 0,
        end_ensemble=cast.get("end_ensemble"),
        start_time=cast.get("start_time"),
        end_time=cast.get("end_time"),
        boxsize_m=vel.get("boxsize_m", 1.0),
        z_max_m=vel.get("z_max_m"),
        motion_correct=vel.get("motion_correct", True),
        attitude=vel.get("attitude", "ahrs"),
        motion=vel.get("motion", "v1"),
        bin_average=vel.get("bin_average", "boxcar"),
        dep_res_m=turb.get("dep_res_m", 3.0),
        max_dep_m=turb.get("max_dep_m", 100.0),
        config_path=p,
    )
