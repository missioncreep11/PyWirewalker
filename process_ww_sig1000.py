#!/usr/bin/env python3
"""Wirewalker Nortek Signature ADCP processing: raw .ad2cp -> gridded product.

Streams a raw .ad2cp (read via mhkit/dolfyn), detects casts, and writes a
(depth x cast) NetCDF. Two products:
  velocity   - motion-corrected ENU currents from the 4 slant beams (L2)
  turbulence - HR beam-5 spectral dissipation eps (ProcessSingleProfile.m method)

Deployment settings come from a JSON config (config_adcp.json by default; see
ww_sig1000/config.py), mirroring the CTD driver. Any CLI flag overrides the config.
Instrument geometry (cell size, blanking, ambiguity velocity, sample rate) is read
from the .ad2cp itself. See ww_sig1000/ and the README.

Examples:
    # config-driven (paths, mooring, params from config_adcp.json)
    python process_ww_sig1000.py --product turbulence
    python process_ww_sig1000.py --product velocity --config /path/to/other.json

    # CLI override of any config value
    python process_ww_sig1000.py --product velocity \
        --file ww_sig1000/test_data/S101913A013_ASTRAL_1_U.ad2cp \
        --out  ww_sig1000/test_data/ASTRAL_1_U_L2.nc --mooring ASTRAL_1_U --boxsize 1.0
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mhkit import dolfyn                       # noqa: E402
from ww_sig1000.config import AmbiguousConfigError, load_adcp_config  # noqa: E402
from ww_sig1000.l2 import build_l2_streaming, save_l2   # noqa: E402
from ww_sig1000.turb_product import build_turbulence_streaming, save_turbulence  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--product", default="velocity", choices=["velocity", "turbulence"],
                    help="which product to build (default velocity)")
    ap.add_argument("--config", default=None,
                    help="path to config_adcp.json (default: ./config_adcp.json or $WW_ADCP_CONFIG)")
    # every flag below overrides the corresponding config value when given
    ap.add_argument("--file", default=None, help="raw .ad2cp path (overrides config ad2cp_file)")
    ap.add_argument("--out", default=None, help="output NetCDF path (default: derived from config)")
    ap.add_argument("--boxsize", type=float, default=None, help="[velocity] depth-bin size (m)")
    ap.add_argument("--z-max", type=float, default=None, help="[velocity] max depth (m); default from data")
    ap.add_argument("--dep-res", type=float, default=None, help="[turbulence] depth resolution (m)")
    ap.add_argument("--max-dep", type=float, default=None, help="[turbulence] max depth (m)")
    ap.add_argument("--chunk", type=int, default=None, help="ensembles per streaming chunk")
    ap.add_argument("--cast-kind", default=None, choices=["both", "up", "down"],
                    help="default both for velocity, up for turbulence")
    ap.add_argument("--min-span-dbar", type=float, default=None, help="min cast pressure span")
    ap.add_argument("--corr-min", type=int, default=None, help="beam correlation threshold")
    ap.add_argument("--no-motion", action="store_true", help="[velocity] disable IMU motion correction")
    ap.add_argument("--mooring", default=None, help="mooring/deployment name for metadata")
    ap.add_argument("--start-ensemble", type=int, default=None,
                    help="first ensemble to process (trims deployment transit)")
    ap.add_argument("--end-ensemble", type=int, default=None,
                    help="stop before this ensemble (default: end of record)")
    ap.add_argument("--start-time", default=None,
                    help="first ensemble at/after this ISO time (overrides --start-ensemble)")
    ap.add_argument("--end-time", default=None,
                    help="stop at this ISO time (overrides --end-ensemble)")
    ap.add_argument("-y", "--yes", action="store_true",
                    help="accept an ambiguous (relative) config path without prompting")
    args = ap.parse_args()

    try:
        cfg = load_adcp_config(args.config, assume_yes=args.yes)
    except AmbiguousConfigError as e:
        ap.error(str(e))
    # overlay CLI overrides (only when explicitly provided)
    if args.file is not None:
        cfg.ad2cp_path = Path(args.file).expanduser()
    if args.mooring is not None:
        cfg.mooring = args.mooring
    if args.boxsize is not None:
        cfg.boxsize_m = args.boxsize
    if args.z_max is not None:
        cfg.z_max_m = args.z_max
    if args.dep_res is not None:
        cfg.dep_res_m = args.dep_res
    if args.max_dep is not None:
        cfg.max_dep_m = args.max_dep
    if args.chunk is not None:
        cfg.chunk = args.chunk
    if args.cast_kind is not None:
        cfg.cast_kind = args.cast_kind
    if args.min_span_dbar is not None:
        cfg.min_span_dbar = args.min_span_dbar
    if args.corr_min is not None:
        cfg.corr_min = args.corr_min
    if args.no_motion:
        cfg.motion_correct = False
    # A CLI bound overrides the config's *other* form of the same bound too —
    # otherwise a `start_time`/`end_time` in the config silently wins over the
    # ensemble flag the user just typed (resolve_trim gives times precedence).
    if args.start_ensemble is not None:
        cfg.start_ensemble, cfg.start_time = args.start_ensemble, None
    if args.end_ensemble is not None:
        cfg.end_ensemble, cfg.end_time = args.end_ensemble, None
    if args.start_time is not None:
        cfg.start_time = args.start_time
    if args.end_time is not None:
        cfg.end_time = args.end_time

    if not str(cfg.ad2cp_path):
        ap.error("no input file: set ad2cp_file in the config or pass --file")
    out = Path(args.out).expanduser() if args.out else cfg.out_path(args.product)
    out.parent.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    print(f"[config] {cfg.mooring or '(no mooring)'} | product={args.product} | "
          f"ad2cp={cfg.ad2cp_path} | out={out}", flush=True)
    if cfg.config_path:
        print(f"[config] loaded {cfg.config_path}", flush=True)

    ens_start, ens_stop = cfg.resolve_trim()
    if ens_start or ens_stop is not None:
        stop_txt = f"{ens_stop:,}" if ens_stop is not None else "end"
        print(f"[config] ensemble range {ens_start:,}:{stop_txt}", flush=True)

    if args.product == "turbulence":
        cast_kind = cfg.cast_kind or "up"
        DS = build_turbulence_streaming(
            str(cfg.ad2cp_path), dolfyn.read, chunk=cfg.chunk, dep_res=cfg.dep_res_m,
            max_dep=cfg.max_dep_m, cast_kind=cast_kind, min_span_dbar=cfg.min_span_dbar,
            corr_min=cfg.corr_min, ens_start=ens_start, total=ens_stop,
            mooring=cfg.mooring, source=cfg.ad2cp_path.name, progress=True)
        _stamp_provenance(DS, cfg)
        save_turbulence(DS, str(out))
        print(f"[done] {DS.sizes['cast']} casts x {DS.sizes['depth']} depths "
              f"(eps spectral) -> {out} in {time.time()-t0:.0f}s")
    else:
        cast_kind = cfg.cast_kind or "both"
        DS = build_l2_streaming(
            str(cfg.ad2cp_path), dolfyn.read, chunk=cfg.chunk, boxsize=cfg.boxsize_m,
            z_max=cfg.z_max_m, cast_kind=cast_kind, min_span_dbar=cfg.min_span_dbar,
            corr_min=cfg.corr_min, ens_start=ens_start, total=ens_stop,
            motion_correct=cfg.motion_correct, mooring=cfg.mooring,
            source=cfg.ad2cp_path.name, progress=True)
        _stamp_provenance(DS, cfg)
        save_l2(DS, str(out))
        print(f"[done] {DS.sizes['cast']} casts x {DS.sizes['depth']} depths "
              f"(look={DS.attrs['instrument_look']}) -> {out} in {time.time()-t0:.0f}s")


def _stamp_provenance(ds, cfg):
    """Record deployment metadata / config source in the product attributes."""
    if cfg.instrument:
        ds.attrs["instrument"] = cfg.instrument
    if cfg.lat is not None:
        ds.attrs["latitude"] = cfg.lat
    if cfg.lon is not None:
        ds.attrs["longitude"] = cfg.lon
    if cfg.config_path:
        ds.attrs["config_file"] = str(cfg.config_path)


if __name__ == "__main__":
    main()
