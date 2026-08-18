#!/usr/bin/env python3
"""Wirewalker Nortek Signature ADCP processing: raw .ad2cp -> gridded product.

Streams a raw .ad2cp (read via mhkit/dolfyn), detects casts, and writes a
(depth x cast) NetCDF. Two products:
  velocity   - motion-corrected ENU currents from the 4 slant beams (L2)
  turbulence - HR beam-5 spectral dissipation eps (ProcessSingleProfile.m method)
Instrument geometry is read from the file. See ww_sig1000/ and WW_Velocity_Processing_SWOT/.

Examples:
    python process_ww_sig1000.py --product velocity \
        --file ww_sig1000/test_data/S101913A013_ASTRAL_1_U.ad2cp \
        --out  ww_sig1000/test_data/ASTRAL_1_U_L2.nc --mooring ASTRAL_1_U --boxsize 1.0

    python process_ww_sig1000.py --product turbulence \
        --file /Users/drew/TLC/2023/Wirewalker/Nortek/S101913A008_TLC_23.ad2cp \
        --out  /Users/drew/TLC/2023/Wirewalker/Nortek/processed/TLC_23_sig1000_turb.nc \
        --mooring TLC_23
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mhkit import dolfyn                       # noqa: E402
from ww_sig1000.l2 import build_l2_streaming, save_l2   # noqa: E402
from ww_sig1000.turb_product import build_turbulence_streaming, save_turbulence  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--product", default="velocity", choices=["velocity", "turbulence"],
                    help="which product to build (default velocity)")
    ap.add_argument("--file", required=True, help="raw .ad2cp path")
    ap.add_argument("--out", required=True, help="output NetCDF path")
    ap.add_argument("--boxsize", type=float, default=1.0, help="[velocity] depth-bin size (m)")
    ap.add_argument("--z-max", type=float, default=None, help="[velocity] max depth (m); default from data")
    ap.add_argument("--dep-res", type=float, default=3.0, help="[turbulence] depth resolution (m)")
    ap.add_argument("--max-dep", type=float, default=100.0, help="[turbulence] max depth (m)")
    ap.add_argument("--chunk", type=int, default=500_000, help="ensembles per streaming chunk")
    ap.add_argument("--cast-kind", default=None, choices=["both", "up", "down"],
                    help="default both for velocity, up for turbulence")
    ap.add_argument("--min-span-dbar", type=float, default=40.0, help="min cast pressure span")
    ap.add_argument("--corr-min", type=int, default=50, help="beam correlation threshold")
    ap.add_argument("--no-motion", action="store_true", help="[velocity] disable IMU motion correction")
    ap.add_argument("--mooring", default="", help="mooring/deployment name for metadata")
    args = ap.parse_args()

    t0 = time.time()
    print(f"[ww_sig1000] streaming {args.product} from {args.file}")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    if args.product == "turbulence":
        cast_kind = args.cast_kind or "up"
        DS = build_turbulence_streaming(
            args.file, dolfyn.read, chunk=args.chunk, dep_res=args.dep_res, max_dep=args.max_dep,
            cast_kind=cast_kind, min_span_dbar=args.min_span_dbar, corr_min=args.corr_min,
            mooring=args.mooring, source=Path(args.file).name, progress=True)
        save_turbulence(DS, args.out)
        print(f"[done] {DS.sizes['cast']} casts x {DS.sizes['depth']} depths "
              f"(eps spectral) -> {args.out} in {time.time()-t0:.0f}s")
    else:
        cast_kind = args.cast_kind or "both"
        DS = build_l2_streaming(
            args.file, dolfyn.read, chunk=args.chunk, boxsize=args.boxsize, z_max=args.z_max,
            cast_kind=cast_kind, min_span_dbar=args.min_span_dbar, corr_min=args.corr_min,
            motion_correct=not args.no_motion, mooring=args.mooring, source=Path(args.file).name,
            progress=True)
        save_l2(DS, args.out)
        print(f"[done] {DS.sizes['cast']} casts x {DS.sizes['depth']} depths "
              f"(look={DS.attrs['instrument_look']}) -> {args.out} in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
