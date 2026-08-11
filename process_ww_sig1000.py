#!/usr/bin/env python3
"""Wirewalker Nortek Signature ADCP processing: raw .ad2cp -> gridded L2 velocity.

Streams a raw .ad2cp (read via mhkit/dolfyn), detects up/down casts, grids each
cast's motion-corrected ENU velocity onto a uniform depth grid, and writes an L2
(depth x cast) NetCDF. Instrument geometry (cell size, blanking, beam angle,
up/down look) is read from the file. See ww_sig1000/ and WW_Velocity_Processing_SWOT/.

Example:
    python process_ww_sig1000.py \
        --file ww_sig1000/test_data/S101913A013_ASTRAL_1_U.ad2cp \
        --out  ww_sig1000/test_data/ASTRAL_1_U_L2.nc \
        --mooring ASTRAL_1_U --boxsize 1.0
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mhkit import dolfyn                       # noqa: E402
from ww_sig1000.l2 import build_l2_streaming, save_l2   # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--file", required=True, help="raw .ad2cp path")
    ap.add_argument("--out", required=True, help="output L2 NetCDF path")
    ap.add_argument("--boxsize", type=float, default=1.0, help="depth-bin size (m)")
    ap.add_argument("--z-max", type=float, default=None, help="max depth (m); default from data")
    ap.add_argument("--chunk", type=int, default=500_000, help="ensembles per streaming chunk")
    ap.add_argument("--cast-kind", default="both", choices=["both", "up", "down"])
    ap.add_argument("--min-span-dbar", type=float, default=40.0, help="min cast pressure span")
    ap.add_argument("--corr-min", type=int, default=50, help="beam correlation threshold")
    ap.add_argument("--no-motion", action="store_true", help="disable IMU motion correction")
    ap.add_argument("--mooring", default="", help="mooring/deployment name for metadata")
    args = ap.parse_args()

    t0 = time.time()
    print(f"[ww_sig1000] streaming {args.file}")
    L2 = build_l2_streaming(
        args.file, dolfyn.read, chunk=args.chunk, boxsize=args.boxsize, z_max=args.z_max,
        cast_kind=args.cast_kind, min_span_dbar=args.min_span_dbar, corr_min=args.corr_min,
        motion_correct=not args.no_motion, mooring=args.mooring, source=Path(args.file).name,
        progress=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    save_l2(L2, args.out)
    print(f"[done] {L2.sizes['cast']} casts x {L2.sizes['depth']} depths "
          f"(look={L2.attrs['instrument_look']}) -> {args.out} in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
