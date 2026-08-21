#!/usr/bin/env python3
"""Wirewalker RBR Concerto CTD processing: raw .rsk -> L1/L2/L3 NetCDF archive.

Thin CLI entry point; all logic lives in the `ww_rbr` package. The Ruskin software
already detected profiles (region / regionCast tables), which we reuse.

Levels:
    L1  full-resolution converted time series (measured channels + depth + flags)
    L2  upcasts bin-averaged onto a depth grid; thermal-mass + TEOS-10 here
    L3  regular depth x time grid (+ gap-filled companion)

Usage:
    python3 process_wirewalker_rbr.py --level all
    python3 process_wirewalker_rbr.py --level 1 --max-casts 50      # quick test
    python3 process_wirewalker_rbr.py --level 2 --config /path/to/other.json
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ww_rbr import load_config, build_L1, build_L2, build_L3   # noqa: E402
from ww_rbr.config import AmbiguousConfigError                 # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--level", choices=["1", "2", "3", "all"], required=True)
    ap.add_argument("--config", default=None,
                    help="path to config_ctd.json (default: ./config_ctd.json or $WW_CONFIG)")
    ap.add_argument("--max-casts", type=int, default=None,
                    help="process only the first N casts (for testing)")
    ap.add_argument("-y", "--yes", action="store_true",
                    help="accept an ambiguous (relative) config path without prompting")
    args = ap.parse_args()

    try:
        cfg = load_config(args.config, assume_yes=args.yes)
    except AmbiguousConfigError as e:
        ap.error(str(e))
    print(f"[config] {cfg.mooring} | rsk={cfg.rsk_path} | out={cfg.output_dir}")
    if args.level in ("1", "all"):
        build_L1(cfg, args.max_casts)
    if args.level in ("2", "all"):
        build_L2(cfg, args.max_casts)
    if args.level in ("3", "all"):
        build_L3(cfg)


if __name__ == "__main__":
    main()
