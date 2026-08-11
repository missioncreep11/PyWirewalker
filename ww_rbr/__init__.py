"""Wirewalker RBR Concerto CTD processing (L1 -> L2 -> L3).

Config-driven, channel-aware, level-based NetCDF archive. Usage:

    from ww_rbr import load_config, build_L1, build_L2, build_L3
    cfg = load_config("config.json")
    build_L1(cfg); build_L2(cfg); build_L3(cfg)      # products at cfg.l1_path, ...
"""
from .config import Config, load_config
from .levels import build_L1, build_L2, build_L3

__all__ = ["Config", "load_config", "build_L1", "build_L2", "build_L3"]
