"""Tests for the ADCP deployment config loader."""
import json
from pathlib import Path

import ww_sig1000
from ww_sig1000.config import AdcpConfig, load_adcp_config


def test_dataclass_defaults():
    # the built-in defaults the loader falls back to when no config file is found,
    # so the driver can still run purely from CLI flags
    cfg = AdcpConfig()
    assert cfg.boxsize_m == 1.0
    assert cfg.dep_res_m == 3.0
    assert cfg.corr_min == 50
    assert cfg.cast_kind is None
    assert cfg.motion_correct is True
    assert cfg.out_path("velocity").name == "adcp_L2_grid1m.nc"


def test_repo_template_config_parses():
    # the tracked config_adcp.example.json template should always load cleanly.
    # (Real per-deployment config_adcp.json files live with their data and are
    # gitignored; only the example is tracked.) Name it explicitly: a bare
    # load_adcp_config() resolves against the working directory, so this would
    # otherwise test whichever deployment pytest ran from.
    template = Path(ww_sig1000.__file__).resolve().parent.parent / "config_adcp.example.json"
    cfg = load_adcp_config(template)
    assert isinstance(cfg, AdcpConfig)
    assert cfg.config_path == template


def test_loads_and_derives_paths(tmp_path, monkeypatch):
    monkeypatch.delenv("WW_AD2CP", raising=False)
    monkeypatch.delenv("WW_OUTPUT_DIR", raising=False)
    p = tmp_path / "config_adcp.json"
    p.write_text(json.dumps({
        "ad2cp_file": "/data/DEPLOY.ad2cp",
        "output_dir": "/out",
        "basename": "DEPLOY_sig1000",
        "mooring": "DEPLOY",
        "instrument": "Nortek Signature1000",
        "velocity": {"boxsize_m": 2.0, "motion_correct": False},
        "turbulence": {"dep_res_m": 3.0, "max_dep_m": 120.0},
        "cast": {"kind": "up", "corr_min": 60, "chunk": 250000},
    }))
    cfg = load_adcp_config(p)
    assert cfg.mooring == "DEPLOY"
    assert cfg.boxsize_m == 2.0 and cfg.motion_correct is False
    assert cfg.max_dep_m == 120.0 and cfg.corr_min == 60 and cfg.chunk == 250000
    assert cfg.cast_kind == "up"
    assert cfg.config_path == p
    # derived, config-tracking output names
    assert cfg.velocity_path.name == "DEPLOY_sig1000_L2_grid2m.nc"
    assert cfg.turbulence_path.name == "DEPLOY_sig1000_turb_dep3m.nc"
    assert cfg.out_path("turbulence") == cfg.turbulence_path


def test_env_overrides_paths(tmp_path, monkeypatch):
    p = tmp_path / "config_adcp.json"
    p.write_text(json.dumps({"ad2cp_file": "/data/orig.ad2cp", "output_dir": "/out",
                             "basename": "b"}))
    monkeypatch.setenv("WW_AD2CP", "/override/new.ad2cp")
    monkeypatch.setenv("WW_OUTPUT_DIR", "/override/out")
    cfg = load_adcp_config(p)
    assert str(cfg.ad2cp_path) == "/override/new.ad2cp"
    assert str(cfg.output_dir) == "/override/out"
