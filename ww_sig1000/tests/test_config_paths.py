"""Config-path resolution must not silently pick another deployment.

A relative `--config config_adcp.json` resolves against the process working
directory. On 2026-08-20 a stale cwd made that load the repo's TLC config instead
of the intended deployment's: the run processed the wrong raw file and reported
success. An ambiguous path is now allowed, but only after the user is shown which
deployment it resolves to and actively agrees.
"""
import json

import pytest

from ww_rbr import config as rbr_config
from ww_sig1000 import config as adcp_config

BOTH = [(adcp_config, "config_adcp.json", "ad2cp_file"),
        (rbr_config, "config_ctd.json", "rsk_file")]


def _write(path, mooring="TEST", raw="/data/raw.bin"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"mooring": mooring, "ad2cp_file": raw, "rsk_file": raw}))
    return path


@pytest.fixture
def tty(monkeypatch):
    """Pretend a terminal is attached, and script the answer to the prompt."""
    def _set(answer):
        monkeypatch.setattr(adcp_config.sys.stdin, "isatty", lambda: True)
        monkeypatch.setattr(rbr_config.sys.stdin, "isatty", lambda: True)
        monkeypatch.setattr("builtins.input", lambda *_: answer)
    return _set


@pytest.mark.parametrize("mod, name, _", BOTH)
def test_relative_path_without_a_terminal_refuses(mod, name, _, tmp_path, monkeypatch):
    """No tty means there is no way to agree, so it must not proceed."""
    monkeypatch.chdir(tmp_path)
    _write(tmp_path / name)
    with pytest.raises(mod.AmbiguousConfigError, match="no terminal"):
        mod._resolve_config(name)


@pytest.mark.parametrize("mod, name, _", BOTH)
def test_relative_path_accepted_when_agreed(mod, name, _, tmp_path, monkeypatch, tty):
    monkeypatch.chdir(tmp_path)
    p = _write(tmp_path / name)
    tty("y")
    assert mod._resolve_config(name) == p.resolve()


@pytest.mark.parametrize("mod, name, _", BOTH)
def test_relative_path_declined_aborts(mod, name, _, tmp_path, monkeypatch, tty):
    monkeypatch.chdir(tmp_path)
    _write(tmp_path / name)
    tty("n")
    with pytest.raises(mod.AmbiguousConfigError, match="Aborted"):
        mod._resolve_config(name)


@pytest.mark.parametrize("mod, name, _", BOTH)
def test_bare_enter_is_not_agreement(mod, name, _, tmp_path, monkeypatch, tty):
    """The prompt is [y/N] — the default must be no."""
    monkeypatch.chdir(tmp_path)
    _write(tmp_path / name)
    tty("")
    with pytest.raises(mod.AmbiguousConfigError, match="Aborted"):
        mod._resolve_config(name)


@pytest.mark.parametrize("mod, name, _", BOTH)
@pytest.mark.parametrize("interrupt", [EOFError, KeyboardInterrupt])
def test_eof_or_ctrl_c_at_the_prompt_aborts_cleanly(mod, name, _, interrupt,
                                                    tmp_path, monkeypatch):
    """Ctrl-D / Ctrl-C must abort, not escape as a traceback."""
    monkeypatch.chdir(tmp_path)
    _write(tmp_path / name)
    monkeypatch.setattr(mod.sys.stdin, "isatty", lambda: True)

    def _boom(*_):
        raise interrupt
    monkeypatch.setattr("builtins.input", _boom)
    with pytest.raises(mod.AmbiguousConfigError, match="Aborted"):
        mod._resolve_config(name)


@pytest.mark.parametrize("mod, name, _", BOTH)
def test_yes_flag_proceeds_without_a_terminal(mod, name, _, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    p = _write(tmp_path / name)
    assert mod._resolve_config(name, assume_yes=True) == p.resolve()


@pytest.mark.parametrize("mod, name, _", BOTH)
def test_absolute_path_needs_no_confirmation(mod, name, _, tmp_path):
    p = _write(tmp_path / name)
    assert mod._resolve_config(str(p)) == p


@pytest.mark.parametrize("mod, env", [(adcp_config, "WW_ADCP_CONFIG"),
                                      (rbr_config, "WW_CONFIG")])
def test_relative_env_var_also_prompts(mod, env, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv(env, "config.json")
    _write(tmp_path / "config.json")
    with pytest.raises(mod.AmbiguousConfigError, match="no terminal"):
        mod._resolve_config()


def test_warning_names_the_deployment_it_would_load(tmp_path, monkeypatch, capsys):
    """The warning must show mooring + raw file — that is what reveals a wrong config."""
    monkeypatch.chdir(tmp_path)
    _write(tmp_path / adcp_config.CONFIG_NAME, mooring="TLC_23", raw="~/TLC/tlc.ad2cp")
    adcp_config._resolve_config(adcp_config.CONFIG_NAME, assume_yes=True)
    err = capsys.readouterr().err
    assert "TLC_23" in err and "tlc.ad2cp" in err and str(tmp_path) in err


def test_two_candidates_require_agreement(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _write(tmp_path / adcp_config.CONFIG_NAME)
    monkeypatch.setattr(adcp_config, "_HERE", tmp_path / "other" / "pkg")
    _write(tmp_path / "other" / adcp_config.CONFIG_NAME)
    with pytest.raises(adcp_config.AmbiguousConfigError, match="no terminal"):
        adcp_config._resolve_config()


def test_two_candidates_choose_the_working_directory_when_agreed(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    p = _write(tmp_path / adcp_config.CONFIG_NAME)
    monkeypatch.setattr(adcp_config, "_HERE", tmp_path / "other" / "pkg")
    _write(tmp_path / "other" / adcp_config.CONFIG_NAME)
    assert adcp_config._resolve_config(assume_yes=True) == p.resolve()


def test_single_candidate_in_cwd_is_used_silently(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    p = _write(tmp_path / adcp_config.CONFIG_NAME)
    monkeypatch.setattr(adcp_config, "_HERE", tmp_path / "nowhere" / "pkg")
    assert adcp_config._resolve_config() == p.resolve()
    assert capsys.readouterr().err == ""


def test_missing_absolute_config_still_raises_filenotfound(tmp_path):
    with pytest.raises(FileNotFoundError):
        adcp_config._resolve_config(str(tmp_path / "absent.json"))


def test_tilde_paths_are_absolute_and_need_no_confirmation(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    p = _write(tmp_path / adcp_config.CONFIG_NAME)
    assert adcp_config._resolve_config(f"~/{adcp_config.CONFIG_NAME}") == p
