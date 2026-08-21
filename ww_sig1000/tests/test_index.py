"""Index decoding and the record-trim config.

The date convention is the fragile part: the `.ad2cp.index` stores the raw Nortek
year byte (year - 1900) and a **1-based** month, which is NOT what the data record
itself uses (`nortek2_lib._calc_time` adds 1 to the data record's month), and is
also not what `dolfyn.time._fullyear` handles. These cases pin the values that were
verified against dolfyn's own decoded timestamps on the NOPP1-California record.
"""
import numpy as np
import pytest
from mhkit.dolfyn.io import nortek2_lib as lib

from ww_sig1000.config import AdcpConfig
from ww_sig1000.index import index_times


def _rows(triples):
    """Build index-shaped records from (year_byte, month, day, h, m, s) tuples."""
    dt = np.dtype(lib._index_dtype[lib._index_version])
    out = np.zeros(len(triples), dtype=dt)
    for i, (y, mo, d, h, mi, s) in enumerate(triples):
        out[i]["year"], out[i]["month"], out[i]["day"] = y, mo, d
        out[i]["hour"], out[i]["minute"], out[i]["second"] = h, mi, s
    return out


def test_index_times_matches_dolfyn_decoded_timestamps():
    """Verified against dolfyn.read on S100601A032_nopp2-001.ad2cp."""
    rows = _rows([(123, 9, 19, 7, 1, 45),      # first ensemble  -> 2023-09-19T07:01:45
                  (123, 12, 8, 21, 13, 40),    # ensemble 15e6   -> 2023-12-08T21:13:40
                  (124, 2, 27, 23, 17, 54)])   # last ensemble   -> 2024-02-27T23:17:54
    got = index_times(rows)
    assert list(got.astype(str)) == ["2023-09-19T07:01:45",
                                     "2023-12-08T21:13:40",
                                     "2024-02-27T23:17:54"]


def test_index_month_is_one_based_not_zero_based():
    """A zero-based reading would make this October; dolfyn says September."""
    assert str(index_times(_rows([(123, 9, 19, 0, 0, 0)]))[0]).startswith("2023-09-19")


def test_index_year_byte_is_offset_from_1900():
    """`dolfyn.time._fullyear` passes 123 through unchanged — it is for 2-digit years."""
    assert str(index_times(_rows([(123, 1, 1, 0, 0, 0)]))[0]).startswith("2023")
    assert str(index_times(_rows([(99, 1, 1, 0, 0, 0)]))[0]).startswith("1999")


def test_index_times_is_monotonic_across_a_year_boundary():
    t = index_times(_rows([(123, 12, 31, 23, 59, 59), (124, 1, 1, 0, 0, 0)]))
    assert t[1] > t[0]


def test_resolve_trim_defaults_to_the_whole_record():
    assert AdcpConfig().resolve_trim() == (0, None)


def test_resolve_trim_passes_through_ensemble_bounds():
    cfg = AdcpConfig(start_ensemble=96_000, end_ensemble=700_000)
    assert cfg.resolve_trim() == (96_000, 700_000)


def test_resolve_trim_gives_times_precedence_over_ensembles():
    """Documents why the driver clears the time form when a CLI ensemble bound is
    given — otherwise a config time silently wins over the flag just typed."""
    cfg = AdcpConfig(start_ensemble=1, end_ensemble=2)
    assert cfg.resolve_trim() == (1, 2)
    cfg.end_time = "2024-10-28T22:30:00"
    with pytest.raises(Exception):        # would consult the (absent) index
        cfg.resolve_trim()


def test_resolve_trim_only_touches_the_index_when_times_are_given():
    """No ad2cp_path is set here, so hitting the index would raise."""
    cfg = AdcpConfig(start_ensemble=5)
    assert cfg.resolve_trim() == (5, None)
    cfg.start_time = "2023-09-19T18:30:00"
    with pytest.raises(Exception):
        cfg.resolve_trim()
