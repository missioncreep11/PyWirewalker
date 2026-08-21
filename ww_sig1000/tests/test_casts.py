"""Cast detection: burst awareness and complete/truncated flagging."""
import numpy as np

from ww_sig1000.casts import Cast, detect_bursts, detect_casts


def _triangle(n_profiles, half, p_max=500.0):
    """Pressure for `n_profiles` half-cycles (down, up, down, ...) of `half` samples."""
    ramp = np.linspace(0.0, p_max, half)
    legs = [ramp if i % 2 == 0 else ramp[::-1] for i in range(n_profiles)]
    return np.concatenate(legs)


def test_detect_bursts_continuous_is_one_block():
    t = np.arange(1000) * 0.125
    assert detect_bursts(t) == [(0, 1000)]


def test_detect_bursts_splits_on_gaps():
    t = np.concatenate([np.arange(500) * 0.125,
                        6000 + np.arange(500) * 0.125,
                        12000 + np.arange(500) * 0.125])
    assert detect_bursts(t) == [(0, 500), (500, 1000), (1000, 1500)]


def test_continuous_record_unchanged_by_burst_logic():
    """One burst -> the filter sees the whole record, as before."""
    p = _triangle(4, 2000)
    t = np.arange(p.size) * 0.125
    casts = detect_casts(p, t, thhold=200)
    assert len(casts) == 4
    assert [c.direction for c in casts] == ["down", "up", "down", "up"]
    assert all(c.burst == 0 for c in casts)


def test_first_and_last_cast_of_a_record_are_truncated():
    p = _triangle(4, 2000)
    t = np.arange(p.size) * 0.125
    casts = detect_casts(p, t, thhold=200)
    assert casts[0].truncated and casts[-1].truncated
    assert all(c.complete for c in casts[1:-1])


def test_casts_detected_within_each_burst():
    """Two bursts of 3 profiles: detection runs per burst, edges flagged truncated."""
    p1, p2 = _triangle(3, 2000), _triangle(3, 2000)
    p = np.concatenate([p1, p2])
    t1 = np.arange(p1.size) * 0.125
    t = np.concatenate([t1, t1[-1] + 6000 + np.arange(p2.size) * 0.125])

    casts = detect_casts(p, t, thhold=200)
    assert {c.burst for c in casts} == {0, 1}
    # every burst's own first and last cast is clipped by the gap / record edge
    for bi in (0, 1):
        inb = [c for c in casts if c.burst == bi]
        assert inb[0].truncated, "cast at a burst start must be flagged truncated"
        assert inb[-1].truncated, "cast at a burst end must be flagged truncated"
        assert all(c.complete for c in inb[1:-1])
    # indices are offset back into the full record
    assert all(c.stop < p.size for c in casts)
    assert all(p[c.start:c.stop + 1].size == c.n for c in casts)


def test_burst_edges_are_not_orphaned():
    """Per-burst detection reaches the burst boundary; whole-record filtering does not."""
    p1 = _triangle(3, 2000)
    p = np.concatenate([p1, p1])
    t1 = np.arange(p1.size) * 0.125
    t = np.concatenate([t1, t1[-1] + 6000 + np.arange(p1.size) * 0.125])

    casts = detect_casts(p, t, thhold=200)
    burst0 = [c for c in casts if c.burst == 0]
    burst1 = [c for c in casts if c.burst == 1]
    assert burst0[0].start == 0
    assert burst0[-1].stop == p1.size - 1        # right up to the gap
    assert burst1[0].start == p1.size            # resumes on the first sample after it


def test_first_is_continuation_suppresses_the_leading_flag():
    """The streaming readers vouch for a carried cast, so index 0 is not auto-flagged."""
    p = _triangle(4, 2000)
    t = np.arange(p.size) * 0.125
    casts = detect_casts(p, t, thhold=200, first_is_continuation=True)
    assert casts[0].complete
    assert casts[-1].truncated                   # buffer end is still unknown territory


def test_cast_defaults_stay_backward_compatible():
    c = Cast(start=0, stop=9, direction="up")
    assert c.n == 10 and c.complete and not c.truncated and c.burst == 0
