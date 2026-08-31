"""Data-free unit tests for the pressure-based cast detector (rsk.detect_casts)."""
import numpy as np

from ww_rbr.rsk import detect_casts, time_continuity


def _sawtooth(n_profiles=5, fs=8.0, depth=100.0, speed=0.5, dwell_s=30.0):
    """Synthetic Wirewalker pressure record: repeated down/up ramps to `depth`
    (dbar) at `speed` (dbar/s), separated by a surface dwell. Returns (t_ms, P)."""
    dt = 1.0 / fs
    leg = np.arange(0, depth / speed, dt)          # one ramp, seconds
    down = speed * leg
    up = depth - speed * leg
    dwell = np.full(int(dwell_s * fs), 0.0)        # sit at the surface
    P = np.concatenate([np.concatenate([down, up, dwell]) for _ in range(n_profiles)])
    t = (np.arange(P.size) * dt * 1000.0).astype(np.int64)
    return t, P


def test_counts_and_directions_alternate():
    t, P = _sawtooth(n_profiles=5)
    casts = detect_casts(t, P)
    ups = [c for c in casts if c.direction == 1]
    downs = [c for c in casts if c.direction == 0]
    assert len(ups) == 5 and len(downs) == 5
    # casts are time-ordered and non-overlapping
    for a, b in zip(casts, casts[1:]):
        assert a.t2 <= b.t1
    # cast_number is sequential
    assert [c.cast_number for c in casts] == list(range(len(casts)))


def test_upcast_pressure_decreases():
    t, P = _sawtooth(n_profiles=3)
    casts = detect_casts(t, P)
    for c in casts:
        seg = P[(t >= c.t1) & (t < c.t2)]
        if c.direction == 1:      # up: net pressure drop
            assert seg[-1] < seg[0]
        else:                     # down: net pressure rise
            assert seg[-1] > seg[0]


def test_surface_dwell_is_not_a_cast():
    # A flat record (never moves) yields no casts: nothing spans min_span_dbar.
    t = (np.arange(8 * 60) * 125).astype(np.int64)   # 60 s at 8 Hz
    P = np.full(t.size, 10.0)
    assert detect_casts(t, P) == []


def test_rate_derived_from_record():
    # fs is read from the record, not passed in: a 4 Hz record is detected fine.
    t, P = _sawtooth(n_profiles=3, fs=4.0)
    fs, dt_med, gaps = time_continuity(t)
    assert abs(fs - 4.0) < 1e-6 and dt_med == 250.0 and gaps == []
    casts = detect_casts(t, P)                        # no fs argument
    assert sum(c.direction == 1 for c in casts) == 3


def test_data_drop_is_flagged_and_does_not_merge_casts():
    # Splice a multi-minute time gap into the middle of a continuous record.
    t, P = _sawtooth(n_profiles=6, fs=8.0)
    cut = t.size // 2
    t2 = t.copy()
    t2[cut:] += 5 * 60 * 1000                         # shift the tail 5 min later
    fs, dt_med, gaps = time_continuity(t2, gap_factor=4.0)
    assert len(gaps) == 1 and gaps[0][0] == cut - 1   # one drop, just before the cut
    # the gap must not fabricate or swallow casts: still 6 up + 6 down
    casts = detect_casts(t2, P, gap_factor=4.0)
    assert sum(c.direction == 1 for c in casts) == 6
    assert sum(c.direction == 0 for c in casts) == 6
    # and no cast interval straddles the gap
    tgap0, tgap1 = t2[cut - 1], t2[cut]
    assert all(not (c.t1 <= tgap0 and c.t2 > tgap1) for c in casts)


def test_continuous_record_has_no_gaps():
    t, P = _sawtooth(n_profiles=4)
    _, _, gaps = time_continuity(t)
    assert gaps == []


def test_profile_pairs_down_with_up():
    t, P = _sawtooth(n_profiles=4)
    casts = detect_casts(t, P)
    # each down->up pair shares a profile_number
    for c in casts:
        if c.direction == 0:
            partner = next((x for x in casts if x.t1 >= c.t2 and x.direction == 1), None)
            if partner is not None:
                assert partner.profile_number == c.profile_number


def test_empty_input():
    assert detect_casts(np.array([], np.int64), np.array([])) == []
