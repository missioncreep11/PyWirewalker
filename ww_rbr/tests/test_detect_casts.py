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


def _short_fast_casts(n_profiles=6, fs=16.0, top=1.4, bottom=18.5,
                      down_speed=0.25, up_speed=1.1, turn_s=2.0):
    """Shallow, fast Wirewalker-like record: slow wave-driven descent, ~15 s buoyant
    ascent, and rounded (cosine) turnarounds of `turn_s` at the apex and nadir, so
    the vehicle decelerates and reverses the way it does at sea. Returns (t_ms, P)."""
    dt = 1.0 / fs
    span = bottom - top
    down = top + down_speed * np.arange(0, span / down_speed, dt)
    up = bottom - up_speed * np.arange(0, span / up_speed, dt)
    tt = np.arange(0, turn_s, dt)
    nadir = bottom + 0.3 * np.sin(np.pi * tt / turn_s)   # overshoot then rise
    apex = top - 0.3 * np.sin(np.pi * tt / turn_s)       # overshoot then sink back
    cyc = np.concatenate([down, nadir, up, apex])
    P = np.concatenate([cyc] * n_profiles)
    t = (np.arange(P.size) * dt * 1000.0).astype(np.int64)
    return t, P


def test_short_casts_end_on_turnarounds():
    # Regression: with the default 5 s slope / 7.5 s debounce windows a ~15 s
    # upcast kept ~1.5 s of post-apex samples (vehicle sinking back through its
    # wake) and pre-nadir samples (still sinking). Upcasts must now span exactly
    # deepest -> shallowest, and be monotone in the sense that no sample is
    # deeper than the first or shallower than the last.
    t, P = _short_fast_casts()
    casts = detect_casts(t, P)
    ups = [c for c in casts if c.direction == 1]
    assert len(ups) == 6
    for c in ups:
        seg = P[(t >= c.t1) & (t < c.t2)]
        assert seg[0] >= seg.max() - 0.04            # starts at the nadir
        assert seg[-1] <= seg.min() + 0.04           # ends at the apex
        # nothing after the apex: the last sample is the shallowest within tol
        i_min = np.argmin(seg)
        assert seg.size - 1 - i_min <= 1 or np.all(seg[i_min:] <= seg.min() + 0.04)
    # casts still tile the record without overlap
    for a, b in zip(casts, casts[1:]):
        assert a.t2 <= b.t1


def test_snapping_does_not_change_long_cast_counts():
    t, P = _sawtooth(n_profiles=5, fs=2.0, depth=400.0, speed=0.3)
    casts = detect_casts(t, P)
    assert sum(c.direction == 1 for c in casts) == 5
    for c in casts:
        if c.direction == 1:
            seg = P[(t >= c.t1) & (t < c.t2)]
            assert seg[0] >= 400.0 - 0.3 and seg[-1] <= 0.3


def test_nan_pressure_does_not_crash_and_keeps_casts():
    # Real .rsk files can hold NULL (-> NaN) pressure samples, singly or in runs;
    # the old detector tolerated them, so the turnaround snapping must too.
    t, P = _short_fast_casts(n_profiles=4)
    P = P.copy()
    P[::97] = np.nan                                  # scattered NaNs
    P[300:340] = np.nan                               # a NaN run inside a cast
    i_apex = int(np.nanargmin(P[:1500]))
    P[i_apex - 5:i_apex + 5] = np.nan                 # NaNs right at a turnaround
    casts = detect_casts(t, P)
    assert sum(c.direction == 1 for c in casts) == 4
    for a, b in zip(casts, casts[1:]):
        assert a.t2 <= b.t1
