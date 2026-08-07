"""Split a Wirewalker pressure record into up/down casts.

Port of the cast-detection logic in ``get_aqd_2G.m`` / ``create_profiles.m``:
low-pass the pressure, find turning points (sign changes of its slope), segment
into casts, drop segments shorter than a threshold, cut on large time gaps (duty
cycling), and label each segment up (deep->shallow, buoyant rise) or down.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import butter, filtfilt


@dataclass
class Cast:
    start: int          # first ping index (inclusive)
    stop: int           # last ping index (inclusive)
    direction: str      # 'up' or 'down'

    @property
    def n(self) -> int:
        return self.stop - self.start + 1


def detect_casts(pressure, time_s, thhold: int = 20, lp_period_samples: float = 200.0,
                 gap_s: float = 30.0, order: int = 3) -> list[Cast]:
    """Segment a pressure record into casts.

    Parameters
    ----------
    pressure : (n,) dbar.
    time_s : (n,) time in **seconds** (monotonic).
    thhold : minimum samples for a segment to count as a cast.
    lp_period_samples : Butterworth low-pass cutoff period, in samples
        (matches MATLAB ``fc = 1/(200*dt)``).
    gap_s : start a new cast when the time step exceeds this (duty-cycle gaps).

    Returns a list of :class:`Cast` in time order.
    """
    p = np.asarray(pressure, float)
    t = np.asarray(time_s, float)
    n = p.size
    if n < max(thhold, 3 * order + 1):
        return []

    dt = float(np.median(np.diff(t)))
    fnb = 1.0 / (2.0 * dt)                       # Nyquist
    fc = 1.0 / (lp_period_samples * dt)          # cutoff
    b, a = butter(order, min(fc / fnb, 0.999), "low")
    pf = filtfilt(b, a, p)

    # turning points: where consecutive slopes have opposite (or zero) sign
    d = np.diff(pf)
    turn = np.flatnonzero(d[:-1] * d[1:] <= 0)   # index i is a local extremum of pf

    # segment boundaries: extrema (+1) and large time-gap cuts
    dtd = np.diff(t, prepend=t[0])
    gaps = np.flatnonzero(dtd > gap_s)
    starts = np.unique(np.concatenate(([0], turn + 1, gaps)))
    # build [start, stop] pairs from consecutive starts
    stops = np.concatenate((starts[1:] - 1, [n - 1]))

    casts: list[Cast] = []
    for s, e in zip(starts, stops):
        if e - s + 1 < thhold:
            continue
        seg = p[s:e + 1]
        # up = deep->shallow: the max occurs before the min (pressure falling)
        i_min = int(np.argmin(seg))
        i_max = int(np.argmax(seg))
        direction = "down" if i_max > i_min else "up"
        casts.append(Cast(start=int(s), stop=int(e), direction=direction))
    return casts


def upcasts(casts: list[Cast]) -> list[Cast]:
    """Filter to upcasts (buoyant rise) — the segments processed into L2 by default."""
    return [c for c in casts if c.direction == "up"]
