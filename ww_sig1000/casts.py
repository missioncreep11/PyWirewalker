"""Split a Wirewalker pressure record into up/down casts.

Port of the cast-detection logic in ``get_aqd_2G.m`` / ``create_profiles.m``:
low-pass the pressure, find turning points (sign changes of its slope), segment
into casts, drop segments shorter than a threshold, and label each segment up
(deep->shallow, buoyant rise) or down.

Burst (duty-cycled) sampling
----------------------------
Many deployments duty-cycle: a burst of continuous sampling, then a long gap.
``detect_bursts`` finds those contiguous blocks from the time base, and
``detect_casts`` runs the low-pass and turning-point detection **within each
burst independently**. Filtering across a multi-hour gap smears the pressure
discontinuity into the burst edges and orphans samples there; on a continuously
sampled record there is exactly one burst, so this reduces to the original
whole-record behaviour.

Complete vs. truncated casts
----------------------------
A cast that runs into a burst boundary is clipped by the duty cycle and covers
only part of the profile. Those are flagged ``truncated`` (``complete`` is its
inverse) so downstream gridding can weight or exclude them; they are still valid
velocity data over the depth range they do cover.
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
    truncated: bool = False   # clipped by a burst boundary or the record edge
    burst: int = 0            # index of the sampling burst this cast came from

    @property
    def n(self) -> int:
        return self.stop - self.start + 1

    @property
    def complete(self) -> bool:
        """True when the cast is bounded by turning points, not by a data gap."""
        return not self.truncated


def detect_bursts(time_s, gap_s: float = 30.0) -> list[tuple[int, int]]:
    """Contiguous sampling blocks in a time base, as ``[(start, stop_exclusive), ...]``.

    A gap longer than `gap_s` starts a new burst. A continuously sampled record
    returns a single block spanning everything.
    """
    t = np.asarray(time_s, float)
    if t.size == 0:
        return []
    cuts = np.flatnonzero(np.diff(t) > gap_s) + 1
    edges = np.concatenate(([0], cuts, [t.size]))
    return [(int(a), int(b)) for a, b in zip(edges[:-1], edges[1:])]


def _casts_in_block(p, t, thhold: int, lp_period_samples: float, order: int) -> list[Cast]:
    """Turning-point cast detection within one contiguous (gap-free) block."""
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

    starts = np.unique(np.concatenate(([0], turn + 1)))
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


def detect_casts(pressure, time_s, thhold: int = 20, lp_period_samples: float = 200.0,
                 gap_s: float = 30.0, order: int = 3,
                 first_is_continuation: bool = False) -> list[Cast]:
    """Segment a pressure record into casts, burst by burst.

    Parameters
    ----------
    pressure : (n,) dbar.
    time_s : (n,) time in **seconds** (monotonic).
    thhold : minimum samples for a segment to count as a cast.
    lp_period_samples : Butterworth low-pass cutoff period, in samples
        (matches MATLAB ``fc = 1/(200*dt)``).
    gap_s : a time step longer than this starts a new burst (duty cycling).
    first_is_continuation : the caller knows sample 0 continues a cast already
        assessed upstream (the streaming readers carry a boundary cast between
        chunks), so a cast starting at index 0 is not flagged truncated on that
        account alone.

    Returns a list of :class:`Cast` in time order, each flagged ``truncated`` when
    a burst boundary or the buffer edge clips it.
    """
    p = np.asarray(pressure, float)
    t = np.asarray(time_s, float)
    bursts = detect_bursts(t, gap_s)
    n_bursts = len(bursts)

    casts: list[Cast] = []
    for bi, (a, b) in enumerate(bursts):
        for c in _casts_in_block(p[a:b], t[a:b], thhold, lp_period_samples, order):
            c.start += a
            c.stop += a
            c.burst = bi
            # clipped at the start: a real gap precedes this burst, or we are at
            # the buffer edge and the caller has not vouched for it
            at_block_start = c.start == a
            at_block_stop = c.stop == b - 1
            clipped_start = at_block_start and (bi > 0 or not first_is_continuation)
            # clipped at the end: a real gap follows, or we ran out of buffer
            clipped_stop = at_block_stop and (bi < n_bursts - 1 or b == t.size)
            c.truncated = bool(clipped_start or clipped_stop)
            casts.append(c)
    return casts


def upcasts(casts: list[Cast]) -> list[Cast]:
    """Filter to upcasts (buoyant rise) — the segments processed into L2 by default."""
    return [c for c in casts if c.direction == "up"]


def complete_casts(casts: list[Cast]) -> list[Cast]:
    """Filter to casts not clipped by a burst boundary or the record edge."""
    return [c for c in casts if c.complete]
