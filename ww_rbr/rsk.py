"""RBR `.rsk` (SQLite) reading, channel discovery, and NetCDF metadata.

The mapping of stored `data` columns (channel01, channel02, ...) to physical
channels is read from the rsk `channels` + `instrumentChannels` tables at runtime
(`discover_channels`), so any extra sensors (turbidity, chlorophyll, CDOM, DO,
PAR, irradiance, extra thermistors, ...) are picked up automatically.
"""
from __future__ import annotations

import re
import sqlite3
import time as _time
from dataclasses import dataclass

import gsw
import numpy as np


def _slug(text: str) -> str:
    """Filesystem/NetCDF-safe lower_snake_case name from a channel long name."""
    s = re.sub(r"[^0-9a-zA-Z]+", "_", (text or "").strip().lower()).strip("_")
    return s or "channel"


def discover_channels(con: sqlite3.Connection) -> list[dict]:
    """Discover the measured channels actually stored in the `data` table.

    Joins `channels` to `instrumentChannels` to map each measured channel to its
    `channelNN` column (NN = channelOrder), keeping only channels whose column is
    physically present in `data`. Conductivity, the C-T cell thermistor (temp14)
    and the measured pressure get canonical names (`conductivity`/`temperature`/
    `pressure`); every other measured channel is carried through under a name
    slugged from its long name (deduped). Returns dicts {col, var, long_name, units, short}.
    """
    datacols = {r[1] for r in con.execute("PRAGMA table_info(data)").fetchall()}
    rows = con.execute(
        "SELECT ic.channelOrder, c.shortName, c.longNamePlainText, c.unitsPlainText "
        "FROM channels c JOIN instrumentChannels ic ON ic.channelID = c.channelID "
        "WHERE c.isMeasured = 1 ORDER BY ic.channelOrder"
    ).fetchall()

    chans: list[dict] = []
    used: set[str] = set()
    have = {"conductivity": False, "temperature": False, "pressure": False}
    for order, short, longname, units in rows:
        col = f"channel{int(order):02d}"
        if col not in datacols:
            continue  # measured channel not stored in this deployment's data table
        short = (short or "").lower()
        if short.startswith("cond") and not have["conductivity"]:
            var = "conductivity"; have["conductivity"] = True
        elif short == "temp14" and not have["temperature"]:
            var = "temperature"; have["temperature"] = True
        elif short.startswith("pres") and not have["pressure"]:
            var = "pressure"; have["pressure"] = True
        else:
            base = _slug(longname or short); var = base; k = 2
            while var in used or var in have:
                var = f"{base}_{k}"; k += 1
        used.add(var)
        chans.append({"col": col, "var": var, "long_name": longname or short,
                      "units": units or "", "short": short})

    if not have["temperature"]:
        for ch in chans:
            if ch["short"].startswith("temp"):
                ch["var"] = "temperature"; have["temperature"] = True
                break

    missing = [k for k, v in have.items() if not v]
    if missing:
        raise ValueError(f"rsk is missing required CTD channel(s) {missing}; "
                         f"found: {[c['short'] for c in chans]}")
    return chans


@dataclass
class Cast:
    cast_number: int       # sequential 0..N-1 over all casts
    profile_number: int    # Ruskin profile index (regionProfileID)
    direction: int         # 0 down, 1 up
    t1: int                # start tstamp (ms)
    t2: int                # end tstamp (ms)


def load_casts(con: sqlite3.Connection) -> list[Cast]:
    """Read cast intervals from the Ruskin region tables, ordered in time."""
    q = """
        SELECT r.tstamp1, r.tstamp2, rc.type, rc.regionProfileID
        FROM regionCast rc
        JOIN region r ON r.regionID = rc.regionID
        ORDER BY r.tstamp1 ASC
    """
    rows = con.execute(q).fetchall()
    casts: list[Cast] = []
    for i, (t1, t2, ctype, prof) in enumerate(rows):
        direction = 1 if ctype.upper() == "UP" else 0
        casts.append(Cast(cast_number=i, profile_number=int(prof),
                           direction=direction, t1=int(t1), t2=int(t2)))
    return casts


def time_continuity(tstamp: np.ndarray, gap_factor: float = 4.0):
    """Check the record for continuity and derive its steady sampling rate.

    An RBR cannot skip samples on its own, but a real-time telemetered file can
    drop stretches of data, so the record is not guaranteed contiguous in time.
    The sampling *rate* is steady within a continuous record, so we take the
    nominal interval to be the median sample spacing and flag every interval
    longer than ``gap_factor`` times that median as a data drop.

    Returns ``(fs_hz, dt_median_ms, gaps)`` where ``fs_hz = 1000 / dt_median_ms``
    and ``gaps`` is a list of ``(index_before, dt_ms)`` — one per drop, at the
    sample just before the gap. An empty ``gaps`` means the record is continuous.
    """
    t = np.asarray(tstamp, np.int64)
    if t.size < 2:
        return float("nan"), float("nan"), []
    dt = np.diff(t).astype(float)                 # ms between consecutive samples
    dt_med = float(np.median(dt))
    fs = 1000.0 / dt_med if dt_med > 0 else float("nan")
    gaps = []
    if dt_med > 0 and gap_factor:
        for i in np.flatnonzero(dt > gap_factor * dt_med):
            gaps.append((int(i), float(dt[i])))
    return fs, dt_med, gaps


def _detect_casts_segment(t, P, half, m, thresh_per_sample, min_span_dbar,
                          cn0, prof0):
    """Detect casts within a single time-continuous segment (no gaps inside).

    Works in sample space at the segment's steady rate: `half`/`m` are sample
    counts and `thresh_per_sample` is dbar/sample. Returns (casts, next_cn,
    next_prof) so the caller can number casts/profiles across segments.
    """
    n = P.size
    slope = np.zeros(n)                            # dbar/sample; edges stay 0 (-> down)
    fl = 2 * half
    if n > fl:
        slope[half:n - half] = (P[fl:] - P[:n - fl]) / fl
    up = (slope < -thresh_per_sample).astype(np.int8)

    csum = np.concatenate(([0], np.cumsum(up, dtype=np.int64)))
    idx = np.arange(n)
    lo = np.clip(idx - m, 0, n)
    hi = np.clip(idx + m + 1, 0, n)
    mean_up = (csum[hi] - csum[lo]) / (hi - lo)
    upf = (mean_up >= 0.5).astype(np.int8)         # majority vote (round-half-up)

    chg = np.flatnonzero(np.diff(upf)) + 1
    starts = np.concatenate(([0], chg))
    ends = np.concatenate((chg, [n]))

    casts, cn, prof = [], cn0, prof0
    for s, e in zip(starts, ends):
        seg = P[s:e]
        if seg.max() - seg.min() < min_span_dbar:
            continue                               # not a real cast
        direction = int(upf[s])                    # 1 up, 0 down
        if direction == 0:
            prof += 1                              # a downcast opens a new profile
        casts.append(Cast(cast_number=cn, profile_number=max(prof, 0),
                          direction=direction, t1=int(t[s]), t2=int(t[e - 1]) + 1))
        cn += 1
    return casts, cn, prof


def detect_casts(tstamp: np.ndarray, pressure: np.ndarray, fs: float | None = None, *,
                 slope_window_s: float = 5.0,
                 debounce_window_s: float = 7.5,
                 min_slope_dbar_per_s: float = 0.04,
                 min_span_dbar: float = 5.0,
                 gap_factor: float = 4.0) -> list[Cast]:
    """Detect up/down casts directly from the CTD pressure record.

    A time-aware, vectorised port of the historical MATLAB ``get_upcastRBR.m``.
    Rather than trusting the Ruskin region tables, it classifies every sample as
    rising (up) or sinking (down) from the sign of the local pressure slope,
    debounces the flag, then splits it into contiguous cast intervals.

    Everything is specified in physical units — profiling **speed** (dbar/s) and
    **time** (s) — and the sampling rate is read from the record itself, so the
    same config applies to instruments logging at different rates:

    - **Continuity first.** The steady rate ``fs`` is the median sample spacing
      (see `time_continuity`); the record is split at every data drop (interval
      > ``gap_factor`` x median) so no analysis window spans a gap. Detection then
      runs independently within each continuous segment. Pass ``fs`` explicitly
      only to override the derived rate (e.g. for a degenerate 1-sample record).
    - **Slope.** Centred difference ``(P[i+h] - P[i-h]) / (2h)`` with
      ``h = round(slope_window_s * fs / 2)`` gives dbar/sample; a sample is *up*
      when the corresponding speed exceeds ``min_slope_dbar_per_s`` downward (an
      ascent lowers pressure). (``sum(diff(window))`` in the MATLAB telescopes to
      this end-to-end difference.)
    - **Debounce.** Majority vote over a centred ``+/- debounce_window_s`` window,
      removing brief flips at the apex/nadir turnarounds and from noise.
    - **Segment.** Each maximal run of the debounced flag is one cast; runs whose
      pressure span is under ``min_span_dbar`` (surface dwell, telemetry stops)
      are dropped, not emitted.

    Casts are numbered sequentially in time; ``profile_number`` pairs each upcast
    with the downcast that precedes it (a new profile opens on each downcast),
    mirroring the Ruskin down+up profile convention.
    """
    t = np.asarray(tstamp, np.int64)
    P = np.asarray(pressure, float)
    n = P.size
    if n == 0:
        return []

    fs_est, _dt_med, gaps = time_continuity(t, gap_factor)
    if not np.isfinite(fs_est) or fs_est <= 0:
        fs_est = fs                                # fall back to caller's rate
    if fs_est is None or not np.isfinite(fs_est) or fs_est <= 0:
        return []                                  # cannot establish a sampling rate

    # window sizes (samples) and slope threshold (dbar/sample) at the steady rate
    half = max(1, int(round(0.5 * slope_window_s * fs_est)))
    m = max(1, int(round(debounce_window_s * fs_est)))
    thresh_per_sample = min_slope_dbar_per_s / fs_est

    # split at data drops so no window straddles a gap; detect within each segment
    brk = np.array([i + 1 for i, _ in gaps], dtype=int)
    seg_starts = np.concatenate(([0], brk))
    seg_ends = np.concatenate((brk, [n]))

    casts: list[Cast] = []
    cn, prof = 0, -1
    for s0, s1 in zip(seg_starts, seg_ends):
        if s1 - s0 < 2 * half + 1:
            continue                               # segment too short to classify
        seg_casts, cn, prof = _detect_casts_segment(
            t[s0:s1], P[s0:s1], half, m, thresh_per_sample, min_span_dbar, cn, prof)
        casts.extend(seg_casts)
    return casts


def casts_from_pressure(con: sqlite3.Connection, chans: list[dict], cfg) -> list[Cast]:
    """Load the full pressure record, report its time continuity, and detect casts.

    Reads only the timestamp + pressure columns (two columns over the whole
    deployment), so its footprint is small even for multi-week records. The
    sampling rate is derived from the record (not `cfg.fs`); see `detect_casts`.
    """
    pcol = next((ch["col"] for ch in chans if ch["var"] == "pressure"), None)
    if pcol is None:
        raise ValueError("no pressure channel discovered; cannot detect casts from pressure")
    rows = con.execute(f"SELECT tstamp, {pcol} FROM data ORDER BY tstamp ASC").fetchall()
    arr = np.asarray(rows, dtype=np.float64)
    if arr.size == 0:
        return []
    t = arr[:, 0].astype(np.int64)
    P = arr[:, 1]

    fs_est, dt_med, gaps = time_continuity(t, cfg.cast_gap_factor)
    print(f"[casts] {t.size:,} samples; steady rate {fs_est:.3f} Hz "
          f"(median dt {dt_med:.0f} ms)")
    if gaps:
        gmax = max(g[1] for g in gaps) / 1000.0
        print(f"[casts] time NOT continuous: {len(gaps)} data drop(s) "
              f"(> {cfg.cast_gap_factor:g}x median), largest {gmax:.1f} s; "
              f"detection runs per continuous segment")
    else:
        print("[casts] time is continuous (no data drops)")

    return detect_casts(
        t, P,
        slope_window_s=cfg.cast_slope_window_s,
        debounce_window_s=cfg.cast_debounce_window_s,
        min_slope_dbar_per_s=cfg.cast_min_slope_dbar_per_s,
        min_span_dbar=cfg.cast_min_span_dbar,
        gap_factor=cfg.cast_gap_factor,
    )


def read_cast_data(con: sqlite3.Connection, cast: Cast,
                   chans: list[dict]) -> dict[str, np.ndarray]:
    """Pull the raw samples for one cast as float64 arrays, keyed by channel var name."""
    cols = ["tstamp"] + [ch["col"] for ch in chans]
    q = (f"SELECT {', '.join(cols)} FROM data "
         f"WHERE tstamp >= ? AND tstamp < ? ORDER BY tstamp ASC")
    rows = con.execute(q, (cast.t1, cast.t2)).fetchall()
    if not rows:
        return {}
    arr = np.asarray(rows, dtype=np.float64)
    out = {"tstamp": arr[:, 0].astype(np.int64)}
    for i, ch in enumerate(chans, start=1):
        out[ch["var"]] = arr[:, i]
    return out


# --------------------------------------------------------------------------- #
# NetCDF metadata
# --------------------------------------------------------------------------- #
VAR_META = {
    "conductivity":             ("S m-1?mS/cm", "sea_water_electrical_conductivity", "mS/cm"),
    "temperature":              ("temp14 (C-T cell thermistor)", "sea_water_temperature", "degC"),
    "pressure":                 ("total pressure", "sea_water_pressure", "dbar"),
    "sea_pressure":             ("pressure - atmosphere", "sea_water_pressure_due_to_sea_water", "dbar"),
    "depth":                    ("depth, positive down", "depth", "m"),
    "practical_salinity":       ("PSS-78", "sea_water_practical_salinity", "1"),
    "absolute_salinity":        ("TEOS-10 SA", "sea_water_absolute_salinity", "g kg-1"),
    "conservative_temperature": ("TEOS-10 CT", "sea_water_conservative_temperature", "degC"),
    "sigma0":                   ("potential density anomaly ref 0 dbar", "sea_water_sigma_theta", "kg m-3"),
    "sound_speed":              ("TEOS-10 sound speed", "speed_of_sound_in_sea_water", "m s-1"),
}


def global_attrs(cfg, level: str) -> dict:
    """Common NetCDF global attributes for a product level, from the config."""
    return {
        "title": f"{cfg.mooring} Wirewalker RBR Concerto CTD - {level}",
        "instrument": cfg.instrument,
        "mooring": cfg.mooring,
        "source_file": cfg.rsk_path.name,
        "processing_level": level,
        "geospatial_lat": cfg.lat,
        "geospatial_lon": cfg.lon,
        "atmospheric_pressure_dbar": cfg.atm_dbar,
        "TEOS10_note": "Salinity uses C-T cell thermistor (temp14). gsw " + gsw.__version__,
        "date_created": _time.strftime("%Y-%m-%dT%H:%M:%S"),
        "Conventions": "CF-1.8",
    }
