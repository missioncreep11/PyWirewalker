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
