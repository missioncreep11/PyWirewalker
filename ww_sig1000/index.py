"""Read the dolfyn `.ad2cp.index` sidecar directly.

dolfyn builds a small index next to every raw `.ad2cp` (one record per ping: byte
offset, ensemble number, record ID, timestamp). Reading it answers two questions
without touching the multi-GB raw file:

- how many ensembles the file holds (`count_ensembles`) — exactly, rather than by
  probing the reader;
- which ensemble a wall-clock time falls on (`ensemble_at_time`), so a deployment
  can be trimmed by date instead of by opaque ensemble number.

`nortek2_lib.get_index` defaults to `eof=2**32` and would silently stop indexing a
file larger than 4 GB; dolfyn's own reader passes the true file length, and so do we.
Building the index is a one-off scan (~4 min for a 17 GiB file); afterwards it is
cached on disk and these calls are fast.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
from mhkit.dolfyn.io import nortek2_lib as _lib

BURST_ID = 21          # 0x15, the standard burst data record

# Index date convention, verified against dolfyn's own decoded timestamps on the
# NOPP1-California record: `year` is the raw Nortek byte (year - 1900) and `month`
# is **1-based**. Note this is NOT the convention of the data record itself, where
# the month is zero-based and `nortek2_lib._calc_time` adds 1 — `_create_index`
# stores an already-incremented month. Do not use `dolfyn.time._fullyear` here
# either; it is for 2-digit years and passes 123 through unchanged.


def read_index(fn, rebuild: bool = False) -> np.ndarray:
    """The dolfyn index for `fn` as a structured array, building it if absent."""
    fn = str(fn)
    index, _ = _lib.get_index(fn, pos=0, eof=os.path.getsize(fn), rebuild=rebuild)
    return index


def _ensemble_starts(index) -> np.ndarray:
    """Index rows that begin a new ensemble (dolfyn's own `_ens_pos` definition)."""
    return np.flatnonzero(_lib._boolarray_firstensemble_ping(index))


def count_ensembles(fn, index=None) -> int:
    """Exact number of ensembles dolfyn will read from `fn`.

    Matches `Ad2cpReader.readfile`'s `nens_total`, so it is a valid upper bound for
    `nens=[start, stop]`. Replaces probing the reader by bisection, which cost many
    reads and rounded *down* to a 5000-ensemble tolerance — silently dropping up to
    ~10 min of data at 8 Hz.
    """
    if index is None:
        index = read_index(fn)
    return int(_ensemble_starts(index).size)


def index_times(rows) -> np.ndarray:
    """datetime64[s] for index records, vectorized (see the convention note above)."""
    y = 1900 + rows["year"].astype("int64")
    years = (y - 1970).astype("timedelta64[Y]") + np.datetime64("1970", "Y")
    days = (years.astype("datetime64[M]") + (rows["month"].astype("int64") - 1)
            ).astype("datetime64[D]") + (rows["day"].astype("int64") - 1)
    return (days.astype("datetime64[s]")
            + rows["hour"].astype("timedelta64[h]")
            + rows["minute"].astype("timedelta64[m]")
            + rows["second"].astype("timedelta64[s]"))


def ensemble_at_time(fn, when, index=None) -> int:
    """First ensemble at or after `when` (anything numpy.datetime64 accepts).

    Returns 0 if `when` precedes the record and `count_ensembles(fn)` if it follows
    it, so the result is always a usable slice bound.
    """
    if index is None:
        index = read_index(fn)
    t = index_times(index[_ensemble_starts(index)])
    return int(np.searchsorted(t, np.datetime64(when, "s"), side="left"))


def summarize(fn) -> dict:
    """Quick description of a raw file from its index alone (no data read)."""
    index = read_index(fn)
    starts = _ensemble_starts(index)
    ids, counts = np.unique(index["ID"], return_counts=True)
    ends = index_times(index[[0, -1]])
    return {"path": str(fn),
            "bytes": os.path.getsize(str(fn)),
            "index_bytes": os.path.getsize(str(fn) + ".index"),
            "records": int(index.size),
            "ensembles": int(starts.size),
            "record_ids": {int(i): int(c) for i, c in zip(ids, counts)},
            "has_beam5": bool(24 in ids),          # 0x18 interleaved burst
            "time_start": str(ends[0]),
            "time_end": str(ends[1])}


def _main():
    import argparse
    import json
    ap = argparse.ArgumentParser(description="Build/inspect a dolfyn .ad2cp index.")
    ap.add_argument("file")
    ap.add_argument("--rebuild", action="store_true", help="discard any existing index")
    ap.add_argument("--at-time", default=None, help="report the ensemble at this ISO time")
    args = ap.parse_args()

    if args.rebuild:
        Path(args.file + ".index").unlink(missing_ok=True)
    info = summarize(args.file)
    print(json.dumps(info, indent=2))
    if args.at_time:
        print(f"\nensemble at {args.at_time}: {ensemble_at_time(args.file, args.at_time):,}")


if __name__ == "__main__":
    _main()
