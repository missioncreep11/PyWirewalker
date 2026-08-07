"""Per-cast velocity processing -> vertically gridded ENU profile.

Port of the core of ``WWvel_upward.m`` (motion correction handled separately in
``motion.py``). For one cast:

1. mask samples where any beam correlation < ``corr_min``;
2. compute each beam's tilt-corrected cell depth (``geometry.cell_depths``);
3. interpolate each beam's velocity onto a per-ping nominal depth grid so the four
   beams are depth-aligned (MATLAB ``interp1(z_coords, beam_vel, dpth_temp)``);
4. rotate the aligned beam velocities to ENU (``transforms.beam2enu``);
5. box-average the ENU samples onto the output depth grid ``0 : boxsize : z_max``.

Returns a dict of (nz,) arrays for one cast.
"""
from __future__ import annotations

import numpy as np

from .transforms import beam2enu
from .geometry import cell_depths
from .motion import beam_motion_correction


def output_grid(boxsize: float, z_max: float) -> np.ndarray:
    """Depth-bin centres (m, positive down): boxsize/2, ..., up to z_max."""
    edges = np.arange(0.0, z_max + boxsize, boxsize)
    return 0.5 * (edges[:-1] + edges[1:])


def _nominal_depth_grid(pressure, n_cells, cell_size, blank_dist, direction):
    """Per-ping nominal (no-tilt) cell depths, z positive up (MATLAB ``dpth_temp``)."""
    cell = np.arange(n_cells)
    off = blank_dist + cell_size * cell            # (n_cells,)
    if direction == "up":
        return -pressure[:, None] + off[None, :]   # cells above instrument
    return -pressure[:, None] - off[None, :]        # cells below (downward-looking)


def process_cast(dsc, *, corr_min=50, boxsize=1.0, z_max=110.0, direction="up",
                 min_bin_samples=10, cell_size=None, blank_dist=None, beam_angle=None,
                 motion_correct=True):
    """Grid one cast (a dolfyn Dataset subset in **beam** coords) to a depth profile.

    Returns dict with keys velE, velN, velU, amp, corr_mean, n_obs (each (nz,)),
    plus scalars time (datetime64), pressure_max, and the depth grid `z`.
    """
    a = dsc.attrs
    cs = a["cell_size"] if cell_size is None else cell_size
    bd = a["blank_dist"] if blank_dist is None else blank_dist
    ba = a.get("beam_angle", 25) if beam_angle is None else beam_angle

    vel = dsc["vel"].values.astype(float)      # (beam, range, time), beam coords
    corr = dsc["corr"].values                  # (beam, range, time)
    amp = dsc["amp"].values.astype(float)      # (beam, range, time)
    pitch = dsc["pitch"].values
    roll = dsc["roll"].values
    head = dsc["heading"].values
    press = dsc["pressure"].values
    nb, nc, npg = vel.shape

    # 1. correlation mask: any beam below threshold at a (cell, ping) -> drop all beams there
    bad = (corr < corr_min).any(axis=0)        # (range, time)
    vel[:, bad] = np.nan

    # 2. tilt-corrected per-beam cell depths (nping, ncell, 4), positive up
    z, ranges, bZ = cell_depths(press, pitch, roll, nc, cs, bd, beam_angle_deg=ba)

    # 3. per-ping nominal depth grid + interpolate each beam onto it
    dpth_nom = _nominal_depth_grid(press, nc, cs, bd, direction)   # (nping, ncell)
    eq_vel = np.full((nb, nc, npg), np.nan)
    for n in range(npg):
        if not np.isfinite(z[n, 0, 0]):        # ping's beams not pointing up -> skip
            continue
        tgt = dpth_nom[n]
        for b in range(nb):
            zc = z[n, :, b]
            v = vel[b, :, n]
            ok = np.isfinite(zc) & np.isfinite(v)
            if ok.sum() < 2:
                continue
            zc, v = zc[ok], v[ok]
            o = np.argsort(zc)
            eq_vel[b, :, n] = np.interp(tgt, zc[o], v[o], left=np.nan, right=np.nan)

    # 3b. platform-motion correction: add the platform's along-beam velocity
    if motion_correct and "accel" in dsc:
        ts = dsc["time"].values.astype("datetime64[ns]").astype("int64") / 1e9
        corr_beam, _ = beam_motion_correction(
            ts, press, dsc["accel"].values, pitch, roll, head, float(a["fs"]), beam_angle=ba)
        eq_vel = eq_vel + corr_beam[:, None, :]                    # broadcast over cells

    # 4. beam -> ENU on the depth-aligned velocities
    enu = beam2enu(eq_vel, head, pitch, roll, theta_deg=ba)        # (3, ncell, nping)

    # 5. box-average ENU samples onto the output depth grid
    zc_grid = output_grid(boxsize, z_max)
    nz = zc_grid.size
    edges = np.arange(0.0, z_max + boxsize, boxsize)
    depth = -dpth_nom.T                          # (ncell, nping) positive down
    ib = np.digitize(depth.ravel(), edges) - 1
    valid = (ib >= 0) & (ib < nz) & np.isfinite(depth.ravel())

    out = {"z": zc_grid}
    comps = {"velE": enu[0], "velN": enu[1], "velU": enu[2],
             "amp": np.nanmean(amp, axis=0)}     # amp: mean over beams (ncell, nping)
    nobs = np.zeros(nz, int)
    for name, arr in comps.items():
        flat = arr.ravel()
        m = valid & np.isfinite(flat)
        ssum = np.bincount(ib[m], weights=flat[m], minlength=nz)
        cnt = np.bincount(ib[m], minlength=nz).astype(float)
        if name == "velE":
            nobs = cnt.astype(int)
        cnt[cnt < min_bin_samples] = np.nan
        out[name] = (ssum / cnt).astype(np.float32)
    out["n_obs"] = nobs

    t = dsc["time"].values
    out["time"] = t[t.size // 2]
    out["pressure_max"] = float(np.nanmax(press))
    return out
