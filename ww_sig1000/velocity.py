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
from .motion import beam_motion_correction, beam_motion_correction_v2


def output_grid(boxsize: float, z_max: float) -> np.ndarray:
    """Depth-bin centres (m, positive down): boxsize/2, ..., up to z_max."""
    edges = np.arange(0.0, z_max + boxsize, boxsize)
    return 0.5 * (edges[:-1] + edges[1:])


# --------------------------------------------------------------------------- #
# per-beam Doppler noise vs beam correlation
# --------------------------------------------------------------------------- #
# Measured on NOPP_d2 (Sig1000, 1 m cells, 8 Hz broadband) from within-ping
# adjacent-cell velocity differences, where platform motion cancels exactly
# (S100430A038_NOPP_d2/corr_weighting.py, 8 casts pooled). Instrument- and
# configuration-specific: remeasure for other cell sizes / bandwidths.
_SIGMA_CORR_MID = np.array([32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5,
                            67.5, 72.5, 77.5, 82.5, 87.5, 92.5, 97.5])
_SIGMA_CORR_VAL = np.array([0.1845, 0.1562, 0.1321, 0.1143, 0.0996, 0.0891,
                            0.0797, 0.0713, 0.0650, 0.0587, 0.0545, 0.0514,
                            0.0482, 0.0430])


def beam_sigma(corr) -> np.ndarray:
    """Per-sample beam-velocity Doppler noise (m/s) from beam correlation."""
    return np.interp(np.asarray(corr, float), _SIGMA_CORR_MID, _SIGMA_CORR_VAL)


# --------------------------------------------------------------------------- #
# depth-gated per-bin notch (opt-in alternative to the boxcar bin mean)
# --------------------------------------------------------------------------- #
# Above this depth the bin mean is replaced by the constant of a ridged
# constant + wave-band-sinusoid fit over the bin's dwell; below it the plain
# mean is kept. Synthetic-truth study on real NOPP_d2 geometry
# (S100430A038_NOPP_d2/wave_inversion.py): 7-17% rms improvement in the top
# 50 m, a comparable penalty where no waves exist - hence the gate.
#
# The gain lives entirely on GAPPY dwells (real ones keep only 17-32% of
# near-surface samples): gaps break the boxcar's wave cancellation, and the fit
# restores it by modelling the wave on the actual sample support. On complete
# uniform sampling the harmonic columns are orthogonal to the constant and the
# notch reduces to the boxcar - it cannot hurt there.
NOTCH_MAX_DEPTH_M = 60.0
NOTCH_BAND = (0.04, 0.5)          # Hz, the surface-wave band
NOTCH_RIDGE = 0.05                # prior scale (m/s) on the wave coefficients
NOTCH_SIGMA = 0.13                # sample-noise scale setting the ridge strength


def _nan_gaussian(x, win=8, alpha=2.5, axis=-1):
    """NaN-aware Gaussian smoothing, matching MATLAB smoothdata(...,'gaussian',win)
    (window of `win` samples, alpha = 2.5)."""
    from scipy.ndimage import convolve1d
    m = np.arange(win) - (win - 1) / 2.0
    k = np.exp(-0.5 * (2 * alpha * m / (win - 1)) ** 2)
    k /= k.sum()
    w = np.isfinite(x).astype(float)
    num = convolve1d(np.nan_to_num(x), k, axis=axis, mode="nearest")
    den = convolve1d(w, k, axis=axis, mode="nearest")
    with np.errstate(invalid="ignore"):
        out = num / den
    out[w == 0] = np.nan
    return out


def _notch_constant(tj, vj, *, band=NOTCH_BAND, ridge=NOTCH_RIDGE, sigma=NOTCH_SIGMA):
    """Constant of a ridged (constant + wave-band sin/cos) fit; None if unfittable.

    The ridge on the wave columns is required: gappy dwells make them
    near-collinear and the unridged fit diverges catastrophically.
    """
    T = float(tj.max() - tj.min())
    if T <= 0:
        return None
    f = np.arange(max(np.ceil(band[0] * T), 1), np.floor(band[1] * T) + 1) / T
    p = 1 + 2 * f.size
    if f.size == 0 or p >= vj.size:
        return None
    X = np.empty((vj.size, p))
    X[:, 0] = 1.0
    ph = 2 * np.pi * np.outer(tj - tj.min(), f)
    X[:, 1::2] = np.cos(ph)
    X[:, 2::2] = np.sin(ph)
    R = np.zeros((p - 1, p))
    R[:, 1:] = np.eye(p - 1) * (sigma / ridge)
    c, *_ = np.linalg.lstsq(np.vstack([X, R]),
                            np.concatenate([vj, np.zeros(p - 1)]), rcond=None)
    return float(c[0])


def _nominal_depth_grid(pressure, n_cells, cell_size, blank_dist, direction):
    """Per-ping nominal (no-tilt) cell depths, z positive up (MATLAB ``dpth_temp``)."""
    cell = np.arange(n_cells)
    off = blank_dist + cell_size * cell            # (n_cells,)
    if direction == "up":
        return -pressure[:, None] + off[None, :]   # cells above instrument
    return -pressure[:, None] - off[None, :]        # cells below (downward-looking)


def process_cast(dsc, *, corr_min=50, boxsize=1.0, z_max=110.0, direction="up",
                 min_bin_samples=10, cell_size=None, blank_dist=None, beam_angle=None,
                 motion_correct=True, motion="v1", bin_average="boxcar"):
    """Grid one cast (a dolfyn Dataset subset in **beam** coords) to a depth profile.

    Returns dict with keys velE, velN, velU, amp, shearE, shearN, n_obs and
    per-component ``*_sem`` standard errors (each (nz,)), plus scalars time
    (datetime64), pressure_max, pressure_min, and the depth grid `z`.

    ``shearE``/``shearN`` are the beam-differenced shear (WWvel_upward's
    ``beamshear``): centred cell differences along each beam before any
    rotation, so per-ping common-mode errors - platform motion, attitude
    leakage, the sail term - cancel exactly and no motion correction applies.

    ``bin_average="notch"`` replaces the bin mean above `NOTCH_MAX_DEPTH_M` with
    the constant of a ridged constant + wave-band fit over the bin's dwell
    (surface-wave residual suppression); the plain mean is kept below the gate
    and everywhere when "boxcar" (the default).

    The ``*_sem`` values are the Doppler-noise standard error of each bin mean,
    from the measured `beam_sigma` correlation-noise relation. They exclude
    correlated errors (attitude systematics, residual wave contamination), so
    they are a floor, not a total error bar.

    ``motion="v2"`` uses the buoyant-ascent correction (`beam_motion_correction_v2`):
    its low-passed accel tilt then replaces the AHRS pitch/roll in *every* rotation
    (cell depths, beam->ENU, and the correction itself), spike pings are excluded
    from the bin averages, and the result records ``tilt_source``. When the v2 tilt
    guard fails the cast falls back to the v1 path (``tilt_source = "ahrs_fallback"``).
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

    # v2 motion model: attitude and correction come from the raw sensors together
    v2 = None
    tilt_source = None
    heading_source = None
    if motion == "v2" and "accel" in dsc:
        ts_v2 = dsc["time"].values.astype("datetime64[ns]").astype("int64") / 1e9
        v2 = beam_motion_correction_v2(ts_v2, press, dsc["accel"].values, head,
                                       float(a["fs"]),
                                       mag=dsc["mag"].values if "mag" in dsc else None,
                                       beam_angle=ba)
        if v2.usable:
            pitch, roll = v2.pitch_deg, v2.roll_deg
            head = v2.heading_deg          # Stage-2 compass when the field is sane
            tilt_source = "lp_accel"
            heading_source = v2.heading_source
        else:
            v2 = None
            tilt_source = "ahrs_fallback"
            heading_source = "ahrs"

    # 1. correlation mask: any beam below threshold at a (cell, ping) -> drop all beams there
    bad = (corr < corr_min).any(axis=0)        # (range, time)
    vel[:, bad] = np.nan

    # 2. tilt-corrected per-beam cell depths (nping, ncell, 4), positive up
    z, ranges, bZ = cell_depths(press, pitch, roll, nc, cs, bd, beam_angle_deg=ba)

    # 2b. beam shear (port of WWvel_upward's beamshear): centered difference
    # (v[c+2]-v[c])/(z[c+2]-z[c]) along each beam on the RAW beam velocities with
    # the tilt-corrected per-beam cell heights. Anything constant across a ping's
    # cells - platform translation, attitude-error leakage of the ascent, the
    # sail term - cancels exactly in the difference, so beam shear is immune to
    # the whole motion/attitude error family and gets no motion correction.
    z_bcn = np.transpose(z, (2, 1, 0))         # (beam, cell, ping)
    with np.errstate(invalid="ignore", divide="ignore"):
        bshear = ((vel[:, 2:, :] - vel[:, :-2, :])
                  / (z_bcn[:, 2:, :] - z_bcn[:, :-2, :]))  # centred at cells 1..nc-2

    # 3. per-ping nominal depth grid + interpolate each beam onto it
    dpth_nom = _nominal_depth_grid(press, nc, cs, bd, direction)   # (nping, ncell)
    eq_vel = np.full((nb, nc, npg), np.nan)
    eq_shear = np.full((nb, max(nc - 2, 1), npg), np.nan)
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
            if nc > 2:
                zs = z[n, 1:-1, b]
                s = bshear[b, :, n]
                oks = np.isfinite(zs) & np.isfinite(s)
                if oks.sum() >= 2:
                    zss, ss = zs[oks], s[oks]
                    osrt = np.argsort(zss)
                    eq_shear[b, :, n] = np.interp(tgt[1:-1], zss[osrt], ss[osrt],
                                                  left=np.nan, right=np.nan)

    # 3b. platform-motion correction: add the platform's along-beam velocity
    if motion_correct and v2 is not None:
        eq_vel = eq_vel + v2.corr_beam[:, None, :]                 # broadcast over cells
        eq_vel[:, :, ~v2.ping_ok] = np.nan                         # spike pings excluded
    elif motion_correct and "accel" in dsc:
        ts = dsc["time"].values.astype("datetime64[ns]").astype("int64") / 1e9
        corr_beam, _ = beam_motion_correction(
            ts, press, dsc["accel"].values, pitch, roll, head, float(a["fs"]), beam_angle=ba)
        eq_vel = eq_vel + corr_beam[:, None, :]                    # broadcast over cells

    # 4. beam -> ENU on the depth-aligned velocities; likewise the beam shears
    # (d/dz is linear, so rotating the four beam shears is rotating the shear).
    # MATLAB then Gaussian-smooths along pings (window 8 ~ 1 s) and masks the two
    # cells nearest the transducer (stagnation-point contamination).
    enu = beam2enu(eq_vel, head, pitch, roll, theta_deg=ba)        # (3, ncell, nping)
    shear_enu = beam2enu(eq_shear, head, pitch, roll, theta_deg=ba)[:2]
    shear_enu[:, :2, :] = np.nan
    shear_enu = _nan_gaussian(shear_enu, axis=-1)

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
    # per-sample Doppler variance from the corr-noise relation: beam-mean beam
    # variance projected through the beam->ENU geometry (E/N amplified by
    # 1/(sqrt(2) sin th), U reduced by 1/(2 cos th)). corr is indexed on the
    # measured cells while comps sit on the per-ping nominal grid - a benign
    # approximation, the two differ by the tilt correction only.
    th = np.deg2rad(ba)
    s2_beam = np.mean(beam_sigma(corr) ** 2, axis=0)          # (ncell, nping)
    var_fac = {"velE": 1.0 / (2 * np.sin(th) ** 2),
               "velN": 1.0 / (2 * np.sin(th) ** 2),
               "velU": 1.0 / (4 * np.cos(th) ** 2)}
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
        if name in var_fac:                       # var(mean) = sum(sigma_i^2)/n^2
            v2 = (s2_beam * var_fac[name]).ravel()
            ssum2 = np.bincount(ib[m], weights=v2[m], minlength=nz)
            out[f"{name}_sem"] = (np.sqrt(ssum2) / cnt).astype(np.float32)
    out["n_obs"] = nobs

    # 5b. bin the beam-differenced shear onto the same grid (interior cells at
    # their nominal depths). The Gaussian ping-smoothing is mean-preserving, so
    # the SEM uses the unsmoothed per-sample variance: 2 sigma_b^2 / dz^2 per
    # beam through the same beam->ENU geometry factor.
    if nc > 2:
        depth_sh = -dpth_nom[:, 1:-1].T                  # (ncell-2, nping), +down
        ib_sh = np.digitize(depth_sh.ravel(), edges) - 1
        valid_sh = (ib_sh >= 0) & (ib_sh < nz) & np.isfinite(depth_sh.ravel())
        ib_sh = np.clip(ib_sh, 0, nz - 1)
        dz_eff = 2.0 * cs * np.cos(th)                   # centred-difference span
        v2s = (s2_beam[1:-1, :] * var_fac["velE"] * 2.0 / dz_eff ** 2).ravel()
        for name, arr in (("shearE", shear_enu[0]), ("shearN", shear_enu[1])):
            flat = arr.ravel()
            m = valid_sh & np.isfinite(flat)
            ssum = np.bincount(ib_sh[m], weights=flat[m], minlength=nz)
            cnt = np.bincount(ib_sh[m], minlength=nz).astype(float)
            cnt[cnt < min_bin_samples] = np.nan
            out[name] = (ssum / cnt).astype(np.float32)
            ssum2 = np.bincount(ib_sh[m], weights=v2s[m], minlength=nz)
            out[f"{name}_sem"] = (np.sqrt(ssum2) / cnt).astype(np.float32)
    else:                                                # too few cells for shear
        for name in ("shearE", "shearN"):
            out[name] = np.full(nz, np.nan, np.float32)
            out[f"{name}_sem"] = np.full(nz, np.nan, np.float32)

    # depth-gated notch: refit the near-surface bins where surface waves live
    if bin_average == "notch":
        ts_s = dsc["time"].values.astype("datetime64[ns]").astype("int64") / 1e9
        tflat = np.broadcast_to(ts_s[None, :], depth.shape).ravel()
        gate = np.flatnonzero(zc_grid < NOTCH_MAX_DEPTH_M)
        for name in ("velE", "velN", "velU"):
            flat = comps[name].ravel()
            for j in gate:
                m = valid & np.isfinite(flat) & (ib == j)
                if m.sum() < max(min_bin_samples, 20):
                    continue
                c = _notch_constant(tflat[m], flat[m])
                if c is not None:
                    out[name][j] = c

    t = dsc["time"].values
    out["time"] = t[t.size // 2]
    out["pressure_max"] = float(np.nanmax(press))
    out["pressure_min"] = float(np.nanmin(press))
    if tilt_source is not None:
        out["tilt_source"] = tilt_source
        out["heading_source"] = heading_source
    return out
