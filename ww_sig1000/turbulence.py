"""HR beam-5 turbulent dissipation (spectral method).

Python port of the spectral path of ``WWturb_upward.m`` (Northcott et al. 2026,
JPO) — the structure-function path is intentionally omitted. For each upcast, the
pulse-coherent HR center-beam (beam 5) velocity is unwrapped/demeaned, the
stagnation zone and spikes removed, detrended in range, then per depth bin the
along-beam velocity-fluctuation wavenumber spectra are averaged and fit to an
inertial-subrange + noise model ``S(k) = N + A k^{-5/3}`` (A = C_k eps^{2/3},
C_k = 0.53) to give eps = (A/C_k)^{3/2}.

Notes
-----
- ``vel_b5`` from dolfyn is the *wrapped* pulse-coherent velocity (dolfyn does not
  unwrap); the angular demean here removes the platform rise, which aliases when it
  exceeds the ambiguity velocity ``v_a`` (= ``ambig_vel_b5``).
- Beam 5 is the instrument-frame vertical beam (phi=90 deg); its earth-frame
  vertical component ``bZ`` (for cell depth) comes from the attitude via
  ``transforms.get_unit_vectors`` with the toolbox's roll/pitch order.
"""
from __future__ import annotations

import numpy as np

from .transforms import get_unit_vectors

C_K = 0.53          # Kolmogorov constant (Sreenivasan 1995), radian-wavenumber spectra


# --------------------------------------------------------------------------- #
# Pure helpers
# --------------------------------------------------------------------------- #
def fit_kolmogorov(k, S):
    """Least-squares fit ``S(k) = N + A k^{-5/3}`` (port of FitKolmogorov.m).

    k, S 1-D (NaNs dropped). Returns (N, A). A<0 is allowed (invalid fit -> eps NaN).
    """
    k = np.asarray(k, float); S = np.asarray(S, float)
    m = np.isfinite(k) & np.isfinite(S)
    k, S = k[m], S[m]
    if k.size < 2:
        return np.nan, np.nan
    G = np.column_stack([np.ones_like(k), k ** (-5.0 / 3.0)])
    coeffs, *_ = np.linalg.lstsq(G, S, rcond=None)
    return float(coeffs[0]), float(coeffs[1])          # N, A


def epsilon_from_A(A, c_k=C_K):
    """eps = (A/C_k)^{3/2}; NaN where A<=0."""
    A = np.asarray(A, float)
    out = np.full(A.shape, np.nan) if A.ndim else np.nan
    with np.errstate(invalid="ignore"):
        val = (A / c_k) ** 1.5
    return np.where(A > 0, val, np.nan) if np.ndim(A) else (val if A > 0 else np.nan)


def _nanmean_complex(z, axis):
    """Mean of a complex array ignoring non-finite entries."""
    m = np.isfinite(z)
    s = np.where(m, z, 0).sum(axis=axis)
    c = m.sum(axis=axis)
    with np.errstate(invalid="ignore"):
        return s / np.where(c > 0, c, np.nan)


def angular_demean(vel, v_a, skip_cells=15):
    """Unwrap/demean the wrapped HR velocity (paper Eq 2-3; WWturb demean block).

    vel : (npings, ncells) wrapped radial velocity (NaN where masked).
    v_a : scalar ambiguity velocity. Per-ping mean phase (over interior cells
    ``skip_cells:-1``) is removed. Returns demeaned velocity (same shape).
    """
    vel = np.asarray(vel, float)
    phase = vel / v_a * 2 * np.pi
    z = np.exp(1j * phase)                       # NaN vel -> non-finite z
    mp = _nanmean_complex(z[:, skip_cells:-1], axis=1)     # per-ping mean complex
    z = z * np.exp(-1j * np.angle(mp))[:, None]
    return np.angle(z) / (2 * np.pi) * v_a


def _detrend_rows_nan(x):
    """Linearly detrend each row vs column index, ignoring NaN (MATLAB detrend omitnan)."""
    x = np.asarray(x, float)
    n, m = x.shape
    idx = np.arange(m, dtype=float)
    out = x.copy()
    for i in range(n):
        row = x[i]; ok = np.isfinite(row)
        if ok.sum() < 2:
            continue
        s, b = np.polyfit(idx[ok], row[ok], 1)
        out[i] = row - (s * idx + b)
    return out


def _raw_spectrum(vel_seg, cellsize):
    """One-sided radian-wavenumber spectrum of a velocity segment (ProcessSingleProfile).

    Demean, zero-fill NaN, Hamming window, |FFT|^2 with variance-preserving
    normalization / radian-wavenumber density. Returns (f, ft) with f in cyc/m.
    """
    v = np.asarray(vel_seg, float)
    v = v - np.nanmean(v)
    v = np.nan_to_num(v, nan=0.0)
    n = v.size
    w = np.hamming(n)
    ft = np.abs(np.fft.fft(v * w)) ** 2 * (2.0 / np.sum(w ** 2)) / n
    f = (1.0 / cellsize) * np.arange(n // 2 + 1) / n          # cyc/m
    df = np.median(2 * np.pi * np.diff(f))
    return f, ft[: n // 2 + 1] / df


def _k_grid(cellsize, ks_len):
    """Wavenumber grid ks (cyc/m) and per-ks band edges (ProcessSingleProfile)."""
    ks = (1.0 / cellsize) * np.arange(ks_len // 2) / ks_len
    dk = (1.0 / cellsize) / ks_len
    mid = (ks[1:] + ks[:-1]) / 2
    lo = np.concatenate([mid - dk, [mid[-1]]])
    hi = np.concatenate([mid, [mid[-1] + dk]])
    return ks, lo, hi


def hr_bins(cellsize, dep_res, max_dep):
    """Spectral depth-bin centres and averaging window (ProcessSingleProfile).

    ``vel_inds = round(dep_res/cellsize/2)``; centres are
    ``vel_inds*cellsize/2 : vel_inds*cellsize : max_dep`` (spacing ~= dep_res/2);
    the depth-averaging window is ``dep_res``. (dep_res=3 m, cellsize 0.06 ->
    67 centres 0.75, 2.25, ... 99.75 m, matching NortekTurbulenceData.nc.)
    """
    vel_inds = int(round(dep_res / cellsize / 2))
    step = vel_inds * cellsize
    centres = np.arange(step / 2, max_dep, step)
    return centres, dep_res


# second-order structure-function constant in D(r) = N + C_SF eps^2/3 r^2/3.
# Literature C_v^2 ~ 2.0-2.1 (Wiles et al. 2006), but that leaves the SF eps ~0.20 dex
# below the spectral eps. C_SF is cross-calibrated to the validated spectral eps on the
# high-scattering TLC deployment (median offset zeroed over high-SNR bins <80 m), so the
# two estimators share one absolute scale. Note: eps-N coupling (the QC use) is invariant
# to C_SF, so this only sets the absolute level.
C_SF = 1.476


def _structure_function_profile(vel_use, z_use, centres, *, cellsize, dep_res,
                                fit_lags=15, min_ref_pings=6):
    """Second-order velocity structure function eps per depth bin (WWturb_upward.m).

    For each grid depth, take each ping's cell nearest that depth as the anchor and
    accumulate D(r) = <(v(anchor) - v(anchor +/- r))^2> over both sides and all pings,
    counting only valid (non-NaN) pairs (no zero-filling — unlike the spectral path).
    Fit D(r) = N + A r^2/3 over the first `fit_lags` lags; eps = (A/C_SF)^3/2.
    Returns (eps, N, A) arrays on `centres`.
    """
    vel_use = np.asarray(vel_use, float)
    npings, ncols = vel_use.shape
    vel_inds = min(int(round(dep_res / cellsize / 2)), ncols - 1)
    nb = centres.size
    eps = np.full(nb, np.nan); Nsf = np.full(nb, np.nan); Asf = np.full(nb, np.nan)
    if vel_inds < 3:
        return eps, Nsf, Asf
    r = cellsize * np.arange(1, vel_inds + 1)
    nfit = int(min(fit_lags, vel_inds))
    G = np.column_stack([np.ones(nfit), r[:nfit] ** (2.0 / 3.0)])
    # squared velocity differences at each lag L (v[c]-v[c+L])^2, reused for both sides
    sq = [(vel_use[:, :ncols - L] - vel_use[:, L:]) ** 2 for L in range(1, vel_inds + 1)]
    pidx = np.arange(npings)
    for ib, zc in enumerate(centres):
        ref = np.argmin(np.abs(z_use - zc), axis=1)                    # nearest cell / ping
        good = (np.abs(z_use[pidx, ref] - zc) <= cellsize / 2) & (ref > 0) & (ref < ncols - 1)
        if good.sum() < min_ref_pings:
            continue
        gp, rc = pidx[good], ref[good]
        D = np.full(vel_inds, np.nan)
        for li, L in enumerate(range(1, vel_inds + 1)):
            vals = []
            up = rc <= ncols - 1 - L
            if up.any():
                vals.append(sq[li][gp[up], rc[up]])                    # (v[c]-v[c+L])^2
            dn = rc - L >= 0
            if dn.any():
                vals.append(sq[li][gp[dn], rc[dn] - L])                # (v[c-L]-v[c])^2
            if vals:
                allv = np.concatenate(vals)
                if np.isfinite(allv).any():
                    D[li] = np.nanmean(allv)
        y = D[:nfit]; m = np.isfinite(y)
        if m.sum() < 3:
            continue
        coeffs, *_ = np.linalg.lstsq(G[m], y[m], rcond=None)
        Nsf[ib], Asf[ib] = float(coeffs[0]), float(coeffs[1])
        if Asf[ib] > 0:
            eps[ib] = (Asf[ib] / C_SF) ** 1.5
    return eps, Nsf, Asf


# --------------------------------------------------------------------------- #
# Per-cast spectral turbulence
# --------------------------------------------------------------------------- #
def process_cast_turbulence(vel, corr, pressure, pitch, roll, v_a, *,
                            cellsize, blockdis, dep_res=3.0, max_dep=100.0,
                            time=None, amp=None, corr_min=50, fit_min_k=0.5, min_specs=6,
                            spec_floor=1e-8, skip_cells=14, deep_thresh=50.0,
                            phi_deg=90.0, azi_deg=0.0, structure_function=True):
    """Spectral eps for one HR beam-5 upcast (port of ProcessSingleProfile.m).

    vel, corr : (npings, ncells) HR beam-5 velocity / correlation.
    pressure, pitch, roll : (npings,) aligned to the HR pings (deg for pitch/roll).
    v_a : scalar ambiguity velocity (median ambig_vel_b5 over the cast).
    dep_res : depth resolution (m); bins are dep_res/2 spaced, dep_res wide (3 m -> 67
        bins 0.75..99.75, matching NortekTurbulenceData.nc).

    Returns dict with (nbin,) arrays eps, N, SNR, A, corr, num and (nbin, nk) spec,
    plus z (bin centres), k (ks), time. None if the stagnation cutoff fails.
    """
    vel = np.array(vel, float)
    corr = np.asarray(corr)
    npings, ncells = vel.shape

    # 1. correlation mask
    vel[corr < corr_min] = np.nan
    # 2. angular demean (unwrap + remove platform motion)
    vel_z = angular_demean(vel, v_a, skip_cells)
    # 3. cell vertical position: z_coord = pressure - bZ*range
    _, _, bZ = get_unit_vectors(np.deg2rad(phi_deg), np.deg2rad(azi_deg),
                                np.deg2rad(roll), np.deg2rad(pitch))
    ranges = blockdis + cellsize * np.arange(ncells)
    z_coord = pressure[:, None] - bZ[:, None] * ranges[None, :]
    # 4. stagnation-zone cutoff: median deep-ping profile, first cell where it drops
    #    within 4 mm/s of its own median; of those crossings pick the one nearest 1 m.
    deep = pressure > deep_thresh
    if deep.sum() < 2:
        deep = pressure > np.nanmedian(pressure)
    avg_deep = np.nanmedian(vel_z[deep], axis=0)
    cutoff_vel = np.nanmedian(avg_deep) + 0.004
    cross = np.flatnonzero(np.diff((avg_deep > cutoff_vel).astype(int)) == -1)
    if cross.size == 0:
        return None
    cutoff = int(cross[np.argmin(np.abs(ranges[cross] - 1.0))])
    # 5. subtract the deep (stagnation) profile from every ping
    vel_z = vel_z - avg_deep[None, :]
    # 6. despike (single point with large dv/dz on both sides)
    filled = np.nan_to_num(vel_z, nan=0.0)
    d1 = np.diff(filled[:, :-1], axis=1)
    d2 = np.diff(filled[:, 1:], axis=1)
    sp = np.zeros((npings, ncells), bool)
    sp[:, 1:-1] = ((d1 < -0.005) & (d2 > 0.005)) | ((d1 > 0.005) & (d2 < -0.005))
    vel_z[sp] = np.nan
    # 7. detrend in range + subset [cutoff : ncells-2]
    sl = slice(cutoff, ncells - 2)
    vel_use = _detrend_rows_nan(vel_z[:, sl])
    z_use = z_coord[:, sl]
    corr_use = corr[:, sl].astype(float)
    amp_use = np.asarray(amp, float)[:, sl] if amp is not None else None
    ncols = vel_use.shape[1]

    # 8. common wavenumber grid + band edges (ProcessSingleProfile: ks_len from data)
    ks_len = int(min(round(dep_res / cellsize), ncols))
    ks, klo, khi = _k_grid(cellsize, ks_len)
    centres, win = hr_bins(cellsize, dep_res, max_dep)
    nb, nk = centres.size, ks.size
    out = {k: np.full(nb, np.nan) for k in
           ("eps", "N", "SNR", "A", "corr", "num", "amp", "eps_sf", "N_sf", "A_sf")}
    spec_out = np.full((nb, nk), np.nan)
    fit_cut = int(np.flatnonzero(ks < fit_min_k)[-1] + 1) if np.any(ks < fit_min_k) else 1

    # 9. spectral eps per depth bin
    for i, zc in enumerate(centres):
        dep_mask = (z_use >= zc - win / 2) & (z_use < zc + win / 2)   # (npings, ncols)
        prof = np.flatnonzero(dep_mask.sum(axis=1) > 4)
        if prof.size == 0:
            continue
        specs = np.full((nk, prof.size), np.nan)
        for jj, j in enumerate(prof):
            f, ft = _raw_spectrum(vel_use[j, dep_mask[j]], cellsize)
            for ik in range(nk):                      # band-average onto ks bins
                km = (f > klo[ik]) & (f < khi[ik])
                if km.any():
                    specs[ik, jj] = np.nanmean(ft[km])
        num_specs = int(np.max(np.sum(~np.isnan(specs), axis=1))) if specs.size else 0
        # cull high-power outlier spectra. MATLAB threshold: per-column power sum over
        # ALL columns (all-NaN columns count as 0), std ddof=1.
        colsum = np.nansum(specs, axis=0)
        allnan = np.all(np.isnan(specs), axis=0)
        if specs.shape[1] > 1:
            bad = allnan | (colsum > colsum.mean() + 2 * colsum.std(ddof=1))
        else:
            bad = allnan
        specs = specs[:, ~bad]
        if specs.size == 0:
            continue
        spec_fit = np.nanmean(specs, axis=1)
        spec_fit[spec_fit < spec_floor] = np.nan

        N, A = fit_kolmogorov(2 * np.pi * ks[fit_cut:], spec_fit[fit_cut:])
        if num_specs > min_specs:
            out["eps"][i] = epsilon_from_A(A)
            out["N"][i] = N
            out["SNR"][i] = A / N if N else np.nan
            out["A"][i] = A
            out["corr"][i] = np.nanmean(corr_use[dep_mask])
            if amp_use is not None:
                out["amp"][i] = np.nanmean(amp_use[dep_mask])
            out["num"][i] = num_specs
            spec_out[i] = spec_fit

    # 10. structure-function eps on the same depth grid (shares all preprocessing;
    #     computed from velocity differences, so gaps are skipped not zero-filled)
    if structure_function:
        out["eps_sf"], out["N_sf"], out["A_sf"] = _structure_function_profile(
            vel_use, z_use, centres, cellsize=cellsize, dep_res=dep_res)

    out["spec"] = spec_out
    out["z"] = centres
    out["k"] = ks
    out["time"] = np.datetime64(time, "ns") if time is not None else np.datetime64("NaT")
    return out
