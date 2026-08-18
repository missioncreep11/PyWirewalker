"""Gold-standard validation: run the ported spectral turbulence method on Devon's
RawProfile2_20231005-1132.nc and compare eps(depth) to the paper's published
NortekTurbulenceData.nc (Fig 3 profile). Mirrors ProcessSingleProfile.m."""
import sys, warnings
from pathlib import Path
import numpy as np
import xarray as xr
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from ww_sig1000.turbulence import process_cast_turbulence
warnings.simplefilter("ignore")

D = ("/Users/drew/TLC/2023/Wirewalker/Nortek/"
     "doi_10_5061_dryad_8sf7m0d44__v20260319/DA_final/DA_final/Data")
raw = xr.open_dataset(f"{D}/RawProfile2_20231005-1132.nc")
turb = xr.open_dataset(f"{D}/NortekTurbulenceData.nc")

# RawProfile2 'time' is MATLAB datenum (label is wrong); convert to unix seconds.
t_datenum = float(raw["time"].values[0])
t_unix = (t_datenum - 719529) * 86400.0
vel = raw["vel"].values.T          # (Range,Time) -> (pings, cells)
corr = raw["corr"].values.T
press = raw["pressure"].values
pitch = raw["pitch"].values
roll = raw["roll"].values
v_a = float(np.nanmedian(raw["v_a"].values))
print(f"RawProfile2: {vel.shape[0]} pings x {vel.shape[1]} cells | v_a={v_a:.4f} | "
      f"press {press.min():.1f}-{press.max():.1f} dbar | unix t0={t_unix:.0f}")

r = process_cast_turbulence(vel, corr, press, pitch, roll, v_a,
                            cellsize=0.06, blockdis=0.1, dep_res=3.0, max_dep=100.0)
z = r["z"]; eps = r["eps"]
print(f"got {np.isfinite(eps).sum()} finite eps bins of {z.size}")

# nearest profile in the published product
tt = turb["time"].values
it = int(np.argmin(np.abs(tt - t_unix)))
print(f"nearest NortekTurbulenceData profile: index {it}, dt={abs(tt[it]-t_unix):.0f}s")
z_ref = turb["depth"].values
eps_ref = turb["Epsilon"].values[it]

# published Fig-3 log10(eps) values printed by ProcessSingleProfile.m
plot_inds = [20, 34, 40, 45]
print("\ndepth-bin match at ProcessSingleProfile plot_inds (0-based i = MATLAB i):")
print(f"  {'bin':>4} {'z(m)':>6} {'log10 ours':>11} {'log10 ref':>10} {'d(dex)':>7}")
for i in plot_inds:
    lo = np.log10(eps[i]) if eps[i] > 0 else np.nan
    lr = np.log10(eps_ref[i]) if eps_ref[i] > 0 else np.nan
    print(f"  {i:>4} {z[i]:>6.2f} {lo:>11.2f} {lr:>10.2f} {lo-lr:>+7.2f}")

# overall offset across matched bins
m = np.isfinite(eps) & np.isfinite(eps_ref) & (eps > 0) & (eps_ref > 0)
d = np.log10(eps[m]) - np.log10(eps_ref[m])
r2 = np.corrcoef(np.log10(eps[m]), np.log10(eps_ref[m]))[0, 1]
print(f"\nALL bins: n={m.sum()} median offset={np.median(d):+.2f} dex "
      f"rms={np.std(d):.2f} corr(log)={r2:.3f}")

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(11, 6))
ax[0].plot(eps, z, "-o", ms=4, label="ours (port)")
ax[0].plot(eps_ref, z_ref, "-s", ms=4, alpha=0.7, label="paper (Dryad)")
ax[0].set_xscale("log"); ax[0].invert_yaxis()
ax[0].set_xlabel("eps (W/kg)"); ax[0].set_ylabel("depth (m)")
ax[0].set_title("RawProfile2 eps(z)"); ax[0].legend()
ax[1].scatter(eps_ref[m], eps[m], s=20)
lims = [1e-10, 1e-4]
ax[1].plot(lims, lims, "k-")
ax[1].set_xscale("log"); ax[1].set_yscale("log")
ax[1].set_xlim(lims); ax[1].set_ylim(lims)
ax[1].set_xlabel("eps paper"); ax[1].set_ylabel("eps ours")
ax[1].set_title(f"n={m.sum()} offset={np.median(d):+.2f} dex corr={r2:.2f}")
fig.tight_layout(); fig.savefig("/Users/drew/PyWirewalker/RawProfile2_validation.png", dpi=120)
print("wrote RawProfile2_validation.png")
