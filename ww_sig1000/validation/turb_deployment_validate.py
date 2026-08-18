"""Validate the full-deployment turbulence product against the paper's published
NortekTurbulenceData.nc: bin-wise offset/correlation + a Fig-4c-style eps(time,depth)
section (ours vs paper)."""
import sys, warnings
from pathlib import Path
import numpy as np
import xarray as xr
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
warnings.simplefilter("ignore")

OURS = "/Users/drew/TLC/2023/Wirewalker/Nortek/processed/TLC_23_sig1000_turb.nc"
REF = ("/Users/drew/TLC/2023/Wirewalker/Nortek/"
       "doi_10_5061_dryad_8sf7m0d44__v20260319/DA_final/DA_final/Data/NortekTurbulenceData.nc")

ours = xr.open_dataset(OURS)
ref = xr.open_dataset(REF)
print(f"ours: {ours.sizes['cast']} casts x {ours.sizes['depth']} depths")
print(f"ref : {ref.sizes['time']} casts x {ref.sizes['depth']} depths")

# match each of our casts to the nearest published profile in time
t_ours = ours["time"].values.astype("datetime64[s]").astype("int64")
t_ref = ref["time"].values
z = ours["depth"].values
e_ours = ours["epsilon"].values                 # (depth, cast)
e_ref_all = ref["Epsilon"].values               # (time, depth)

idx = np.array([int(np.argmin(np.abs(t_ref - t))) for t in t_ours])
dt = np.abs(t_ref[idx] - t_ours)
print(f"time match: median dt={np.median(dt):.0f}s  max dt={dt.max():.0f}s")
e_ref = e_ref_all[idx].T                          # (depth, cast) aligned to ours

m = np.isfinite(e_ours) & np.isfinite(e_ref) & (e_ours > 0) & (e_ref > 0)
d = np.log10(e_ours[m]) - np.log10(e_ref[m])
r2 = np.corrcoef(np.log10(e_ours[m]), np.log10(e_ref[m]))[0, 1]
print(f"\nALL matched bins: n={m.sum()} median offset={np.median(d):+.3f} dex "
      f"rms={np.std(d):.3f} corr(log)={r2:.3f}")
# depth-mean profiles
with np.errstate(invalid="ignore"):
    lo = np.nanmean(np.where(e_ours > 0, np.log10(e_ours), np.nan), axis=1)
    lr = np.nanmean(np.where(e_ref > 0, np.log10(e_ref), np.nan), axis=1)
print("\n depth  <log10 eps> ours / ref:")
for k in range(0, len(z), 6):
    print(f"  {z[k]:5.1f}   {lo[k]:+.2f}   {lr[k]:+.2f}")

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
td = ours["time"].values
fig, ax = plt.subplots(3, 1, figsize=(13, 11), sharex=True)
vmin, vmax = -9, -5
pc = ax[0].pcolormesh(td, z, np.log10(e_ours), vmin=vmin, vmax=vmax, cmap="inferno", shading="nearest")
ax[0].set_title("ours (ww_sig1000 port)  log10 eps"); ax[0].invert_yaxis(); ax[0].set_ylabel("depth (m)")
fig.colorbar(pc, ax=ax[0], label="log10 eps")
pc = ax[1].pcolormesh(td, z, np.log10(e_ref), vmin=vmin, vmax=vmax, cmap="inferno", shading="nearest")
ax[1].set_title("paper (NortekTurbulenceData.nc)  log10 eps"); ax[1].invert_yaxis(); ax[1].set_ylabel("depth (m)")
fig.colorbar(pc, ax=ax[1], label="log10 eps")
diff = np.log10(e_ours) - np.log10(e_ref)
pc = ax[2].pcolormesh(td, z, diff, vmin=-1, vmax=1, cmap="RdBu_r", shading="nearest")
ax[2].set_title(f"difference (ours - paper), dex   [median {np.median(d):+.2f}, corr {r2:.2f}]")
ax[2].invert_yaxis(); ax[2].set_ylabel("depth (m)"); ax[2].set_xlabel("time")
fig.colorbar(pc, ax=ax[2], label="dex")
fig.tight_layout()
fig.savefig("/Users/drew/PyWirewalker/TLC_turb_deployment.png", dpi=110)
print("\nwrote TLC_turb_deployment.png")
