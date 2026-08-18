"""Wire HR beam-5 turbulence to real TLC data and validate eps vs the paper's
published product NortekTurbulenceData.nc (2116 profiles x 67 depth bins, from the
Dryad DA_final set — the same grid ProcessSingleProfile.m produces).
(dev script — reads a profiling window, runs process_cast_turbulence per upcast.)"""
import sys, warnings
from pathlib import Path
import numpy as np
import xarray as xr
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from mhkit import dolfyn
from ww_sig1000.casts import detect_casts, upcasts
from ww_sig1000.turbulence import process_cast_turbulence
warnings.simplefilter("ignore")

FN = "/Users/drew/TLC/2023/Wirewalker/Nortek/S101913A008_TLC_23.ad2cp"
REF = ("/Users/drew/TLC/2023/Wirewalker/Nortek/"
       "doi_10_5061_dryad_8sf7m0d44__v20260319/DA_final/DA_final/Data/NortekTurbulenceData.nc")
START = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000
N = int(sys.argv[2]) if len(sys.argv) > 2 else 300_000

ds = dolfyn.read(FN, nens=[START, START + N])
a = ds.attrs
va = float(a["ambig_vel_b5"])
fs = float(a["fs"])
ts = ds["time"].values.astype("datetime64[ns]").astype("int64") / 1e9
press = ds["pressure"].values
casts = detect_casts(press, ts, thhold=int(0.5 * fs * 60))
ups = [c for c in upcasts(casts) if np.ptp(press[c.start:c.stop + 1]) > 0.2 * 93]
print(f"{len(ups)} substantial upcasts in window")

pitch = ds["pitch"].values; roll = ds["roll"].values
velb5 = ds["vel_b5"].values; corrb5 = ds["corr_b5"].values          # (range_b5, time_b5)
tarr = ds["time"].values

results = []
for c in ups:
    sl = slice(c.start, c.stop + 1)
    r = process_cast_turbulence(velb5[:, sl].T, corrb5[:, sl].T, press[sl],
                                pitch[sl], roll[sl], va,
                                cellsize=0.06, blockdis=0.1, dep_res=3.0, max_dep=100.0,
                                time=tarr[sl][(c.stop - c.start) // 2])
    if r is not None:
        results.append(r)
print(f"{len(results)} casts produced turbulence\n")

# --- compare to the published product ---
turb = xr.open_dataset(REF)
tt = turb["time"].values                         # unix seconds
z_ref = turb["depth"].values
eps_ref = turb["Epsilon"].values                 # (2116, 67)

OUR, REF_E = [], []
for r in results:
    t_unix = r["time"].astype("datetime64[s]").astype("int64")
    it = int(np.argmin(np.abs(tt - t_unix)))
    OUR.append(r["eps"]); REF_E.append(eps_ref[it])
OUR = np.concatenate(OUR); REF_E = np.concatenate(REF_E)
mm = np.isfinite(OUR) & np.isfinite(REF_E) & (OUR > 0) & (REF_E > 0)
off = np.median(np.log10(OUR[mm]) - np.log10(REF_E[mm]))
r2 = np.corrcoef(np.log10(OUR[mm]), np.log10(REF_E[mm]))[0, 1]
print(f"ALL matched bins: n={mm.sum()} median_offset={off:+.2f} dex corr(log)={r2:.3f}")
