"""Cross-check ww_adcp.transforms.beam2enu (ported MATLAB Beam2ENU) against
mhkit/dolfyn's rotate2(beam->earth) on a real .ad2cp sample."""
import sys, warnings
import numpy as np
sys.path.insert(0, "/Users/drew/PyWirewalker")
from ww_adcp.transforms import beam2enu
from mhkit import dolfyn

FN = sys.argv[1] if len(sys.argv) > 1 else "/Users/drew/Downloads/S100430A002_M3_d2-004.ad2cp"
NENS = int(sys.argv[2]) if len(sys.argv) > 2 else 1500
warnings.simplefilter("ignore")

ds = dolfyn.read(FN, nens=NENS)
assert ds.attrs["coord_sys"] == "beam", f"expected beam coords, got {ds.attrs['coord_sys']}"

# --- ported toolbox path: beam velocities + attitude -> ENU ---
beam = ds["vel"].values                     # (dir=beam 1..4, range, time)
head = ds["heading"].values
pitch = ds["pitch"].values
roll = ds["roll"].values
print(f"pings={beam.shape[2]} cells={beam.shape[1]} | "
      f"|pitch|<={np.nanmax(np.abs(pitch)):.1f} deg  |roll|<={np.nanmax(np.abs(roll)):.1f} deg")
enu = beam2enu(beam, head, pitch, roll)     # (3, range, time): E, N, U

# --- dolfyn path: rotate2 beam -> earth ---
de = dolfyn.rotate2(ds, "earth", inplace=False)
dE = de["vel"].sel(dir="E").values
dN = de["vel"].sel(dir="N").values
dU = 0.5 * (de["vel"].sel(dir="U1").values + de["vel"].sel(dir="U2").values)

def stats(a, b, name):
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    d = a - b
    rms = np.sqrt(np.mean(d**2))
    r = np.corrcoef(a, b)[0, 1] if a.size > 2 else np.nan
    print(f"  {name}: n={a.size:>7d}  rms_diff={rms*100:7.3f} cm/s  "
          f"max|d|={np.max(np.abs(d))*100:7.3f} cm/s  corr={r:.5f}  "
          f"rms_signal={np.sqrt(np.mean(b**2))*100:.2f} cm/s")
    return rms, r

print("\nported beam2enu  vs  dolfyn rotate2('earth'):")
re = stats(enu[0], dE, "E")
rn = stats(enu[1], dN, "N")
ru = stats(enu[2], dU, "U")

# cross-terms to detect an E/N swap or 90-deg heading-offset mismatch
print("\ncross-correlations (detect axis swap / heading convention):")
def corr(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    return np.corrcoef(a[m], b[m])[0, 1]
print(f"  corr(portE, dolfynE)={corr(enu[0],dE):+.4f}   corr(portE, dolfynN)={corr(enu[0],dN):+.4f}")
print(f"  corr(portN, dolfynE)={corr(enu[1],dE):+.4f}   corr(portN, dolfynN)={corr(enu[1],dN):+.4f}")

ok = re[1] > 0.99 and rn[1] > 0.99 and ru[1] > 0.99
print("\nVERDICT:", "AGREE (corr>0.99 all comps)" if ok
      else "DIVERGENCE — investigate convention (see cross-correlations)")
