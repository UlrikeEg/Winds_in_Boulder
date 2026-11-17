"""HRRR by lead hour (hourly), and strong-wind event detection vs. other deterministic models."""
import warnings; warnings.filterwarnings("ignore")
import glob
import numpy as np, pandas as pd, dynamical_catalog as dc
from verify import load_obs

LAT, LON = 39.9106, -105.2347
obs = load_obs().obs_hub
ds = dc.open("noaa-hrrr-forecast-48-hour")
lat2, lon2 = ds.latitude.values, ds.longitude.values
iy, ix = np.unravel_index(np.argmin((lat2 - LAT) ** 2 + ((lon2 - LON) * 0.77) ** 2), lat2.shape)
ydim, xdim = ds.latitude.dims
p = ds.isel({ydim: iy, xdim: ix}).sel(init_time=pd.date_range("2025-12-01", "2025-12-31"),
                                      lead_time=pd.to_timedelta(np.arange(1, 49), unit="h"))
p = p[["wind_u_80m", "wind_v_80m"]].load()
h = np.hypot(p.wind_u_80m, p.wind_v_80m).to_series().rename("fc").reset_index()
h["valid_time"] = h.init_time + h.lead_time
h["lead_h"] = h.lead_time / pd.Timedelta(hours=1)
h["obs"] = obs.reindex(h.valid_time).to_numpy()
h = h.dropna()
h["bin"] = pd.cut(h.lead_h, [0, 6, 12, 24, 36, 48])
print("HRRR hourly, by lead bin:")
print(h.groupby("bin").apply(lambda g: pd.Series({"MAE": (g.fc - g.obs).abs().mean(), "bias": (g.fc - g.obs).mean(),
                                                   "corr": g.fc.corr(g.obs), "n": len(g)})).round(2).to_string())

# event detection at 6-hourly leads <= 48 h, threshold 12 m/s, all deterministic models + ensemble probabilities
fc = pd.concat([pd.read_parquet(f) for f in glob.glob("data/processed/fc_*_20251201_20251231.parquet")])
fc = fc[fc.lead_h <= 48]
fc["obs"] = obs.reindex(fc.valid_time).to_numpy()
fc = fc.dropna(subset=["obs"])
thr = 12.0
rows = []
for m, g in fc.groupby("model"):
    gg = g.groupby(["init_time", "lead_h"]).agg(p=("ws_hub", lambda s: (s > thr).mean()), obs=("obs", "first"))
    ev = gg.obs > thr
    yes = gg.p >= 0.5 if g.member.nunique() == 1 else gg.p >= 0.2      # ensembles: >= 20 % of members
    rows.append(dict(model=m, events=int(ev.sum()), hits=int((yes & ev).sum()), false_alarms=int((yes & ~ev).sum()),
                     brier=((gg.p - ev) ** 2).mean()))
print(f"\nStrong wind (> {thr} m/s at hub height), 6-hourly valid times, leads 6-48 h:")
print(pd.DataFrame(rows).round(3).to_string(index=False))
