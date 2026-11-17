"""HRRR diagnostics: hourly leads, 7x7 neighbourhood around the tower, MAE per grid point."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, dynamical_catalog as dc
from verify import load_obs

LAT, LON = 39.9106, -105.2347
ds = dc.open("noaa-hrrr-forecast-48-hour")
lat2, lon2 = ds.latitude.values, ds.longitude.values
d2 = (lat2 - LAT) ** 2 + ((lon2 - LON) * np.cos(np.deg2rad(LAT))) ** 2
iy, ix = np.unravel_index(np.argmin(d2), d2.shape)
ydim, xdim = ds.latitude.dims
sub = ds.isel({ydim: slice(iy - 3, iy + 4), xdim: slice(ix - 3, ix + 4)})
inits = pd.date_range("2025-12-01", "2025-12-31", freq="1D")
sub = sub.sel(init_time=inits, lead_time=pd.to_timedelta(np.arange(1, 49), unit="h"))
sub = sub[["wind_u_80m", "wind_v_80m"]].load()
ws = np.hypot(sub.wind_u_80m, sub.wind_v_80m)
obs = load_obs().obs_hub
valid = (sub.init_time + sub.lead_time).values            # (init, lead)
o = obs.reindex(pd.DatetimeIndex(valid.ravel())).to_numpy().reshape(valid.shape)
ok = ~np.isnan(o)
mae = np.abs(ws.values - o[..., None, None])               # dims: init, lead, y, x
mae_map = np.nanmean(np.where(ok[..., None, None], mae, np.nan), axis=(0, 1))
bias_map = np.nanmean(np.where(ok[..., None, None], ws.values - o[..., None, None], np.nan), axis=(0, 1))
dist = np.sqrt(d2[iy - 3: iy + 4, ix - 3: ix + 4]) * 111
print("grid point MAE [m/s], nearest point in the centre (3,3); rows = y (south->north), cols = x (west->east)")
print(np.round(mae_map, 2))
print("bias [m/s]"); print(np.round(bias_map, 2))
print("distance to tower [km]"); print(np.round(dist, 1))
print("terrain? lat/lon of centre:", lat2[iy, ix], lon2[iy, ix])
c = ws.values[:, :, 3, 3]
print("centre point, hourly leads: MAE %.2f, bias %.2f, n=%d" % (np.nanmean(np.abs(c - o)[ok]), np.nanmean((c - o)[ok]), ok.sum()))
print("centre point, 6-hourly leads only: MAE %.2f" % np.nanmean(np.abs(c - o)[:, 5::6][ok[:, 5::6]]))
best = np.unravel_index(np.nanargmin(mae_map), mae_map.shape)
print("best point", best, "MAE %.2f" % mae_map[best], "dist %.1f km" % dist[best])
