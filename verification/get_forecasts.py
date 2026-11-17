"""
Extract point forecasts at the M2 tower from dynamical.org cloud archives (Icechunk Zarr on AWS).

One parquet file per model, long format:
    init_time, lead_h, valid_time, member, ws10, ws_hub
ws_hub = 100 m wind for the global models, 80 m wind for HRRR (tower top is 80 m).
Deterministic models get member = 0.

Usage: python get_forecasts.py 2025-12-01 2025-12-31 [model ...]
"""
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import dynamical_catalog as dc

warnings.filterwarnings("ignore")

LAT, LON = 39.9106, -105.2347           # M2 tower
OUT = Path(__file__).parent / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)

LEADS_H = np.arange(6, 169, 6)          # 6 h ... 7 days, 6-hourly (common to all models)

MODELS = {   # name: (dataset id, hub-height wind level, ensemble?, max lead h)
    "AIFS_ENS":    ("ecmwf-aifs-ens-forecast",                   "100m", True,  168),
    "IFS_ENS":     ("ecmwf-ifs-ens-forecast-15-day-0-25-degree", "100m", True,  168),
    "WN2":         ("earth-engine: weathernext_2_0_0",           "100m", True,  168),   # see extract_wn2_ee
    "AIFS_Single": ("ecmwf-aifs-single-forecast",                "100m", False, 168),
    "GFS":         ("noaa-gfs-forecast",                         "100m", False, 168),
    "HRRR":        ("noaa-hrrr-forecast-48-hour",                "80m",  False, 48),
}


def nearest_point(ds):
    """Select the grid point nearest to the tower (regular lat/lon grid or 2-D projected grid)."""
    if "latitude" not in ds.coords and "y" in ds.dims and ds["y"].ndim == 1:   # regular grid stored as y/x (WeatherNext 2)
        ds = ds.rename(y="latitude", x="longitude")
    if "latitude" in ds.dims:
        lon = LON % 360 if float(ds.longitude.max()) > 180 else LON
        return ds.sel(latitude=LAT, longitude=lon, method="nearest")
    lat2, lon2 = ds.latitude.values, ds.longitude.values
    lon2 = np.where(lon2 > 180, lon2 - 360, lon2)
    d2 = (lat2 - LAT) ** 2 + ((lon2 - LON) * np.cos(np.deg2rad(LAT))) ** 2
    iy, ix = np.unravel_index(np.argmin(d2), d2.shape)
    ydim, xdim = ds.latitude.dims
    return ds.isel({ydim: iy, xdim: ix})


def extract_wn2_ee(start, end):
    """WeatherNext 2 ensemble via Earth Engine (server-side point extraction).
    The dynamical.org copy stores one global field per chunk, which makes point reads extremely slow."""
    import sys as _sys
    import ee
    _sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from access_keys import ee_project
    ee.Initialize(project=ee_project)
    coll = ee.ImageCollection("projects/gcp-public-data-weathernext/assets/weathernext_2_0_0")
    bands = ["10m_u_component_of_wind", "10m_v_component_of_wind", "100m_u_component_of_wind", "100m_v_component_of_wind"]
    pt = ee.Geometry.Point([LON, LAT])
    frames = []
    for init in pd.date_range(start, end, freq="1D"):
        sub = coll.filter(ee.Filter.date(init.strftime("%Y-%m-%dT%H:%M:%S"))) \
                  .filter(ee.Filter.lte("forecast_hour", int(LEADS_H.max()))).select(bands)
        r = sub.getRegion(pt, 27830).getInfo()
        d = pd.DataFrame(r[1:], columns=r[0])
        d["member"] = d["id"].str.split("_").str[-1].astype(int)
        d["valid_time"] = pd.to_datetime(d["id"].str.split("_").str[1], format="%Y%m%d%H%M")
        d["init_time"] = init
        frames.append(d)
        print(f"  WN2 {init.date()}: {len(d)} rows", flush=True)
    df = pd.concat(frames, ignore_index=True)
    df["lead_h"] = ((df.valid_time - df.init_time) / pd.Timedelta(hours=1)).astype(int)
    df = df[df.lead_h.isin(LEADS_H)]
    df["ws10"] = np.hypot(df[bands[0]], df[bands[1]])
    df["ws_hub"] = np.hypot(df[bands[2]], df[bands[3]])
    df["model"] = "WN2"
    df["grid_lat"], df["grid_lon"] = df.latitude.iloc[0], df.longitude.iloc[0]
    return df[["model", "init_time", "lead_h", "valid_time", "member", "ws10", "ws_hub", "grid_lat", "grid_lon"]]


def extract(model, start, end):
    if model == "WN2":
        return extract_wn2_ee(start, end)
    dsid, hub, is_ens, max_lead = MODELS[model]
    ds = dc.open(dsid)
    pt = nearest_point(ds)
    print(f"{model}: grid point {float(pt.latitude):.3f}, {float(pt.longitude):.3f}", flush=True)
    inits = pd.date_range(start, end, freq="1D") + pd.Timedelta(hours=0)       # 00 UTC runs
    leads = pd.to_timedelta(LEADS_H[LEADS_H <= max_lead], unit="h")
    pt = pt.sel(init_time=inits, lead_time=leads)
    pt = pt[[f"wind_u_10m", "wind_v_10m", f"wind_u_{hub}", f"wind_v_{hub}"]].load()
    df = pd.DataFrame({
        "ws10":   np.hypot(pt["wind_u_10m"], pt["wind_v_10m"]).to_series(),
        "ws_hub": np.hypot(pt[f"wind_u_{hub}"], pt[f"wind_v_{hub}"]).to_series(),
    }).reset_index()
    if "ensemble_member" in df:
        df = df.rename(columns={"ensemble_member": "member"})
    else:
        df["member"] = 0
    df["lead_h"] = (df["lead_time"] / pd.Timedelta(hours=1)).astype(int)
    df["valid_time"] = df["init_time"] + df["lead_time"]
    df["model"] = model
    df["grid_lat"], df["grid_lon"] = float(pt.latitude), float(pt.longitude)
    return df[["model", "init_time", "lead_h", "valid_time", "member", "ws10", "ws_hub", "grid_lat", "grid_lon"]]


if __name__ == "__main__":
    start, end = pd.Timestamp(sys.argv[1]), pd.Timestamp(sys.argv[2])
    models = sys.argv[3:] or list(MODELS)
    for m in models:
        t0 = time.time()
        try:
            df = extract(m, start, end)
            f = OUT / f"fc_{m}_{start:%Y%m%d}_{end:%Y%m%d}.parquet"
            df.to_parquet(f)
            print(f"  saved {f.name}: {len(df)} rows, {time.time()-t0:.0f} s", flush=True)
        except Exception as e:
            print(f"  {m} FAILED: {type(e).__name__}: {str(e)[:300]}", flush=True)
