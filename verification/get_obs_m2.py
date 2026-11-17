"""
Download NREL/NLR Flatirons Campus M2 tower data (1-min) and store hourly means.

Tower: 39.9106 N, -105.2347 W, 1855 m MSL, ~8 km south of Boulder. Wind at 2, 5, 10, 20, 50, 80 m.
Data API: https://midcdmz.nlr.gov/apps/data_api_doc.pl
Time stamps in the raw data are MST (UTC-7, no daylight saving); converted to UTC here.

Usage: python get_obs_m2.py 2025-11-01 2026-03-31
"""
import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import requests

OUT = Path(__file__).parent / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)

URL = "https://midcdmz.nlr.gov/apps/data_api.pl?site=NWTC&begin={b}&end={e}"
HEIGHTS = [80, 50, 20, 10]


def fetch(begin, end):
    r = requests.get(URL.format(b=begin.strftime("%Y%m%d"), e=end.strftime("%Y%m%d")), timeout=180)
    r.raise_for_status()
    df = pd.read_csv(io.StringIO(r.text))
    t = (pd.to_datetime(df["Year"].astype(str), format="%Y")
         + pd.to_timedelta(df["DOY"] - 1, unit="D")
         + pd.to_timedelta(df["MST"] // 100, unit="h")
         + pd.to_timedelta(df["MST"] % 100, unit="min"))
    out = pd.DataFrame(index=pd.DatetimeIndex(t + pd.Timedelta(hours=7), name="time"))   # MST -> UTC
    for h in HEIGHTS:
        out[f"ws{h}"] = df[f"Avg Wind Speed @ {h}m [m/s]"].values
        out[f"wd{h}"] = df[f"Avg Wind Direction @ {h}m [deg]"].values
    peak = [c for c in df.columns if c.startswith("Peak Wind Speed @ 80m")]
    out["ws80_peak"] = df[peak[0]].values if peak else np.nan
    return out


def to_hourly(df):
    """Hourly mean speed; direction as vector mean; peak = max 1-min peak."""
    res = {}
    for h in HEIGHTS:
        res[f"ws{h}"] = df[f"ws{h}"].resample("1h").mean()
        rad = np.deg2rad(df[f"wd{h}"])
        u = (-df[f"ws{h}"] * np.sin(rad)).resample("1h").mean()
        v = (-df[f"ws{h}"] * np.cos(rad)).resample("1h").mean()
        res[f"wd{h}"] = (np.rad2deg(np.arctan2(-u, -v)) + 360) % 360
    res["ws80_peak"] = df["ws80_peak"].resample("1h").max()
    res["n_min"] = df["ws80"].resample("1h").count()      # data completeness (minutes per hour)
    return pd.DataFrame(res)


if __name__ == "__main__":
    start, end = pd.Timestamp(sys.argv[1]), pd.Timestamp(sys.argv[2])
    chunks = []
    for b in pd.date_range(start, end, freq="7D"):
        e = min(b + pd.Timedelta(days=6), end)
        print(f"  {b.date()} - {e.date()}", flush=True)
        chunks.append(fetch(b, e))
    raw = pd.concat(chunks).sort_index()
    raw = raw[~raw.index.duplicated()]
    raw = raw.mask(raw < -900)                            # missing-value flags
    hourly = to_hourly(raw)
    f = OUT / f"m2_hourly_{start:%Y%m%d}_{end:%Y%m%d}.parquet"
    hourly.to_parquet(f)
    print("saved", f, hourly.shape)
