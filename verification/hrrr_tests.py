"""
Extra tests: does HRRR (3 km) really lose against 25 km global models at the M2 tower?

A. Latest HRRR run: every valid hour, the most recent HRRR cycle (00/06/12/18 UTC) with lead 1-6 h.
B. Daily maximum hub-height wind (MST day, day-1 forecast from the 00 UTC run): windstorm strength,
   independent of the exact hour. CRPS for ensembles, MAE for deterministic models.
C. HRRR neighbourhood ensemble: 5x5 grid points (about +-6 km) x +-1 h = 75 "members" -> CRPS,
   the standard way to verify convection-allowing models without the double penalty.

Usage: python hrrr_tests.py 2025-12-01 2025-12-31
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import dynamical_catalog as dc

from verify import fair_crps, load_obs, load_forecasts

warnings.filterwarnings("ignore")
LAT, LON = 39.9106, -105.2347
D = Path(__file__).parent / "data" / "processed"


def hrrr_box(start, end, half=2):
    """HRRR 80 m wind, all 4 daily cycles, leads 1-48 h, (2*half+1)^2 grid points around the tower."""
    f = D / f"hrrr_box_{start:%Y%m%d}_{end:%Y%m%d}.parquet"
    if f.exists():
        return pd.read_parquet(f)
    ds = dc.open("noaa-hrrr-forecast-48-hour")
    lat2, lon2 = ds.latitude.values, ds.longitude.values
    iy, ix = np.unravel_index(np.argmin((lat2 - LAT) ** 2 + ((lon2 - LON) * np.cos(np.deg2rad(LAT))) ** 2), lat2.shape)
    ydim, xdim = ds.latitude.dims
    sub = ds.isel({ydim: slice(iy - half, iy + half + 1), xdim: slice(ix - half, ix + half + 1)})
    inits = pd.date_range(start, end + pd.Timedelta(hours=18), freq="6h")
    sub = sub.sel(init_time=inits, lead_time=pd.to_timedelta(np.arange(1, 49), unit="h"))
    sub = sub[["wind_u_80m", "wind_v_80m"]].load()
    ws = np.hypot(sub.wind_u_80m, sub.wind_v_80m).rename("ws")
    df = ws.to_dataframe().reset_index()[["init_time", "lead_time", ydim, xdim, "ws"]]
    df["dy"] = df[ydim].rank(method="dense").astype(int) - 1 - half      # grid index offsets from the centre
    df["dx"] = df[xdim].rank(method="dense").astype(int) - 1 - half
    df["lead_h"] = (df.lead_time / pd.Timedelta(hours=1)).astype(int)
    df["valid_time"] = df.init_time + df.lead_time
    df = df[["init_time", "lead_h", "valid_time", "dy", "dx", "ws"]]
    df.to_parquet(f)
    return df


def gfs_hourly(start, end):
    """GFS 80 m wind, 00 UTC runs, hourly leads 1-48 h (nearest point)."""
    f = D / f"gfs_hourly_{start:%Y%m%d}_{end:%Y%m%d}.parquet"
    if f.exists():
        return pd.read_parquet(f)
    ds = dc.open("noaa-gfs-forecast")
    p = ds.sel(latitude=LAT, longitude=LON, method="nearest")
    p = p.sel(init_time=pd.date_range(start, end), lead_time=pd.to_timedelta(np.arange(1, 49), unit="h"))
    p = p[["wind_u_80m", "wind_v_80m"]].load()
    df = np.hypot(p.wind_u_80m, p.wind_v_80m).rename("ws").to_series().reset_index()
    df["lead_h"] = (df.lead_time / pd.Timedelta(hours=1)).astype(int)
    df["valid_time"] = df.init_time + df.lead_time
    df = df[["init_time", "lead_h", "valid_time", "ws"]]
    df.to_parquet(f)
    return df


def stats(fc, ob):
    e = fc - ob
    return pd.Series({"MAE": np.abs(e).mean(), "bias": e.mean(), "corr": np.corrcoef(fc, ob)[0, 1], "n": len(e)})


if __name__ == "__main__":
    start, end = pd.Timestamp(sys.argv[1]), pd.Timestamp(sys.argv[2])
    obs = load_obs().obs_hub
    box = hrrr_box(start, end)
    centre = box[(box.dy == 0) & (box.dx == 0)].copy()
    centre["obs"] = obs.reindex(centre.valid_time).to_numpy()
    centre = centre.dropna()
    month = (centre.valid_time >= start) & (centre.valid_time < end + pd.Timedelta(days=1))

    # ---- A: latest run vs 00 UTC run, hourly valid times
    print("\nA. HRRR at the nearest grid point, every hour of the month")
    latest = centre[month & centre.lead_h.between(1, 6)]
    run00 = centre[month & (centre.init_time.dt.hour == 0)]
    rows = {"latest run (lead 1-6 h)": stats(latest.ws, latest.obs)}
    for a, b in [(1, 6), (7, 12), (13, 24), (25, 48)]:
        g = run00[run00.lead_h.between(a, b)]
        rows[f"00 UTC run, lead {a}-{b} h"] = stats(g.ws, g.obs)
    gfs = gfs_hourly(start, end)
    gfs["obs"] = obs.reindex(gfs.valid_time).to_numpy()
    gfs = gfs.dropna()
    for a, b in [(1, 6), (7, 12), (13, 24), (25, 48)]:
        g = gfs[gfs.lead_h.between(a, b)]
        rows[f"GFS 00 UTC run, lead {a}-{b} h"] = stats(g.ws, g.obs)
    print(pd.DataFrame(rows).T.round(2).to_string())


    # ---- B: daily maximum hub-height wind, MST day D, from the 00 UTC run of day D (leads 7-30 h)
    print("\nB. Daily maximum 80 m wind (MST day, 00 UTC run of the same day, leads 7-30 h)")
    mst = lambda t: (t - pd.Timedelta(hours=7)).dt.floor("D")
    o = obs[start + pd.Timedelta(hours=7): end + pd.Timedelta(hours=31)]
    obs_max = o.groupby((o.index - pd.Timedelta(hours=7)).floor("D")).max()
    res = {}
    day1 = lambda df: df[df.lead_h.between(7, 30)].assign(day=lambda d: mst(d.valid_time))
    h = day1(run00).groupby("day").ws.max()
    res["HRRR (hourly)"] = h
    h6 = day1(run00[run00.lead_h % 6 == 0]).groupby("day").ws.max()
    res["HRRR (6-hourly samples)"] = h6
    res["GFS (hourly)"] = day1(gfs).groupby("day").ws.max()
    fc = load_forecasts(start, end)
    out = {}
    for name, s in res.items():
        ob = obs_max.reindex(s.index)
        ok = ob.notna()
        out[name] = pd.Series({"CRPS/MAE": np.abs(s - ob)[ok].mean(), "bias": (s - ob)[ok].mean(), "n": ok.sum()})
    for m in ["AIFS_Single", "IFS_Control", "IFS_ENS", "AIFS_ENS", "WN2"]:
        g = day1(fc[fc.model == m])
        mx = g.groupby(["day", "member"]).ws_hub.max().unstack()          # daily max per member
        ob = obs_max.reindex(mx.index)
        ok = ob.notna().to_numpy()
        crps = fair_crps(mx.to_numpy()[ok], ob.to_numpy()[ok])
        steps = sorted(g.lead_h.unique())
        out[f"{m} ({len(steps)} steps/day)"] = pd.Series({"CRPS/MAE": crps.mean(),
                                                          "bias": (mx.mean(axis=1) - ob)[ok].mean(), "n": ok.sum()})
    print(pd.DataFrame(out).T.round(2).sort_values("CRPS/MAE").to_string())


    # ---- C: HRRR neighbourhood ensemble (5x5 points x +-1 h) vs. global ensembles, same valid times
    print("\nC. CRPS at 6-hourly valid times, 00 UTC runs, leads 6-48 h")
    r00 = box[box.init_time.dt.hour == 0]
    piv = r00.pivot_table(index=["init_time", "lead_h"], columns=["dy", "dx"], values="ws")
    rows = []
    for (init, lead) in piv.index:
        if lead % 6 or lead < 6:
            continue
        mem = np.concatenate([piv.loc[(init, l)].to_numpy() for l in (lead - 1, lead, lead + 1) if (init, l) in piv.index])
        y = obs.get(init + pd.Timedelta(hours=lead), np.nan)
        if np.isnan(y):
            continue
        rows.append(dict(init_time=init, lead_h=lead, crps_nbh=fair_crps(mem[None, :], np.array([y]))[0],
                         ae_centre=abs(piv.loc[(init, lead), (0, 0)] - y)))
    nb = pd.DataFrame(rows)
    sc = pd.read_parquet(D / f"scores_hub_{start:%Y%m%d}_{end:%Y%m%d}.parquet")
    sc = sc[sc.lead_h <= 48]
    keys = nb[["init_time", "lead_h"]]
    tab = {"HRRR neighbourhood ens. (75)": nb.crps_nbh.mean(), "HRRR nearest point": nb.ae_centre.mean()}
    for m in ["WN2", "IFS_ENS", "AIFS_ENS", "GFS", "AIFS_Single", "IFS_Control"]:
        tab[m] = keys.merge(sc[sc.model == m], on=["init_time", "lead_h"]).crps.mean()
    print(pd.Series(tab, name="CRPS / MAE").round(2).sort_values().to_string())
