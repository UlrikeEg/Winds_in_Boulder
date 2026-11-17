"""
Verify point wind forecasts against the M2 tower: CRPS (ensembles) and MAE (deterministic) by lead time.

CRPS uses the "fair" ensemble estimator, so ensembles of different size (51 vs 64 members) compare fairly.
For a single-value forecast CRPS equals the absolute error, so MAE and CRPS share one axis.

Usage: python verify.py 2025-12-01 2025-12-31
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

D = Path(__file__).parent / "data" / "processed"


def fair_crps(members, obs):
    """members: (n, M) array, obs: (n,) -> (n,) fair CRPS."""
    x = np.sort(members, axis=1)
    M = x.shape[1]
    term1 = np.abs(x - obs[:, None]).mean(axis=1)
    if M == 1:
        return term1
    k = np.arange(1, M + 1)
    pair_sum = 2 * ((2 * k - M - 1) * x).sum(axis=1)          # sum_{i,j} |x_i - x_j|
    return term1 - pair_sum / (2 * M * (M - 1))


def load_obs():
    obs = pd.concat([pd.read_parquet(f) for f in sorted(D.glob("m2_hourly_*.parquet"))]).sort_index()
    obs = obs[~obs.index.duplicated()]
    ok = obs.n_min >= 45
    # hourly means are labelled by hour start; the mean of the hours starting at t-1h and t is centred on t
    ws80 = obs.ws80.where(ok)
    ws10 = obs.ws10.where(ok)
    centred = lambda s: (s + s.shift(1)) / 2
    return pd.DataFrame({"obs_hub": centred(ws80), "obs10": centred(ws10)})


def score(fc, obs, var="hub"):
    """One row per (model, init_time, lead_h): CRPS, ensemble-mean error, and the forecast mean."""
    o = obs[f"obs_{var}" if var == "hub" else "obs10"]
    col = "ws_hub" if var == "hub" else "ws10"
    rows = []
    for (model, init, lead), g in fc.groupby(["model", "init_time", "lead_h"]):
        y = o.get(g["valid_time"].iloc[0], np.nan)
        if np.isnan(y):
            continue
        m = g[col].to_numpy()
        rows.append(dict(model=model, init_time=init, lead_h=lead, obs=y, fc_mean=m.mean(),
                         crps=fair_crps(m[None, :], np.array([y]))[0],
                         ae_mean=abs(m.mean() - y), spread=m.std(ddof=1) if len(m) > 1 else 0.0))
    return pd.DataFrame(rows)


def load_forecasts(start, end):
    fc = pd.concat([pd.read_parquet(f) for f in sorted(D.glob(f"fc_*_{start:%Y%m%d}_{end:%Y%m%d}.parquet"))])
    # Global models give 100 m wind; interpolate to the 80 m tower height with a power law fitted between
    # 10 m and 100 m for each forecast (alpha = ln(ws100/ws10)/ln(10)). HRRR and GFS provide 80 m directly.
    glob_ = ~fc.model.isin(["HRRR", "GFS"])                 # these two provide 80 m directly
    alpha = (np.log(fc.ws_hub.clip(lower=0.1) / fc.ws10.clip(lower=0.1)) / np.log(10)).clip(0, 0.6)
    fc["ws_100m_raw"] = fc.ws_hub
    fc.loc[glob_, "ws_hub"] = fc.ws10[glob_].clip(lower=0.1) * 8 ** alpha[glob_]
    # IFS ENS control member as the physics-based deterministic counterpart of AIFS Single
    ctrl = fc[(fc.model == "IFS_ENS") & (fc.member == 0)].assign(model="IFS_Control")
    return pd.concat([fc, ctrl], ignore_index=True)


if __name__ == "__main__":
    start, end = pd.Timestamp(sys.argv[1]), pd.Timestamp(sys.argv[2])
    obs, fc = load_obs(), load_forecasts(start, end)
    for var in ["hub", "10m"]:
        s = score(fc, obs, var)
        s.to_parquet(D / f"scores_{var}_{start:%Y%m%d}_{end:%Y%m%d}.parquet")
        s["bias"] = s.fc_mean - s.obs
        tab = s.groupby(["model", "lead_h"]).agg(crps=("crps", "mean"), ae_mean=("ae_mean", "mean"),
                                                bias=("bias", "mean"), spread=("spread", "mean"), n=("crps", "size"))
        print(f"\n=== {var}: mean CRPS (= MAE for deterministic), selected lead times ===")
        print(tab["crps"].unstack("model").loc[[6, 24, 48, 72, 120, 168]].round(2).to_string())
        print(f"\n=== {var}: mean bias (forecast mean - obs) ===")
        print(tab["bias"].unstack("model").loc[[6, 24, 48, 72, 120, 168]].round(2).to_string())
