"""
Figures for the December 2025 verification at the M2 tower.

Color = model family, line style = ensemble (solid) vs deterministic (dashed).
Usage: python plot_scores.py 2025-12-01 2025-12-31
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

D = Path(__file__).parent / "data" / "processed"
FIG = Path(__file__).parent / "figures"
FIG.mkdir(exist_ok=True)

STYLE = {  # model: (label, color, linestyle, marker)
    "IFS_ENS":     ("ECMWF IFS ENS (physics, ens.)",    "#2a78d6", "-",  "o"),
    "IFS_Control": ("ECMWF IFS control (physics, det.)", "#2a78d6", "--", None),
    "AIFS_ENS":    ("ECMWF AIFS ENS (AI, ens.)",        "#eb6834", "-",  "o"),
    "AIFS_Single": ("ECMWF AIFS Single (AI, det.)",     "#eb6834", "--", None),
    "WN2":         ("Google WeatherNext 2 (AI, ens.)",  "#1baf7a", "-",  "o"),
    "GFS":         ("NOAA GFS (physics, det.)",         "#8a8a83", "--", None),
    "HRRR":        ("NOAA HRRR 3 km (physics, det.)",   "#e87ba4", "--", "s"),
}
INK, MUTED, GRID = "#1a1a19", "#6b6a63", "#e4e3dd"


def style_ax(ax):
    ax.grid(color=GRID, lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=9)


def skill_figure(start, end):
    s = pd.read_parquet(D / f"scores_hub_{start:%Y%m%d}_{end:%Y%m%d}.parquet")
    s["bias"] = s.fc_mean - s.obs
    s["lead_day"] = ((s.lead_h - 1) // 24 + 1).astype(int)              # day 1 = leads 6-24 h, ...
    s["hour_mst"] = ((s.init_time + pd.to_timedelta(s.lead_h, unit="h")).dt.hour - 7) % 24
    by_day = s.groupby(["model", "lead_day"]).crps.mean()
    by_hour = s.groupby(["model", "hour_mst"]).bias.mean()
    fig, (a, b) = plt.subplots(1, 2, figsize=(14, 5.2))
    for m, (lab, c, ls, mk) in STYLE.items():
        if m not in by_day.index.get_level_values(0):
            continue
        a.plot(by_day.loc[m].index, by_day.loc[m], ls=ls, color=c, lw=2, marker=mk or "o", ms=5, label=lab)
        h = by_hour.loc[m]
        b.plot(h.index, h, ls=ls, color=c, lw=2, marker=mk or "o", ms=5, label=lab)
    a.set_ylabel("CRPS (ensembles) / MAE (deterministic)  [m/s]", color=MUTED)
    a.set_title("Error by lead time, 80 m (models interpolated from 100 m where needed)", color=INK, fontsize=11)
    a.set_xlabel("Lead time (day)", color=MUTED)
    b.axhline(0, color=INK, lw=0.8)
    b.set_ylabel("Bias: forecast (ens. mean) − observed  [m/s]", color=MUTED)
    b.set_title("Bias by time of day, hub height (all lead times)", color=INK, fontsize=11)
    b.set_xlabel("Time of day (MST)", color=MUTED)
    b.set_xticks([5, 11, 17, 23])
    for ax in (a, b):
        style_ax(ax)
    a.set_ylim(bottom=0)
    a.legend(fontsize=8.5, frameon=False, loc="lower left", handlelength=3.5)
    n = s.init_time.nunique()
    fig.suptitle(f"Wind forecasts at the NREL M2 tower near Boulder, {start:%B %Y} ({n} runs at 00 UTC); lower is better",
                 color=INK, fontsize=12)
    fig.tight_layout()
    f = FIG / f"skill_vs_lead_{start:%Y%m}.png"
    fig.savefig(f, dpi=150)
    print("saved", f)


def timeseries_figure(start, end):
    """Observed 80 m wind vs. day-1 forecasts. Every valid time comes from exactly one run:
    the 00 UTC run of the same UTC day (leads 6-24 h for the global ensembles, 1-24 h hourly for HRRR)."""
    from verify import load_obs, load_forecasts
    obs = load_obs()
    fc = load_forecasts(start, end)                               # 80 m (interpolated from 100 m where needed)
    fc = fc[fc.lead_h.between(6, 24)]
    t0, t1 = start, end + pd.Timedelta(days=1)
    fig, ax = plt.subplots(figsize=(15, 4.8))
    o = obs.obs_hub[t0:t1]
    ax.plot(o.index, o, color=INK, lw=1.2, label="Observed, M2 tower 80 m (hourly)")
    for m in ["IFS_ENS", "WN2"]:
        lab, c, _, _ = STYLE[m]
        for init, g in fc[fc.model == m].groupby("init_time"):
            q = g.groupby("valid_time").ws_hub.quantile([0.1, 0.5, 0.9]).unstack()
            ax.fill_between(q.index, q[0.1], q[0.9], color=c, alpha=0.2, lw=0)
            ax.plot(q.index, q[0.5], color=c, lw=1.6, marker="o", ms=3,
                    label=f"{lab.split(' (')[0]}: median and 10–90 %" if init == fc.init_time.min() else None)
    hb = D / f"hrrr_box_{start:%Y%m%d}_{end:%Y%m%d}.parquet"
    if hb.exists():
        h = pd.read_parquet(hb)
        h = h[(h.dy == 0) & (h.dx == 0) & (h.init_time.dt.hour == 0) & h.lead_h.between(1, 24)]
        lab, c, _, _ = STYLE["HRRR"]
        for init, g in h.groupby("init_time"):
            ax.plot(g.valid_time, g.ws, color=c, lw=1.4,
                    label="NOAA HRRR 3 km (hourly)" if init == h.init_time.min() else None)
    ax.set_xlim(t0, t1)
    ax.set_ylabel("Wind speed at 80 m [m/s]", color=MUTED)
    ax.set_title(f"Day-1 forecasts from each day's 00 UTC run (one run per day, no averaging) vs. observations, {start:%B %Y}",
                 color=INK, fontsize=11)
    style_ax(ax)
    ax.legend(fontsize=8.5, frameon=False, loc="upper left", ncol=2)
    fig.tight_layout()
    f = FIG / f"timeseries_day1_{start:%Y%m}.png"
    fig.savefig(f, dpi=150)
    print("saved", f)


if __name__ == "__main__":
    start, end = pd.Timestamp(sys.argv[1]), pd.Timestamp(sys.argv[2])
    skill_figure(start, end)
    timeseries_figure(start, end)
