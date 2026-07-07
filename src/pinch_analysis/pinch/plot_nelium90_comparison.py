#!/usr/bin/env python
"""
plot_nelium90_comparison.py
---------------------------
Two-panel composite curve comparison:
  (a) Baseline He-Ne (80 % He / 20 % Ne, 3 bar)
  (b) Nelium-90     (10 % He / 90 % Ne, 10 bar, Wilhelmsen et al.)

Reads:
    nelium90_composite_baseline.csv
    nelium90_composite_nelium90.csv
    nelium90_summary.csv

Output: nelium90_comparison.pdf / .png
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

from PLOT_SETTINGS import (
    graphic_font,
    math_font,
    spine_width,
    tick_width,
    tick_length,
    minor_tick_width,
    minor_tick_length,
    tick_labelsize,
    style_legend,
    save_figure,
)

here = Path(__file__).parent

# ---- Load data ----------------------------------------------------------
def load_curves(csv_path):
    df = pd.read_csv(csv_path)
    hm = ~df["Q_hot_kW"].isna()
    cm = ~df["Q_cold_kW"].isna()
    return (df.loc[hm, "Q_hot_kW"].values,  df.loc[hm, "T_hot_K"].values,
            df.loc[cm, "Q_cold_kW"].values, df.loc[cm, "T_cold_K"].values)

Q_hot_b, T_hot_b, Q_cold_b, T_cold_b = load_curves(
    here / "nelium90_composite_baseline.csv")
Q_hot_n, T_hot_n, Q_cold_n, T_cold_n = load_curves(
    here / "nelium90_composite_nelium90.csv")

# Force baseline cold composite to end at the same Q as the hot stream
Q_hot_end_b = Q_hot_b[-1]
if Q_cold_b[-1] > Q_hot_end_b:
    # Truncate cold composite at Q_hot_end
    mask = Q_cold_b <= Q_hot_end_b
    T_end = np.interp(Q_hot_end_b, Q_cold_b, T_cold_b)
    Q_cold_b = np.append(Q_cold_b[mask], Q_hot_end_b)
    T_cold_b = np.append(T_cold_b[mask], T_end)
elif Q_cold_b[-1] < Q_hot_end_b:
    # Extend cold composite to Q_hot_end (extrapolate)
    T_end = np.interp(Q_hot_end_b, Q_cold_b, T_cold_b)
    Q_cold_b = np.append(Q_cold_b, Q_hot_end_b)
    T_cold_b = np.append(T_cold_b, T_end)

# Load summary for VLE info
summary = pd.read_csv(here / "nelium90_summary.csv")
nel_row = summary[summary["case"] == "nelium90"].iloc[0]
T_bub = nel_row["T_bub_K"]
T_dew = nel_row["T_dew_K"]

# ---- Pinch computation --------------------------------------------------
def compute_pinch(Q_hot, T_hot, Q_cold, T_cold, n=4001):
    Q_max = min(Q_hot[-1], Q_cold[-1])
    Qg = np.linspace(0, Q_max, n)
    Th = np.interp(Qg, Q_hot,  T_hot)
    Tc = np.interp(Qg, Q_cold, T_cold)
    dT = Th - Tc
    i = int(np.argmin(dT))
    return dT[i], Qg[i], Th[i], Tc[i]

pinch_b = compute_pinch(Q_hot_b, T_hot_b, Q_cold_b, T_cold_b)
pinch_n = compute_pinch(Q_hot_n, T_hot_n, Q_cold_n, T_cold_n)

print(f"Baseline:  dT_min = {pinch_b[0]:.3f} K  @ Q = {pinch_b[1]:.1f} kW")
print(f"Nelium-90: dT_min = {pinch_n[0]:.3f} K  @ Q = {pinch_n[1]:.1f} kW")

# ---- Axis styling -------------------------------------------------------
def style_axis(ax):
    for sp in ax.spines.values():
        sp.set_linewidth(spine_width)
    ax.tick_params(axis="both", which="major", direction="in",
                   width=tick_width, length=tick_length,
                   labelsize=tick_labelsize,
                   bottom=True, top=True, left=True, right=True)
    ax.tick_params(axis="both", which="minor", direction="in",
                   width=minor_tick_width, length=minor_tick_length,
                   bottom=True, top=True, left=True, right=True)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))


def draw_arrow(ax, Q, T, frac, color, direction="forward"):
    n = len(Q)
    if n < 4:
        return
    i = int(np.clip(frac * n, 2, n - 3))
    di = max(1, n // 200)
    j_tail = i
    j_head = min(n - 1, i + di) if direction == "forward" else max(0, i - di)
    ax.annotate("", xy=(Q[j_head], T[j_head]),
                xytext=(Q[j_tail], T[j_tail]),
                arrowprops=dict(arrowstyle="->", color=color,
                                lw=1.2, mutation_scale=12),
                zorder=5)


# ---- Colors -------------------------------------------------------------
hot_color  = (0.75, 0.10, 0.10)
cold_color = (0.10, 0.25, 0.75)

# ---- Common limits ------------------------------------------------------
T_lo = min(T_hot_b.min(), T_cold_b.min(), T_hot_n.min(), T_cold_n.min())
T_hi = max(T_hot_b.max(), T_cold_b.max(), T_hot_n.max(), T_cold_n.max())
Q_max = max(Q_hot_b.max(), Q_cold_b.max(), Q_hot_n.max(), Q_cold_n.max())
y_lim = (np.floor(T_lo) - 4, np.ceil(T_hi) + 6)
x_lim = (0.0, 1.06 * Q_max)

# ---- Figure -------------------------------------------------------------
with plt.style.context(["ieee"]):
    plt.rcParams["font.family"]      = graphic_font
    plt.rcParams["mathtext.fontset"] = math_font
    plt.rcParams["text.usetex"]      = True

    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(8.5, 3.6), sharey=False,
        gridspec_kw={"wspace": 0.28},
    )

style_axis(ax_a)
style_axis(ax_b)

# =========================================================================
# Panel (a): Baseline
# =========================================================================
ax_a.plot(Q_hot_b, T_hot_b, color=hot_color, lw=1.4, ls="-",
          label=r"Hot H$_2$ stream")
ax_a.plot(Q_cold_b, T_cold_b, color=cold_color, lw=1.4, ls="-",
          label=r"Cold composite (He-Ne 80/20)")

# Pinch marker
dT, Qp, Thp, Tcp = pinch_b
ax_a.plot([Qp, Qp], [Tcp, Thp], "k-", lw=1.2)
ax_a.text(Qp + 0.03* x_lim[1], Tcp - 2,
          rf"$\Delta T_{{\min}}={dT:.1f}\;$K",
          ha="left", va="center", fontsize=9)

# Arrows
draw_arrow(ax_a, Q_hot_b,  T_hot_b,  0.55, hot_color,  "backward")
draw_arrow(ax_a, Q_cold_b, T_cold_b, 0.40, cold_color, "forward")

# Endpoint labels
ax_a.annotate(r"$T_{h,o}$", xy=(Q_hot_b[0], T_hot_b[0]),
              xytext=(18, 19), textcoords="offset points",
              fontsize=8, color=hot_color, ha="right", va="bottom")
ax_a.annotate(r"$T_{h,i}$", xy=(Q_hot_b[-1], T_hot_b[-1]),
              xytext=(-16, 3), textcoords="offset points",
              fontsize=8, color=hot_color, ha="left", va="bottom")
ax_a.annotate(r"$T_{c,i}$", xy=(Q_cold_b[0], T_cold_b[0]),
              xytext=(18, 1), textcoords="offset points",
              fontsize=8, color=cold_color, ha="right", va="top")
ax_a.annotate(r"$T_{c,o}$", xy=(Q_cold_b[-1], T_cold_b[-1]),
              xytext=(-15, -15), textcoords="offset points",
              fontsize=8, color=cold_color, ha="left", va="top")

# Grey dashed endpoint lines
Q_end_b = Q_hot_b[-1]
ax_a.axvline(Q_end_b, color="grey", ls="--", lw=0.8, zorder=1)

# =========================================================================
# Panel (b): Nelium-90
# =========================================================================
ax_b.plot(Q_hot_n, T_hot_n, color=hot_color, lw=1.4, ls="-",
          label=r"Hot H$_2$ stream")
ax_b.plot(Q_cold_n, T_cold_n, color=cold_color, lw=1.4, ls="-",
          label=r"Cold composite (He-Ne 10/90)")

# Shade two-phase region
if not (np.isnan(T_bub) or np.isnan(T_dew)):
    Q_bub = np.interp(T_bub, T_cold_n, Q_cold_n)
    Q_dew = np.interp(T_dew, T_cold_n, Q_cold_n)
    mask_2ph = (Q_cold_n >= Q_bub) & (Q_cold_n <= Q_dew)
    if np.any(mask_2ph):
        ax_b.fill_between(
            Q_cold_n[mask_2ph], T_cold_n[mask_2ph],
            y2=y_lim[0],
            color=cold_color, alpha=0.08, zorder=1,
        )
        # Bubble and dew point markers
        ax_b.plot(Q_bub, T_bub, "o", ms=5,
                  mfc="white", mec=cold_color, mew=1.2, zorder=6)
        ax_b.plot(Q_dew, T_dew, "o", ms=5,
                  mfc="white", mec=cold_color, mew=1.2, zorder=6)
        ax_b.annotate(rf"$T_{{\mathrm{{bub}}}}={T_bub:.1f}\;$K",
                      xy=(Q_bub, T_bub),
                      xytext=(10, -14), textcoords="offset points",
                      fontsize=7.5, color=cold_color, ha="left", va="top")
        ax_b.annotate(rf"$T_{{\mathrm{{dew}}}}={T_dew:.1f}\;$K",
                      xy=(Q_dew, T_dew),
                      xytext=(10, 4), textcoords="offset points",
                      fontsize=7.5, color=cold_color, ha="left", va="bottom")

# Pinch marker
dT, Qp, Thp, Tcp = pinch_n
ax_b.plot([Qp, Qp], [Tcp, Thp], "k-", lw=1.2)
ax_b.text(Qp + 0.03* x_lim[1], Tcp - 2,
          rf"$\Delta T_{{\min}}={dT:.1f}\;$K",
          ha="left", va="center", fontsize=9)

# Arrows
draw_arrow(ax_b, Q_hot_n,  T_hot_n,  0.55, hot_color,  "backward")
draw_arrow(ax_b, Q_cold_n, T_cold_n, 0.40, cold_color, "forward")

# Endpoint labels
ax_b.annotate(r"$T_{h,o}$", xy=(Q_hot_n[0], T_hot_n[0]),
              xytext=(18, 19), textcoords="offset points",
              fontsize=8, color=hot_color, ha="right", va="bottom")
ax_b.annotate(r"$T_{h,i}$", xy=(Q_hot_n[-1], T_hot_n[-1]),
              xytext=(-16, 3), textcoords="offset points",
              fontsize=8, color=hot_color, ha="left", va="bottom")
ax_b.annotate(r"$T_{c,i}$", xy=(Q_cold_n[0], T_cold_n[0]),
              xytext=(18, 1), textcoords="offset points",
              fontsize=8, color=cold_color, ha="right", va="top")
ax_b.annotate(r"$T_{c,o}$", xy=(Q_cold_n[-1], T_cold_n[-1]),
              xytext=(-15, 10), textcoords="offset points",
              fontsize=8, color=cold_color, ha="left", va="top")

# Grey dashed endpoint lines
Q_end_n = Q_hot_n[-1]
ax_b.axvline(Q_end_n, color="grey", ls="--", lw=0.8, zorder=1)

# ---- Formatting ---------------------------------------------------------
for ax in (ax_a, ax_b):
    ax.set_xlim(*x_lim)
    ax.set_ylim(*y_lim)
    ax.set_xlabel(r"$\dot{Q}\;/\;[\mathrm{kW}]$")
    ax.set_ylabel(r"$T\;/\;[\mathrm{K}]$")

# Panel labels
ax_a.text(-0.17, 1.02, r"\textbf{(a)}", transform=ax_a.transAxes,
          ha="left", va="bottom", fontsize=11)
ax_b.text(-0.17, 1.02, r"\textbf{(b)}", transform=ax_b.transAxes,
          ha="left", va="bottom", fontsize=11)

# Legends
style_legend(ax_a, loc="upper left", frame=False)
style_legend(ax_b, loc="upper left", frame=False)

# ---- Save ---------------------------------------------------------------
save_figure(fig, str(here / "nelium90_comparison.png"))
save_figure(fig, str(here / "nelium90_comparison.pdf"))
print("Saved nelium90_comparison.pdf and .png")