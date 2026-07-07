#!/usr/bin/env python
"""
plot_nelium_2x1.py
------------------
1x2 panel composite curves for Nelium-90 (He-Ne 10/90, 3 bar) at
two hydrogen pressures:
  (a) H2 at 75 bar (baseline)
  (b) H2 at 21 bar (ejector, Pp=40)

Reads: nelium_case_B.csv, nelium_case_D.csv, nelium_4case_summary.csv
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

from PLOT_SETTINGS import (
    graphic_font, math_font, spine_width,
    tick_width, tick_length, minor_tick_width, minor_tick_length,
    tick_labelsize, style_legend, save_figure,
)

here = Path(__file__).parent
hot_color  = (0.75, 0.10, 0.10)
cold_color = (0.10, 0.25, 0.75)

# ---- Load ---------------------------------------------------------------
def load_curves(tag):
    df = pd.read_csv(here / f"nelium_case_{tag}.csv")
    hm = ~df["Q_hot_kW"].isna();  cm = ~df["Q_cold_kW"].isna()
    return (df.loc[hm, "Q_hot_kW"].values, df.loc[hm, "T_hot_K"].values,
            df.loc[cm, "Q_cold_kW"].values, df.loc[cm, "T_cold_K"].values)

data_B = load_curves("B")
data_D = load_curves("D")
summary = pd.read_csv(here / "nelium_4case_summary.csv")
sum_B = summary[summary["case"] == "B"].iloc[0]
sum_D = summary[summary["case"] == "D"].iloc[0]
T_bub = sum_B["T_bub_K"]
T_dew = sum_B["T_dew_K"]

# ---- Helpers -------------------------------------------------------------
def compute_pinch(Qh, Th, Qc, Tc, n=4001):
    Qm = min(Qh[-1], Qc[-1])
    Qg = np.linspace(0, Qm, n)
    dT = np.interp(Qg, Qh, Th) - np.interp(Qg, Qc, Tc)
    i = int(np.argmin(dT))
    return dT[i], Qg[i], np.interp(Qg[i], Qh, Th), np.interp(Qg[i], Qc, Tc)

def style_ax(ax):
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

def draw_arrow(ax, Q, T, frac, color, backward=False):
    n = len(Q)
    i = int(np.clip(frac * n, 2, n - 3))
    di = max(1, n // 200)
    i0, i1 = i, min(n - 1, i + di)
    if backward:
        i0, i1 = i1, i0
    ax.annotate("", xy=(Q[i1], T[i1]),
                xytext=(Q[i0], T[i0]),
                arrowprops=dict(arrowstyle="->", color=color,
                                lw=1.2, mutation_scale=12), zorder=5)

def force_endpoint(Qh, Qc, Tc):
    """Force cold composite to end at same Q as hot stream."""
    Q_end = Qh[-1]
    if Qc[-1] > Q_end:
        mask = Qc <= Q_end
        T_end = np.interp(Q_end, Qc, Tc)
        return np.append(Qc[mask], Q_end), np.append(Tc[mask], T_end)
    elif Qc[-1] < Q_end:
        T_end = np.interp(Q_end, Qc, Tc)
        return np.append(Qc, Q_end), np.append(Tc, T_end)
    return Qc, Tc

# ---- Limits --------------------------------------------------------------
all_Q = np.concatenate([data_B[0], data_B[2], data_D[0], data_D[2]])
all_T = np.concatenate([data_B[1], data_B[3], data_D[1], data_D[3]])
x_lim = (0, 1.06 * all_Q.max())
y_lim = (np.floor(all_T.min()) - 4, np.ceil(all_T.max()) + 6)

# ---- Figure --------------------------------------------------------------
with plt.style.context(["ieee"]):
    plt.rcParams["font.family"]      = graphic_font
    plt.rcParams["mathtext.fontset"] = math_font
    plt.rcParams["text.usetex"]      = True

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(8.5, 3.6),
                                      gridspec_kw={"wspace": 0.28})

panels = [
    ("B", ax_a, data_B, sum_B, r"$P_{\mathrm{H_2}} = 75$ bar, He-Ne 10/90 (10 bar)"),
    ("D", ax_b, data_D, sum_D, r"$P_{\mathrm{H_2}} = 21$ bar, He-Ne 10/90 (10 bar)"),
]

for tag, ax, dat, s, title in panels:
    style_ax(ax)
    Qh, Th, Qc, Tc = dat
    Qc, Tc = force_endpoint(Qh, Qc, Tc)

    ax.plot(Qh, Th, color=hot_color, lw=1.4, ls="-", label=r"Hot H$_2$ stream")
    ax.plot(Qc, Tc, color=cold_color, lw=1.4, ls="-", label=r"Cold composite (H$_2$ and He-Ne)")

    # Grey endpoint line
    ax.axvline(Qh[-1], color="grey", ls="--", lw=0.8, zorder=1)

    # T_R,o marker: open circle at the R-stream outlet on the cold composite
    T_Ro = 46.65
    Q_Ro = np.interp(T_Ro, Tc, Qc)
    ax.plot(Q_Ro, T_Ro, "o", ms=6, mfc="white", mec=cold_color, mew=1.3, zorder=6)
    ax.annotate(r"$T_{R,o}$", xy=(Q_Ro, T_Ro),
                xytext=(12, -14), textcoords="offset points",
                fontsize=8, color=cold_color, ha="left", va="top",
                arrowprops=dict(arrowstyle="-", color=cold_color, lw=0.7))

    # Shade two-phase region
    if not (np.isnan(T_bub) or np.isnan(T_dew)):
        try:
            Q_bub = np.interp(T_bub, Tc, Qc)
            Q_dew = np.interp(T_dew, Tc, Qc)
            mask_2ph = (Qc >= Q_bub) & (Qc <= Q_dew)
            if np.any(mask_2ph):
                ax.fill_between(Qc[mask_2ph], Tc[mask_2ph],
                                y2=y_lim[0], color=cold_color,
                                alpha=0.08, zorder=1)
                ax.plot(Q_bub, T_bub, "o", ms=5,
                        mfc="white", mec=cold_color, mew=1.2, zorder=6)
                ax.plot(Q_dew, T_dew, "o", ms=5,
                        mfc="white", mec=cold_color, mew=1.2, zorder=6)
                ax.annotate(rf"$T_{{\mathrm{{bub}}}}={T_bub:.1f}\;$K",
                            xy=(Q_bub, T_bub),
                            xytext=(10, -14), textcoords="offset points",
                            fontsize=7, color=cold_color, ha="left")
                ax.annotate(rf"$T_{{\mathrm{{dew}}}}={T_dew:.1f}\;$K",
                            xy=(Q_dew, T_dew),
                            xytext=(10, 4), textcoords="offset points",
                            fontsize=7, color=cold_color, ha="left")
        except Exception:
            pass

    # Pinch
    dT, Qp, Thp, Tcp_ = compute_pinch(Qh, Th, Qc, Tc)
    ax.plot([Qp, Qp], [Tcp_, Thp], "k-", lw=1.2)
    ax.text(Qp + 0.04 * x_lim[1], (Tcp_ + Thp) / 2,
            rf"$\Delta T_{{\min}}={dT:.1f}\;$K",
            ha="left", va="center", fontsize=9)

    # Arrows
    draw_arrow(ax, Qh, Th, 0.55, hot_color, backward=True)
    draw_arrow(ax, Qc, Tc, 0.40, cold_color)

    # Labels
    ax.annotate(r"$T_{h,o}$", xy=(Qh[0], Th[0]),
                xytext=(-6, 10), textcoords="offset points",
                fontsize=8, color=hot_color, ha="right", va="bottom")
    ax.annotate(r"$T_{h,i}$", xy=(Qh[-1], Th[-1]),
                xytext=(5, 4), textcoords="offset points",
                fontsize=8, color=hot_color, ha="left", va="bottom")
    ax.annotate(r"$T_{c,i}$", xy=(Qc[0], Tc[0]),
                xytext=(-6, -6), textcoords="offset points",
                fontsize=8, color=cold_color, ha="right", va="top")
    ax.annotate(r"$T_{c,o}$", xy=(Qc[-1], Tc[-1]),
                xytext=(5, -6), textcoords="offset points",
                fontsize=8, color=cold_color, ha="left", va="top")

    ax.set_xlim(*x_lim)
    ax.set_ylim(*y_lim)
    ax.set_xlabel(r"$\dot{Q}\;/\;[\mathrm{kW}]$")
    ax.set_ylabel(r"$T\;/\;[\mathrm{K}]$")
    style_legend(ax, loc="upper left", frame=False, fontsize=7)

# Panel labels
ax_a.text(-0.17, 1.02, r"\textbf{(a)}", transform=ax_a.transAxes,
          ha="left", va="bottom", fontsize=11)
ax_b.text(-0.17, 1.02, r"\textbf{(b)}", transform=ax_b.transAxes,
          ha="left", va="bottom", fontsize=11)

save_figure(fig, str(here / "nelium_2x1.png"))
save_figure(fig, str(here / "nelium_2x1.pdf"))
print("Saved nelium_2x1.pdf and .png")