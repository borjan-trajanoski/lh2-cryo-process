#!/usr/bin/env python
"""
plot_nelium_2x2.py
------------------
2x2 panel composite curve comparison:
  Row 1: H2 at 75 bar  |  (a) He-Ne 80/20  |  (b) He-Ne 10/90
  Row 2: H2 at 21 bar  |  (c) He-Ne 80/20  |  (d) He-Ne 10/90

Reads: nelium_case_A.csv .. nelium_case_D.csv, nelium_4case_summary.csv
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

data = {t: load_curves(t) for t in "ABCD"}
summary = pd.read_csv(here / "nelium_4case_summary.csv")

def get_summary(tag):
    return summary[summary["case"] == tag].iloc[0]

# ---- Pinch ---------------------------------------------------------------
def compute_pinch(Qh, Th, Qc, Tc, n=4001):
    Qm = min(Qh[-1], Qc[-1])
    Qg = np.linspace(0, Qm, n)
    dT = np.interp(Qg, Qh, Th) - np.interp(Qg, Qc, Tc)
    i = int(np.argmin(dT))
    return dT[i], Qg[i], np.interp(Qg[i], Qh, Th), np.interp(Qg[i], Qc, Tc)

# ---- Style ---------------------------------------------------------------
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

# ---- Titles for panels ---------------------------------------------------
panel_titles = {
    "A": r"$P_{\mathrm{H_2}} = 75$ bar, He-Ne 80/20 (3 bar)",
    "B": r"$P_{\mathrm{H_2}} = 75$ bar, He-Ne 10/90 (10 bar)",
    "C": r"$P_{\mathrm{H_2}} = 21$ bar, He-Ne 80/20 (3 bar)",
    "D": r"$P_{\mathrm{H_2}} = 21$ bar, He-Ne 10/90 (10 bar)",
}
panel_labels = {"A": r"\textbf{(a)}", "B": r"\textbf{(b)}",
                "C": r"\textbf{(c)}", "D": r"\textbf{(d)}"}

# ---- Limits --------------------------------------------------------------
all_Q = np.concatenate([np.concatenate([d[0], d[2]]) for d in data.values()])
all_T = np.concatenate([np.concatenate([d[1], d[3]]) for d in data.values()])
x_lim = (0, 1.06 * all_Q.max())
y_lim = (np.floor(all_T.min()) - 4, np.ceil(all_T.max()) + 6)

# ---- Figure --------------------------------------------------------------
with plt.style.context(["ieee"]):
    plt.rcParams["font.family"]      = graphic_font
    plt.rcParams["mathtext.fontset"] = math_font
    plt.rcParams["text.usetex"]      = True

    fig, axes = plt.subplots(2, 2, figsize=(8.5, 6.8),
                             gridspec_kw={"wspace": 0.28, "hspace": 0.32})

layout = [("A", axes[0, 0]), ("B", axes[0, 1]),
          ("C", axes[1, 0]), ("D", axes[1, 1])]

for tag, ax in layout:
    style_ax(ax)
    Qh, Th, Qc, Tc = data[tag]
    s = get_summary(tag)

    # Force cold composite to end at same Q as hot
    Q_end = Qh[-1]
    if Qc[-1] > Q_end:
        mask = Qc <= Q_end
        T_end = np.interp(Q_end, Qc, Tc)
        Qc = np.append(Qc[mask], Q_end)
        Tc = np.append(Tc[mask], T_end)
    elif Qc[-1] < Q_end:
        T_end = np.interp(Q_end, Qc, Tc)
        Qc = np.append(Qc, Q_end)
        Tc = np.append(Tc, T_end)

    # Cold composite label (matching reference style)
    cold_label = r"Cold composite (H$_2$ and He-Ne)"

    # Plot
    ax.plot(Qh, Th, color=hot_color, lw=1.4, ls="-", label=r"Hot H$_2$ stream")
    ax.plot(Qc, Tc, color=cold_color, lw=1.4, ls="-", label=cold_label)

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

    # Pinch
    dT, Qp, Thp, Tcp = compute_pinch(Qh, Th, Qc, Tc)
    ax.plot([Qp, Qp], [Tcp, Thp], "k-", lw=1.2)
    ax.text(Qp + 0.04 * x_lim[1], (Tcp + Thp) / 2,
            rf"$\Delta T_{{\min}}={dT:.1f}\;$K",
            ha="left", va="center", fontsize=8)

    # Arrows
    draw_arrow(ax, Qh, Th, 0.55, hot_color, backward=True)
    draw_arrow(ax, Qc, Tc, 0.40, cold_color)

    # Panel label (no title)
    ax.text(-0.17, 1.02, panel_labels[tag], transform=ax.transAxes,
            ha="left", va="bottom", fontsize=11)

    ax.set_xlim(*x_lim)
    ax.set_ylim(*y_lim)
    ax.set_xlabel(r"$\dot{Q}\;/\;[\mathrm{kW}]$")
    ax.set_ylabel(r"$T\;/\;[\mathrm{K}]$")

    style_legend(ax, loc="upper left", frame=False, fontsize=7)

save_figure(fig, str(here / "nelium_2x2.png"))
save_figure(fig, str(here / "nelium_2x2.pdf"))
print("Saved nelium_2x2.pdf and .png")