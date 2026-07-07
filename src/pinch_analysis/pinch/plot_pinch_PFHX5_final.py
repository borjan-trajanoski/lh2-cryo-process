#!/usr/bin/env python
"""
Plot PFHX-5 pinch analysis: baseline (panel a), ejector low-eta (panel b),
ejector high-eta (panel c) side by side.

Reads three CSVs from the same directory:
    pinch_curves_baseline.csv
    pinch_curves_ejector_low_eta.csv     # f = 0.20, eta = 0.20
    pinch_curves_ejector_high_eta.csv    # f = 0.10, eta = 0.30

Pinch metrics (dT_min, Q_pinch, T_hot_pinch, T_cold_pinch) are computed
on the fly from the curve data.

Annotations match the v9 single-ejector-panel version one-for-one:
- Endpoint symbol labels: T_h,i, T_h,o, T_c,i, T_c,o (no numbers).
- Direction arrows on each curve (hot stream cools right-to-left, cold
  composite warms left-to-right).
- T_R,o kink marker on the cold composite at the temperature where the
  R-stream finishes warming and stops contributing. Same value (46.65 K)
  in both ejector panels because T_R_out is fixed in the model.
- Pinch marker (vertical black segment between curves) and dT_min text.

Output: pinch_PFHX5_combined.png and .pdf
"""

from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt

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


# --------- Stream-end / kink temperatures (constants from the model;
#           used for annotation placement only, sanity-checked at runtime
#           against each panel's cold-composite range) -------------------
T_R_OUT_BASELINE  = 46.65   # K, R-stream outlet, baseline case
T_R_OUT_EJECTOR   = 46.65   # K, R-stream outlet, both ejector cases


# ---------- Helpers ----------------------------------------------------
def load_curves(csv_path):
    """Load a pinch_curves CSV; return Q_hot, T_hot, Q_cold, T_cold."""
    df = pd.read_csv(csv_path)
    hot_mask  = ~df['Q_hot_kW'].isna()
    cold_mask = ~df['Q_cold_kW'].isna()
    Q_hot  = df.loc[hot_mask,  'Q_hot_kW'].values
    T_hot  = df.loc[hot_mask,  'T_hot_K'].values
    Q_cold = df.loc[cold_mask, 'Q_cold_kW'].values
    T_cold = df.loc[cold_mask, 'T_cold_K'].values
    return Q_hot, T_hot, Q_cold, T_cold


def compute_pinch(Q_hot, T_hot, Q_cold, T_cold, n_grid=4001):
    """
    Interpolate both composites onto a common Q grid over
    [0, min(Q_hot_max, Q_cold_max)] and return
    (dT_min, Q_pinch, T_hot_pinch, T_cold_pinch).
    """
    Q_max_common = min(Q_hot[-1], Q_cold[-1])
    Q_grid = np.linspace(0, Q_max_common, n_grid)
    T_h = np.interp(Q_grid, Q_hot,  T_hot)
    T_c = np.interp(Q_grid, Q_cold, T_cold)
    dT = T_h - T_c
    i = int(np.argmin(dT))
    return dT[i], Q_grid[i], T_h[i], T_c[i]


def style_axis(ax):
    """Apply IEEE-style spines and tick params to a given axis."""
    for spine in ax.spines.values():
        spine.set_linewidth(spine_width)
    ax.tick_params(axis='both', which='major', direction='in',
                   width=tick_width, length=tick_length,
                   labelsize=tick_labelsize,
                   bottom=True, top=True, left=True, right=True)
    ax.tick_params(axis='both', which='minor', direction='in',
                   width=minor_tick_width, length=minor_tick_length,
                   bottom=True, top=True, left=True, right=True)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))


def draw_arrow_along(ax, Q, T, frac, color, direction='forward'):
    """
    Draw a direction arrow on the curve (Q, T) at fractional position frac.
    direction = 'forward' -> increasing-Q; 'backward' -> decreasing-Q.
    """
    n = len(Q)
    if n < 4:
        return
    i = int(np.clip(frac * n, 2, n - 3))
    di = max(1, n // 200)
    if direction == 'forward':
        j_tail, j_head = i, min(n - 1, i + di)
    elif direction == 'backward':
        j_tail, j_head = i, max(0, i - di)
    else:
        raise ValueError(f"direction must be 'forward' or 'backward', "
                         f"got {direction!r}")
    ax.annotate(
        '', xy=(Q[j_head], T[j_head]), xytext=(Q[j_tail], T[j_tail]),
        arrowprops=dict(arrowstyle='->', color=color,
                        lw=1.2, mutation_scale=12),
        zorder=5,
    )


def annotate_endpoints(ax, Q, T, label_in, label_out, color,
                       offset_in=(8, -6), offset_out=(-8, 6)):
    """Place text labels at (Q[0], T[0]) and (Q[-1], T[-1])."""
    def _ha_va(off):
        ha = 'left'   if off[0] >= 0 else 'right'
        va = 'bottom' if off[1] >= 0 else 'top'
        return ha, va

    ha_i, va_i = _ha_va(offset_in)
    ha_o, va_o = _ha_va(offset_out)
    ax.annotate(
        label_in, xy=(Q[0], T[0]), xytext=offset_in,
        textcoords='offset points',
        ha=ha_i, va=va_i, fontsize=8, color=color,
    )
    ax.annotate(
        label_out, xy=(Q[-1], T[-1]), xytext=offset_out,
        textcoords='offset points',
        ha=ha_o, va=va_o, fontsize=8, color=color,
    )


def draw_panel(ax, Q_hot, T_hot, Q_cold, T_cold,
               hot_color, cold_color, pinch_info,
               Q_max_common, annotate_offset,
               panel_id, T_R_out,
               offset_T_co=(-4, -10),
               offset_T_ci=(8, 4),
               offset_T_ho=(4, 24),
               offset_T_hi=(5, 0),
               cold_arrow_frac=0.30,
               R_kink_label_offset=(-30, 18)):
    """Draw composites + annotations on ax. Same logic as v9."""
    dT_min, Q_p, T_h_p, T_c_p = pinch_info

    ax.plot(Q_hot,  T_hot,
            color=hot_color, linewidth=1.4, linestyle='-',
            label='Hot H$_2$ stream')
    ax.plot(Q_cold, T_cold,
            color=cold_color, linewidth=1.4, linestyle='-',
            label='Cold composite (H$_2$ and He-Ne)')

    grey = (0.40, 0.40, 0.40)
    ax.axvline(Q_hot[-1],  color=grey, linewidth=0.9, linestyle='--',
               alpha=0.9, label=r'End of hot $\mathrm{H_2}$ stream')
    ax.axvline(Q_cold[-1], color=grey, linewidth=0.9, linestyle=':',
               alpha=0.9, label='End of cold composite')

    # Pinch marker
    ax.plot([Q_p, Q_p], [T_c_p, T_h_p],
            color='black', linewidth=1.2, linestyle='-')
    dx_frac, dy = annotate_offset
    # Round dT_min to 1 decimal for display: 2.48 K shows as 2.5 K.
    # The exact value is reported in the console summary.
    ax.text(Q_p + dx_frac * Q_max_common,
            T_c_p + dy,
            rf'$\Delta T_{{\min}}={dT_min:.1f}\;$K',
            ha='left', va='top', fontsize=9)

    # Direction arrows
    draw_arrow_along(ax, Q_hot,  T_hot,  frac=0.55,
                     color=hot_color,  direction='backward')
    draw_arrow_along(ax, Q_cold, T_cold, frac=cold_arrow_frac,
                     color=cold_color, direction='forward')

    # Endpoint symbol labels (per-panel offsets)
    annotate_endpoints(
        ax, Q_hot, T_hot,
        label_in =r'$T_{h,o}$',
        label_out=r'$T_{h,i}$',
        color=hot_color,
        offset_in=offset_T_ho,
        offset_out=offset_T_hi,
    )
    annotate_endpoints(
        ax, Q_cold, T_cold,
        label_in =r'$T_{c,i}$',
        label_out=r'$T_{c,o}$',
        color=cold_color,
        offset_in=offset_T_ci,
        offset_out=offset_T_co,
    )

    # R-stream kink marker (T_R,o)
    Q_kink = np.interp(T_R_out, T_cold, Q_cold)
    ax.plot(Q_kink, T_R_out, marker='o', markersize=5,
            markerfacecolor='white',
            markeredgecolor=cold_color,
            markeredgewidth=1.2, zorder=6)
    ax.annotate(
        r'$T_{R,o}$',
        xy=(Q_kink, T_R_out),
        xytext=R_kink_label_offset, textcoords='offset points',
        ha='center', va='bottom', fontsize=9, color=cold_color,
        arrowprops=dict(arrowstyle='-', color=cold_color, lw=0.6),
    )


# ---------- Main -------------------------------------------------------
here = Path(__file__).parent

# Load all three cases
Q_hot_b, T_hot_b, Q_cold_b, T_cold_b = load_curves(
    here / 'pinch_curves_baseline.csv')
Q_hot_l, T_hot_l, Q_cold_l, T_cold_l = load_curves(
    here / 'pinch_curves_ejector_low_eta.csv')
Q_hot_h, T_hot_h, Q_cold_h, T_cold_h = load_curves(
    here / 'pinch_curves_ejector_high_eta.csv')

# Compute pinch for each
pinch_b = compute_pinch(Q_hot_b, T_hot_b, Q_cold_b, T_cold_b)
pinch_l = compute_pinch(Q_hot_l, T_hot_l, Q_cold_l, T_cold_l)
pinch_h = compute_pinch(Q_hot_h, T_hot_h, Q_cold_h, T_cold_h)

print(f'Baseline (a):           dT_min = {pinch_b[0]:.3f} K @ Q = {pinch_b[1]:.1f} kW '
      f'(T_hot = {pinch_b[2]:.2f} K, T_cold = {pinch_b[3]:.2f} K)')
print(f'Ejector low-eta (b):    dT_min = {pinch_l[0]:.3f} K @ Q = {pinch_l[1]:.1f} kW '
      f'(T_hot = {pinch_l[2]:.2f} K, T_cold = {pinch_l[3]:.2f} K)')
print(f'Ejector high-eta (c):   dT_min = {pinch_h[0]:.3f} K @ Q = {pinch_h[1]:.1f} kW '
      f'(T_hot = {pinch_h[2]:.2f} K, T_cold = {pinch_h[3]:.2f} K)')

# Cross-check: T_R_out must lie inside each panel's cold-composite range.
for label, T_R, T_cold in [('baseline', T_R_OUT_BASELINE, T_cold_b),
                            ('low_eta',  T_R_OUT_EJECTOR,  T_cold_l),
                            ('high_eta', T_R_OUT_EJECTOR,  T_cold_h)]:
    if not (T_cold[0] <= T_R <= T_cold[-1]):
        raise ValueError(
            f'T_R_out = {T_R} K is outside the {label} '
            f'cold-composite range [{T_cold[0]}, {T_cold[-1]}] K. '
            f'Update the constant or check the input CSV.'
        )

# Common axis limits across all three panels so visual comparison is fair
T_lo = min(T_hot_b.min(), T_cold_b.min(),
           T_hot_l.min(), T_cold_l.min(),
           T_hot_h.min(), T_cold_h.min())
T_hi = max(T_hot_b.max(), T_cold_b.max(),
           T_hot_l.max(), T_cold_l.max(),
           T_hot_h.max(), T_cold_h.max())
Q_max = max(Q_hot_b.max(), Q_cold_b.max(),
            Q_hot_l.max(), Q_cold_l.max(),
            Q_hot_h.max(), Q_cold_h.max())

y_lim = (np.floor(T_lo) - 4, np.ceil(T_hi) + 6)
x_lim = (0.0, 1.06 * Q_max)

# Figure: 1x3 panels, slightly wider than the v9 1x2 layout
with plt.style.context(['ieee']):
    plt.rcParams['font.family']      = graphic_font
    plt.rcParams['mathtext.fontset'] = math_font
    plt.rcParams['text.usetex']      = True

    fig, (ax_a, ax_b, ax_c) = plt.subplots(
        1, 3, figsize=(12.5, 3.6), sharey=False,
        gridspec_kw={'wspace': 0.28},
    )

style_axis(ax_a)
style_axis(ax_b)
style_axis(ax_c)

hot_color  = (0.75, 0.10, 0.10)
cold_color = (0.10, 0.25, 0.75)

# Panel (a): baseline
draw_panel(ax_a, Q_hot_b, T_hot_b, Q_cold_b, T_cold_b,
           hot_color, cold_color, pinch_b,
           Q_max_common=x_lim[1], annotate_offset=(0.07, -1.5),
           panel_id='baseline', T_R_out=T_R_OUT_BASELINE,
           offset_T_co=(-4, -16),
           offset_T_ci=(8, 4),
           offset_T_ho=(4, 24),
           offset_T_hi=(5, 0))

# Panel (b): ejector, low_eta (f=0.20, eta=0.20)
draw_panel(ax_b, Q_hot_l, T_hot_l, Q_cold_l, T_cold_l,
           hot_color, cold_color, pinch_l,
           Q_max_common=x_lim[1], annotate_offset=(0.05, -1.5),
           panel_id='ejector_low_eta', T_R_out=T_R_OUT_EJECTOR,
           offset_T_co=(6, 0),
           offset_T_ci=(8, 4),
           offset_T_ho=(4, 24),
           offset_T_hi=(5, 0))

# Panel (c): ejector, high_eta (f=0.10, eta=0.30)
draw_panel(ax_c, Q_hot_h, T_hot_h, Q_cold_h, T_cold_h,
           hot_color, cold_color, pinch_h,
           Q_max_common=x_lim[1], annotate_offset=(-0.115, -7.0),
           panel_id='ejector_high_eta', T_R_out=T_R_OUT_EJECTOR,
           offset_T_co=(6, 0),
           offset_T_ci=(8, 4),
           offset_T_ho=(4, 24),
           offset_T_hi=(5, 0))

# Apply identical limits to all three panels
for ax in (ax_a, ax_b, ax_c):
    ax.set_xlim(*x_lim)
    ax.set_ylim(*y_lim)
    ax.set_xlabel(r'$\dot{Q}\;/\;[\mathrm{kW}]$')
    ax.set_ylabel(r'$T\;/\;[\mathrm{K}]$')

# Panel labels (a) / (b) / (c) outside axes top-left, LaTeX bold
ax_a.text(-0.17, 1.02, r'\textbf{(a)}',
          transform=ax_a.transAxes, ha='left', va='bottom', fontsize=11)
ax_b.text(-0.17, 1.02, r'\textbf{(b)}',
          transform=ax_b.transAxes, ha='left', va='bottom', fontsize=11)
ax_c.text(-0.17, 1.02, r'\textbf{(c)}',
          transform=ax_c.transAxes, ha='left', va='bottom', fontsize=11)

# Legends (frameless)
style_legend(ax_a, loc='upper left', frame=False)
style_legend(ax_b, loc='upper left', frame=False)
style_legend(ax_c, loc='upper left', frame=False)

# Save
save_figure(fig, str(here / 'pinch_PFHX5_combined.png'))
save_figure(fig, str(here / 'pinch_PFHX5_combined.pdf'))
