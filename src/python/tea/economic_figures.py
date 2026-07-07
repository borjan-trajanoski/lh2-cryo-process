"""
Economic figures for cryogenic-section TEA
==========================================
Generates two PDF figures for the thesis chapter:

  1. capex_breakdown_cryo.pdf  - stacked CAPEX bar chart, one bar per
     scenario (baseline, adapted), segments per equipment item from the
     largest (compressors, bottom) to the smallest (ejector, top).
     Totals are annotated above the bars; no segment value labels, the
     exact values live in the companion table in the chapter.

  2. tornado_slc_cryo.pdf      - SLC sensitivity tornado chart for the
     baseline scenario. Each parameter is varied by +/-50% from its
     nominal value; the resulting SLC range is plotted as two bars
     (one for the low-value input, one for the high-value input)
     emanating from the baseline SLC centerline.

All styling is taken from the thesis-wide PLOT_SETTINGS.py, which must
sit in the same directory: IEEE style context, Arial font, LaTeX
rendering, spine and tick geometry, legend styling (frameless), color
palette (ps.colors in palette order), and 1200 DPI saving via
PLOT_SETTINGS.save_figure(). The only local additions are
AutoMinorLocator(2) on numerical axes and the bracket-unit axis-label
convention ("X / [unit]"), per the thesis plot conventions.

Cost basis: 2024 USD, textbook correlations as in
cryo_textbook_tea.py. All numerical values match the script and
the chapter tables.
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

import PLOT_SETTINGS as ps

# Colors from the PLOT_SETTINGS palette
# ps.colors = [red, green, blue, orange, purple, brown, pink]
# Segment order: Compressors, LH2 tank, Heat exchanger, Turbines, Ejector
SEG_COLORS  = [ps.colors[2],   # blue   - compressors
               ps.colors[1],   # green  - LH2 tank
               ps.colors[0],   # red    - heat exchanger (PFHX-5)
               ps.colors[3],   # orange - turbines
               ps.colors[4]]   # purple - ejector
COLOR_NEG   = ps.colors[4]   # purple #984ea3 - parameter at -50% (tornado)
COLOR_POS   = ps.colors[3]   # orange #ff7f00 - parameter at +50% (tornado)


# ===========================================================================
# Styled multi-axes initializer (PLOT_SETTINGS.plot_init generalized)
# ===========================================================================
def init_axes(nrows=1, ncols=1, figsize=None):
    """Create a styled figure with an arbitrary subplot grid.

    Applies exactly the same style stack as PLOT_SETTINGS.plot_init()
    (ieee context, Arial, mathtext font, usetex, spine widths, in-facing
    major/minor ticks on all four sides) but supports multi-panel
    figures and a custom figure size, which plot_init() does not.
    """
    if figsize is None:
        figsize = ps.plot_size

    with plt.style.context(['ieee']):
        plt.rcParams['font.family'] = ps.graphic_font
        plt.rcParams['mathtext.fontset'] = ps.math_font
        plt.rcParams['text.usetex'] = True

        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        axes_list = np.atleast_1d(axes).ravel()

        for ax in axes_list:
            for spine in ax.spines.values():
                spine.set_linewidth(ps.spine_width)
            ax.tick_params(axis='both', which='major', direction='in',
                           width=ps.tick_width, length=ps.tick_length,
                           labelsize=ps.tick_labelsize,
                           bottom=True, top=True, left=True, right=True)
            ax.tick_params(axis='both', which='minor', direction='in',
                           width=ps.minor_tick_width, length=ps.minor_tick_length,
                           bottom=True, top=True, left=True, right=True)

        return fig, axes


# ===========================================================================
# Cost correlations (identical to cryo_textbook_tea.py)
# ===========================================================================
CEPCI = {1997: 386.5, 1998: 389.5, 2001: 394.3, 2010: 550.8, 2024: 800.0}
GBP_TO_USD_1997 = 1.64

# Amos (1998) LH2 dewar correlation, NREL/TP-570-25106 Table 10
RHO_LH2_KG_M3        = 71.0      # LH2 density at NBP (~20 K, 1 atm)
AMOS_BASE_KG         = 45.0      # base dewar capacity, kg H2
AMOS_BASE_USD_PER_KG = 441.0     # specific cost AT the base size, USD 1998
AMOS_EXP             = 0.70      # sizing exponent


def compressor_usd_2024(S_kw):
    """Towler & Sinnott (2010) Table 7.2"""
    return (580_000 + 20_000 * S_kw**0.6) * CEPCI[2024] / CEPCI[2010]


def turbine_usd_2024(S_kw):
    """Turton et al. (2018) Table A.1 radial expander"""
    L = math.log10(S_kw)
    return 10**(2.2476 + 1.4965 * L - 0.1618 * L * L) * CEPCI[2024] / CEPCI[2001]


def pfhx_usd_2024(V_m3):
    """ESDU 97006 (1997) plate-fin 4-to-6 streams"""
    return 81_186.396 * V_m3**0.35 * GBP_TO_USD_1997 * CEPCI[2024] / CEPCI[1997]


def tank_usd_2024(V_m3):
    """Amos (1998) NREL/TP-570-25106 Table 10 - liquid hydrogen dewar.
    C_base = 441 USD/kg * 45 kg = 19_845 USD (1998);
    Ce_1998 = C_base * (M/45)^0.70 with M = 71 kg/m^3 * V.
    441 USD/kg is the unit cost AT the base size, not the power-law
    coefficient. Extrapolated ~2_761x above the 45 kg base capacity."""
    M_kg = RHO_LH2_KG_M3 * V_m3
    c_base_1998 = AMOS_BASE_USD_PER_KG * AMOS_BASE_KG
    return c_base_1998 * (M_kg / AMOS_BASE_KG)**AMOS_EXP * CEPCI[2024] / CEPCI[1998]


def crf_fn(i, n):
    """Capital recovery factor"""
    return i * (1 + i)**n / ((1 + i)**n - 1)


# ===========================================================================
# Base TEA parameters
# ===========================================================================
INTEREST_BASE   = 0.09
LIFETIME_BASE   = 20
CRF_BASE        = crf_fn(INTEREST_BASE, LIFETIME_BASE)
FIXED_OPEX_BASE = 0.05
ELEC_BASE       = 0.110
HOURS_PER_YEAR  = 8_322
MTR_EFF         = 0.96
GEN_EFF         = 0.80
LH2_BASE_KG_DAY = 86_000
EJECTOR_USD     = 300_000

# Baseline equipment sizes
C_BASE_KW       = [2973.0, 3064.0, 2210.0, 2157.0]
T_BASE_KW       = [579.2, 352.9]
PFHX5_BASE_M3   = 21.0
TANK_M3         = 1750.0

# Adapted equipment sizes
C_ADAPT_KW      = [5741.2, 5788.4, 4092.7, 4140.8]
T_ADAPT_KW      = [1197.8, 620.2]
PFHX5_ADAPT_M3  = 28.7

# Baseline equipment costs (USD)
c_cost_base    = [compressor_usd_2024(p) for p in C_BASE_KW]
t_cost_base    = [turbine_usd_2024(p)    for p in T_BASE_KW]
pfhx_cost_base = pfhx_usd_2024(PFHX5_BASE_M3)
tank_cost      = tank_usd_2024(TANK_M3)

# Adapted equipment costs (USD)
c_cost_adapt    = [compressor_usd_2024(p) for p in C_ADAPT_KW]
t_cost_adapt    = [turbine_usd_2024(p)    for p in T_ADAPT_KW]
pfhx_cost_adapt = pfhx_usd_2024(PFHX5_ADAPT_M3)

# Totals
capex_base  = sum(c_cost_base)  + sum(t_cost_base)  + pfhx_cost_base  + tank_cost
capex_adapt = sum(c_cost_adapt) + sum(t_cost_adapt) + pfhx_cost_adapt + tank_cost + EJECTOR_USD


# ===========================================================================
# SLC computation (parameterized for sensitivity analysis)
# ===========================================================================
def slc_cryo(elec_price=ELEC_BASE,
             interest=INTEREST_BASE, lifetime=LIFETIME_BASE,
             fixed_opex=FIXED_OPEX_BASE,
             comp_mult=1.0, turb_mult=1.0,
             pfhx_mult=1.0, tank_mult=1.0,
             capex_mult=1.0):
    """Cryogenic SLC with optional perturbations to each input."""
    capex = (sum(c_cost_base) * comp_mult
             + sum(t_cost_base) * turb_mult
             + pfhx_cost_base   * pfhx_mult
             + tank_cost        * tank_mult)
    capex *= capex_mult
    crf_v = crf_fn(interest, lifetime)
    elec_kw = sum(C_BASE_KW) / MTR_EFF - sum(T_BASE_KW) * GEN_EFF
    annual_cap   = capex * crf_v
    annual_fixed = capex * fixed_opex
    annual_elec  = elec_kw * HOURS_PER_YEAR * elec_price
    lh2_annual   = LH2_BASE_KG_DAY * 365 * 0.95
    return (annual_cap + annual_fixed + annual_elec) / lh2_annual


SLC_BASE = slc_cryo()
print(f'Baseline SLC = {SLC_BASE:.4f} USD/kg')


# ===========================================================================
# FIGURE 1: CAPEX breakdown (stacked bars, one per scenario)
# ===========================================================================
def plot_capex_breakdown():
    # Segments bottom-to-top: largest (compressors) to smallest (ejector),
    # matching the row order of the chapter table. The ejector does not
    # exist in the baseline (0).
    segments = [
        ('Compressors',  sum(c_cost_base)  / 1e6,   sum(c_cost_adapt) / 1e6),
        (r'LH$_2$ tank', tank_cost         / 1e6,   tank_cost         / 1e6),
        ('Heat exchanger', pfhx_cost_base  / 1e6,   pfhx_cost_adapt   / 1e6),
        ('Turbines',     sum(t_cost_base)  / 1e6,   sum(t_cost_adapt) / 1e6),
        ('Ejector',      0.0,                       EJECTOR_USD       / 1e6),
    ]

    fig, axes = init_axes(1, 1, figsize=(5.2, 3.4))
    ax = np.atleast_1d(axes).ravel()[0]

    x = np.array([0, 1])
    width = 0.55
    bottoms = np.zeros(2)

    for (name, v_base, v_adapt), color in zip(segments, SEG_COLORS):
        vals = np.array([v_base, v_adapt])
        ax.bar(x, vals, width, bottom=bottoms, color=color,
               edgecolor='black', linewidth=ps.minor_tick_width, label=name)
        bottoms += vals

    # Bold totals above each bar (the only numbers in the plot)
    for xi, tot in zip(x, [capex_base / 1e6, capex_adapt / 1e6]):
        ax.text(xi, tot + 0.6, f'${tot:.2f}$', ha='center', va='bottom',
                fontsize=ps.tick_labelsize, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(['Baseline', 'Adapted'], fontsize=ps.tick_labelsize)
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylabel(r'Direct cost / [M USD]', fontsize=ps.tick_labelsize)
    ax.set_ylim(0, capex_adapt / 1e6 * 1.12)
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    # Legend outside the axes on the right, frameless, segment order
    # top-to-bottom matching the visual stack top-to-bottom.
    handles, labels = ax.get_legend_handles_labels()
    ps.style_legend(ax, loc='center left', frame=False,
                    bbox_to_anchor=(1.02, 0.5),
                    handles=handles[::-1], labels=labels[::-1])

    fig.tight_layout()
    ps.save_figure(fig, 'capex_breakdown_cryo.pdf')
    plt.close(fig)
    print('Wrote capex_breakdown_cryo.pdf')


# ===========================================================================
# FIGURE 2: Tornado sensitivity of SLC
# ===========================================================================
def plot_tornado():
    # Each entry: (display name, perturbation fn at -50%, at +50%)
    perturbations = [
        ('Electricity price',
         lambda: slc_cryo(elec_price=ELEC_BASE   * 0.5),
         lambda: slc_cryo(elec_price=ELEC_BASE   * 1.5)),
        (r'Fixed CAPEX',
         lambda: slc_cryo(capex_mult=0.5),
         lambda: slc_cryo(capex_mult=1.5)),
        ('Fixed OPEX',
         lambda: slc_cryo(fixed_opex=FIXED_OPEX_BASE * 0.5),
         lambda: slc_cryo(fixed_opex=FIXED_OPEX_BASE * 1.5)),
        ('Interest rate',
         lambda: slc_cryo(interest=INTEREST_BASE  * 0.5),
         lambda: slc_cryo(interest=INTEREST_BASE  * 1.5)),
        ('Project lifetime',
         lambda: slc_cryo(lifetime=LIFETIME_BASE  * 0.5),
         lambda: slc_cryo(lifetime=LIFETIME_BASE  * 1.5)),
        (r'LH$_2$ tank cost',
         lambda: slc_cryo(tank_mult=0.5),
         lambda: slc_cryo(tank_mult=1.5)),
        ('Compressor cost',
         lambda: slc_cryo(comp_mult=0.5),
         lambda: slc_cryo(comp_mult=1.5)),
    ]
    # The heat exchanger and turbine cost levers are omitted from the
    # figure: their combined impact is below one cent per kilogram (the
    # chapter prose quantifies this), and at the plot scale the bars are
    # invisible slivers.

    # Compute the SLC at each perturbation
    rows = []
    for name, fn_low, fn_high in perturbations:
        slc_lo = fn_low()    # SLC at -50% input
        slc_hi = fn_high()   # SLC at +50% input
        # Total impact = |slc_lo - base| + |slc_hi - base|
        impact = abs(slc_lo - SLC_BASE) + abs(slc_hi - SLC_BASE)
        rows.append((name, slc_lo, slc_hi, impact))

    # Sort by impact descending (largest at top)
    rows.sort(key=lambda r: r[3], reverse=True)

    names    = [r[0] for r in rows]
    slc_lows = np.array([r[1] for r in rows])
    slc_his  = np.array([r[2] for r in rows])

    # Bar geometry: from baseline SLC to slc_low (one bar), and to slc_high
    fig, axes = init_axes(1, 1, figsize=(6.4, 4.5))
    ax = np.atleast_1d(axes).ravel()[0]

    y = np.arange(len(names))
    height = 0.4
    for i in range(len(names)):
        ax.barh(y[i], slc_lows[i] - SLC_BASE, left=SLC_BASE, height=height,
                color=COLOR_NEG, edgecolor='black',
                linewidth=ps.minor_tick_width,
                label=r'$-50\%$' if i == 0 else None)
        ax.barh(y[i], slc_his[i]  - SLC_BASE, left=SLC_BASE, height=height,
                color=COLOR_POS, edgecolor='black',
                linewidth=ps.minor_tick_width,
                label=r'$+50\%$' if i == 0 else None)

    # Centerline at base SLC
    ax.axvline(SLC_BASE, color='black', linewidth=ps.spine_width * 0.8, zorder=10)

    ax.set_yticks(y)
    ax.set_yticklabels(names)
    ax.invert_yaxis()
    ax.set_xlabel(r'Cryogenic SLC / [USD/kg LH$_2$]', fontsize=ps.tick_labelsize)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))

    # Symmetric x-limits around base
    x_extent = max(abs(slc_lows - SLC_BASE).max(),
                   abs(slc_his  - SLC_BASE).max()) * 1.18

    ax.set_xlim(SLC_BASE - x_extent, SLC_BASE + x_extent)

    ps.style_legend(ax, loc='lower right', frame=False,
                    fontsize=ps.tick_labelsize, markerscale=1.2)

    fig.tight_layout()
    ps.save_figure(fig, 'tornado_slc_cryo.pdf')
    plt.close(fig)
    print('Wrote tornado_slc_cryo.pdf')

    # Also print the numerical sensitivity table to stdout
    print('\nTornado sensitivity (baseline scenario):')
    print(f'  Baseline SLC = {SLC_BASE:.4f} USD/kg')
    print(f'  {"Parameter":<25s} {"SLC @ -50%":>12s} {"SLC @ +50%":>12s} {"Impact":>10s}')
    print('-' * 70)
    for name, slc_lo, slc_hi, impact in rows:
        print(f'  {name:<25s} {slc_lo:>10.4f}   {slc_hi:>10.4f}   {impact:>8.4f}')


# ===========================================================================
# Run both
# ===========================================================================
if __name__ == '__main__':
    plot_capex_breakdown()
    plot_tornado()
    print('\nDone.')
