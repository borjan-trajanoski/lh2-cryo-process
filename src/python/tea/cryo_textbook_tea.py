"""
Simplified Cryogenic SLC - all textbook correlations
====================================================
SCOPE: cryogenic section only
  C-1..C-4, T-1, T-2, PFHX-5, LH2 tank, ejector (retrofit only)

CAPEX        = sum of purchased equipment costs (no Hand factor, no stack)
Annual CAPEX = CAPEX * CRF                 [pure annuity, i=9%, N=20]
Annual fixed OPEX = 5% of CAPEX            [mid-range of 3-6% maintenance
                                            + overhead per Towler & Sinnott]
Annual elec OPEX  = elec_kw * 8322 h/yr * price
SLC          = (annual CAPEX + annual fixed + annual elec) / annual LH2

ALL COST CORRELATIONS ARE TEXTBOOK / PUBLISHED REFERENCES
---------------------------------------------------------
COMPRESSOR (Towler & Sinnott 2010, "Chemical Engineering Design" 2nd ed,
            Table 7.2 - single-stage centrifugal compressor):

    Ce_2010 = 580_000 + 20_000 * S^0.6        [USD 2010, S = shaft kW]
    Validity: 75 <= S <= 30_000 kW

TURBINE (Turton et al. 2018, "Analysis, Synthesis and Design of Chemical
         Processes" 5th ed, Table A.1 - Radial gas/liquid expander):

    log10(Ce_2001) = K1 + K2*log10(S) + K3*(log10(S))^2
    K1 = 2.2476,  K2 = 1.4965,  K3 = -0.1618    [USD 2001, S = shaft kW]
    Validity: 100 <= S <= 1500 kW

PLATE-FIN HEAT EXCHANGER (ESDU 97006, 1997, plate & fin 4-to-6 streams):

    Ce_1997_GBP = 81_186.396 * V^0.35           [GBP 1997, V = volume m^3]
    Validity: 0.01 <= V <= 2 m^3    *** PFHX-5 IS BEYOND VALIDITY ***
    Best available textbook correlation for plate-fin HEX. Underestimates
    real cryogenic plate-fin HEX cost at 21-29 m^3; vendor quotes from
    Linde/Sumitomo/Chart are typically 5-10x higher at this size range
    due to material premiums and small production volumes.

LH2 STORAGE TANK (Amos 1998, "Costs of Storing and Transporting Hydrogen",
                  NREL/TP-570-25106, Table 10 - liquid hydrogen dewar):

    Table 10 lists, for LH2 dewars:
        base size           M_base = 45 kg H2 capacity
        base unit cost      441 USD per kg H2 capacity   [USD 1998]
        sizing exponent     0.70

    441 USD/kg is the SPECIFIC cost AT the base size, not the power-law
    coefficient. The total base cost is computed first and then scaled:

        C_base  = 441 * 45 = 19_845 USD                  [USD 1998]
        Ce_1998 = C_base * (M / M_base)^0.70             [USD 1998, M in kg]

    Volume-to-mass conversion at the LH2 normal boiling point:
        M = rho_LH2 * V,  rho_LH2 = 71 kg/m^3
        1750 m^3  ->  M = 124_250 kg

    *** EXTRAPOLATION: M/M_base ~ 2_761x above the 45 kg dewar base. ***
    The exponent 0.70 is in the standard 0.6-0.7 vessel-scaling range,
    and the resulting specific cost is of the same order as recent
    estimates for large LH2 tanks (30-50 USD/kg-H2, tank level).
    Treated as INSTALLED cost (NREL storage capital basis), entering
    the summation without further mark-up, as before.

EJECTOR (vendor flat estimate; literature range 100-500k USD):

    Ce_2024_USD = 300_000

INFLATION & CURRENCY
--------------------
CEPCI: 1997=386.5, 1998=389.5, 2001=394.3, 2010=550.8, 2024=800.0
GBP to USD (1997): 1.64 USD/GBP  (BoE historical average)
(The Amos correlation is denominated in USD; no currency conversion.)

POWER VALUES verified against thesis Table tab:tea-hene-comparison.
"""

# ===========================================================================
# Constants
# ===========================================================================
CEPCI = {1997: 386.5, 1998: 389.5, 2001: 394.3, 2010: 550.8, 2024: 800.0}
GBP_TO_USD_1997 = 1.64

CRF              = 0.10955     # 9% interest, 20-yr lifetime
HOURS_PER_YEAR   = 8322        # 95% utilization
MTR_EFF          = 0.96
GEN_EFF          = 0.80
FIXED_OPEX_FRAC  = 0.05        # 5% of CAPEX/yr (mid-range maintenance+overhead)
LH2_BASE_KG_DAY  = 86_000
LH2_TANK_M3      = 1750.0
EJECTOR_USD_2024 = 300_000
ELEC_BASE_USD_KWH = 0.110

# Amos (1998) LH2 dewar correlation, NREL/TP-570-25106 Table 10
RHO_LH2_KG_M3        = 71.0      # LH2 density at NBP (~20 K, 1 atm)
AMOS_BASE_KG         = 45.0      # base dewar capacity, kg H2
AMOS_BASE_USD_PER_KG = 441.0     # specific cost AT the base size, USD 1998
AMOS_EXP             = 0.70      # sizing exponent

# Electricity prices to sweep
ELEC_PRICES = [0.04, 0.06, 0.08, 0.11, 0.15, 0.20]


# ===========================================================================
# Textbook cost correlations (each returns 2024 USD)
# ===========================================================================
def compressor_purchased_usd_2024(shaft_power_kw, warn=False):
    """Towler & Sinnott (2010) Table 7.2 - single-stage centrifugal compressor.
    Ce_2010 = 580_000 + 20_000 * S^0.6; validity 75 to 30_000 kW."""
    if warn and not (75 <= shaft_power_kw <= 30_000):
        print(f"  WARNING: compressor at {shaft_power_kw} kW outside TS-2010 validity 75-30000 kW")
    ce_2010 = 580_000 + 20_000 * shaft_power_kw**0.6
    return ce_2010 * CEPCI[2024] / CEPCI[2010]


def turbine_purchased_usd_2024(shaft_power_kw, warn=False):
    """Turton et al. (2018) Table A.1 - radial gas/liquid expander.
    log10(Ce_2001) = 2.2476 + 1.4965*log10(S) - 0.1618*(log10(S))^2;
    validity 100 to 1500 kW."""
    import math
    if warn and not (100 <= shaft_power_kw <= 1500):
        print(f"  WARNING: turbine at {shaft_power_kw} kW outside Turton-2018 validity 100-1500 kW")
    L = math.log10(shaft_power_kw)
    log_ce = 2.2476 + 1.4965 * L - 0.1618 * L * L
    ce_2001 = 10**log_ce
    return ce_2001 * CEPCI[2024] / CEPCI[2001]


def pfhx_purchased_usd_2024(volume_m3, warn=False):
    """ESDU 97006 (1997) - plate & fin heat exchanger, 4-to-6 streams.
    Ce_1997_GBP = 81_186.396 * V^0.35; validity 0.01 to 2 m^3.
    PFHX-5 at 21-29 m^3 is beyond validity; vendor quotes for cryogenic
    plate-fin HEX at this size are typically 5-10x higher due to material
    premiums and small production volumes."""
    if warn and volume_m3 > 2.0:
        print(f"  WARNING: PFHX at {volume_m3} m3 BEYOND ESDU 97006 validity (max 2 m3)")
        print(f"           Result extrapolated; likely underestimates real cryogenic PFHX cost.")
    ce_1997_gbp = 81_186.396 * volume_m3**0.35
    ce_1997_usd = ce_1997_gbp * GBP_TO_USD_1997
    return ce_1997_usd * CEPCI[2024] / CEPCI[1997]


def lh2_tank_installed_usd_2024(volume_m3, warn=False):
    """Amos (1998) NREL/TP-570-25106 Table 10 - liquid hydrogen dewar.
    C_base  = 441 USD/kg * 45 kg = 19_845 USD (1998)
    Ce_1998 = C_base * (M / 45)^0.70,  M = 71 kg/m^3 * V
    NOTE: 441 USD/kg is the unit cost AT the base size, not the power-law
    coefficient. Extrapolated ~2_761x above the 45 kg base capacity."""
    mass_kg = RHO_LH2_KG_M3 * volume_m3
    if warn and mass_kg > AMOS_BASE_KG:
        print(f"  WARNING: LH2 tank at {mass_kg:,.0f} kg is {mass_kg/AMOS_BASE_KG:,.0f}x "
              f"above the Amos 45 kg dewar base (extrapolated)")
    c_base_1998 = AMOS_BASE_USD_PER_KG * AMOS_BASE_KG
    ce_1998 = c_base_1998 * (mass_kg / AMOS_BASE_KG)**AMOS_EXP
    return ce_1998 * CEPCI[2024] / CEPCI[1998]


# ===========================================================================
# Scenarios (values verified against thesis Table tab:tea-hene-comparison)
# ===========================================================================
SCENARIOS = {
    'Baseline': {
        'c_kw':            [2973.0, 3064.0, 2210.0, 2157.0],
        't_kw':            [579.2, 352.9],
        'pfhx5_m3':        21.0,
        'include_ejector': False,
        'prod_factor':     1.0,
    },
    'Adapted (Pp=40 bar)': {
        'c_kw':            [5741.2, 5788.4, 4092.7, 4140.8],
        't_kw':            [1197.8, 620.2],
        'pfhx5_m3':        28.7,
        'include_ejector': True,
        'prod_factor':     1.214,
    },
}


# ===========================================================================
# Compute one scenario at one electricity price
# ===========================================================================
def compute(s, elec_price, quiet=True):
    if not quiet:
        print(f"\n  -- {s.get('label', '?')} --")
    c_cost     = [compressor_purchased_usd_2024(p) for p in s['c_kw']]
    t_cost     = [turbine_purchased_usd_2024(p)    for p in s['t_kw']]
    pfhx5_cost = pfhx_purchased_usd_2024(s['pfhx5_m3'])
    tank_cost  = lh2_tank_installed_usd_2024(LH2_TANK_M3)
    ej_cost    = EJECTOR_USD_2024 if s['include_ejector'] else 0.0

    capex = sum(c_cost) + sum(t_cost) + pfhx5_cost + tank_cost + ej_cost

    elec_kw      = sum(s['c_kw']) / MTR_EFF - sum(s['t_kw']) * GEN_EFF
    annual_cap   = capex * CRF
    annual_fixed = FIXED_OPEX_FRAC * capex
    annual_elec  = elec_kw * HOURS_PER_YEAR * elec_price

    lh2_kg = LH2_BASE_KG_DAY * s['prod_factor'] * 365 * 0.95
    slc    = (annual_cap + annual_fixed + annual_elec) / lh2_kg

    return {
        'c_cost': c_cost, 't_cost': t_cost,
        'pfhx5_cost': pfhx5_cost, 'tank_cost': tank_cost, 'ej_cost': ej_cost,
        'capex': capex, 'elec_kw': elec_kw,
        'annual_cap': annual_cap,
        'annual_fixed': annual_fixed,
        'annual_elec': annual_elec,
        'annual_total': annual_cap + annual_fixed + annual_elec,
        'lh2_kg': lh2_kg, 'slc': slc,
    }


# ===========================================================================
# Run
# ===========================================================================
print('='*78)
print('TEXTBOOK COST CORRELATIONS - sanity check on each component')
print('='*78)

# Pre-flight sanity check, baseline values
print('\nBaseline equipment costs (each item, 2024 USD, no Hand factor)')
print('-'*78)
for i, kw in enumerate(SCENARIOS['Baseline']['c_kw'], 1):
    c = compressor_purchased_usd_2024(kw, warn=True)
    print(f"  C-{i} ({kw:>7.1f} kW)   TS-2010 compressor:    {c/1e3:>8.0f} k USD")
for i, kw in enumerate(SCENARIOS['Baseline']['t_kw'], 1):
    t = turbine_purchased_usd_2024(kw, warn=True)
    print(f"  T-{i} ({kw:>7.1f} kW)   Turton-2018 expander:  {t/1e3:>8.0f} k USD")
v = SCENARIOS['Baseline']['pfhx5_m3']
p = pfhx_purchased_usd_2024(v, warn=True)
print(f"  PFHX-5 ({v:>5.1f} m3)   ESDU 97006:                {p/1e3:>8.0f} k USD")
t = lh2_tank_installed_usd_2024(LH2_TANK_M3, warn=True)
print(f"  LH2 tank ({LH2_TANK_M3:>5.0f} m3) Amos 1998 (inst.):     {t/1e3:>8.0f} k USD")

print('\nAdapted equipment costs (each item, 2024 USD, no Hand factor)')
print('-'*78)
for i, kw in enumerate(SCENARIOS['Adapted (Pp=40 bar)']['c_kw'], 1):
    c = compressor_purchased_usd_2024(kw, warn=True)
    print(f"  C-{i} ({kw:>7.1f} kW)   TS-2010 compressor:    {c/1e3:>8.0f} k USD")
for i, kw in enumerate(SCENARIOS['Adapted (Pp=40 bar)']['t_kw'], 1):
    t = turbine_purchased_usd_2024(kw, warn=True)
    print(f"  T-{i} ({kw:>7.1f} kW)   Turton-2018 expander:  {t/1e3:>8.0f} k USD")
v = SCENARIOS['Adapted (Pp=40 bar)']['pfhx5_m3']
p = pfhx_purchased_usd_2024(v, warn=True)
print(f"  PFHX-5 ({v:>5.1f} m3)   ESDU 97006:                {p/1e3:>8.0f} k USD")

b = compute(SCENARIOS['Baseline'], ELEC_BASE_USD_KWH)
a = compute(SCENARIOS['Adapted (Pp=40 bar)'], ELEC_BASE_USD_KWH)


def line(label, vb, va, unit='M USD', fmt='.2f'):
    d = (va - vb) / vb * 100 if vb != 0 else float('nan')
    if unit == 'M USD':
        print(f'{label:<25s}{vb/1e6:>12{fmt}} M$ {va/1e6:>12{fmt}} M$ {d:>+9.1f} %')
    elif unit == 'kg':
        print(f'{label:<25s}{vb:>13,.0f} kg{va:>13,.0f} kg{d:>+9.1f} %')
    elif unit == 'kW':
        print(f'{label:<25s}{vb:>13,.0f} kW{va:>13,.0f} kW{d:>+9.1f} %')
    else:
        print(f'{label:<25s}{vb:>13.3f} {unit}{va:>13.3f} {unit}{d:>+9.1f} %')


print()
print(f"{'='*78}")
print(f"COST BREAKDOWN at electricity price = {ELEC_BASE_USD_KWH:.3f} USD/kWh")
print(f"{'='*78}")
print(f"{'Item':<25s}{'Baseline':>16s}{'Adapted':>16s}{'Delta':>11s}")
print('-' * 78)
print('--- PURCHASE COSTS (2024 USD) ---')
line('C-train (sum 4)', sum(b['c_cost']), sum(a['c_cost']))
line('T-train (sum 2)', sum(b['t_cost']), sum(a['t_cost']))
line('PFHX-5',          b['pfhx5_cost'],  a['pfhx5_cost'])
line('LH2 tank',        b['tank_cost'],   a['tank_cost'])
line('Ejector',         b['ej_cost'],     a['ej_cost'])
print('-' * 78)
line('TOTAL CAPEX',     b['capex'],       a['capex'])
print()
print('--- ANNUAL COSTS (2024 USD/yr) ---')
line('Annual CAPEX',         b['annual_cap'],   a['annual_cap'])
line('Annual fixed OPEX',    b['annual_fixed'], a['annual_fixed'])
line('Annual elec OPEX',     b['annual_elec'],  a['annual_elec'])
line('Annual total',         b['annual_total'], a['annual_total'])
print()
print('--- KPIs ---')
line('Electricity',     b['elec_kw'],     a['elec_kw'],            unit='kW')
line('LH2 production',  b['lh2_kg'],      a['lh2_kg'],             unit='kg')
line('SLC_cryo',        b['slc'],         a['slc'],                unit='/kg')


# ---- Electricity price sensitivity ----------------------------------------
print()
print(f"{'='*78}")
print(f"ELECTRICITY PRICE SENSITIVITY")
print(f"{'='*78}")
print(f"{'Elec price':<15s}{'Baseline SLC':>15s}{'Adapted SLC':>15s}"
      f"{'Delta':>15s}{'Delta %':>10s}")
print(f"{'(USD/kWh)':<15s}{'(USD/kg)':>15s}{'(USD/kg)':>15s}{'(USD/kg)':>15s}")
print('-' * 78)
for p in ELEC_PRICES:
    rb = compute(SCENARIOS['Baseline'], p)
    ra = compute(SCENARIOS['Adapted (Pp=40 bar)'], p)
    d_abs = ra['slc'] - rb['slc']
    d_pct = d_abs / rb['slc'] * 100
    marker = ' <- base' if abs(p - ELEC_BASE_USD_KWH) < 1e-6 else ''
    print(f"{p:<15.3f}{rb['slc']:>15.3f}{ra['slc']:>15.3f}"
          f"{d_abs:>+15.3f}{d_pct:>+10.1f}{marker}")


# ---- Crossover analysis ---------------------------------------------------
def intercept(r):
    return (r['annual_cap'] + r['annual_fixed']) / r['lh2_kg']
def slope(r):
    return r['elec_kw'] * HOURS_PER_YEAR / r['lh2_kg']

rb0 = compute(SCENARIOS['Baseline'], 0.0)
ra0 = compute(SCENARIOS['Adapted (Pp=40 bar)'], 0.0)
int_b = intercept(rb0); int_a = intercept(ra0)
slp_b = slope(rb0);     slp_a = slope(ra0)

print()
print(f"{'='*78}")
print(f"CROSSOVER ANALYSIS")
print(f"{'='*78}")
print(f"  SLC_baseline(price) = {int_b:.5f} + {slp_b:.4f} * price")
print(f"  SLC_adapted(price)  = {int_a:.5f} + {slp_a:.4f} * price")
print(f"  At free electricity (intercepts):")
print(f"    Baseline: {int_b:.5f} USD/kg")
print(f"    Adapted:  {int_a:.5f} USD/kg  ({(int_a-int_b)/int_b*100:+.1f}% vs baseline)")
print(f"  Specific electricity (slopes):")
print(f"    Baseline: {slp_b:.4f} kWh/kg LH2")
print(f"    Adapted:  {slp_a:.4f} kWh/kg LH2  ({(slp_a-slp_b)/slp_b*100:+.1f}% vs baseline)")

if slp_a != slp_b:
    p_cross = (int_b - int_a) / (slp_a - slp_b)
    print(f"  Crossover electricity price (adapted = baseline):")
    print(f"    p* = {p_cross:.5f} USD/kWh  ({p_cross*1000:.1f} USD/MWh)")
    if p_cross < 0:
        print(f"    NEGATIVE - adapted never wins at any realistic price")
    elif p_cross < 0.02:
        print(f"    Below industrial range - adapted only wins at nearly-free elec")
    elif p_cross < 0.20:
        print(f"    Within industrial range - adapted wins if elec < {p_cross*1000:.0f} USD/MWh")
    else:
        print(f"    Above industrial range - adapted wins at most realistic prices")

print()
print(f"{'='*78}")
print(f"BOTTOM LINE at base price {ELEC_BASE_USD_KWH:.3f} USD/kWh:")
print(f"  Baseline SLC_cryo = {b['slc']:.3f}  USD/kg LH2")
print(f"  Adapted  SLC_cryo = {a['slc']:.3f}  USD/kg LH2")
print(f"  Retrofit delta    = {a['slc']-b['slc']:+.3f}  USD/kg "
      f"({(a['slc']-b['slc'])/b['slc']*100:+.1f}%)")
print(f"{'='*78}")
