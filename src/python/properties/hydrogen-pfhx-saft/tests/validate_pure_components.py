"""
Pure component transport property validation: O'Neill vs SAFT-RES vs CoolProp.

CoolProp uses the same reference correlations as REFPROP for pure He and Ne
(Arp et al. for He, Rabinovich et al. for Ne), so it serves as a reliable
reference. By comparing pure component values we isolate whether differences
come from the dilute gas correlations or from the mixing rule.

Run:
    python tests/validate_pure_components.py
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from hydrogen_pfhx import saft_hene


# ============================================================================
# O'Neill's original correlations (exact copy from fluids.py)
# ============================================================================

_T0 = 298.15
_MU_0 = np.array([19.8253, 31.7088]) * 1e-6   # Pa*s  [He, Ne]
_LAM_0 = np.array([155.0008, 49.1732]) * 1e-3  # W/(m*K)  [He, Ne]

_A_COEFFS = np.array([
    [6.8257552e-1,  6.75404e-1],
    [1.4496203e-2, -2.03522e-2],
    [1.1987706e-3,  1.61102e-2],
    [-6.7722412e-5, -4.88074e-3],
    [4.9875650e-5,  5.32334e-4],
    [-6.1456994e-6,  2.93695e-4],
    [1.3189407e-6, -1.55155e-4],
    [-3.7245774e-7,  3.10797e-5],
    [1.3671981e-8, -2.50504e-6],
    [5.0354149e-8,  2.74563e-8],
    [-1.5714379e-8, 0.0],
    [1.4720785e-9,  0.0],
])

_B_COEFFS = np.array([
    [6.8192175e-1,  6.76478e-1],
    [1.4441872e-2, -2.13734e-2],
    [1.2138429e-3,  1.63523e-2],
    [-7.4912205e-5, -4.79402e-3],
    [5.2123986e-5,  4.23959e-4],
    [-6.0795764e-6,  3.35013e-4],
    [8.6135146e-7, -1.53714e-4],
    [-2.6311453e-7,  2.45053e-5],
    [6.8367328e-8, -5.03701e-7],
    [1.6608814e-8, -1.67000e-7],
    [-9.0525341e-9, 0.0],
    [9.9986305e-10, 0.0],
])


def oneill_pure_viscosity(T, species_idx):
    """O'Neill polynomial: species_idx 0=He, 1=Ne. Returns Pa*s."""
    inds = np.arange(1, _A_COEFFS.shape[0] + 1)
    terms = _A_COEFFS[:, species_idx] * np.log(T / _T0) ** inds
    return np.exp(np.sum(terms)) * _MU_0[species_idx]


def oneill_pure_conductivity(T, species_idx):
    """O'Neill polynomial: species_idx 0=He, 1=Ne. Returns W/(m*K)."""
    inds = np.arange(1, _B_COEFFS.shape[0] + 1)
    terms = _B_COEFFS[:, species_idx] * np.log(T / _T0) ** inds
    return np.exp(np.sum(terms)) * _LAM_0[species_idx]


# ============================================================================
# SAFT-RES pure component (Chapman-Enskog dilute gas, no residual at low P)
# ============================================================================

def saft_res_pure_viscosity(T, species):
    """SAFT-RES: Chapman-Enskog dilute + residual via entropy scaling."""
    if species == 'He':
        params = saft_hene.HELIUM_VISC_PARAMS
    else:
        params = saft_hene.NEON_VISC_PARAMS
    return saft_hene.chapman_enskog_viscosity(
        T, params['M'], params['sigma'], params['epsilon_k'])


def saft_res_pure_conductivity(T, species):
    """SAFT-RES: polynomial dilute gas thermal conductivity."""
    return saft_hene.dilute_gas_thermal_conductivity(T, species)


# ============================================================================
# CoolProp reference (same correlations as REFPROP)
# ============================================================================

def coolprop_viscosity(T, P_Pa, fluid):
    """Pure component viscosity via CoolProp or REFPROP backend. Returns Pa*s."""
    # Try CoolProp native first, then REFPROP backend
    for backend_fluid in [fluid, f'REFPROP::{fluid}']:
        try:
            return CP.PropsSI('V', 'T', T, 'P', P_Pa, backend_fluid)
        except (ValueError, RuntimeError):
            continue
    return np.nan


def coolprop_conductivity(T, P_Pa, fluid):
    """Pure component thermal conductivity via CoolProp or REFPROP backend. Returns W/(m*K)."""
    for backend_fluid in [fluid, f'REFPROP::{fluid}']:
        try:
            return CP.PropsSI('L', 'T', T, 'P', P_Pa, backend_fluid)
        except (ValueError, RuntimeError):
            continue
    return np.nan


# ============================================================================
# Main validation
# ============================================================================

def validate_and_plot():
    P_Pa = 500e3  # 5 bar (same as PFHX conditions)

    # Check which reference backends are available
    print("Checking reference backends...")
    for fluid in ['Helium', 'Neon']:
        for backend in [fluid, f'REFPROP::{fluid}']:
            try:
                CP.PropsSI('V', 'T', 40, 'P', P_Pa, backend)
                print(f"  {fluid} viscosity:      {backend} OK")
                break
            except Exception:
                if 'REFPROP' in backend:
                    print(f"  {fluid} viscosity:      NOT AVAILABLE (CoolProp + REFPROP both failed)")
        for backend in [fluid, f'REFPROP::{fluid}']:
            try:
                CP.PropsSI('L', 'T', 40, 'P', P_Pa, backend)
                print(f"  {fluid} conductivity:   {backend} OK")
                break
            except Exception:
                if 'REFPROP' in backend:
                    print(f"  {fluid} conductivity:   NOT AVAILABLE (CoolProp + REFPROP both failed)")
    print()

    # Neon melts at ~24.56 K, so CoolProp refuses below that.
    # Use species-specific temperature ranges for CoolProp, but a common
    # range for O'Neill and SAFT-RES (which are dilute gas correlations
    # and don't care about the melting point).
    T_common = np.linspace(25, 80, 50)  # safe for both species

    species_list = [
        {'name': 'Helium', 'idx': 0, 'short': 'He', 'coolprop': 'Helium'},
        {'name': 'Neon',   'idx': 1, 'short': 'Ne', 'coolprop': 'Neon'},
    ]

    # Collect data
    data = {}
    for sp in species_list:
        key = sp['short']
        data[key] = {
            'T': T_common,
            'mu_oneill': [], 'mu_saft': [], 'mu_coolprop': [],
            'lam_oneill': [], 'lam_saft': [], 'lam_coolprop': [],
        }
        for T in T_common:
            data[key]['mu_oneill'].append(oneill_pure_viscosity(T, sp['idx']))
            data[key]['mu_saft'].append(saft_res_pure_viscosity(T, sp['short']))
            data[key]['mu_coolprop'].append(coolprop_viscosity(T, P_Pa, sp['coolprop']))
            data[key]['lam_oneill'].append(oneill_pure_conductivity(T, sp['idx']))
            data[key]['lam_saft'].append(saft_res_pure_conductivity(T, sp['short']))
            data[key]['lam_coolprop'].append(coolprop_conductivity(T, P_Pa, sp['coolprop']))

        for k in data[key]:
            data[key][k] = np.array(data[key][k])

    # ---- Print table at key temperatures ----
    print("=" * 90)
    print("PURE COMPONENT TRANSPORT PROPERTY VALIDATION")
    print("Reference: CoolProp (same correlations as REFPROP)")
    print(f"Conditions: P = {P_Pa/1e5:.0f} bar")
    print("=" * 90)

    for sp in species_list:
        key = sp['short']
        print(f"\n--- {sp['name']} ---")
        print(f"  {'T(K)':>5}  {'O\'Neill':>10} {'SAFT-RES':>10} {'CoolProp':>10}"
              f"  {'ON err%':>8} {'SAFT err%':>9}")

        # Viscosity
        print(f"\n  Viscosity (uPa*s):")
        for T_idx, T in enumerate([25, 30, 40, 50, 60, 70, 80]):
            idx = np.argmin(np.abs(T_common - T))
            mu_on = data[key]['mu_oneill'][idx] * 1e6
            mu_sf = data[key]['mu_saft'][idx] * 1e6
            mu_cp = data[key]['mu_coolprop'][idx] * 1e6
            if np.isnan(mu_cp):
                print(f"  {T:5.0f}  {mu_on:10.4f} {mu_sf:10.4f} {'N/A':>10}"
                      f"  {'N/A':>8} {'N/A':>9}")
            else:
                err_on = (mu_on - mu_cp) / mu_cp * 100
                err_sf = (mu_sf - mu_cp) / mu_cp * 100
                print(f"  {T:5.0f}  {mu_on:10.4f} {mu_sf:10.4f} {mu_cp:10.4f}"
                      f"  {err_on:+8.2f}% {err_sf:+9.2f}%")

        # Thermal conductivity
        print(f"\n  Thermal conductivity (mW/(m*K)):")
        for T_idx, T in enumerate([25, 30, 40, 50, 60, 70, 80]):
            idx = np.argmin(np.abs(T_common - T))
            lam_on = data[key]['lam_oneill'][idx] * 1e3
            lam_sf = data[key]['lam_saft'][idx] * 1e3
            lam_cp = data[key]['lam_coolprop'][idx] * 1e3
            if np.isnan(lam_cp):
                print(f"  {T:5.0f}  {lam_on:10.4f} {lam_sf:10.4f} {'N/A':>10}"
                      f"  {'N/A':>8} {'N/A':>9}")
            else:
                err_on = (lam_on - lam_cp) / lam_cp * 100
                err_sf = (lam_sf - lam_cp) / lam_cp * 100
                print(f"  {T:5.0f}  {lam_on:10.4f} {lam_sf:10.4f} {lam_cp:10.4f}"
                      f"  {err_on:+8.2f}% {err_sf:+9.2f}%")

    # ---- Plots ----
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.suptitle('Pure Component Transport Properties: O\'Neill vs SAFT-RES vs CoolProp',
                 fontsize=14, fontweight='bold')

    colors = {'oneill': '#d62728', 'saft': '#1f77b4', 'coolprop': '#2ca02c'}
    styles = {'oneill': '--', 'saft': '-.', 'coolprop': '-'}
    labels = {'oneill': "O'Neill (Xiao/Huber)", 'saft': 'SAFT-RES (Chapman-Enskog)',
              'coolprop': 'CoolProp (reference)'}

    for col, sp in enumerate(species_list):
        key = sp['short']

        # Viscosity (top row)
        ax = axes[0, col]
        ax.plot(T_common, data[key]['mu_coolprop'] * 1e6, styles['coolprop'],
                color=colors['coolprop'], linewidth=2.2, label=labels['coolprop'])
        ax.plot(T_common, data[key]['mu_oneill'] * 1e6, styles['oneill'],
                color=colors['oneill'], linewidth=1.8, label=labels['oneill'])
        ax.plot(T_common, data[key]['mu_saft'] * 1e6, styles['saft'],
                color=colors['saft'], linewidth=1.8, label=labels['saft'])
        ax.set_xlabel('Temperature (K)', fontsize=11)
        ax.set_ylabel(r'Viscosity ($\mu$Pa$\cdot$s)', fontsize=11)
        ax.set_title(f'{sp["name"]} - Viscosity', fontsize=12)
        ax.legend(fontsize=9, frameon=False)
        ax.set_xlim(25, 80)
        ax.set_ylim(bottom=0)

        # Thermal conductivity (bottom row)
        ax = axes[1, col]
        ax.plot(T_common, data[key]['lam_coolprop'] * 1e3, styles['coolprop'],
                color=colors['coolprop'], linewidth=2.2, label=labels['coolprop'])
        ax.plot(T_common, data[key]['lam_oneill'] * 1e3, styles['oneill'],
                color=colors['oneill'], linewidth=1.8, label=labels['oneill'])
        ax.plot(T_common, data[key]['lam_saft'] * 1e3, styles['saft'],
                color=colors['saft'], linewidth=1.8, label=labels['saft'])
        ax.set_xlabel('Temperature (K)', fontsize=11)
        ax.set_ylabel(r'$\lambda$ (mW$\cdot$m$^{-1}\cdot$K$^{-1}$)', fontsize=11)
        ax.set_title(f'{sp["name"]} - Thermal conductivity', fontsize=12)
        ax.legend(fontsize=9, frameon=False)
        ax.set_xlim(25, 80)
        ax.set_ylim(bottom=0)

    fig.tight_layout()
    fig.savefig('pure_component_validation.png', dpi=200, bbox_inches='tight')

    # ---- Relative error plots ----
    fig2, axes2 = plt.subplots(2, 2, figsize=(13, 10))
    fig2.suptitle('Pure Component Errors vs CoolProp Reference',
                  fontsize=14, fontweight='bold')

    for col, sp in enumerate(species_list):
        key = sp['short']

        # Viscosity error (top row)
        ax = axes2[0, col]
        err_on = (data[key]['mu_oneill'] - data[key]['mu_coolprop']) / data[key]['mu_coolprop'] * 100
        err_sf = (data[key]['mu_saft'] - data[key]['mu_coolprop']) / data[key]['mu_coolprop'] * 100
        ax.axhline(0, color=colors['coolprop'], linewidth=1.5, label='CoolProp (reference)')
        ax.plot(T_common, err_on, styles['oneill'], color=colors['oneill'],
                linewidth=2, label=labels['oneill'])
        ax.plot(T_common, err_sf, styles['saft'], color=colors['saft'],
                linewidth=2, label=labels['saft'])
        ax.set_xlabel('Temperature (K)', fontsize=11)
        ax.set_ylabel('Relative error vs CoolProp (%)', fontsize=11)
        ax.set_title(f'{sp["name"]} - Viscosity error', fontsize=12)
        ax.legend(fontsize=9, frameon=False)
        ax.set_xlim(25, 80)

        # Thermal conductivity error (bottom row)
        ax = axes2[1, col]
        err_on = (data[key]['lam_oneill'] - data[key]['lam_coolprop']) / data[key]['lam_coolprop'] * 100
        err_sf = (data[key]['lam_saft'] - data[key]['lam_coolprop']) / data[key]['lam_coolprop'] * 100
        ax.axhline(0, color=colors['coolprop'], linewidth=1.5, label='CoolProp (reference)')
        ax.plot(T_common, err_on, styles['oneill'], color=colors['oneill'],
                linewidth=2, label=labels['oneill'])
        ax.plot(T_common, err_sf, styles['saft'], color=colors['saft'],
                linewidth=2, label=labels['saft'])
        ax.set_xlabel('Temperature (K)', fontsize=11)
        ax.set_ylabel('Relative error vs CoolProp (%)', fontsize=11)
        ax.set_title(f'{sp["name"]} - Thermal conductivity error', fontsize=12)
        ax.legend(fontsize=9, frameon=False)
        ax.set_xlim(25, 80)

    fig2.tight_layout()
    fig2.savefig('pure_component_errors.png', dpi=200, bbox_inches='tight')

    plt.show()
    print("\nFigures saved: pure_component_validation.png, pure_component_errors.png")


if __name__ == '__main__':
    validate_and_plot()
