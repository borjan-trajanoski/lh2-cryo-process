"""
compare_models.py
=================
Run the hydrogen PFHX model with both coolant property frameworks and compare
reactor-level results:

    1. O'Neill baseline:  CoolProp HEOS + Tkaczuk, Xiao/Huber polynomials,
                          linear mole-fraction mixing
    2. SAFT-RES:          feos SAFT-VRQ-Mie, Chapman-Enskog + Li et al. (2024)
                          residual entropy scaling, Wilke/Wassiljewa mixing

Usage:
    python compare_models.py

Output:
    - Console summary table with key reactor-level metrics
    - output/comparison_temperature.png
    - output/comparison_pressure.png
    - output/comparison_velocity.png
    - output/comparison_conversion.png
    - output/comparison_summary.png   (combined 4-panel figure)
    - output/oneill_results.csv
    - output/saft_results.csv
"""

import os
import sys
import time
import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Ensure the package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from hydrogen_pfhx import model, utils


# ============================================================================
# Configuration
# ============================================================================

CONFIG_DIR = os.path.join(os.path.dirname(__file__), 'tests')
ONEILL_CONFIG = os.path.join(CONFIG_DIR, 'helium_neon_oneill.yaml')
SAFT_CONFIG = os.path.join(CONFIG_DIR, 'helium_neon_configuration.yaml')
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), 'output')


# ============================================================================
# Plot styling (matching O'Neill thesis style)
# ============================================================================

def setup_plot_style():
    """Configure matplotlib for publication-quality figures."""
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman'],
        'font.size': 11,
        'axes.titlesize': 12,
        'axes.labelsize': 11,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'mathtext.fontset': 'custom',
        'mathtext.it': 'Times New Roman:italic',
    })


# ============================================================================
# Main comparison logic
# ============================================================================

def run_simulation(config_path, label):
    """Run a single simulation and return results with timing."""
    print(f"\n{'='*60}")
    print(f"Running: {label}")
    print(f"Config:  {os.path.basename(config_path)}")
    print(f"{'='*60}")

    config = utils.load_config(config_path)
    t_start = time.time()
    results = model.model(config)
    elapsed = time.time() - t_start

    converged = results.attrs.get('converged', 'unknown')
    solver_msg = results.attrs.get('solver_message', '')
    coolant_model = results.attrs.get('coolant_model', label)

    print(f"\nCompleted in {elapsed:.1f}s")
    print(f"Converged: {converged}")
    print(f"Solver:    {solver_msg}")
    print(f"Coolant:   {coolant_model}")
    print(f"Nodes:     {len(results)}")

    return results, elapsed


def extract_metrics(results, label):
    """Extract key reactor-level metrics from a results DataFrame."""
    z = results['Z (m)'].values
    T_r = results['Reactant temperature (K)'].values
    T_c = results['Coolant temperature (K)'].values
    P_r = results['Reactant pressure (kPa)'].values
    P_c = results['Coolant pressure (kPa)'].values
    u_r = results['Reactant velocity (m/s)'].values
    u_c = results['Coolant velocity (m/s)'].values
    xp = results['Actual para-hydrogen fraction (mol/mol)'].values
    xp_eq = results['Equilibrium para-hydrogen fraction (mol/mol)'].values

    metrics = {
        'Model': label,
        # Temperature
        'T_reactant_in (K)': T_r[0],
        'T_reactant_out (K)': T_r[-1],
        'T_coolant_in (K)': T_c[-1],   # counter-current: coolant enters at z=L
        'T_coolant_out (K)': T_c[0],
        'dT_approach (K)': T_r[-1] - T_c[-1],
        'LMTD (K)': _lmtd(T_r[0] - T_c[0], T_r[-1] - T_c[-1]),
        # Pressure
        'dP_reactant (kPa)': P_r[-1] - P_r[0],
        'dP_coolant (kPa)': P_c[0] - P_c[-1],
        # Conversion
        'x_para_in': xp[0],
        'x_para_out': xp[-1],
        'x_para_equil_out': xp_eq[-1],
        'conversion_efficiency': (xp[-1] - xp[0]) / (xp_eq[-1] - xp[0]) if (xp_eq[-1] - xp[0]) > 0 else 0,
        # Velocity
        'u_reactant_avg (m/s)': np.mean(u_r),
        'u_coolant_avg (m/s)': np.mean(u_c),
    }
    return metrics


def _lmtd(dT1, dT2):
    """Log-mean temperature difference, handling dT1 ~= dT2."""
    if abs(dT1 - dT2) < 1e-6:
        return (dT1 + dT2) / 2.0
    return (dT1 - dT2) / np.log(dT1 / dT2)


# ============================================================================
# Plotting
# ============================================================================

ONEILL_COLOR = (0.8, 0.3, 0.2)
SAFT_COLOR = (0.2, 0.3, 0.8)
EQUIL_COLOR = 'k'


def plot_temperature_comparison(res_on, res_sf, save_path=None):
    """Temperature profiles: reactant and coolant for both models."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    # Reactant temperature
    ax1.plot(res_on['Z (m)'], res_on['Reactant temperature (K)'],
             '-', color=ONEILL_COLOR, linewidth=1.5, label="O'Neill")
    ax1.plot(res_sf['Z (m)'], res_sf['Reactant temperature (K)'],
             '--', color=SAFT_COLOR, linewidth=1.5, label='SAFT-RES')
    ax1.set_xlabel('Length along reactor (m)')
    ax1.set_ylabel('Temperature (K)')
    ax1.set_title('Reactant (H$_2$) temperature')
    ax1.legend()
    _style_axis(ax1)

    # Coolant temperature
    ax2.plot(res_on['Z (m)'], res_on['Coolant temperature (K)'],
             '-', color=ONEILL_COLOR, linewidth=1.5, label="O'Neill")
    ax2.plot(res_sf['Z (m)'], res_sf['Coolant temperature (K)'],
             '--', color=SAFT_COLOR, linewidth=1.5, label='SAFT-RES')
    ax2.set_xlabel('Length along reactor (m)')
    ax2.set_ylabel('Temperature (K)')
    ax2.set_title('Coolant (He-Ne) temperature')
    ax2.legend()
    _style_axis(ax2)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path)
        print(f"  Saved: {save_path}")
    return fig


def plot_pressure_comparison(res_on, res_sf, save_path=None):
    """Pressure profiles for both models."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    ax1.plot(res_on['Z (m)'], res_on['Reactant pressure (kPa)'],
             '-', color=ONEILL_COLOR, linewidth=1.5, label="O'Neill")
    ax1.plot(res_sf['Z (m)'], res_sf['Reactant pressure (kPa)'],
             '--', color=SAFT_COLOR, linewidth=1.5, label='SAFT-RES')
    ax1.set_xlabel('Length along reactor (m)')
    ax1.set_ylabel('Pressure (kPa)')
    ax1.set_title('Reactant pressure')
    ax1.legend()
    _style_axis(ax1)

    ax2.plot(res_on['Z (m)'], res_on['Coolant pressure (kPa)'],
             '-', color=ONEILL_COLOR, linewidth=1.5, label="O'Neill")
    ax2.plot(res_sf['Z (m)'], res_sf['Coolant pressure (kPa)'],
             '--', color=SAFT_COLOR, linewidth=1.5, label='SAFT-RES')
    ax2.set_xlabel('Length along reactor (m)')
    ax2.set_ylabel('Pressure (kPa)')
    ax2.set_title('Coolant pressure')
    ax2.legend()
    _style_axis(ax2)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path)
        print(f"  Saved: {save_path}")
    return fig


def plot_velocity_comparison(res_on, res_sf, save_path=None):
    """Velocity profiles for both models."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    ax1.plot(res_on['Z (m)'], res_on['Reactant velocity (m/s)'],
             '-', color=ONEILL_COLOR, linewidth=1.5, label="O'Neill")
    ax1.plot(res_sf['Z (m)'], res_sf['Reactant velocity (m/s)'],
             '--', color=SAFT_COLOR, linewidth=1.5, label='SAFT-RES')
    ax1.set_xlabel('Length along reactor (m)')
    ax1.set_ylabel('Velocity (m/s)')
    ax1.set_title('Reactant velocity')
    ax1.legend()
    _style_axis(ax1)

    ax2.plot(res_on['Z (m)'], res_on['Coolant velocity (m/s)'],
             '-', color=ONEILL_COLOR, linewidth=1.5, label="O'Neill")
    ax2.plot(res_sf['Z (m)'], res_sf['Coolant velocity (m/s)'],
             '--', color=SAFT_COLOR, linewidth=1.5, label='SAFT-RES')
    ax2.set_xlabel('Length along reactor (m)')
    ax2.set_ylabel('Velocity (m/s)')
    ax2.set_title('Coolant velocity')
    ax2.legend()
    _style_axis(ax2)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path)
        print(f"  Saved: {save_path}")
    return fig


def plot_conversion_comparison(res_on, res_sf, save_path=None):
    """Para-hydrogen conversion profiles for both models."""
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))

    ax.plot(res_on['Z (m)'], res_on['Actual para-hydrogen fraction (mol/mol)'],
            '-', color=ONEILL_COLOR, linewidth=1.5, label="O'Neill")
    ax.plot(res_sf['Z (m)'], res_sf['Actual para-hydrogen fraction (mol/mol)'],
            '--', color=SAFT_COLOR, linewidth=1.5, label='SAFT-RES')
    ax.plot(res_on['Z (m)'], res_on['Equilibrium para-hydrogen fraction (mol/mol)'],
            ':', color=EQUIL_COLOR, linewidth=1.0, label="Equilibrium (O'Neill)")
    ax.plot(res_sf['Z (m)'], res_sf['Equilibrium para-hydrogen fraction (mol/mol)'],
            '-.', color='gray', linewidth=1.0, label='Equilibrium (SAFT-RES)')

    ax.set_xlabel('Length along reactor (m)')
    ax.set_ylabel('Para-H$_2$ fraction (mol/mol)')
    ax.set_title('Ortho-para conversion')
    ax.set_ylim(0, 1)
    ax.legend(loc='lower right')
    _style_axis(ax)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path)
        print(f"  Saved: {save_path}")
    return fig


def plot_summary(res_on, res_sf, save_path=None):
    """Combined 4-panel comparison figure for the presentation."""
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    # (a) Temperature
    ax = axes[0, 0]
    ax.plot(res_on['Z (m)'], res_on['Reactant temperature (K)'],
            '-', color=ONEILL_COLOR, linewidth=1.5, label="Reactant (O'Neill)")
    ax.plot(res_sf['Z (m)'], res_sf['Reactant temperature (K)'],
            '--', color=SAFT_COLOR, linewidth=1.5, label='Reactant (SAFT-RES)')
    ax.plot(res_on['Z (m)'], res_on['Coolant temperature (K)'],
            '-', color=ONEILL_COLOR, linewidth=1.0, alpha=0.5)
    ax.plot(res_sf['Z (m)'], res_sf['Coolant temperature (K)'],
            '--', color=SAFT_COLOR, linewidth=1.0, alpha=0.5)
    ax.set_xlabel('Length (m)')
    ax.set_ylabel('Temperature (K)')
    ax.set_title('(a) Temperature profiles')
    ax.legend(fontsize=8)
    _style_axis(ax)

    # (b) Conversion
    ax = axes[0, 1]
    ax.plot(res_on['Z (m)'], res_on['Actual para-hydrogen fraction (mol/mol)'],
            '-', color=ONEILL_COLOR, linewidth=1.5, label="O'Neill")
    ax.plot(res_sf['Z (m)'], res_sf['Actual para-hydrogen fraction (mol/mol)'],
            '--', color=SAFT_COLOR, linewidth=1.5, label='SAFT-RES')
    ax.plot(res_on['Z (m)'], res_on['Equilibrium para-hydrogen fraction (mol/mol)'],
            ':', color=EQUIL_COLOR, linewidth=0.8, label='Equilibrium')
    ax.set_xlabel('Length (m)')
    ax.set_ylabel('Para-H$_2$ fraction')
    ax.set_title('(b) Ortho-para conversion')
    ax.set_ylim(0, 1)
    ax.legend(fontsize=8, loc='lower right')
    _style_axis(ax)

    # (c) Coolant velocity
    ax = axes[1, 0]
    ax.plot(res_on['Z (m)'], res_on['Coolant velocity (m/s)'],
            '-', color=ONEILL_COLOR, linewidth=1.5, label="O'Neill")
    ax.plot(res_sf['Z (m)'], res_sf['Coolant velocity (m/s)'],
            '--', color=SAFT_COLOR, linewidth=1.5, label='SAFT-RES')
    ax.set_xlabel('Length (m)')
    ax.set_ylabel('Velocity (m/s)')
    ax.set_title('(c) Coolant velocity')
    ax.legend(fontsize=8)
    _style_axis(ax)

    # (d) Pressure drop
    ax = axes[1, 1]
    ax.plot(res_on['Z (m)'], res_on['Coolant pressure (kPa)'],
            '-', color=ONEILL_COLOR, linewidth=1.5, label="O'Neill")
    ax.plot(res_sf['Z (m)'], res_sf['Coolant pressure (kPa)'],
            '--', color=SAFT_COLOR, linewidth=1.5, label='SAFT-RES')
    ax.set_xlabel('Length (m)')
    ax.set_ylabel('Pressure (kPa)')
    ax.set_title('(d) Coolant pressure')
    ax.legend(fontsize=8)
    _style_axis(ax)

    fig.suptitle(
        "PFHX Reactor Comparison: O'Neill vs SAFT-VRQ-Mie + RES\n"
        "He-Ne coolant (x$_{He}$ = 0.80), 100 tpd H$_2$",
        fontsize=13, fontweight='bold', y=1.02)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path)
        print(f"  Saved: {save_path}")
    return fig


def _style_axis(ax):
    """Apply consistent axis styling."""
    ax.minorticks_on()
    ax.tick_params(direction='in', which='minor', length=2,
                   bottom=True, top=True, left=True, right=True)
    ax.tick_params(direction='in', which='major', length=4,
                   bottom=True, top=True, left=True, right=True)


# ============================================================================
# Summary table
# ============================================================================

def print_comparison_table(metrics_on, metrics_sf):
    """Print a side-by-side comparison table to console."""
    print("\n")
    print("=" * 72)
    print("REACTOR-LEVEL COMPARISON: O'Neill vs SAFT-RES")
    print("=" * 72)
    print(f"{'Metric':<35} {'O Neill':>15} {'SAFT-RES':>15}")
    print("-" * 72)

    # Pair up metrics for side-by-side display
    skip_keys = {'Model'}
    for key in metrics_on:
        if key in skip_keys:
            continue
        v_on = metrics_on[key]
        v_sf = metrics_sf[key]

        # Format based on magnitude
        if isinstance(v_on, (int, float)):
            if abs(v_on) < 1:
                fmt = f"{v_on:>15.4f}"
                fmt_sf = f"{v_sf:>15.4f}"
            else:
                fmt = f"{v_on:>15.2f}"
                fmt_sf = f"{v_sf:>15.2f}"
        else:
            fmt = f"{v_on!s:>15}"
            fmt_sf = f"{v_sf!s:>15}"

        print(f"{key:<35} {fmt} {fmt_sf}")

    # Relative differences for key quantities
    print("-" * 72)
    print("RELATIVE DIFFERENCES (SAFT-RES vs O'Neill):")
    diff_keys = [
        ('T_reactant_out (K)', 'Reactant outlet T'),
        ('T_coolant_out (K)', 'Coolant outlet T'),
        ('LMTD (K)', 'LMTD'),
        ('dP_reactant (kPa)', 'Reactant dP'),
        ('dP_coolant (kPa)', 'Coolant dP'),
        ('x_para_out', 'Para-H2 outlet fraction'),
        ('conversion_efficiency', 'Conversion efficiency'),
        ('u_coolant_avg (m/s)', 'Coolant avg velocity'),
    ]
    for key, label in diff_keys:
        v_on = metrics_on[key]
        v_sf = metrics_sf[key]
        if abs(v_on) > 1e-10:
            rel_diff = (v_sf - v_on) / abs(v_on) * 100
            print(f"  {label:<30} {rel_diff:+.2f}%")
        else:
            print(f"  {label:<30} N/A (baseline ~0)")

    print("=" * 72)


# ============================================================================
# Entry point
# ============================================================================

def main():
    setup_plot_style()
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Run both simulations
    print("\n" + "#" * 60)
    print("# PFHX MODEL COMPARISON: O'Neill vs SAFT-VRQ-Mie + RES")
    print("#" * 60)

    results_oneill, t_on = run_simulation(ONEILL_CONFIG, "O'Neill (baseline)")
    results_saft, t_sf = run_simulation(SAFT_CONFIG, "SAFT-VRQ-Mie + RES")

    # Save raw results
    results_oneill.to_csv(os.path.join(OUTPUT_DIR, 'oneill_results.csv'), index=False)
    results_saft.to_csv(os.path.join(OUTPUT_DIR, 'saft_results.csv'), index=False)
    print(f"\nSaved CSV results to {OUTPUT_DIR}/")

    # Extract metrics
    metrics_on = extract_metrics(results_oneill, "O'Neill")
    metrics_sf = extract_metrics(results_saft, "SAFT-RES")

    # Print comparison
    print_comparison_table(metrics_on, metrics_sf)
    print(f"\nSolver time: O'Neill = {t_on:.1f}s, SAFT-RES = {t_sf:.1f}s")

    # Generate plots
    print("\nGenerating comparison plots...")
    plot_temperature_comparison(
        results_oneill, results_saft,
        os.path.join(OUTPUT_DIR, 'comparison_temperature.png'))
    plot_pressure_comparison(
        results_oneill, results_saft,
        os.path.join(OUTPUT_DIR, 'comparison_pressure.png'))
    plot_velocity_comparison(
        results_oneill, results_saft,
        os.path.join(OUTPUT_DIR, 'comparison_velocity.png'))
    plot_conversion_comparison(
        results_oneill, results_saft,
        os.path.join(OUTPUT_DIR, 'comparison_conversion.png'))
    plot_summary(
        results_oneill, results_saft,
        os.path.join(OUTPUT_DIR, 'comparison_summary.png'))

    print("\nDone. All results in output/")
    plt.show()


if __name__ == '__main__':
    main()
