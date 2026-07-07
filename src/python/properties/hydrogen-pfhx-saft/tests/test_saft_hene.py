"""
Smoke test for the SAFT-VRQ-Mie He-Ne coolant integration.

Validates that:
  1. The saft_hene module loads and computes properties at a single state point.
  2. The HeliumNeonSAFT class works as a drop-in replacement with the same
     public interface as the original HeliumNeon class.
  3. Transport properties are physically reasonable.
  4. Comparison with O'Neill's actual implementation (Xiao/Huber polynomials
     + linear mole-fraction mixing).

Run:
    python tests/test_saft_hene.py
"""

import sys
import os
import numpy as np

# Ensure the source is on the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from hydrogen_pfhx import saft_hene
from hydrogen_pfhx.fluids import HeliumNeonSAFT


# ============================================================================
# O'Neill's original correlations (reproduced exactly from fluids.py)
# ============================================================================
# These are the Xiao et al. (He) and Huber (Ne) polynomial correlations
# that O'Neill used in the original HeliumNeon class.

_T0 = 298.15
_MU_0 = np.array([19.8253, 31.7088]) * 1e-6   # Pa*s
_LAM_0 = np.array([155.0008, 49.1732]) * 1e-3  # W/(m*K)

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


def oneill_pure_viscosities(T):
    """Compute pure He and Ne viscosities using O'Neill's Xiao/Huber polynomials."""
    inds = np.arange(1, _A_COEFFS.shape[0] + 1)
    mu_vals = np.zeros(2)
    for f_idx in (0, 1):
        terms = _A_COEFFS[:, f_idx] * np.log(T / _T0) ** inds
        mu_vals[f_idx] = np.exp(np.sum(terms)) * _MU_0[f_idx]
    return mu_vals  # [mu_He, mu_Ne] in Pa*s


def oneill_pure_conductivities(T):
    """Compute pure He and Ne thermal conductivities using O'Neill's polynomials."""
    inds = np.arange(1, _B_COEFFS.shape[0] + 1)
    lam_vals = np.zeros(2)
    for f_idx in (0, 1):
        terms = _B_COEFFS[:, f_idx] * np.log(T / _T0) ** inds
        lam_vals[f_idx] = np.exp(np.sum(terms)) * _LAM_0[f_idx]
    return lam_vals  # [lam_He, lam_Ne] in W/(m*K)


def oneill_mixture_viscosity(T, x_He):
    """O'Neill's full mixture viscosity: polynomial pure + linear mixing."""
    mu_vals = oneill_pure_viscosities(T)
    return x_He * mu_vals[0] + (1 - x_He) * mu_vals[1]


def oneill_mixture_conductivity(T, x_He):
    """O'Neill's full mixture thermal conductivity: polynomial pure + linear mixing."""
    lam_vals = oneill_pure_conductivities(T)
    return x_He * lam_vals[0] + (1 - x_He) * lam_vals[1]


# ============================================================================
# Tests
# ============================================================================

def test_saft_module_single_point():
    """Test that the SAFT module can evaluate properties at a representative point."""
    print("=" * 70)
    print("TEST 1: saft_hene module -- single state point evaluation")
    print("=" * 70)

    T_K = 40.0
    P_Pa = 500e3
    x_He = 0.80

    print(f"  Conditions: T = {T_K} K, P = {P_Pa / 1e5} bar, x_He = {x_He}")

    state = saft_hene.create_state(T_K, P_Pa, x_He)
    print(f"  State created successfully")

    props = saft_hene.extract_thermodynamic_properties(state, x_He)
    print(f"  Molar density:     {props['molar_density']:.2f} mol/m^3")
    print(f"  Mass density:      {props['mass_density']:.4f} kg/m^3")
    print(f"  Molecular mass:    {props['molecular_mass'] * 1000:.4f} g/mol")
    print(f"  Cp (mass):         {props['specific_heat_capacity']:.2f} J/(kg*K)")
    print(f"  Enthalpy (mass):   {props['enthalpy']:.2f} J/kg")
    print(f"  s_res (molar):     {props['residual_entropy']:.6f} J/(mol*K)")

    x = np.array([x_He, 1.0 - x_He])
    mu = saft_hene.mixture_viscosity(state, x)
    lam = saft_hene.mixture_thermal_conductivity(state, x)

    print(f"  Viscosity:         {mu * 1e6:.4f} uPa*s")
    print(f"  Thermal cond:      {lam * 1e3:.4f} mW/(m*K)")

    assert props['mass_density'] > 0, "Density must be positive"
    assert props['specific_heat_capacity'] > 0, "Cp must be positive"
    assert mu > 0, "Viscosity must be positive"
    assert lam > 0, "Thermal conductivity must be positive"
    assert 3.0e-6 < mu < 20.0e-6, f"Viscosity {mu*1e6:.2f} uPa*s out of expected range"
    assert 5.0e-3 < lam < 50.0e-3, f"Conductivity {lam*1e3:.2f} mW/(m*K) out of expected range"

    print("\n  PASSED\n")


def test_kij_temperature_dependence():
    """Test that k_ij(T) varies correctly with temperature."""
    print("=" * 70)
    print("TEST 2: k_ij(T) -- linear temperature dependence")
    print("=" * 70)

    temperatures = [20, 27, 30, 35, 40, 50, 80]
    print(f"  {'T (K)':>8}  {'k_ij':>10}")
    print("  " + "-" * 22)
    kij_prev = None
    for T in temperatures:
        kij = saft_hene.kij_of_T(T)
        print(f"  {T:8.1f}  {kij:+10.6f}")
        if kij_prev is not None and 27 <= T <= 42:
            assert kij <= kij_prev + 0.001, "k_ij should decrease with T in fitted range"
        kij_prev = kij

    kij_20 = saft_hene.kij_of_T(20.0)
    kij_26_95 = saft_hene.kij_of_T(26.95)
    assert abs(kij_20 - kij_26_95) < 1e-10, "k_ij should be clamped below fitted range"

    kij_80 = saft_hene.kij_of_T(80.0)
    kij_41_9 = saft_hene.kij_of_T(41.9)
    assert abs(kij_80 - kij_41_9) < 1e-10, "k_ij should be clamped above fitted range"

    print("\n  PASSED\n")


def test_helium_neon_saft_class():
    """Test that HeliumNeonSAFT has the full FluidStream-compatible interface."""
    print("=" * 70)
    print("TEST 3: HeliumNeonSAFT class -- interface compatibility")
    print("=" * 70)

    mass_flow_rate = 6.944
    helium_fraction = 0.80
    coolant = HeliumNeonSAFT(mass_flow_rate, helium_fraction)

    assert coolant.mass_flow_rate == mass_flow_rate
    assert coolant.temperature is None
    print("  Instantiation: OK")

    T_inlet = 20.0
    P_inlet = 500e3
    coolant.update_conditions(T_inlet, P_inlet)
    assert coolant.temperature == T_inlet
    assert coolant.pressure == P_inlet
    print(f"  update_conditions({T_inlet} K, {P_inlet/1e3} kPa): OK")

    coolant.set_properties()
    print(f"  set_properties(): OK")
    print(f"    density:        {coolant.mass_density:.4f} kg/m^3")
    print(f"    viscosity:      {coolant.viscosity * 1e6:.4f} uPa*s")
    print(f"    conductivity:   {coolant.thermal_conductivity * 1e3:.4f} mW/(m*K)")
    print(f"    Cp:             {coolant.specific_heat_capacity:.2f} J/(kg*K)")
    print(f"    Prandtl:        {coolant.prandtl_number:.4f}")
    print(f"    enthalpy:       {coolant.enthalpy:.2f} J/kg")
    print(f"    molecular mass: {coolant.molecular_mass * 1000:.4f} g/mol")

    required_attrs = [
        'molecular_mass', 'molar_density', 'mass_density', 'viscosity',
        'specific_heat_capacity', 'thermal_conductivity', 'prandtl_number',
        'enthalpy', 'entropy', 'molar_flow_rate', 'critical_temperature',
        'critical_pressure',
    ]
    for attr in required_attrs:
        val = getattr(coolant, attr)
        assert val is not None, f"Attribute '{attr}' is None after set_properties()"

    print("  All required attributes populated: OK")

    cross_section_area = 0.01
    coolant.calculate_velocity(cross_section_area)
    assert coolant.velocity > 0
    print(f"  calculate_velocity(A={cross_section_area}): v = {coolant.velocity:.3f} m/s")

    char_length = 0.005
    coolant.calculate_reynolds_number(char_length)
    assert coolant.reynolds_number > 0
    print(f"  calculate_reynolds_number(L={char_length}): Re = {coolant.reynolds_number:.1f}")

    T_mid = 50.0
    P_mid = 499e3
    coolant.update_conditions(T_mid, P_mid)
    coolant.set_properties()
    print(f"\n  At T = {T_mid} K:")
    print(f"    viscosity:      {coolant.viscosity * 1e6:.4f} uPa*s")
    print(f"    conductivity:   {coolant.thermal_conductivity * 1e3:.4f} mW/(m*K)")
    print(f"    Cp:             {coolant.specific_heat_capacity:.2f} J/(kg*K)")

    assert coolant.viscosity > 0

    print("\n  PASSED\n")


def test_comparison_with_oneill():
    """
    Compare SAFT-RES transport properties with O'Neill's actual implementation:
    Xiao et al. (He) / Huber (Ne) polynomial correlations + linear mole-fraction
    mixing. This is exactly what the original HeliumNeon class computes.
    """
    print("=" * 70)
    print("TEST 4: Comparison -- SAFT-RES vs O'Neill original")
    print("=" * 70)

    x_He = 0.80
    x = np.array([x_He, 1.0 - x_He])

    print(f"\n  {'T (K)':>6}  {'O\'Neill mu':>12}  {'SAFT mu':>12}  {'diff':>8}"
          f"  {'O\'Neill lam':>12}  {'SAFT lam':>12}  {'diff':>8}")
    print(f"  {'':>6}  {'(uPa*s)':>12}  {'(uPa*s)':>12}  {'(%)':>8}"
          f"  {'(mW/m/K)':>12}  {'(mW/m/K)':>12}  {'(%)':>8}")
    print("  " + "-" * 82)

    for T_K in [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]:
        P_Pa = 500e3

        # O'Neill original (exact reproduction)
        mu_oneill = oneill_mixture_viscosity(T_K, x_He)
        lam_oneill = oneill_mixture_conductivity(T_K, x_He)

        # SAFT-RES (new)
        state = saft_hene.create_state(T_K, P_Pa, x_He)
        mu_saft = saft_hene.mixture_viscosity(state, x)
        lam_saft = saft_hene.mixture_thermal_conductivity(state, x)

        mu_diff = (mu_saft - mu_oneill) / mu_oneill * 100
        lam_diff = (lam_saft - lam_oneill) / lam_oneill * 100

        print(f"  {T_K:6.1f}  {mu_oneill*1e6:12.4f}  {mu_saft*1e6:12.4f}  {mu_diff:+7.2f}%"
              f"  {lam_oneill*1e3:12.4f}  {lam_saft*1e3:12.4f}  {lam_diff:+7.2f}%")

    # Also show pure component values at one temperature for transparency
    T_check = 40.0
    mu_pure = oneill_pure_viscosities(T_check)
    lam_pure = oneill_pure_conductivities(T_check)
    print(f"\n  Pure component values at T = {T_check} K (O'Neill polynomials):")
    print(f"    He viscosity:      {mu_pure[0]*1e6:.4f} uPa*s")
    print(f"    Ne viscosity:      {mu_pure[1]*1e6:.4f} uPa*s")
    print(f"    He conductivity:   {lam_pure[0]*1e3:.4f} mW/(m*K)")
    print(f"    Ne conductivity:   {lam_pure[1]*1e3:.4f} mW/(m*K)")
    print(f"    Linear mix mu:     {x_He*mu_pure[0]*1e6 + (1-x_He)*mu_pure[1]*1e6:.4f} uPa*s")
    print(f"    Linear mix lam:    {x_He*lam_pure[0]*1e3 + (1-x_He)*lam_pure[1]*1e3:.4f} mW/(m*K)")

    print("\n  DONE\n")


if __name__ == '__main__':
    test_saft_module_single_point()
    test_kij_temperature_dependence()
    test_helium_neon_saft_class()
    test_comparison_with_oneill()
    print("=" * 70)
    print("ALL TESTS PASSED")
    print("=" * 70)
