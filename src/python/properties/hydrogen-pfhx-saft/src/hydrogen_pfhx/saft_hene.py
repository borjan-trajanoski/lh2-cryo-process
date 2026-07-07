"""
SAFT-VRQ-Mie equation of state and Residual Entropy Scaling (RES) transport
properties for Helium-Neon mixtures.

Based on:
  - Aasen et al. (2020) SAFT-VRQ-Mie pure-component parameters
  - Li et al. (2024) residual entropy scaling for viscosity and thermal
    conductivity, with Chapman-Enskog dilute gas, Wilke/Wassiljewa mixing
  - VLE-optimised temperature-dependent k_ij(T) for He-Ne

References
----------
Aasen A, Hammer M, Erber S, et al. J Chem Phys 2020;152:074507.
Li Z, Aasen A, Wilhelmsen O, et al. J Chem Eng Data 2024.
"""

import json
import os
import tempfile

import numpy as np

try:
    import feos
    import si_units as si
    FEOS_AVAILABLE = True
except ImportError:
    FEOS_AVAILABLE = False


# ============================================================================
# Physical constants (CODATA 2018)
# ============================================================================
K_B = 1.380649e-23       # Boltzmann constant [J/K]
N_A = 6.02214076e23      # Avogadro number [1/mol]
R_GAS = 8.314462618      # Gas constant [J/(mol*K)]


# ============================================================================
# Li et al. (2024) entropy scaling parameters
# ============================================================================

# -- Viscosity --
# Helium: Group 1 (light gases with quantum effects)
HELIUM_VISC_PARAMS = {
    'has_specific': False,
    'M': 4.002602,        # g/mol
    'sigma': 2.7443,      # Angstrom
    'epsilon_k': 5.4195,  # K
    'xi': 1.6627,
    'ng1': -0.449854,
    'ng2': 3.219854,
    'ng3': -5.298638,
    'ng4': 2.975827,
}

# Neon: Group 2 (noble gases) -- fluid-specific parameters
NEON_VISC_PARAMS = {
    'has_specific': True,
    'M': 20.1797,
    'sigma': 2.7778,
    'epsilon_k': 37.501,
    'xi': 0.8513,
    'n1': -0.053490,
    'n2': 0.212844,
    'n3': 0.249905,
    'n4': -0.065619,
}

# -- Thermal conductivity --
# Helium: Group 1
HELIUM_THERMAL_PARAMS = {
    'has_specific': False,
    'M': 4.002602,
    'sigma': 2.7443,
    'epsilon_k': 5.4195,
    'xi': 1.2966,
    'ng1': 2.391631,
    'ng2': -8.1473,
    'ng3': 12.52226,
    'ng4': -4.38311,
}

# Neon: Group 2
NEON_THERMAL_PARAMS = {
    'has_specific': False,
    'M': 20.1797,
    'sigma': 2.7778,
    'epsilon_k': 37.501,
    'xi': 1.0,
    'ng1': 2.173335,
    'ng2': -4.8767,
    'ng3': 5.754321,
    'ng4': -1.18193,
}


# ============================================================================
# Dilute gas thermal conductivity polynomial coefficients
# from Li et al. (2024) supporting information (Dilute_gas_TC.txt)
# lambda_0 / [W/(m*K)] = n0*T^4 + n1*T^3 + n2*T^2 + n3*T + n4
# ============================================================================
_DILUTE_TC_COEFFS = {
    'He': (-3.378269e-14, 1.709281e-10, -3.347972e-07, 5.381719e-04, 1.925738e-02),
    'Ne': (-3.046444e-13, 6.061412e-10, -4.717322e-07, 2.600919e-04, -1.659178e-04),
}


# ============================================================================
# EoS construction
# ============================================================================

# Module-level cache to avoid rebuilding EoS for repeated k_ij values
_eos_cache = {}


def _data_dir():
    """Return the path to the bundled data directory."""
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')


def _load_kij_fit():
    """Load k_ij(T) polynomial coefficients from JSON."""
    kij_path = os.path.join(_data_dir(), 'kij_fit_coefficients.json')
    with open(kij_path, 'r') as f:
        data = json.load(f)

    order = data['active_order']
    coeffs = np.array(data[f'poly{order}_coeffs'])
    fit_range = data['fitted_range_K']
    return np.poly1d(coeffs), fit_range


# Load once at import time
_kij_poly, _kij_fit_range = _load_kij_fit()


def kij_of_T(T_K):
    """
    Temperature-dependent binary interaction parameter k_ij(T).

    Uses the polynomial within its fitted range. Beyond the fitted
    range the value is clamped to the boundary value to prevent
    extrapolation blow-up.

    Parameters
    ----------
    T_K : float
        Temperature in Kelvin.

    Returns
    -------
    float
        k_ij value at the given temperature.
    """
    T_clamped = np.clip(T_K, _kij_fit_range[0], _kij_fit_range[1])
    return float(_kij_poly(T_clamped))


def build_eos(kij, params_json=None, binary_json=None):
    """
    Build a SAFT-VRQ-Mie mixture EoS for He-Ne with a specific k_ij.

    Uses a temporary JSON file to pass k_ij to feos, working around
    the API limitation where Parameters.new_binary stores k_ij but the
    SAFT-VRQ-Mie model ignores it unless loaded from a JSON file.

    Parameters
    ----------
    kij : float
        Binary interaction parameter.
    params_json : str, optional
        Path to the pure-component parameter file. Defaults to bundled file.
    binary_json : str, optional
        Unused; kept for API compatibility.

    Returns
    -------
    feos.EquationOfState
        SAFT-VRQ-Mie EoS instance.
    """
    if params_json is None:
        params_json = os.path.join(_data_dir(), 'parameters.json')

    binary_record = [{
        "id1": {"cas": "7440-59-7", "name": "helium"},
        "id2": {"cas": "7440-01-9", "name": "neon"},
        "k_ij": kij,
        "l_ij": 0.0,
    }]
    tmp_path = os.path.join(tempfile.gettempdir(), f"binary_kij_{kij:.6f}.json")
    with open(tmp_path, 'w') as f:
        json.dump(binary_record, f)

    params = feos.Parameters.from_json(
        ["helium", "neon"], params_json, binary_path=tmp_path
    )
    return feos.EquationOfState.saftvrqmie(params)


def get_eos_at_T(T_K):
    """
    Build or retrieve from cache a mixture EoS using k_ij(T).

    Parameters
    ----------
    T_K : float
        Temperature in Kelvin.

    Returns
    -------
    feos.EquationOfState
        Cached SAFT-VRQ-Mie EoS instance.
    """
    kij = kij_of_T(T_K)
    cache_key = round(kij, 6)
    if cache_key not in _eos_cache:
        _eos_cache[cache_key] = build_eos(kij)
    return _eos_cache[cache_key]


# ============================================================================
# Thermodynamic state evaluation
# ============================================================================

def create_state(T_K, P_Pa, x_He):
    """
    Create a feos thermodynamic State for a He-Ne mixture.

    Parameters
    ----------
    T_K : float
        Temperature [K].
    P_Pa : float
        Pressure [Pa].
    x_He : float
        Mole fraction of helium (0 to 1).

    Returns
    -------
    feos.State
        Thermodynamic state object.
    """
    eos = get_eos_at_T(T_K)
    molefracs = np.array([x_He, 1.0 - x_He])
    return feos.State(
        eos,
        temperature=T_K * si.KELVIN,
        pressure=P_Pa * si.PASCAL,
        molefracs=molefracs,
    )


def extract_thermodynamic_properties(state, x_He):
    """
    Extract all thermodynamic properties needed by the PFHX model from a
    feos State, returned as plain floats in SI units.

    SAFT-VRQ-Mie in feos does not ship an ideal gas model, so calling
    total-property methods (specific_enthalpy, specific_isobaric_heat_capacity,
    etc.) raises a PanicException.  We work around this by fetching the
    *residual* contributions from feos and adding the ideal gas part
    analytically.  He and Ne are both monatomic, so:

        Cp_ig = 5/2 R   (exact, no temperature dependence)
        H_ig(T) = 5/2 R T   (relative to 0 K)

    Parameters
    ----------
    state : feos.State
        Thermodynamic state.
    x_He : float
        Mole fraction of helium.

    Returns
    -------
    dict
        Keys: molar_density [mol/m^3], mass_density [kg/m^3],
        molecular_mass [kg/mol], specific_heat_capacity [J/(kg*K)],
        enthalpy [J/kg], entropy [J/(kg*K)], residual_entropy [J/(mol*K)].
    """
    x_Ne = 1.0 - x_He
    M_He = HELIUM_VISC_PARAMS['M'] / 1000.0  # kg/mol
    M_Ne = NEON_VISC_PARAMS['M'] / 1000.0    # kg/mol
    M_mix = x_He * M_He + x_Ne * M_Ne        # kg/mol

    # SI unit helpers
    mol_per_m3 = si.MOL / (si.METER ** 3)
    J_per_mol = si.JOULE / si.MOL
    J_per_mol_K = si.JOULE / si.MOL / si.KELVIN
    kg_per_m3 = si.KILOGRAM / (si.METER ** 3)

    T = float(state.temperature / si.KELVIN)

    # --- Density (purely from EoS, no ideal gas needed) ---
    rho_molar = float(state.density / mol_per_m3)               # mol/m^3
    rho_mass = float(state.mass_density() / kg_per_m3)           # kg/m^3

    # --- Residual contributions from feos ---

    # Residual Cp: not all feos versions expose Contributions for heat capacity,
    # so we fall back to zero (= pure ideal gas Cp) if it fails.  For noble
    # gases far from their critical points the residual Cp is typically < 2%.
    try:
        cp_res_molar = float(
            state.molar_isobaric_heat_capacity(feos.Contributions.Residual) / J_per_mol_K
        )  # J/(mol*K)
    except Exception:
        cp_res_molar = 0.0

    h_res_molar = float(
        state.molar_enthalpy(feos.Contributions.Residual) / J_per_mol
    )  # J/mol
    s_res_molar = float(
        state.molar_entropy(feos.Contributions.Residual) / J_per_mol_K
    )  # J/(mol*K)

    # --- Ideal gas for monatomic species (exact) ---
    CP_IG_MOLAR = 2.5 * R_GAS   # 20.786 J/(mol*K), same for He and Ne

    # Total molar properties = ideal gas + residual
    cp_molar = CP_IG_MOLAR + cp_res_molar        # J/(mol*K)
    h_molar = CP_IG_MOLAR * T + h_res_molar      # J/mol  (relative to 0 K)

    # For entropy, the ideal gas part includes a pressure-dependent term:
    #   S_ig(T,P) = Cp_ig * ln(T/T_ref) - R * ln(P/P_ref) + S_ref
    # In the PFHX model only entropy *differences* matter, so the reference
    # cancels out.  We use T_ref = 298.15 K, P_ref = 101325 Pa.
    P = float(state.pressure() / si.PASCAL)
    T_REF = 298.15
    P_REF = 101325.0
    s_ig_molar = CP_IG_MOLAR * np.log(T / T_REF) - R_GAS * np.log(P / P_REF)
    s_molar = s_ig_molar + s_res_molar            # J/(mol*K)

    # Convert molar -> mass-specific
    cp_mass = cp_molar / M_mix                     # J/(kg*K)
    h_mass = h_molar / M_mix                       # J/kg
    s_mass = s_molar / M_mix                       # J/(kg*K)

    return {
        'molar_density': rho_molar,
        'mass_density': rho_mass,
        'molecular_mass': M_mix,
        'specific_heat_capacity': cp_mass,
        'enthalpy': h_mass,
        'entropy': s_mass,
        'residual_entropy': s_res_molar,
    }


# ============================================================================
# Dilute gas viscosity: Chapman-Enskog theory
# ============================================================================

def chapman_enskog_viscosity(T, M, sigma, epsilon_k):
    """
    Dilute gas viscosity from Chapman-Enskog theory (Eq. 2, Li et al. 2024).

    Parameters
    ----------
    T : float
        Temperature [K].
    M : float
        Molar mass [g/mol].
    sigma : float
        Lennard-Jones diameter [Angstrom].
    epsilon_k : float
        Lennard-Jones well depth / k_B [K].

    Returns
    -------
    float
        Dilute gas viscosity [Pa*s].
    """
    m_kg = (M / 1000.0) / N_A          # mass per molecule [kg]
    sigma_m = sigma * 1e-10             # diameter [m]
    T_star = T / epsilon_k              # reduced temperature

    # Neufeld collision integral (Eq. 3)
    omega_22 = (1.16145 * T_star ** (-0.14874)
                + 0.52487 * np.exp(-0.77320 * T_star)
                + 2.16178 * np.exp(-2.43787 * T_star))

    eta_0 = (5.0 / 16.0) * np.sqrt(m_kg * K_B * T / np.pi) / (sigma_m ** 2 * omega_22)
    return eta_0


# ============================================================================
# Dilute gas thermal conductivity: polynomial from Li et al. (2024) SI
# ============================================================================

def dilute_gas_thermal_conductivity(T, species):
    """
    Dilute gas thermal conductivity from polynomial fit.

    Parameters
    ----------
    T : float
        Temperature [K].
    species : str
        'He' or 'Ne'.

    Returns
    -------
    float
        Dilute gas thermal conductivity [W/(m*K)].
    """
    n0, n1, n2, n3, n4 = _DILUTE_TC_COEFFS[species]
    return n0 * T**4 + n1 * T**3 + n2 * T**2 + n3 * T + n4


# ============================================================================
# Mixing rules for dilute gas transport
# ============================================================================

def wilke_mixing_viscosity(x, eta_0_pure, m_pure):
    """
    Wilke's mixing rule for dilute gas viscosity.

    Parameters
    ----------
    x : array_like
        Mole fractions [x_He, x_Ne].
    eta_0_pure : array_like
        Pure component dilute gas viscosities [Pa*s].
    m_pure : array_like
        Pure component molecular masses per molecule [kg].

    Returns
    -------
    float
        Mixture dilute gas viscosity [Pa*s].
    """
    N = len(x)
    phi = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            ratio_eta = (eta_0_pure[i] / eta_0_pure[j]) ** 0.5
            ratio_m = (m_pure[j] / m_pure[i]) ** 0.25
            phi[i, j] = (1.0 + ratio_eta * ratio_m) ** 2 / (8.0 * (1.0 + m_pure[i] / m_pure[j])) ** 0.5

    eta_mix = 0.0
    for i in range(N):
        denom = sum(x[j] * phi[i, j] for j in range(N))
        eta_mix += x[i] * eta_0_pure[i] / denom

    return eta_mix


def wassiljewa_mixing_conductivity(x, lambda_0_pure, m_pure):
    """
    Wassiljewa equation (Mason-Saxena form) for dilute gas thermal
    conductivity mixing (Eqs. 11-12, Li et al. 2024).

    Same functional form as Wilke but with conductivities instead of
    viscosities.

    Parameters
    ----------
    x : array_like
        Mole fractions [x_He, x_Ne].
    lambda_0_pure : array_like
        Pure component dilute gas thermal conductivities [W/(m*K)].
    m_pure : array_like
        Pure component molecular masses per molecule [kg].

    Returns
    -------
    float
        Mixture dilute gas thermal conductivity [W/(m*K)].
    """
    N = len(x)
    phi = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            ratio_lam = (lambda_0_pure[i] / lambda_0_pure[j]) ** 0.5
            ratio_m = (m_pure[j] / m_pure[i]) ** 0.25
            phi[i, j] = (1.0 + ratio_lam * ratio_m) ** 2 / (8.0 * (1.0 + m_pure[i] / m_pure[j])) ** 0.5

    lam_mix = 0.0
    for i in range(N):
        denom = sum(x[j] * phi[i, j] for j in range(N))
        lam_mix += x[i] * lambda_0_pure[i] / denom

    return lam_mix


# ============================================================================
# Entropy scaling: residual contributions
# ============================================================================

def _get_nk_coeffs(params):
    """
    Extract the n_k polynomial coefficients from a parameter dictionary.

    For fluid-specific parameters (has_specific=True), uses n1..n4 directly.
    For group parameters (has_specific=False), converts ng1..ng4 via xi.

    Returns
    -------
    ndarray, shape (4,)
    """
    if params['has_specific']:
        return np.array([params['n1'], params['n2'], params['n3'], params['n4']])
    else:
        xi = params['xi']
        return np.array([
            params['ng1'] / xi,
            params['ng2'] / xi ** 1.5,
            params['ng3'] / xi ** 2,
            params['ng4'] / xi ** 2.5,
        ])


def _viscosity_residual_plus(s_res, params):
    """
    Dimensionless residual viscosity from residual entropy (Eqs. 6-7).

    Parameters
    ----------
    s_res : float
        Residual molar entropy [J/(mol*K)] (negative value).
    params : dict
        Li et al. parameter dictionary.

    Returns
    -------
    eta_res_plus : float
        Dimensionless residual viscosity.
    s_plus : float
        Plus-scaled residual entropy (-s_res / R).
    """
    s_plus = -s_res / R_GAS

    if s_plus <= 0:
        return 0.0, s_plus

    n_k = _get_nk_coeffs(params)
    ln_eta_res_plus_1 = (n_k[0] * s_plus
                         + n_k[1] * s_plus ** 1.5
                         + n_k[2] * s_plus ** 2
                         + n_k[3] * s_plus ** 2.5)
    eta_res_plus = np.exp(ln_eta_res_plus_1) - 1.0

    return eta_res_plus, s_plus


def _thermal_conductivity_residual_plus(s_res, params):
    """
    Dimensionless residual thermal conductivity from residual entropy.

    Note: thermal conductivity uses s+/xi scaling, different from viscosity.
    (Li et al. 2024 SI code line 212-213)

    Parameters
    ----------
    s_res : float
        Residual molar entropy [J/(mol*K)].
    params : dict
        Li et al. parameter dictionary.

    Returns
    -------
    lambda_res_plus : float
        Dimensionless residual thermal conductivity.
    s_plus : float
        Plus-scaled residual entropy.
    """
    s_plus = -s_res / R_GAS

    if s_plus <= 0:
        return 0.0, s_plus

    s_scaled = s_plus / params['xi']
    lambda_res_plus = (params['ng1'] * s_scaled
                       + params['ng2'] * s_scaled ** 1.5
                       + params['ng3'] * s_scaled ** 2
                       + params['ng4'] * s_scaled ** 2.5)

    return lambda_res_plus, s_plus


# ============================================================================
# Full mixture transport properties
# ============================================================================

def mixture_viscosity(state, x):
    """
    Calculate He-Ne mixture viscosity using entropy scaling.

    Total viscosity: eta = eta_0_mix + eta_res_mix

    Dilute gas: Chapman-Enskog for each component, mixed via Wilke's rule.
    Residual: Li et al. (2024) entropy scaling with linear n_k mixing.

    Parameters
    ----------
    state : feos.State
        Thermodynamic state of the mixture.
    x : array_like
        Mole fractions [x_He, x_Ne].

    Returns
    -------
    float
        Mixture viscosity [Pa*s].
    """
    T = float(state.temperature / si.KELVIN)
    s_res = float(state.molar_entropy(feos.Contributions.Residual) / (si.JOULE / si.MOL / si.KELVIN))
    rho_N = float(state.density / (si.MOL / si.METER ** 3)) * N_A  # number density [1/m^3]

    # Dilute gas viscosities
    eta_0_He = chapman_enskog_viscosity(T, HELIUM_VISC_PARAMS['M'],
                                        HELIUM_VISC_PARAMS['sigma'],
                                        HELIUM_VISC_PARAMS['epsilon_k'])
    eta_0_Ne = chapman_enskog_viscosity(T, NEON_VISC_PARAMS['M'],
                                        NEON_VISC_PARAMS['sigma'],
                                        NEON_VISC_PARAMS['epsilon_k'])

    m_He = HELIUM_VISC_PARAMS['M'] / 1000.0 / N_A  # kg per molecule
    m_Ne = NEON_VISC_PARAMS['M'] / 1000.0 / N_A

    eta_0_mix = wilke_mixing_viscosity(
        x, np.array([eta_0_He, eta_0_Ne]), np.array([m_He, m_Ne])
    )

    # Residual contribution via entropy scaling
    s_plus = -s_res / R_GAS

    if s_plus <= 0:
        return eta_0_mix

    # Linear mixing of n_k coefficients
    nk_He = _get_nk_coeffs(HELIUM_VISC_PARAMS)
    nk_Ne = _get_nk_coeffs(NEON_VISC_PARAMS)
    nk_mix = x[0] * nk_He + x[1] * nk_Ne

    # Mixture molecular mass per molecule
    m_mix = x[0] * m_He + x[1] * m_Ne

    ln_eta_res_plus_1 = (nk_mix[0] * s_plus
                         + nk_mix[1] * s_plus ** 1.5
                         + nk_mix[2] * s_plus ** 2
                         + nk_mix[3] * s_plus ** 2.5)
    eta_res_plus = np.exp(ln_eta_res_plus_1) - 1.0

    # Dimensional residual viscosity (Eq. 4)
    eta_res = eta_res_plus * rho_N ** (2.0 / 3.0) * np.sqrt(m_mix * K_B * T) / s_plus ** (2.0 / 3.0)

    return eta_0_mix + eta_res


def mixture_thermal_conductivity(state, x):
    """
    Calculate He-Ne mixture thermal conductivity using entropy scaling.

    Total conductivity: lambda = lambda_0_mix + lambda_res_mix
    (Critical enhancement omitted -- conditions are far from critical.)

    Dilute gas: polynomial for each component, mixed via Wassiljewa.
    Residual: Li et al. (2024) entropy scaling with linear ng_k mixing.

    Parameters
    ----------
    state : feos.State
        Thermodynamic state of the mixture.
    x : array_like
        Mole fractions [x_He, x_Ne].

    Returns
    -------
    float
        Mixture thermal conductivity [W/(m*K)].
    """
    T = float(state.temperature / si.KELVIN)
    s_res = float(state.molar_entropy(feos.Contributions.Residual) / (si.JOULE / si.MOL / si.KELVIN))
    rho_N = float(state.density / (si.MOL / si.METER ** 3)) * N_A

    # Dilute gas thermal conductivities
    lam_0_He = dilute_gas_thermal_conductivity(T, 'He')
    lam_0_Ne = dilute_gas_thermal_conductivity(T, 'Ne')

    m_He = HELIUM_VISC_PARAMS['M'] / 1000.0 / N_A
    m_Ne = NEON_VISC_PARAMS['M'] / 1000.0 / N_A

    lam_0_mix = wassiljewa_mixing_conductivity(
        x, np.array([lam_0_He, lam_0_Ne]), np.array([m_He, m_Ne])
    )

    # Residual contribution
    s_plus = -s_res / R_GAS

    if s_plus <= 0:
        return lam_0_mix

    # Linear mole-fraction mixing of ng parameters, then apply xi scaling
    # For thermal conductivity, each component's contribution is:
    #   ng_k_i / xi_i  applied to s+  (different from viscosity!)
    # But the mixture approach (Eq. 14 of Li et al.) uses linear mixing of
    # the already-scaled coefficients:
    params_list = [HELIUM_THERMAL_PARAMS, NEON_THERMAL_PARAMS]
    ngk_mix = np.zeros(4)
    for i in range(2):
        xi_i = params_list[i]['xi']
        ngk_i = np.array([
            params_list[i]['ng1'] / xi_i,
            params_list[i]['ng2'] / xi_i ** 1.5,
            params_list[i]['ng3'] / xi_i ** 2,
            params_list[i]['ng4'] / xi_i ** 2.5,
        ])
        ngk_mix += x[i] * ngk_i

    lambda_res_plus = (ngk_mix[0] * s_plus
                       + ngk_mix[1] * s_plus ** 1.5
                       + ngk_mix[2] * s_plus ** 2
                       + ngk_mix[3] * s_plus ** 2.5)

    # Mixture molecular mass for thermal conductivity (Eq. 13):
    #   m_mix = (sum_i y_i * sqrt(m_i))^2   where y_i are mass fractions
    M_arr = np.array([HELIUM_VISC_PARAMS['M'], NEON_VISC_PARAMS['M']])  # g/mol
    M_avg = np.sum(x * M_arr)
    y = x * M_arr / M_avg  # mass fractions
    m_pure = np.array([m_He, m_Ne])
    m_mix = np.sum(y * np.sqrt(m_pure)) ** 2

    # Dimensional residual thermal conductivity (Eq. 3, adapted for lambda)
    lambda_res = (lambda_res_plus / s_plus ** (2.0 / 3.0)) * (
        rho_N ** (2.0 / 3.0) * K_B * np.sqrt(K_B * T / m_mix)
    )

    # Clamp: for pure He at low pressure, residual can go slightly negative
    # due to Group 1 parameter limitations
    if lambda_res < 0 and x[0] > 0.99:
        lambda_res = 0.0

    return lam_0_mix + lambda_res
