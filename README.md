# Large-Scale Hydrogen Liquefaction: Cryogenic Cooling and Boil-off Gas Recovery

Process Modeling and Techno-Economic Analysis of Catalytic Plate-Fin Heat
Exchangers and Ejector-Driven Cycles

MSc Thesis, Mechanical Engineering (Energy, Flow and Process Technology)
Faculty of Mechanical Engineering, Department of Process & Energy
Delft University of Technology

Author: Borjan Trajanoski

---

## Project Overview

This repository contains the source code, simulation models, data, and results
for an MSc thesis on large-scale hydrogen liquefaction. The work takes an 86
tonnes-per-day mixed-refrigerant pre-cooled Joule-Brayton reference process and
quantifies two idealizations common in conceptual large-scale designs.

The first is the property and kinetic modeling of the cryogenic catalytic
plate-fin heat exchanger (PFHX). The helium-neon refrigerant is modeled with
SAFT-VRQ-Mie and residual entropy scaling in place of the dilute-gas
correlations used in earlier studies, which changes the predicted thermal
conductivity by a factor of two near 30 K and drives the required heat
exchanger length toward the single-unit manufacturing limit.

The second is the assumption of full liquid yield, which ignores the boil-off
gas (BOG) generated during LH2 storage and truck loading. A two-vessel storage
model quantifies holding and loading BOG across a range of tank sizes,
geometries, and orientations, and an ejector-driven recovery cycle is proposed
and assessed for that BOG.

A techno-economic analysis of the isolated cryogenic system then compares the
baseline and adapted configurations.

The codebase couples several tools: Aspen HYSYS for the process flowsheet,
REFPROP and FeOs for fluid properties, MATLAB for the ejector and two-vessel
storage models, and Python for data processing, TEA calculations, VLE/transport
property fitting, and all publication figures.

---

## Repository Structure

````text
.
├── README.md
├── .gitignore
│
├── data/                        # Simulation exports and validation datasets (flat, CSV)
│   ├── 1_3200m3_Spherical.csv, 2_3200m3_horizontal.csv, 2a_3200m3_vertical.csv,
│   │   3_1500m3_spherical.csv, 4_1500m3_horizontal.csv, 4a_1500m3_vertical.csv,
│   │   5_500m3_spherical.csv, 6_500m3_horizontal.csv, 6a_500m3_vertical.csv,
│   │   7_100m3_spherical.csv, 8_100m3_horizontal.csv, 8a_100m3_vertical.csv
│   │                            # Storage-tank BOG cases by volume/geometry/orientation
│   ├── 3200m3_NASA_KSC_Case.csv, 3200m3_NASA_KSC_Case_v2.csv
│   │                            # NASA KSC reference tank case
│   ├── Boil-off_Simulation_Results.csv
│   ├── Case_1_self_pressurization.csv, Case_2_DDL_Dewar_100days.csv
│   │                            # Self-pressurization and long-duration dewar validation
│   ├── ksite_selfpress_3.5Wm2.csv, ksite_temp_139cm_interface_3.5Wm2.csv,
│   │   ksite_temp_162cm_vapor_3.5Wm2.csv
│   │                            # K-Site cryogenic test facility validation data
│   ├── llnl_boiloff_summer.csv  # LLNL boil-off validation data
│   ├── assael_1981_he_ne_thermal_conductivity.csv, helium_thermal_conductivity_degroot1978.csv,
│   │   neon_thermal_conductivity_degroot1978.csv
│   │                            # Literature transport-property reference data
│   └── validation_*assael1981*.csv
│                                # SAFT/REFPROP transport-property validation vs. Assael (1981)
│
├── models/
│   ├── LH2_aspen_model.hsc      # Aspen HYSYS process flowsheet
│   └── LH2_aspen_model.bk0      # HYSYS backup
│
├── src/
│   ├── python/
│   │   ├── properties/
│   │   │   ├── hydrogen-pfhx-saft/       # O'Neill PFHX model + SAFT-VRQ-Mie EOS (standalone package)
│   │   │   └── hene_transport_properties/
│   │   │       ├── HeNe_VLE_PR.ipynb     # He-Ne VLE via Peng-Robinson, vs. Heck et al. data
│   │   │       ├── HeNe_VLE_SAFT.ipynb   # He-Ne VLE via SAFT-VRQ-Mie (feos), vs. Heck et al. data
│   │   │       └── HeNe_transport.ipynb  # He-Ne transport property fitting
│   │   └── tea/                          # Techno-economic assessment
│   │
│   ├── matlab/
│   │   └── storage/              # Two-vessel holding/loading BOG model, steady-state loading
│   │
│   └── pinch_analysis/           # PFHX-5 pinch + area pipeline and the ejector BOG-recovery model
│       ├── ejector/              # Ejector solver (adapted from Moro et al.), MATLAB
│       └── pinch/                # Pinch/area calculations (MATLAB) and figure scripts (Python)
│
└── results/
    ├── figures/                 # Publication figures: validation plots, Ts diagrams, TEA/sensitivity
    │                            # figures, pinch curves, transport-property fits, etc.
    └── pfd/                     # Process flow diagrams (baseline, adapted, isolated configurations)
````

The `ejector/` and `pinch/` subfolders under `pinch_analysis/` are coupled by a
fixed workflow (run the ejector step, then `cd ../pinch` and run the pinch
step), so they are kept as siblings rather than split across the
`python/`/`matlab/` boundary.

The thesis is organized by process stage. Chapter 4 covers the cryogenic PFHX
(Stage 3), Chapter 5 covers LH2 storage and BOG management (Stage 5), and
Chapter 6 covers the techno-economic analysis. The `src/` layout mirrors this
split.

`.gitignore` excludes Python build artifacts (`__pycache__/`, `*.pyc`) and
`.mat` files (large, regenerable MATLAB binaries). Everything else — figures,
PDFs, HYSYS case files, and processed data — is tracked directly.

---

## Getting Started / Software Prerequisites

The following software and versions were used. Later versions are likely
compatible but were not tested. Confirm the version rows against your actual
setup before publishing; some are carried over from working notes rather than
read off the manuscript.

| Tool                  | Version         | Purpose                                        |
|-----------------------|-----------------|------------------------------------------------|
| Python                | 3.10 or later   | Data processing, TEA, figures, VLE notebooks   |
| Jupyter               | latest          | He-Ne VLE / transport-property notebooks       |
| MATLAB                | R2023b or later | Ejector and two-vessel storage models          |
| Aspen HYSYS           | V14             | Process flowsheet and simulation               |
| REFPROP               | v9.1 and v10    | Reference fluid properties                     |
| FeOs                  | latest          | SAFT-VRQ-Mie property calculations             |
| BoilFAST              | latest          | Boil-off gas cross-verification                |

### Python environment

```bash
git clone git@github.com:borjan-trajanoski/lh2-cryo-process.git
cd lh2-cryo-process
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
```

Populate `requirements.txt` with the packages actually imported (numpy, scipy,
pandas, matplotlib, feos, jupyter, and the REFPROP Python bindings).

### Licensed and proprietary tools

Aspen HYSYS, REFPROP, MATLAB, and BoilFAST require separate licenses and are not
distributed in this repository. Case files (`.hsc`) require a matching HYSYS
installation to open. REFPROP calls assume the fluid files ship with your local
install and are not committed here.
