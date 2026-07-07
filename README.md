# Large-Scale Hydrogen Liquefaction: Cryogenic Cooling and Boil-off Gas Recovery

Process Modeling and Techno-Economic Analysis of Catalytic Plate-Fin Heat
Exchangers and Ejector-Driven Cycles

MSc Thesis, Mechanical Engineering (Energy, Flow and Process Technology)
Faculty of Mechanical Engineering, Department of Process & Energy
Delft University of Technology

Author: Borjan Trajanoski (6304168)

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

```text
.
├── README.md
├── .gitignore
│
├── data/            # Storage-tank BOG cases, facility/literature validation datasets (CSV)
├── flowsheet/       # Aspen HYSYS process flowsheet
│
├── src/
│   ├── python/
│   │   ├── properties/    # PFHX property models (O'Neill + SAFT-VRQ-Mie) and He-Ne VLE/transport notebooks
│   │   └── tea/            # Techno-economic assessment
│   ├── matlab/
│   │   └── bog/            # Two-vessel holding/loading BOG model
│   └── pinch_analysis/     # PFHX-5 pinch/area pipeline and the coupled ejector BOG-recovery model
│
└── results/
    ├── figures/     # Publication figures
    └── pfd/         # Process flow diagrams
```

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
| REFPROP               | v10             | Reference fluid properties                     |
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
