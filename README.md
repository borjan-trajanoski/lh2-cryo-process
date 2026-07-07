# Large-Scale Hydrogen Liquefaction: Cryogenic Cooling and Boil-off Gas Recovery

Process Modeling and Techno-Economic Analysis of Catalytic Plate-Fin Heat
Exchangers and Ejector-Driven Cycles

MSc Thesis, Mechanical Engineering (Energy, Flow and Process Technology)
Faculty of Mechanical Engineering, Department of Process & Energy
Delft University of Technology

Author: Borjan Trajanoski (6304168)
Supervisors: Dr.ir. Mahinder Ramdin, Dr. Chiara Falsetti
Committee: Dr. Emanuele Zanetti, Prof.dr. Arvind Gangoli Rao
Project duration: November 17, 2025 to July 14, 2026

---

## Project Overview

This repository contains the source code, simulation models, data, and written
documentation for an MSc thesis on large-scale hydrogen liquefaction. The work
takes an 86 tonnes-per-day mixed-refrigerant pre-cooled Joule-Brayton reference
process and quantifies two idealizations common in conceptual large-scale
designs.

The first is the property and kinetic modeling of the cryogenic catalytic
plate-fin heat exchanger (PFHX). The helium-neon refrigerant is modeled with
SAFT-VRQ-Mie and residual entropy scaling in place of the dilute-gas
correlations used in earlier studies, which changes the predicted thermal
conductivity by a factor of two near 30 K and drives the required heat
exchanger length toward the single-unit manufacturing limit.

The second is the assumption of full liquid yield, which ignores the boil-off
gas (BOG) generated during LH2 storage and truck loading. A two-vessel storage
model quantifies holding and loading BOG, and an ejector-driven recovery cycle
is proposed and assessed for that BOG.

A techno-economic analysis of the isolated cryogenic system then compares the
baseline and adapted configurations.

The codebase couples several tools: Aspen HYSYS for the process flowsheet,
REFPROP and FeOs for fluid properties, MATLAB for the ejector and two-vessel
storage models, and Python for data processing, TEA calculations, and all
publication figures.

---

## Repository Structure

```text
.
├── README.md
├── .gitignore
├── LICENSE
│
├── data/
│   ├── raw/                 # Unprocessed simulation exports, property tables
│   ├── processed/           # Cleaned, analysis-ready datasets
│   └── external/            # Third-party reference data (Petitpas, Al Ghafri, etc.)
│
├── src/
│   ├── python/
│   │   ├── properties/
│   │   │   └── hydrogen-pfhx-saft/  # O'Neill PFHX model + SAFT-VRQ-Mie EOS (standalone package)
│   │   └── tea/              # Techno-economic assessment (cryo_textbook_tea.py, economic_figures.py)
│   │
│   ├── matlab/
│   │   └── storage/          # Two-vessel holding/loading BOG model (Petitpas-based), steady-state loading
│   │
│   └── pinch_analysis/       # PFHX-5 pinch + area pipeline and the ejector BOG-recovery model
│       ├── ejector/          # Ejector solver (adapted from Moro et al.), MATLAB
│       └── pinch/            # Pinch/area calculations (MATLAB) and figure scripts (Python)
│
├── models/
│   └── hysys/               # Aspen HYSYS case files (.hsc) and flowsheet notes
│
├── results/
│   ├── figures/             # Final figures (versioned; large exports ignored)
```

The `ejector/` and `pinch/` subfolders under `pinch_analysis/` are coupled by a
fixed workflow (run the ejector step, then `cd ../pinch` and run the pinch
step), so they are kept as siblings rather than split across the
`python/`/`matlab/` boundary.

The thesis is organized by process stage. Chapter 4 covers the cryogenic PFHX
(Stage 3), Chapter 5 covers LH2 storage and BOG management (Stage 5), and
Chapter 6 covers the techno-economic analysis. The `src/` layout mirrors this
split.

Large binary outputs (raw HYSYS exports, `.mat` files, high-resolution figure
exports, compiled PDFs) are intentionally excluded from version control via
`.gitignore`. Only source files and analysis-ready processed data are tracked.

---

## Getting Started / Software Prerequisites

The following software and versions were used. Later versions are likely
compatible but were not tested. Confirm the version rows against your actual
setup before publishing; some are carried over from working notes rather than
read off the manuscript.

| Tool                  | Version         | Purpose                                        |
|-----------------------|-----------------|------------------------------------------------|
| Python                | 3.10 or later   | Data processing, TEA, figures                  |
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
pandas, matplotlib, feos, and the REFPROP Python bindings).

### Licensed and proprietary tools

Aspen HYSYS, REFPROP, MATLAB, and BoilFAST require separate licenses and are not
distributed in this repository. Case files (`.hsc`) require a matching HYSYS
installation to open. REFPROP calls assume the fluid files ship with your local
install and are not committed here.
