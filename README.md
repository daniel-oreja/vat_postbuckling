# Postbuckling Analysis and Design Space Exploration of Tow-Steered VAT Composite Plates

This repository contains the MATLAB and OptiStruct scripts developed for my Bachelor's Thesis (TFG) on the buckling and postbuckling behaviour of Variable Angle Tow (VAT) composite plates with cubic Bézier fibre paths, manufactured through Automated Fibre Placement (AFP).

**Author:** Daniel Oreja González
**Institution:** Universidad Carlos III de Madrid
**Supervisor:** Alberto Racionero Sánchez
**Year:** 2026

## Overview

The repository contains two complementary pieces of code:

1. **Single-case validation code** (`01_validation/`) — runs *one* design at a time (a single laminate / fibre-path definition) through a modal, linear-buckling, or non-linear postbuckling analysis in Altair OptiStruct. It was used to validate the finite element model against literature benchmarks.

2. **Design-space exploration code** (`02_analysis/`) — the main deliverable of the thesis. A Latin Hypercube Sample (LHS) of 160 VAT laminate designs is generated and each design is run through the same buckling + postbuckling.

> **Reproducing the thesis figures and statistics does NOT require running the OptiStruct simulations again.** The two processed result tables are already included in `02_doe_and_postprocessing/` — see "How to reproduce the results" below.

## Repository structure

```
├── 01_validation/
│   └── validation_case.m                       # Run ONE design through OptiStruct (modal / buckling / postbuckling)
│
├── 02_analysis/
│   ├── main_lhs.m                           # [OPTIONAL] Generates the 160-sample LHS and runs OptiStruct for all of them
│   ├── analyze_results.m                    # MAIN SCRIPT: statistics + all figures
│   ├── pareto.m                             # MAIN SCRIPT: Pareto dominance analysis
│   ├── results_buckling.csv                 # Included: Kcr, eigenvalues and design variables for all 160 samples
│   └── results_postbuckling.csv             # Included: postbuckling curve points for all samples
│
└── README.md
```

## Requirements

- MATLAB (R2021a or later recommended) 
- Statistics and Machine Learning Toolbox 
- Altair OptiStruct (student or full license) — **only** required for `01_validation/validation_case.m` and `02_analysis/main_lhs.m`. **Not** required to reproduce the thesis figures and statistics from the included CSVs.

## How to reproduce the results

### Option A — Reproduce the thesis figures and statistics (recommended)

No OptiStruct license, no simulations, no waiting. Just:

1. Open MATLAB with `02_analysis/` as the working folder.
2. Run **`analyze_results.m`**. It reads `results_buckling.csv` and `results_postbuckling.csv` and reproduce the thesis figures and statistics.
3. Run **`pareto_analysis.m`**. It reads the same two CSVs and produces the pareto analysis statistics.

### Option B — Reproduce the full analysis from scratch (requires Altair OptiStruct)

1. **Generate and run the full LHS design space**: run `02_analysis/main_lhs.m`. This generates the 160-sample Latin Hypercube design, builds and runs the corresponding OptiStruct buckling + postbuckling jobs for every sample, and writes `results_buckling.csv` and `results_postbuckling.csv`. **This step requires an OptiStruct license and can take a long time** (hours, depending on the machine and the number of samples).
2. **Analyze the results**: run `analyze_results.m` and `pareto_analysis.m` as described in Option A.

## Data availability

The raw per-sample simulation files (`.fem`, `.out`, `.h3d`, `.spcf` — one set per sample) are **not** included in this repository due to their size (160 samples × several files each). The two processed result tables that all of the analysis is based on are included:

- `results_buckling.csv`
- `results_postbuckling.csv`

## Citation

If you use this code or refer to this work, please cite:

> Daniel Oreja González. "Post-buckling Analysis of Tow-Steered Composites Using Finite Elements Method." Bachelor's Thesis, Universidad Carlos III de Madrid, 2026.

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

