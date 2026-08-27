# Postbuckling Analysis and Design Space Exploration of Tow-Steered VAT Composite Plates

This repository contains the MATLAB and OptiStruct scripts developed for my Bachelor's Thesis (TFG) on the buckling and postbuckling behaviour of Variable Angle Tow (VAT) composite plates with cubic Bézier fibre paths, manufactured through Automated Fibre Placement (AFP).

**Author:** Daniel Oreja González
**Institution:** Universidad Carlos III de Madrid
**Supervisor:** Alberto Racionero Sánchez
**Year:** 2026

## Overview

The repository contains two complementary pieces of code:

1. **Single-case validation code** (`01_single_case_validation/`) — runs *one* design at a time (a single laminate / fibre-path definition) through a modal, linear-buckling, or geometrically non-linear postbuckling analysis in Altair OptiStruct. It was used to validate the finite element model against literature benchmarks (the quasi-isotropic buckling load, and the postbuckling reference curve of Wu et al., 2013) before launching the full design-space sweep.

2. **Design-space exploration code** (`02_doe_and_postprocessing/`) — the main deliverable of the thesis. A Latin Hypercube Sample (LHS) of 160 VAT laminate designs (cubic Bézier fibre paths, 2 unique curves per laminate) is generated and each design is run through the same buckling + postbuckling pipeline. The results are then post-processed to:
   - Compute the buckling load factor ($K_{cr}$) and the postbuckling strain reserve ($\varepsilon_{2N}$) for every design.
   - Run a sensitivity analysis (Pearson/Spearman correlations) of both objectives against the 16 Bézier control-point coordinates.
   - Assess the manufacturability of each design (minimum fibre steering radius) against three representative AFP process scenarios (Standard AFP, Narrow-tow AFP, Continuous Tow Shearing).
   - Run a Pareto dominance analysis between the two objectives, both unconstrained and under each manufacturability scenario.

> **Reproducing the thesis figures and statistics does NOT require running the OptiStruct simulations again.** The two processed result tables are already included in `02_doe_and_postprocessing/` — see "How to reproduce the results" below.

## Repository structure

```
├── 01_single_case_validation/
│   └── run_single_case.m                       # Run ONE design through OptiStruct (modal / buckling / postbuckling)
│
├── 02_doe_and_postprocessing/
│   ├── run_lhs_doe.m                           # [OPTIONAL] Generates the 160-sample LHS and runs OptiStruct for all of them
│   ├── extract_w_from_nlout.m                  # [OPTIONAL] Re-extracts out-of-plane deflection (w) from raw solver output
│   ├── analyze_results.m                       # MAIN SCRIPT: statistics + all figures (reads the CSVs only)
│   ├── pareto_analysis.m                       # MAIN SCRIPT: Pareto dominance analysis + fibre-path maps (reads the CSVs only)
│   ├── results_buckling.csv                    # Included: Kcr, eigenvalues and design variables for all 160 samples
│   └── results_postbuckling_curves_with_w.csv  # Included: postbuckling curve points (eps, Nx, w) for all samples
│
└── README.md
```

## Requirements

- MATLAB (R2021a or later recommended; the scripts use `readtable`, `contains`, `corr`, `polyfit`, `sgtitle`)
- Statistics and Machine Learning Toolbox (used by `analyze_results.m` for the Pearson/Spearman `corr` sensitivity analysis)
- Altair OptiStruct (student or full license) — **only** required for `01_single_case_validation/run_single_case.m` and `02_doe_and_postprocessing/run_lhs_doe.m`. **Not** required to reproduce the thesis figures/statistics from the included CSVs.

## How to reproduce the results

### Option A — Reproduce the thesis figures and statistics (recommended; this is what the evaluation committee should do)

No OptiStruct license, no simulations, no waiting. Just:

1. Open MATLAB with `02_doe_and_postprocessing/` as the working folder.
2. Run **`analyze_results.m`**. It reads `results_buckling.csv` and `results_postbuckling_curves_with_w.csv` directly and produces:
   - Convergence summary and $K_{cr}$ / $\varepsilon_{2N}$ statistics (printed to the command window).
   - Sensitivity analysis: Pearson/Spearman correlations of $K_{cr}$ and $\varepsilon_{2N}$ with the 16 design variables.
   - Manufacturability check under the three AFP scenarios.
   - All 13 summary figures: w-extraction verification, dataset-wide overviews, $K_{cr}$ histogram, postbuckling curves, sensitivity bar charts, manufacturability scatter plots, best/worst postbuckling design comparison, and fibre-path maps.
3. Run **`pareto_analysis.m`**. It reads the same two CSVs and produces:
   - The unconstrained Pareto front between $K_{cr}$ and $\varepsilon_{2N}$, with the best-compromise ("knee point") design highlighted.
   - The manufacturability-constrained Pareto fronts under the three AFP scenarios.
   - Fibre-path maps of the knee-point design for the unconstrained case and for each AFP scenario.

### Option B — Reproduce the full pipeline from scratch (requires Altair OptiStruct)

1. **Validate the FE model** (optional, but this is how the model was checked before running the full sweep): in `01_single_case_validation/run_single_case.m`, set `type_variation`, `type_analysis`, `transverse_edges`, and the other options at the top of the script to the case you want to check, adjust the `altair_path` variable to your local OptiStruct installation, and run it.
2. **Generate and run the full LHS design space**: run `02_doe_and_postprocessing/run_lhs_doe.m`. This generates the 160-sample Latin Hypercube design, builds and runs the corresponding OptiStruct buckling + postbuckling jobs for every sample, and writes `results_buckling.csv` and `results_postbuckling_curves.csv`. **This step requires an OptiStruct license and can take a long time** (hours, depending on the machine and the number of samples).
3. **(Optional) Re-extract the out-of-plane deflection**: run `02_doe_and_postprocessing/extract_w_from_nlout.m` if you want to regenerate `results_postbuckling_curves_with_w.csv` from the raw per-sample solver output files (`sample_XXX_pb_nl.out`). This is only needed if the raw solver files are regenerated; the CSV is already included in this repository.
4. **Analyze the results**: run `analyze_results.m` and `pareto_analysis.m` as described in Option A.

## Data availability

The raw per-sample simulation files (`.fem`, `.out`, `.h3d`, `.spcf` — one set per DOE sample) are **not** included in this repository due to their size (160 samples × several files each). The two processed result tables that all of the analysis is based on are included:

- `results_buckling.csv`
- `results_postbuckling_curves_with_w.csv`

## Citation

If you use this code or refer to this work, please cite:

> Daniel Oreja González. "Post-buckling Analysis of Tow-Steered Composites Using Finite Elements Method." Bachelor's Thesis, Universidad Carlos III de Madrid, 2026.

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

