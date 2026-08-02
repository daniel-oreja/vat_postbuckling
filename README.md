# Postbuckling Analysis and Design Space Exploration of Tow-Steered VAT Composite Plates

This repository contains the MATLAB and OptiStruct scripts developed for my Bachelor's/Master's thesis (TFG) on the buckling and postbuckling behaviour of Variable Angle Tow (VAT) composite plates with cubic Bézier fibre paths, manufactured through Automated Fibre Placement (AFP).

**Author:** Daniel Oreja González
**Institution:** Universidad Carlos III de Madrid
**Supervisor:** Alberto Racionero Sánchez
**Year:** 2026

## Overview

The code implements a fully automated pipeline that:
1. Generates tow-steered laminate designs defined by cubic Bézier fibre paths, sampled via Latin Hypercube Sampling (LHS).
2. Maps the fibre orientations onto a finite element mesh (HyperMesh/OptiStruct format).
3. Runs linear buckling and geometrically non-linear postbuckling analyses in Altair OptiStruct.
4. Post-processes the results to compute the buckling load factor ($K_{cr}$) and the postbuckling stiffness metric ($\varepsilon_{2N}$).
5. Performs a variance-based sensitivity analysis (Polynomial Chaos Expansion / PC-Kriging, via UQLab) and a Pareto dominance analysis on the resulting dataset.
6. Assesses the manufacturability of each design against representative AFP steering-radius constraints.

## Repository structure

## Requirements

- MATLAB R20XX or later
- [UQLab](https://www.uqlab.com/) v2.2 (for PCE/PC-Kriging and Sobol sensitivity analysis)
- Altair OptiStruct (student or full license) for the finite element analyses
- Global Optimization Toolbox (for `gamultiobj`, used in the NSGA-II surrogate-based optimisation)

## How to reproduce the results

1. Run `lhs_pipeline/generate_lhs_designs.m` to generate the 160-sample LHS design matrix and the corresponding `.fem` input files.
2. Run `fem_automation/run_all_samples.m` to launch the OptiStruct analyses (linear buckling + NLSTAT postbuckling) for all samples. **Note:** this step requires an OptiStruct license and can take a long time depending on your machine.
3. Run `fem_automation/parse_results.m` to collect the results into `results_buckling.csv` and `results_postbuckling_curves.csv`.
4. Run `pce_sobol/build_surrogates.m` to fit the PCE and PC-Kriging models and compute the Sobol indices.
5. Run `pareto/pareto_analysis.m` for the Pareto dominance analysis, and `pareto/nsga2_surrogate.m` for the surrogate-based multi-objective optimisation.
6. Run `manufacturability/manufacturability_assessment.m` to compute the minimum steering radius for each design and evaluate the three AFP scenarios.

## Data availability

The raw simulation output files (`.fem`, `.out`, `.h3d`) are not included in this repository due to their size. The processed datasets (`results_buckling.csv`, `results_postbuckling_curves.csv`) are included.

## Citation

If you use this code or refer to this work, please cite:

> [Tu nombre]. "[Título completo de tu TFG]." Bachelor's/Master's Thesis, [Universidad], 2026.

## License

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
