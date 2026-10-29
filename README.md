# Faster, cheaper or more accurate? Valuing smartphone phenotyping in East African bean breeding

Code and data for the article

> Tang, J.-S., Wu, M., Liang, W., Gong, J., Chen, J. Faster, cheaper or more accurate? Valuing smartphone phenotyping in East African bean breeding with a simulation calibrated to public trial data. *Agricultural Systems* (submitted).

The model values a generic smartphone computer-vision phenotyping tool in a four-stage East African common-bean breeding pipeline. A tool can act through three channels: less measurement error, cheaper plots under a fixed budget, and a shorter breeding cycle. The study has three parts.

- **Genetic parameters from public trials** (`empirical/`): within-environment plot heritability and the ratio of genotype-by-environment to genetic variance, estimated with mixed models from public East African multi-environment bean trials, and the error reduction from adjusting plot yield for stand.
- **Layer 1** (`layer1_genetic_sim/`): an AlphaSimR stochastic simulation that compares manual and tool phenotyping under one budget across a factorial grid of tool profiles (measurement-error reduction, per-plot cost reduction, cycle compression, fixed cost, budget). Each profile is run with 100 paired replicates under common random numbers. A comparison with the multi-stage breeder's equation is included.
- **Layer 2** (`layer2_economics/`): converts simulated genetic gain into logistic adoption and 25-year net present value, then runs global sensitivity analysis (exact Sobol and Shapley indices on the grid, PAWN, active subspaces, horizon-dependent Sobol indices), value-of-information analysis (EVPPI) and robustness scenarios.

Three meta-analyses (`paper/lit/`) set the ranges of the tool inputs and the calibration target: realized genetic gain in common bean (MA1), agreement between image-based tools and manual references (MA2), and time and cost savings of digital phenotyping (MA3).

Every number in the manuscript and the supplementary material is a LaTeX macro generated from the outputs (`paper/numbers.tex`; `paper/numbers.csv` lists the source file of each number).

## Versions

| Release | Article version | Genetic parameters | Layer 1 / Layer 2 outputs |
|---|---|---|---|
| 2.0.0 | *Agricultural Systems* submission | estimated from public trials: plot h² 0.10, GxE-to-genetic variance ratio 1.5 | `sweep_v3*`, `outputs_v3` |
| 1.0.0 | earlier submission | assumed: plot h² 0.30, between-environment GxE parameter 0.60 | `sweep_v2*`, `outputs_v2` |

The 1.0.0 outputs are kept. The v2 sweep configurations pin their genetic values, so both versions can be rerun from this release. `layer2_economics/R/00_paths.R` explains how to switch between them with environment variables.

## Repository layout

| Path | Contents |
|---|---|
| `config/` | Genetic, breeding and economic parameters; sweep grids (`config/sweeps/`) |
| `empirical/R/`, `empirical/outputs/` | Cleaning and variance-component analyses of the public trial data; estimates, summaries and cleaning logs |
| `layer1_genetic_sim/R/` | Founder population, breeding pipeline, budget allocation (with an optional fixed cost per environment), sweep driver (`07_param_sweep.R`), breeder's-equation comparison (`08_breeders_equation_check.R`) |
| `layer1_genetic_sim/outputs/sweep_v3*/` | Layer 1 results: per-replicate files (`cell_*.csv.gz`), per-cell summaries and `cell_summaries.csv` for the main grid, two extension blocks, the heritability × GxE block, the environment-cost block and the stress block |
| `layer2_economics/R/` | Economic model, sensitivity and decision analyses |
| `layer2_economics/outputs_v3/` | Layer 2 results (lookup table, sensitivity indices, EVPPI, scenarios) |
| `paper/lit/` | Meta-analysis script, extraction sheets with page or table locators, search logs and PRISMA counts |
| `paper/R/`, `paper/py/` | Number and table export (R); figures, graphical abstract and supplementary figures (Python) |
| `paper/` | Manuscript source (Elsevier `elsarticle`), supplementary material, figures and tables |
| `slurm/` | Scripts used to run Layer 1 as Slurm array jobs |

## Public trial data

The raw trial data are not redistributed here. Download them from their repositories, place them under `data/external/bean_public_2026-10-10/datasets/` (the paths are set at the top of each script in `empirical/R/`), and run the scripts in numeric order:

- TSEZVG: Harvard Dataverse, https://doi.org/10.7910/DVN/TSEZVG (CC BY 4.0; guestbook)
- ZUDPEP: Harvard Dataverse, https://doi.org/10.7910/DVN/ZUDPEP (CC0)
- Dryad dz08kpsbz: https://doi.org/10.5061/dryad.dz08kpsbz

The derived estimates that the simulation uses are included in `empirical/outputs/`, so the steps below do not need the raw data.

## Requirements

- R 4.3 or later with `AlphaSimR` (2.1.0 was used for Layer 1), `data.table`, `yaml`, `arrow`, `mvtnorm`, `metafor` (5.2.1), `lme4` and `glmmTMB`. `layer1_genetic_sim/R/00_install_deps.R` installs the simulation packages.
- Python 3.10 or later with `matplotlib`, `numpy`, `pandas`, `pyyaml` and `pillow`.
- A LaTeX distribution with `latexmk` to compile the manuscript and the supplementary material.

## Reproducing the results

Run all commands from the repository root.

```bash
make numbers      # regenerate paper/numbers.tex and numbers.csv from the stored outputs
make analyses     # rebuild out-of-date Layer 2 analyses (slow)
make figures      # figures and graphical abstract
make supplement   # supplementary figures, tables and PDF
make paper        # numbers, supplementary material and the manuscript PDF
Rscript layer1_genetic_sim/R/08_breeders_equation_check.R   # breeder's equation vs simulation
Rscript paper/lit/meta_analysis.R                            # re-run the three meta-analyses
```

The Layer 1 results are included, so the commands above do not need a cluster. To re-run Layer 1, simulate one cell with

```bash
TOOL_SWEEP_GRID=config/sweeps/sweep_v3_main.yaml TOOL_SWEEP_SUBDIR=sweep_v3 \
  Rscript layer1_genetic_sim/R/07_param_sweep.R <cell_id> 100
```

or submit a whole block with the scripts in `slurm/` (set the allocation for your cluster), then rebuild its cell summaries with

```bash
TOOL_SWEEP_GRID=config/sweeps/sweep_v3_main.yaml TOOL_SWEEP_SUBDIR=sweep_v3 \
  Rscript -e 'source("layer1_genetic_sim/R/07_param_sweep.R"); aggregate_sweep()'
```

## Citation

See `CITATION.cff`. Please cite the article and the archived version of this repository (https://doi.org/10.5281/zenodo.23197205, which resolves to the latest release).
