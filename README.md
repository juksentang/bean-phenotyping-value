# Cycle Time versus Accuracy: the simulated economic value of AI-based phenotyping in smallholder bean breeding

Code and data for the article

> Tang, J.-S.; Chen, J. Cycle Time versus Accuracy: The Simulated Economic Value of AI-Based Phenotyping in Smallholder Bean Breeding. *Agronomy* (submitted).

The model values a generic smartphone computer-vision phenotyping tool in a four-stage East African common-bean breeding pipeline. It couples two layers:

- **Layer 1** (`layer1_genetic_sim/`): an AlphaSimR stochastic simulation that compares manual and tool phenotyping under one budget, across a factorial grid of tool profiles (measurement-error reduction, per-plot cost reduction, cycle compression, fixed cost, budget). Each profile is run with paired replicates under common random numbers.
- **Layer 2** (`layer2_economics/`): converts simulated genetic gain into logistic adoption and 25-year net present value, then runs global sensitivity analysis (exact Sobol and Shapley indices on the grid, PAWN, active subspaces, dynamic Sobol), value-of-information analysis (EVPPI) and robustness scenarios.

Three meta-analyses (`paper/lit/`) anchor the model inputs: realized genetic gain in common bean (MA1), agreement between image-based tools and manual references (MA2), and time and cost savings of digital phenotyping (MA3).

Every number in the manuscript is a LaTeX macro generated from the outputs (`paper/numbers.tex`, with its source file listed in `paper/numbers.csv`).

## Repository layout

| Path | Contents |
|---|---|
| `config/` | Genetic, breeding and economic parameters; sweep grids (`config/sweeps/`) |
| `layer1_genetic_sim/R/` | Founder population, breeding pipeline, budget allocation, sweep driver (`07_param_sweep.R`) |
| `layer1_genetic_sim/outputs/sweep_v2*/` | Layer 1 results: per-replicate files (`cell_*.csv.gz`), per-cell summaries and `cell_summaries.csv` for the main grid, the two extension blocks, the robustness block and the stress block |
| `layer2_economics/R/` | Economic model, sensitivity and decision analyses |
| `layer2_economics/outputs_v2/` | Layer 2 results (lookup table, sensitivity indices, EVPPI, scenarios) |
| `paper/lit/` | Meta-analysis script, extraction sheets with page or table locators, Crossref verification, search logs and PRISMA counts |
| `paper/R/`, `paper/py/` | Number and table export (R); figures, graphical abstract and supplementary material (Python) |
| `paper/` | Manuscript source (MDPI LaTeX template), figures, tables and supplementary material |
| `slurm/` | Scripts used to run Layer 1 as Slurm array jobs |

## Requirements

- R 4.3 or later with `AlphaSimR` (2.1.0 was used for Layer 1), `data.table`, `yaml`, `arrow`, `mvtnorm` and `metafor` (5.2.1). `layer1_genetic_sim/R/00_install_deps.R` installs the R packages.
- Python 3.10 or later with `matplotlib`, `numpy`, `pandas`, `pyyaml` and `pillow`.
- A LaTeX distribution with `latexmk` to compile the manuscript.

## Reproducing the results

Run all commands from the repository root.

```bash
make numbers      # regenerate paper/numbers.tex and numbers.csv from the stored outputs
make analyses     # rebuild out-of-date Layer 2 analyses (slow)
make figures      # figures and graphical abstract
make supplement   # supplementary Figure S1 and Table S1
make paper        # compile paper/main.tex
Rscript paper/lit/meta_analysis.R   # re-run the three meta-analyses
```

The Layer 1 results are included, so the commands above do not need a cluster. To re-run Layer 1, simulate one cell with

```bash
Rscript layer1_genetic_sim/R/07_param_sweep.R <cell_id> 100
```

or submit the whole grid with the scripts in `slurm/` (set the account and paths for your cluster), then run `make aggregate` to rebuild the cell summaries. The grid is selected with the environment variables `TOOL_SWEEP_GRID` (a file in `config/sweeps/`) and `TOOL_SWEEP_SUBDIR` (the output folder).

## Citation

See `CITATION.cff`. Please cite the article and the archived version of this repository.

## License

Code is released under the MIT License (`LICENSE`). Simulation outputs, extraction sheets and other data files are released under the Creative Commons Attribution 4.0 International license (CC BY 4.0). The MDPI LaTeX class files in `paper/Definitions/` remain under their original terms.
