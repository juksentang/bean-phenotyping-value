# ============================================================================
# Reproducible pipeline: sweep summaries → L2 analyses → paper numbers
# ============================================================================
# Every number in the manuscript comes from paper/numbers.tex, generated here.
#
#   make numbers     regenerate paper/numbers.tex (+ numbers.csv) — cheap
#   make analyses    rebuild every out-of-date L2 analysis output   — expensive
#   make paper       numbers + compile paper/main.tex
#   make -n numbers  dry run: show what would be rebuilt
#
# The Layer-1 sweep (AlphaSimR, 100 paired replicates per cell) runs on a Slurm
# cluster via slurm/; the cell summaries are treated as sources here. After
# adding cells, `make aggregate` rebuilds them.
#
# Paths are resolved in layer2_economics/R/00_paths.R (TOOL_L1_SUBDIR,
# TOOL_L2_OUT; defaults: the sweep_v2 blocks and layer2_economics/outputs_v2).

SHELL  := /bin/bash
R      := Rscript
R2     := layer2_economics/R
OUT    := layer2_economics/outputs_v2
SWEEP  := layer1_genetic_sim/outputs/sweep_v2
# Layer 2 stacks the main grid and the extension blocks (00_paths.R)
L1     := $(SWEEP)/cell_summaries.csv $(SWEEP)_extA/cell_summaries.csv $(SWEEP)_extB/cell_summaries.csv

CFG    := config/economic_params.yaml config/breeding_params.yaml config/econ_sweep_grid.yaml config/revision_scenarios.yaml
# layer1 utils.R is deliberately NOT a prerequisite: Layer 2 only sources it for
# load_config(), which the v2 edits left untouched (they change phenotyping and
# the gain baseline, i.e. Layer 1 only). Listing it made every Layer-2 target
# stale, so a bare `make` would have rebuilt (overwritten) the v1 Layer-2
# results. If load_config() changes, force the rebuild explicitly (make -B).
CORE   := $(R2)/00_paths.R $(R2)/01_gain_to_value.R $(R2)/03_npv_irr.R $(CFG)
SOBOLR := $(CORE) $(R2)/07_sweep_economics.R $(R2)/08_break_even.R $(R2)/09_sobol.R
LOOKUP := $(OUT)/lookup_table.parquet

.PHONY: all numbers tables analyses figures paper aggregate clean-numbers
all: numbers

# --- Layer 1 aggregation (sweep cells pulled from cluster) -------------------
aggregate:
	$(R) -e 'source("layer1_genetic_sim/R/07_param_sweep.R"); aggregate_sweep()'

# --- Layer 2 -----------------------------------------------------------------
$(LOOKUP): $(L1) $(CORE) $(R2)/07_sweep_economics.R
	$(R) $(R2)/07_sweep_economics.R 2>&1 | tee $(OUT)/sweep_build.log

$(OUT)/break_even_marginals.csv $(OUT)/break_even_marginals_stressed.csv $(OUT)/profit_surface.csv &: \
		$(LOOKUP) $(CORE) $(R2)/08_break_even.R
	$(R) $(R2)/08_break_even.R

# Attribution and decision analyses are computed EXACTLY on the simulated grid (17_exact_gsa.R):
# Sobol, Shapley (independent and copula), PAWN, horizon-dependent indices, EVPPI, Layer-1 noise check.
# The Monte Carlo estimators of the first submission (09, 11, 12, 13, 15) are kept in the repository;
# their outputs are archived in outputs_v2/archive_mc_nearest_grid.
$(OUT)/sobol_indices_full.csv $(OUT)/sobol_indices_cyc0.csv $(OUT)/shapley_indices_indep.csv \
$(OUT)/shapley_indices_copula.csv $(OUT)/pawn_indices.csv $(OUT)/dynamic_sensitivity_ST_wide.csv \
$(OUT)/evppi.csv $(OUT)/sobol_bootstrap_layer1.csv $(OUT)/exact_gsa_run.log &: $(L1) $(LOOKUP) $(CORE) $(R2)/17_exact_gsa.R
	$(R) -e 'source("$(R2)/17_exact_gsa.R"); run_exact_gsa(B = 120)' 2>&1 | tee $(OUT)/exact_gsa_console.log

$(OUT)/active_subspace_eigenvalues.csv: $(L1) $(SOBOLR) $(R2)/14_active_subspaces.R
	$(R) $(R2)/14_active_subspaces.R 2>&1 | tee $(OUT)/active_sub_run.log

$(OUT)/scenario_summary.csv $(OUT)/stress_block_summary.csv &: $(L1) $(CORE) $(R2)/17_exact_gsa.R $(R2)/18_revision_scenarios.R \
		config/revision_scenarios.yaml layer1_genetic_sim/outputs/sweep_v2_stress/cell_summaries.csv
	$(R) $(R2)/18_revision_scenarios.R 2>&1 | tee $(OUT)/scenario_run.log

$(OUT)/scenario_decomposition.csv: $(LOOKUP) $(R2)/15_decision_voi.R
	$(R) $(R2)/15_decision_voi.R

$(OUT)/figures/method_comparison.png: $(OUT)/sobol_indices_full.csv $(OUT)/shapley_indices_indep.csv \
		$(OUT)/shapley_indices_copula.csv $(OUT)/pawn_indices.csv $(OUT)/evppi.csv $(R2)/16_compare_methods.R
	$(R) $(R2)/16_compare_methods.R

ANALYSES := $(LOOKUP) $(OUT)/break_even_marginals_stressed.csv $(OUT)/sobol_indices_full.csv \
            $(OUT)/sobol_indices_cyc0.csv $(OUT)/shapley_indices_copula.csv \
            $(OUT)/dynamic_sensitivity_ST_wide.csv $(OUT)/pawn_indices.csv \
            $(OUT)/active_subspace_eigenvalues.csv $(OUT)/evppi.csv \
            $(OUT)/scenario_decomposition.csv $(OUT)/scenario_summary.csv $(OUT)/stress_block_summary.csv \
            $(OUT)/sobol_bootstrap_layer1.csv $(OUT)/figures/method_comparison.png
analyses: $(ANALYSES)

# --- Paper numbers -----------------------------------------------------------
paper/numbers.tex paper/numbers.csv &: paper/lit/extract/cost_anchors.csv config/software_versions_hpc.yaml paper/R/export_numbers.R $(L1) $(CFG) config/genetic_params.yaml \
		$(filter-out $(OUT)/figures/method_comparison.png,$(ANALYSES))
	$(R) paper/R/export_numbers.R

numbers: paper/numbers.tex

tables: paper/tables/table2_meta.tex paper/tables/tableA1_ma2_appraisal.tex

paper/tables/table2_meta.tex paper/tables/tableA1_ma2_appraisal.tex &: paper/R/export_tables.R paper/lit/meta_results.csv paper/lit/meta_inputs_used.csv \
		paper/lit/extract/synthesis/bibkeys.tsv paper/lit/refs.bib
	$(R) paper/R/export_tables.R

paper: paper/numbers.tex paper/tables/table2_meta.tex paper/tables/tableA1_ma2_appraisal.tex
	cd paper && latexmk -pdf -interaction=nonstopmode main.tex

clean-numbers:
	rm -f paper/numbers.tex paper/numbers.csv

# --- Figures --------------------------------------------------------------------
FIG_SRC := paper/py/make_figures.py paper/py/figstyle.py
figures: paper/figures/captions.md
paper/figures/captions.md: $(FIG_SRC) $(L1) $(ANALYSES) \
		layer1_genetic_sim/outputs/sweep_v2_robust/cell_summaries.csv \
		paper/lit/meta_inputs_used.csv paper/lit/meta_results.csv
	python3 paper/py/make_figures.py

figures: paper/figures/graphical_abstract.png
paper/figures/graphical_abstract.png: paper/py/make_graphical_abstract.py paper/py/figstyle.py \
		paper/numbers.csv layer1_genetic_sim/outputs/sweep_v2_robust/cell_summaries.csv
	python3 paper/py/make_graphical_abstract.py

# --- Supplementary materials (Table S1, Figure S1) --------------------------------
supplement: paper/supplementary/supplementary.pdf
paper/supplementary/figS1_prisma.pdf paper/supplementary/tableS1_queries.tex &: paper/py/make_supplement.py \
		paper/py/figstyle.py paper/lit/search/prisma/prisma_flow.csv paper/lit/search/prisma/search_queries.csv
	python3 paper/py/make_supplement.py
paper/supplementary/supplementary.pdf: paper/supplementary/supplementary.tex paper/supplementary/figS1_prisma.pdf \
		paper/supplementary/tableS1_queries.tex
	cd paper/supplementary && latexmk -pdf -interaction=nonstopmode supplementary.tex
