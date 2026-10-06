# ============================================================================
# Export every number quoted in the manuscript → paper/numbers.tex (+ .csv)
# ============================================================================
# The manuscript must never hard-code a result. Each value is emitted as a
# LaTeX macro (\NumXxx) and also logged to numbers.csv for auditing/diffing.
# Run via `make numbers` from the project root.

suppressMessages({
  library(data.table)
  library(arrow)
})
source(file.path(getwd(), "layer2_economics/R/08_break_even.R"))  # load_config, baseline_slice, stressed_slice

out_dir  <- l2_out_dir()        # 00_paths.R: TOOL_L2_OUT (default outputs_v2)
l1_path  <- l1_summary_path()   # 00_paths.R: TOOL_L1_SUBDIR (default sweep_v2)
l1_files <- l1_path             # kept: l1_path is overwritten with a label below
tex_path <- "paper/numbers.tex"
csv_path <- "paper/numbers.csv"

# --- Registry ---------------------------------------------------------------
reg <- new.env()
reg$rows <- list()
put <- function(name, value, fmt = "%.2f", src = "") {
  stopifnot(grepl("^[A-Za-z]+$", name), length(value) == 1, !is.na(value))
  if (!is.null(reg$rows[[name]])) stop("duplicate number: ", name)
  txt <- if (is.character(value)) value else sprintf(fmt, value)
  reg$rows[[name]] <- data.table(name = name, value = txt,
                                 raw = as.character(value),
                                 source = sub(paste0("^", getwd(), "/"), "", src, fixed = FALSE))
}
money_M <- function(x) x / 1e6
pct     <- function(x) 100 * x
int     <- function(x) formatC(x, format = "d", big.mark = ",")
# Number-to-word suffixes for macro names (levels of the explored inputs)
word_of <- function(x) {
  w <- c(`0` = "Zero", `5` = "Five", `10` = "Ten", `20` = "Twenty", `30` = "Thirty",
         `40` = "Forty", `50` = "Fifty", `60` = "Sixty", `70` = "Seventy")
  unname(w[as.character(round(x))])
}
# Range text "lo to hi" with a thousands separator
rng <- function(x, scale = 1, digits = 0)
  paste(formatC(range(x) * scale, format = "f", digits = digits, big.mark = ","), collapse = " to ")
# Level list "a, b and c"
lvls <- function(x, digits = 2) {
  v <- formatC(sort(unique(x)), format = "f", digits = digits)
  paste0(paste(v[-length(v)], collapse = ", "), " and ", v[length(v)])
}

read_src <- function(f) paste(readLines(f, warn = FALSE), collapse = "\n")
grab    <- function(txt, pat, i = 1) regmatches(txt, regexec(pat, txt))[[1]][i + 1]
# Parameter name → macro stem
stem <- c(cycle_compression = "Cyc",  discount_rate   = "Disc",
          bean_price        = "Price", adoption_ceiling = "Adopt",
          tool_fixed_cost    = "FixCost", cost_reduction = "CostRed",
          total_budget      = "Budget",  error_reduction = "ErrRed")
num_word <- c("Zero", "One", "Two", "Three")

# --- Configuration (design inputs quoted in Methods) ------------------------
econ  <- load_config("economic_params")
breed <- load_config("breeding_params")
gen   <- load_config("genetic_params")
grid  <- yaml::read_yaml("config/econ_sweep_grid.yaml")

put("CfgCycleYears",   breed$cycle_years, "%d", "breeding_params.yaml")
put("CfgMeanYield",    econ$mean_yield_kgha, "%d", "economic_params.yaml")
put("CfgSigmaA",       econ$sigma_a_kgha, "%d", "economic_params.yaml")
put("CfgHtwoYield",    gen$h2_yield, "%.2f", "genetic_params.yaml")
put("CfgHorizon",      econ$evaluation_horizon, "%d", "economic_params.yaml")

# --- Layer 1: per-cycle genetic gain ----------------------------------------
l1cal <- load_l1_summaries(l1_path)
l1 <- rbindlist(lapply(l1_path, fread), use.names = TRUE, fill = TRUE)  # raw (uncalibrated) gains
l1_path <- paste(basename(dirname(l1_path)), collapse = "+")
put("GainTargetPct", if (is.null(econ$gain_target_pct_per_yr)) "none" else econ$gain_target_pct_per_yr,
    "%.1f", "economic_params.yaml")
put("GainScale",     attr(l1cal, "gain_scale"), "%.2f", "load_l1_summaries()")
put("DGtradPerYearCal", mean(l1cal$mean_dg_trad) / breed$cycle_years, "%.1f", "load_l1_summaries()")
put("NcellsOK",       int(nrow(l1)), src = l1_path)
put("NrepsPerCell",   int(unique(l1$n_reps)[1]), src = l1_path)
put("NcellsNegDiff",  int(sum(l1$mean_dg_diff < 0)), src = l1_path)
put("DGtradPerCycle", mean(l1$mean_dg_trad), "%.0f", l1_path)
put("DGtradPerYear",  mean(l1$mean_dg_trad) / breed$cycle_years, "%.1f", l1_path)
put("DGdiffPerCycle", mean(l1$mean_dg_diff), "%.1f", l1_path)
put("DGdiffPct",      pct(mean(l1$mean_dg_diff) / mean(l1$mean_dg_trad)), "%.1f", l1_path)
if (!is.null(l1$mean_dg_vs_base_trad))
  put("DGvsBaseTradPerCycle", mean(l1$mean_dg_vs_base_trad), "%.0f", l1_path)
put("DGdiffSE",       mean(l1$sd_dg_diff / sqrt(l1$n_reps)), "%.1f", l1_path)
for (e in sort(unique(l1$error_reduction))) {
  put(paste0("DGdiffErr", c("Zero","Ten","Twenty","Thirty","Forty","Fifty")[round(e * 10) + 1]),
      mean(l1[error_reduction == e]$mean_dg_diff), "%.1f", l1_path)
}

# --- Layer 2: lookup table --------------------------------------------------
lk <- as.data.table(read_parquet(file.path(out_dir, "lookup_table.parquet")))
put("NlookupRows",  int(nrow(lk)), src = "lookup_table.parquet")
put("PrNPVpos",     pct(mean(lk$npv > 0)), "%.1f", "lookup_table.parquet")

base <- baseline_slice(lk, econ)
put("BasePrice", base$bean_price[1],       "%.2f", "baseline_slice()")
put("BaseDisc",  pct(base$discount_rate[1]), "%.0f", "baseline_slice()")
put("BaseAdopt", pct(base$adoption_ceiling[1]), "%.0f", "baseline_slice()")
put("BaseArea",  base$total_area_ha[1] / 1e6, "%.0f", "baseline_slice()")
for (k in 0:3) {
  s <- base[cycle_compression == k]
  put(paste0("EAAcyc",    num_word[k + 1]), median(s$eaa_per_ha_yr), "%.2f", "baseline_slice()")
  put(paste0("NPVcyc",    num_word[k + 1]), money_M(median(s$npv)), "%.1f", "baseline_slice()")
  put(paste0("PrPosCyc",  num_word[k + 1]), pct(mean(s$npv > 0)), "%.1f", "baseline_slice()")
}
put("EAAperYearCompressed",
    (median(base[cycle_compression == 3]$eaa_per_ha_yr) -
     median(base[cycle_compression == 0]$eaa_per_ha_yr)) / 3, "%.2f", "baseline_slice()")

str_ <- stressed_slice(lk)
put("StressMedianNPV", money_M(median(str_$npv)), "%.1f", "stressed_slice()")

be_s <- fread(file.path(out_dir, "break_even_marginals_stressed.csv"))
put("StressMinMarginalNPV", money_M(min(be_s$npv)), "%.3f", "break_even_marginals_stressed.csv")

cons <- base[cycle_compression == 0]
put("ConsEAAmean", mean(cons$eaa_per_ha_yr), "%.3f", "baseline_slice() cyc=0")
put("ConsEAAsd",   sd(cons$eaa_per_ha_yr),   "%.3f", "baseline_slice() cyc=0")

# --- Global sensitivity -----------------------------------------------------
sob <- fread(file.path(out_dir, "sobol_indices_full.csv"))
stopifnot(all(sob$S1 >= -1e-12))            # exact indices on the grid: no negative first-order estimates
f3 <- function(x) if (x < 0.0005) "$<$0.001" else if (x < 0.01) sprintf("%.3f", x) else sprintf("%.2f", x)  # small indices: not "0.00"
for (i in seq_len(nrow(sob))) {
  p <- stem[[sob$param[i]]]
  put(paste0("SobolSone", p), f3(max(sob$S1[i], 0)), src = "sobol_indices_full.csv (exact)")
  put(paste0("SobolST",   p), f3(sob$ST[i]), src = "sobol_indices_full.csv (exact)")
}
put("SobolSumSone", sum(pmax(sob$S1, 0)), "%.2f", "sobol_indices_full.csv")
put("SobolSumST",   sum(sob$ST), "%.2f", "sobol_indices_full.csv")
put("SobolEconST",  sum(sob[group == "L2 economic"]$ST), "%.2f", "sobol_indices_full.csv")
put("SobolMinorSTmax", max(sob[param %in% c("error_reduction", "cost_reduction", "tool_fixed_cost", "total_budget")]$ST), "%.3f", "sobol_indices_full.csv")
put("SobolMinorSTsum", sum(sob[param %in% c("error_reduction", "cost_reduction", "tool_fixed_cost", "total_budget")]$ST), "%.3f", "sobol_indices_full.csv")
bsl <- fread(file.path(out_dir, "sobol_bootstrap_layer1.csv"))
bsl <- merge(bsl, sob[, .(param, ST)], by = "param")
put("SobolBootMaxShift", max(pmax(abs(bsl$ST_lo - bsl$ST), abs(bsl$ST_hi - bsl$ST))), "%.3f", "sobol_bootstrap_layer1.csv (largest shift of any total index)")
elog <- read_src(file.path(out_dir, "exact_gsa_run.log"))
put("BootB", as.numeric(grab(elog, "layer1_bootstrap B=([0-9]+)")), "%d", "exact_gsa_run.log")
put("NexactGrid", int(as.numeric(grab(elog, "n = ([0-9]+)"))), src = "exact_gsa_run.log")
sob0 <- fread(file.path(out_dir, "sobol_indices_cyc0.csv"))
for (i in seq_len(nrow(sob0))) {
  put(paste0("SobolCycZeroST", stem[[sob0$param[i]]]), sob0$ST[i], "%.2f", "sobol_indices_cyc0.csv (exact)")
}

shc <- fread(file.path(out_dir, "shapley_indices_copula.csv"))
for (i in seq_len(nrow(shc))) {
  put(paste0("Shapley", stem[[shc$param[i]]]), f3(shc$shapley_effect[i]), src = "shapley_indices_copula.csv (exact)")
}

pawn <- fread(file.path(out_dir, "pawn_indices.csv"))
for (i in seq_len(nrow(pawn))) {
  put(paste0("Pawn", stem[[pawn$param[i]]]), pawn$PAWN_median[i], "%.3f", "pawn_indices.csv (exact)")
}

as_ev <- fread(file.path(out_dir, "active_subspace_eigenvalues.csv"))
put("ASshareOne", pct(as_ev[k == 1]$var_share), "%.1f", "active_subspace_eigenvalues.csv")
put("AScumTwo",   pct(as_ev[k == 2]$cum_share), "%.1f", "active_subspace_eigenvalues.csv")

dyn <- fread(file.path(out_dir, "dynamic_sensitivity_ST_wide.csv"))
put("DynSTCycYearFive", dyn[year == 5]$cycle_compression, "%.2f", "dynamic_sensitivity_ST_wide.csv")
put("DynSTCycFinal",    dyn[year == max(year)]$cycle_compression, "%.2f", "dynamic_sensitivity_ST_wide.csv")

# --- Decision analysis ------------------------------------------------------
ev <- fread(file.path(out_dir, "evppi.csv"))
put("EVprior", money_M(mean(lk$npv)), "%.0f", "lookup_table.parquet")  # = E[NPV] = indifference hurdle
ind <- ev[at_indifference == TRUE]
stopifnot(nrow(ind) > 0)
for (k in seq_len(nrow(ind))) {
  put(paste0("EVPPIIndiff", stem[[ind$param[k]]]), money_M(max(ind$EVPPI[k], 0)), "%.1f", "evppi.csv")
}
scn <- fread(file.path(out_dir, "scenario_decomposition.csv"))
for (s in c("pessimistic", "other", "optimistic")) {
  r <- scn[scen == s]
  S <- tools::toTitleCase(s)
  put(paste0("Scn", S, "Prob"),   pct(r$prob), "%.1f", "scenario_decomposition.csv")
  put(paste0("Scn", S, "Mean"),   money_M(r$mean_npv), "%.1f", "scenario_decomposition.csv")
  put(paste0("Scn", S, "Median"), money_M(r$median_npv), "%.1f", "scenario_decomposition.csv")
  put(paste0("Scn", S, "PrNeg"),  pct(r$pr_neg), "%.1f", "scenario_decomposition.csv")
}


# ============================================================================
# Additional design constants, simulation-settings and result numbers
# ============================================================================
read_src <- function(f) paste(readLines(f, warn = FALSE), collapse = "\n")
grab    <- function(txt, pat, i = 1) regmatches(txt, regexec(pat, txt))[[1]][i + 1]
sweep_dir <- file.path("layer1_genetic_sim/outputs")
rb_path   <- file.path(sweep_dir, "sweep_v2_robust/cell_summaries.csv")
rb <- fread(rb_path)
master <- load_config("master_config")

# --- Design constants: pipeline, genome, costs (config files) ---------------
put("CfgNcrosses",  breed$n_crosses_per_cycle,  "%d", "breeding_params.yaml")
put("CfgNprogeny",  breed$n_progeny_per_cross,  "%d", "breeding_params.yaml")
put("CfgNssd",      breed$n_ssd_generations,    "%d", "breeding_params.yaml")
for (st in c("ON", "PYT", "AYT", "NPT")) {
  sc <- breed$stages[[st]]
  put(paste0("CfgStage", st, "Lines"),   int(sc$n_genotypes),   src = "breeding_params.yaml")
  put(paste0("CfgStage", st, "Reps"),    sc$reps,               "%d", "breeding_params.yaml")
  put(paste0("CfgStage", st, "Envs"),    sc$envs,               "%d", "breeding_params.yaml")
  put(paste0("CfgStage", st, "SelFrac"), pct(sc$selection_intensity), "%.0f", "breeding_params.yaml")
}
# Released lines: rule in 03_budget_optimizer.R, max(3, round(NPT lines x NPT selection intensity))
put("CfgNreleased", max(3, round(breed$stages$NPT$n_genotypes * breed$stages$NPT$selection_intensity)),
    "%d", "03_budget_optimizer.R rule applied to breeding_params.yaml")
put("CfgNfounders", gen$n_founders,         "%d", "genetic_params.yaml")
put("CfgSegSites",  int(gen$seg_sites_per_chr), src = "genetic_params.yaml")
put("CfgQTLperChr", gen$n_qtl_per_chr,      "%d", "genetic_params.yaml")
put("CfgGxEOnstation", gen$gxe_cor_onstation, "%.2f", "genetic_params.yaml")
put("CfgCostPerPlot",  econ$cost_per_plot_traditional, "%d", "economic_params.yaml")
put("CfgDiffK",        econ$diffusion_rate_k, "%.2f", "economic_params.yaml")
put("CfgDiffMid",      econ$diffusion_midpoint_years, "%d", "economic_params.yaml")

# Budget-allocation search settings, read from the optimizer source
opt <- read_src("layer1_genetic_sim/R/03_budget_optimizer.R")
ng  <- as.numeric(strsplit(grab(opt, "n_grid <- seq\\(([^)]*)\\)"), "[ ,]+|by = ")[[1]])
ng  <- ng[!is.na(ng)]
put("CfgPYTmin",  int(ng[1]), src = "03_budget_optimizer.R n_grid")
put("CfgPYTmax",  int(ng[2]), src = "03_budget_optimizer.R n_grid")
put("CfgPYTstep", int(ng[3]), src = "03_budget_optimizer.R n_grid")
put("CfgPYTenvMax", as.numeric(grab(opt, "e_grid <- seq\\(1, ([0-9]+)")), "%d", "03_budget_optimizer.R e_grid")
put("CfgONmult", as.numeric(grab(opt, "n_on <- min\\(n_pyt \\* ([0-9]+)")), "%d", "03_budget_optimizer.R")
put("CfgONcap",  int(as.numeric(grab(opt, "n_on <- min\\(n_pyt \\* [0-9]+, ([0-9]+)"))), src = "03_budget_optimizer.R")

# --- Tool inputs: central values, explored ranges (percent for the two reductions) ---
put("CfgErrRedCentral",  pct(econ$tool_error_reduction_frac), "%.0f", "economic_params.yaml")
put("CfgCostRedCentral", pct(econ$tool_cost_reduction_frac),  "%.0f", "economic_params.yaml")
put("CfgCycCompCentral", econ$tool_cycle_compression,         "%d",   "economic_params.yaml")
put("CfgFixCostCentral", int(econ$tool_fixed_cost),                   src = "economic_params.yaml")
put("CfgBudgetCentral",  int(master$total_budget),                   src = "master_config.yaml")
put("CfgErrRedMin",  pct(min(l1$error_reduction)),  "%.0f", l1_path)
put("CfgErrRedMax",  pct(max(l1$error_reduction)),  "%.0f", l1_path)
put("CfgCostRedMin", pct(min(l1$cost_reduction)),   "%.0f", l1_path)
put("CfgCostRedMax", pct(max(l1$cost_reduction)),   "%.0f", l1_path)
put("CfgErrRedRange",  rng(l1$error_reduction, 100), src = l1_path)
put("CfgCostRedRange", rng(l1$cost_reduction, 100),  src = l1_path)
put("CfgCycRange",     rng(l1$cycle_compression),    src = l1_path)
put("CfgFixCostRange", rng(l1$tool_fixed_cost),       src = l1_path)
put("CfgBudgetRange",  rng(l1$total_budget),         src = l1_path)
put("CfgRobHtwoLevels",   lvls(rb$h2_yield),           src = rb_path)
put("CfgRobGxELevels",  lvls(rb$gxe_cor_onstation),  src = rb_path)
put("CfgRobHtwoLow",    min(rb$h2_yield), "%.2f", rb_path)
put("CfgRobHtwoHigh",   max(rb$h2_yield), "%.2f", rb_path)

# --- Economic grid ranges ---------------------------------------------------
put("CfgPriceRange", rng(grid$econ_grid$bean_price_usd_per_kg, digits = 2),     src = "econ_sweep_grid.yaml")
put("CfgDiscRange",  rng(grid$econ_grid$discount_rate, 100),                    src = "econ_sweep_grid.yaml")
put("CfgAdoptRange", rng(grid$econ_grid$adoption_ceiling, 100),                 src = "econ_sweep_grid.yaml")
put("CfgAreaRange",  rng(grid$econ_grid$total_area_ha, 1e-6),                   src = "econ_sweep_grid.yaml")

# --- Block sizes --------------------------------------------------------------
blk <- vapply(l1_files, function(f) nrow(fread(f)), numeric(1))
put("NcellsMain",   int(blk[1]), src = l1_files[1])
put("NcellsExtA",   int(blk[2]), src = l1_files[2])
put("NcellsExtB",   int(blk[3]), src = l1_files[3])
put("NcellsRobust", int(nrow(rb)), src = rb_path)
put("NeconStates",  int(nrow(unique(lk[, .(bean_price, discount_rate, adoption_ceiling, total_area_ha)]))),
    src = "lookup_table.parquet")

# --- Layer 1: uplift summaries (uplift = 100 x paired difference / manual gain; scale-free) ---
L <- breed$cycle_years
l1[, uplift := 100 * mean_dg_diff / mean_dg_trad]
l1[, annual := 100 * ((1 + uplift / 100) * L / (L - cycle_compression) - 1)]
is_meta <- with(l1, abs(error_reduction - econ$tool_error_reduction_frac) < 1e-9 &
                    abs(cost_reduction  - econ$tool_cost_reduction_frac)  < 1e-9)
src_l1 <- "cell_summaries.csv (main grid and extensions)"
put("CycOneYearPct", 100 * (L / (L - 1) - 1), "%.1f", "closed form, 100 (L/(L-1) - 1)")
put("DGtradPctPerYearRaw", 100 * (mean(l1$mean_dg_trad) / L) / econ$mean_yield_kgha, "%.1f", src_l1)
put("NcellsNegDiffErrZero", int(sum(l1$mean_dg_diff < 0 & l1$error_reduction == 0)), src = src_l1)
for (k in 0:3) {
  put(paste0("AnnUpliftCyc", num_word[k + 1]), median(l1[cycle_compression == k]$annual), "%.0f", src_l1)
  put(paste0("AnnUpliftMetaCyc", num_word[k + 1]), median(l1[cycle_compression == k & is_meta]$annual), "%.0f", src_l1)
}
put("UpliftMeta", median(l1[is_meta]$uplift), "%.1f", src_l1)
for (e in sort(unique(l1$error_reduction)))
  put(paste0("UpliftErr", word_of(100 * e)), median(l1[error_reduction == e]$uplift), "%.1f", src_l1)
for (cr in sort(unique(l1$cost_reduction)))
  put(paste0("UpliftCost", word_of(100 * cr)), median(l1[cost_reduction == cr]$uplift), "%.1f", src_l1)
put("UpliftBudgetLow",  median(l1[total_budget == min(total_budget)]$uplift), "%.1f", src_l1)
put("UpliftBudgetHigh", median(l1[total_budget == max(total_budget)]$uplift), "%.1f", src_l1)

# Robustness block: uplift per (h2, GxE, error reduction) averaged over cost reduction and compression,
# then averaged over GxE; months of compression with the same annual-gain effect (Figure 2a)
rb[, uplift := 100 * mean_dg_diff / mean_dg_trad]
rbt <- rb[, .(uplift = mean(uplift)), by = .(h2_yield, gxe_cor_onstation, error_reduction)]
rbt[, years := L - L / (1 + uplift / 100)]
h2_lv <- sort(unique(rb$h2_yield)); h2_nm <- c("Low", "Mid", "High")
e_lo <- min(rb$error_reduction); e_hi <- max(rb$error_reduction)
e_ev <- econ$tool_error_reduction_frac
for (i in seq_along(h2_lv)) {
  d <- rbt[abs(h2_yield - h2_lv[i]) < 1e-9]
  put(paste0("RobUpliftHtwo", h2_nm[i], "ErrZero"),  mean(d[abs(error_reduction - e_lo) < 1e-9]$uplift), "%.1f", rb_path)
  put(paste0("RobUpliftHtwo", h2_nm[i], "ErrFifty"), mean(d[abs(error_reduction - e_hi) < 1e-9]$uplift), "%.1f", rb_path)
  put(paste0("ExchMoHtwo", h2_nm[i], "ErrTen"),   12 * mean(d[abs(error_reduction - e_ev) < 1e-9]$years), "%.1f", rb_path)
  put(paste0("ExchMoHtwo", h2_nm[i], "ErrFifty"), 12 * mean(d[abs(error_reduction - e_hi) < 1e-9]$years), "%.1f", rb_path)
}

# --- Layer 2: NPV by compression (all profiles and states), evidence-based profile, price slice ---
lk[, meta := abs(error_reduction - econ$tool_error_reduction_frac) < 1e-9 &
             abs(cost_reduction  - econ$tool_cost_reduction_frac)  < 1e-9]
for (k in 0:3) {
  put(paste0("NPVallCyc",  num_word[k + 1]), money_M(median(lk[cycle_compression == k]$npv)), "%.1f", "lookup_table.parquet")
  put(paste0("NPVmetaCyc", num_word[k + 1]), money_M(median(lk[cycle_compression == k & meta]$npv)), "%.1f", "lookup_table.parquet")
}
put("PrPosAllCycZero", pct(mean(lk[cycle_compression == 0]$npv > 0)), "%.1f", "lookup_table.parquet")
bslice <- lk[discount_rate == base$discount_rate[1] & adoption_ceiling == base$adoption_ceiling[1] &
             total_area_ha == base$total_area_ha[1]]
pr <- sort(unique(bslice$bean_price))
put("EAAcycThreePriceLow",  median(bslice[cycle_compression == 3 & bean_price == min(pr)]$eaa_per_ha_yr), "%.2f", "lookup_table.parquet")
put("EAAcycThreePriceHigh", median(bslice[cycle_compression == 3 & bean_price == max(pr)]$eaa_per_ha_yr), "%.2f", "lookup_table.parquet")

dyn_st <- dyn[, setdiff(names(dyn), "year"), with = FALSE]
put("DynSTCycPeak",     max(dyn$cycle_compression), "%.2f", "dynamic_sensitivity_ST_wide.csv")
put("DynSTCycPeakYear", dyn[which.max(cycle_compression)]$year, "%d", "dynamic_sensitivity_ST_wide.csv")
put("DynSTBudgetYearFive", dyn[year == 5]$total_budget, "%.2f", "dynamic_sensitivity_ST_wide.csv")

shi <- fread(file.path(out_dir, "shapley_indices_indep.csv"))
for (i in seq_len(nrow(shi)))
  put(paste0("ShapleyIndep", stem[[shi$param[i]]]), f3(shi$shapley_effect[i]), src = "shapley_indices_indep.csv (exact)")
stopifnot(all(shi$bracket_ok))              # S1 <= Shapley <= ST holds for every input under independence
put("ShapleyBracketN", sum(shi$bracket_ok), "%d", "shapley_indices_indep.csv")
asw <- fread(file.path(out_dir, "active_subspace_W.csv"))
put("ASloadCycOne", abs(asw[param == "cycle_compression"]$w1), "%.2f", "active_subspace_W.csv")  # eigenvector sign is arbitrary

# --- Settings of the sensitivity analyses (run logs and correlation matrix) ---
rho <- as.matrix(fread(file.path(out_dir, "shapley_correlation_matrix.csv"), drop = 1))
rownames(rho) <- colnames(rho)
put("CfgCopRhoPriceDisc",   rho["bean_price", "discount_rate"],      "%.1f", "shapley_correlation_matrix.csv")
put("CfgCopRhoAdoptBudget", rho["adoption_ceiling", "total_budget"], "%.1f", "shapley_correlation_matrix.csv")
put("CfgCopRhoCycErr",      rho["cycle_compression", "error_reduction"], "%.1f", "shapley_correlation_matrix.csv")
put("CfgCopRhoCostFix",     rho["cost_reduction", "tool_fixed_cost"], "%.1f", "shapley_correlation_matrix.csv")
alog <- read_src(file.path(out_dir, "active_sub_run.log"))
put("CfgASN",    as.numeric(grab(alog, "N=([0-9]+) gradients")), "%d", "active_sub_run.log")
put("CfgASstep", as.numeric(grab(alog, "h=([0-9.]+)")),          "%.2f", "active_sub_run.log")
voi <- read_src("layer2_economics/R/15_decision_voi.R")
hm  <- as.numeric(strsplit(grab(voi, "hurdle_mult <- c\\(([^)]*)\\)"), "[ ,]+")[[1]])
put("CfgHurdleMaxMult", max(hm, na.rm = TRUE), "%.0f", "15_decision_voi.R hurdle_mult")

# --- Meta-analyses (paper/lit/meta_results.csv, meta_leave_one_out.csv) ---
mr  <- fread("paper/lit/meta_results.csv")
loo <- fread("paper/lit/meta_leave_one_out.csv")
m   <- function(id) mr[analysis_id == id]
put("MAoneEst",  m("MA1.CB.primary")$estimate, "%.2f", "meta_results.csv MA1.CB.primary")
put("MAoneCIlo", m("MA1.CB.primary")$ci_lb,    "%.2f", "meta_results.csv MA1.CB.primary")
put("MAoneCIhi", m("MA1.CB.primary")$ci_ub,    "%.2f", "meta_results.csv MA1.CB.primary")
put("MAonePIlo", m("MA1.CB.primary")$pi_lb,    "%.2f", "meta_results.csv MA1.CB.primary")
put("MAonePIhi", m("MA1.CB.primary")$pi_ub,    "%.2f", "meta_results.csv MA1.CB.primary")
put("MAoneItwo", m("MA1.CB.primary")$I2_pct,   "%.0f", "meta_results.csv MA1.CB.primary")
put("MAoneK",     m("MA1.CB.primary")$k_effects,  "%d", "meta_results.csv MA1.CB.primary")
put("MAoneKclus", m("MA1.CB.primary")$k_clusters, "%d", "meta_results.csv MA1.CB.primary")
sens1 <- mr[grepl("^MA1\\.CB\\.", analysis_id) & role == "sensitivity"]
put("MAoneSensLo", min(sens1$estimate), "%.2f", "meta_results.csv MA1.CB sensitivity rows")
put("MAoneSensHi", max(sens1$estimate), "%.2f", "meta_results.csv MA1.CB sensitivity rows")
l1o <- loo[analysis_id == "MA1.CB.primary"]
put("MAoneLOOmin", min(l1o$estimate), "%.2f", "meta_leave_one_out.csv")
put("MAoneLOOmax", max(l1o$estimate), "%.2f", "meta_leave_one_out.csv")
put("MAtwoR",     m("MA2.primary")$estimate, "%.2f", "meta_results.csv MA2.primary")
put("MAtwoCIlo",  m("MA2.primary")$ci_lb,    "%.2f", "meta_results.csv MA2.primary")
put("MAtwoCIhi",  m("MA2.primary")$ci_ub,    "%.2f", "meta_results.csv MA2.primary")
put("MAtwoPIlo",  m("MA2.primary")$pi_lb,    "%.2f", "meta_results.csv MA2.primary")
put("MAtwoPIhi",  m("MA2.primary")$pi_ub,    "%.2f", "meta_results.csv MA2.primary")
put("MAtwoItwo",  m("MA2.primary")$I2_pct,   "%.0f", "meta_results.csv MA2.primary")
put("MAtwoK",     m("MA2.primary")$k_effects,  "%d", "meta_results.csv MA2.primary")
put("MAtwoKclus", m("MA2.primary")$k_clusters, "%d", "meta_results.csv MA2.primary")
put("MAtwoNoPre",  m("MA2.nopreprint")$estimate, "%.2f", "meta_results.csv MA2.nopreprint")
put("MAtwoHeight", m("MA2.trait.height")$estimate, "%.2f", "meta_results.csv MA2.trait.height")
put("MAtwoCounts", m("MA2.trait.plant_count")$estimate, "%.2f", "meta_results.csv MA2.trait.plant_count")
put("MAtwoDisease", m("MA2.trait.disease_severity")$estimate, "%.2f", "meta_results.csv MA2.trait.disease_severity")
put("MAtwoManualRel", m("MA2.manualrel")$estimate, "%.2f", "meta_results.csv MA2.manualrel")
put("ERtOne",   pct(m("MA2.ER.T1")$estimate), "%.0f", "meta_results.csv MA2.ER.T1 (median of Monte Carlo draws)")
put("ERtTwo",   pct(m("MA2.ER.T2")$estimate), "%.0f", "meta_results.csv MA2.ER.T2 (median of Monte Carlo draws)")
put("ERhtwoMedian", pct(m("MA2.ER.h2.median")$estimate), "%.0f", "meta_results.csv MA2.ER.h2.median")
put("MAthreeK",      m("MA3.time.primary")$k_effects, "%d", "meta_results.csv MA3.time.primary")
put("MAthreeTime",   m("MA3.time.primary")$estimate, "%.0f", "meta_results.csv MA3.time.primary")
put("MAthreeTimeCIlo", m("MA3.time.primary")$ci_lb, "%.0f", "meta_results.csv MA3.time.primary")
put("MAthreeTimeCIhi", m("MA3.time.primary")$ci_ub, "%.0f", "meta_results.csv MA3.time.primary")
put("MAthreeHandheldLo", m("MA3.time.fieldhandheld_withC")$ci_lb, "%.0f", "meta_results.csv MA3.time.fieldhandheld_withC (min)")
put("MAthreeHandheldHi", m("MA3.time.fieldhandheld_withC")$ci_ub, "%.0f", "meta_results.csv MA3.time.fieldhandheld_withC (max)")
put("MAthreeHandheldK",  m("MA3.time.fieldhandheld_withC")$k_effects, "%d", "meta_results.csv MA3.time.fieldhandheld_withC")
put("MAthreeCostRed",   pct(m("MA3.cost.implied_handheld")$estimate), "%.0f", "meta_results.csv MA3.cost.implied_handheld")
put("MAthreeCostRedLo", pct(m("MA3.cost.implied_handheld")$ci_lb),    "%.0f", "meta_results.csv MA3.cost.implied_handheld")
put("MAthreeCostRedHi", pct(m("MA3.cost.implied_handheld")$ci_ub),    "%.0f", "meta_results.csv MA3.cost.implied_handheld")


# --- Further design ranges, interval summaries and meta-analysis rows (reported for completeness) ---
put("CfgNchr",      gen$n_chr, "%d", "genetic_params.yaml")
put("CfgCycMax",    max(l1$cycle_compression), "%d", l1_path)
put("CfgFixCostMin", int(min(l1$tool_fixed_cost)), src = l1_path)
put("CfgFixCostMax", int(max(l1$tool_fixed_cost)), src = l1_path)
put("CfgBudgetMin",  int(min(l1$total_budget)),   src = l1_path)
put("CfgBudgetMax",  int(max(l1$total_budget)),   src = l1_path)
put("CfgPriceMin", min(grid$econ_grid$bean_price_usd_per_kg), "%.2f", "econ_sweep_grid.yaml")
put("CfgPriceMax", max(grid$econ_grid$bean_price_usd_per_kg), "%.2f", "econ_sweep_grid.yaml")
put("CfgDiscMin",  pct(min(grid$econ_grid$discount_rate)),    "%.0f", "econ_sweep_grid.yaml")
put("CfgDiscMax",  pct(max(grid$econ_grid$discount_rate)),    "%.0f", "econ_sweep_grid.yaml")
put("CfgAdoptMin", pct(min(grid$econ_grid$adoption_ceiling)), "%.0f", "econ_sweep_grid.yaml")
put("CfgAdoptMax", pct(max(grid$econ_grid$adoption_ceiling)), "%.0f", "econ_sweep_grid.yaml")
put("CfgAreaMin",  min(grid$econ_grid$total_area_ha) / 1e6,   "%.0f", "econ_sweep_grid.yaml")
put("CfgAreaMax",  max(grid$econ_grid$total_area_ha) / 1e6,   "%.0f", "econ_sweep_grid.yaml")
put("CfgRobGxELow",  min(rb$gxe_cor_onstation), "%.2f", rb_path)
put("CfgRobGxEHigh", max(rb$gxe_cor_onstation), "%.2f", rb_path)
for (k in 0:3) {
  put(paste0("AnnUpliftCyc", num_word[k + 1], "Qone"),   quantile(l1[cycle_compression == k]$annual, 0.25), "%.0f", src_l1)
  put(paste0("AnnUpliftCyc", num_word[k + 1], "Qthree"), quantile(l1[cycle_compression == k]$annual, 0.75), "%.0f", src_l1)
}
put("ASloadDiscTwo", abs(asw[param == "discount_rate"]$w2), "%.2f", "active_subspace_W.csv")
put("EVPPIzeroHurdleCyc", money_M(max(ev[param == "cycle_compression" & abs(hurdle_mult - 0) < 1e-9]$EVPPI, 0)), "%.1f", "evppi.csv")
put("EVPPItwoXCyc",       money_M(max(ev[param == "cycle_compression" & abs(hurdle_mult - 2) < 1e-9]$EVPPI, 0)), "%.2f", "evppi.csv")
put("MAoneOneper",  m("MA1.CB.clustercollapsed")$estimate, "%.2f", "meta_results.csv MA1.CB.clustercollapsed")
put("MAoneEmbrapa", m("MA1.CB.embrapaB08")$estimate,       "%.2f", "meta_results.csv MA1.CB.embrapaB08")
put("MAoneExIAC",   m("MA1.CB.exB02")$estimate,            "%.2f", "meta_results.csv MA1.CB.exB02")
put("MAoneImputed", m("MA1.CB.imputedA")$estimate,         "%.2f", "meta_results.csv MA1.CB.imputedA")
put("MAoneOther",   m("MA1.OC.primary")$estimate,          "%.2f", "meta_results.csv MA1.OC.primary")
put("MAtwoLegume",    m("MA2.legume")$estimate,    "%.2f", "meta_results.csv MA2.legume")
put("MAtwoNonLegume", m("MA2.nonlegume")$estimate, "%.2f", "meta_results.csv MA2.nonlegume")
put("MAtwoEggerP", m("MA2.egger")$Q_p,    "%.2f", "meta_results.csv MA2.egger (p of the asymmetry test)")
put("MAtwoRankP",  m("MA2.ranktest")$Q_p, "%.2f", "meta_results.csv MA2.ranktest (p of the rank test)")
put("MAthreeWalter", m("MA3.time.fieldhandheld")$estimate, "%.0f", "meta_results.csv MA3.time.fieldhandheld")
put("MAthreeNonHandShareLo", m("MA3.cost.share_nonhandling")$ci_lb, "%.0f", "meta_results.csv MA3.cost.share_nonhandling (min)")
put("MAthreeNonHandShareHi", m("MA3.cost.share_nonhandling")$ci_ub, "%.0f", "meta_results.csv MA3.cost.share_nonhandling (max)")

# --- Search and screening counts (PRISMA flow, reconciled from the search records) ---
pf <- fread("paper/lit/search/prisma/prisma_flow.csv")
sq <- fread("paper/lit/search/prisma/search_queries.csv")
pcount <- function(ma_, box_, source_) {
  v <- pf[analysis == ma_ & box == box_ & source == source_]$count
  if (length(v) != 1 || is.na(suppressWarnings(as.numeric(v)))) stop("PRISMA count missing: ", ma_, " / ", box_)
  as.numeric(v)
}
for (k in 1:3) {
  ma <- paste0("MA", k); w <- c("One", "Two", "Three")[k]
  put(paste0("Prisma", w, "Queries"),  int(nrow(sq[analysis == ma])), src = "search_queries.csv")
  put(paste0("Prisma", w, "Records"),  int(pcount(ma, "Records identified from databases/registers", "all")), src = "prisma_flow.csv")
  put(paste0("Prisma", w, "Unique"),   int(pcount(ma, "Records after duplicates removed", "database search")), src = "prisma_flow.csv")
  put(paste0("Prisma", w, "Other"),    int(pcount(ma, "Records identified via other methods", "other methods")), src = "prisma_flow.csv")
  put(paste0("Prisma", w, "Sought"),   int(pcount(ma, "Reports sought for retrieval", "all")), src = "prisma_flow.csv")
  put(paste0("Prisma", w, "Assessed"), int(pcount(ma, "Reports assessed for eligibility (full text read)", "all")), src = "prisma_flow.csv")
  put(paste0("Prisma", w, "Included"), int(pcount(ma, "Studies included in review", "all")), src = "prisma_flow.csv")
}
# --- FAOSTAT dry-bean area harvested (QCL item 176, element 5312), 2020-2023 ---
fao <- fread("paper/lit/faostat/QCL_beans_dry_area_harvested_extract.csv")[`Element Code` == 5312 & Year %in% 2020:2023]
six <- c(226, 114, 215, 184, 29, 238)   # Uganda, Kenya, Tanzania, Rwanda, Burundi, Ethiopia
fao_six <- fao[`Area Code` %in% six, .(ha = sum(Value)), by = Year]$ha / 1e6
fao_ea  <- fao[`Area Code` == 5101]$Value / 1e6   # FAOSTAT region "Eastern Africa"
put("FaoYearFrom", "2020", src = "FAOSTAT QCL extract"); put("FaoYearTo", "2023", src = "FAOSTAT QCL extract")
put("FaoSixMin", min(fao_six), "%.1f", "FAOSTAT QCL extract (six countries)")
put("FaoSixMax", max(fao_six), "%.1f", "FAOSTAT QCL extract (six countries)")
put("FaoEAMin",  min(fao_ea),  "%.1f", "FAOSTAT QCL extract (area 5101)")
put("FaoEAMax",  max(fao_ea),  "%.1f", "FAOSTAT QCL extract (area 5101)")
nov <- fread("paper/lit/search/novelty_rerun/novelty_search_log.csv")
put("NoveltyDate", format(as.Date(unique(nov$date)), "%-d %B %Y"), src = "novelty_search_log.csv")
put("NoveltyQueries", int(uniqueN(nov[grepl("^Crossref API", source)]$query)), src = "novelty_search_log.csv (queries run in both the federated search and Crossref)")
# --- Operating cost per cycle of single national programs (cost_anchors.csv) ---
ca <- fread("paper/lit/extract/cost_anchors.csv")
bud <- range(l1$total_budget)
put("NatCostMin", int(round(min(ca$usd_per_cycle))), src = "cost_anchors.csv")
put("NatCostMax", int(round(max(ca$usd_per_cycle))), src = "cost_anchors.csv")
put("NetProgMin", bud[1] / max(ca$usd_per_cycle), "%.0f", "budget grid / largest national cost")
put("NetProgMax", bud[2] / min(ca$usd_per_cycle), "%.0f", "budget grid / smallest national cost")
put("SearchDate", format(as.Date(unique(sq$date)), "%-d %B %Y"), src = "search_queries.csv")
mac <- read_src("paper/lit/meta_analysis.R")
put("MAtwoMCdraws", int(as.numeric(grab(mac, "\nN <- ([0-9e.+]+)"))), src = "meta_analysis.R")
sinfo <- read_src("paper/lit/extract/synthesis/sessionInfo.txt")
put("RVersionMA",     grab(sinfo, "R version ([0-9.]+)"),   src = "sessionInfo.txt")
put("metaforVersion", sub("-", ".", grab(sinfo, "metafor_([0-9.-]+)")), src = "sessionInfo.txt")
put("RVersionLocal",  paste(R.version$major, R.version$minor, sep = "."), src = "R.version of the session that ran export_numbers.R")
put("RVersionHPC",    grab(read_src("slurm/sweep_v2_array.sbatch"), "module load [^\n]* r/([0-9.]+)"), src = "slurm/sweep_v2_array.sbatch")
hpc_versions <- yaml::read_yaml("config/software_versions_hpc.yaml")
put("AlphaSimRVersionHPC", hpc_versions$AlphaSimR, src = "config/software_versions_hpc.yaml")

# ============================================================================
# Revision additions: rank checks, evidence routes, robustness scenarios, stress block, break-even costs
# ============================================================================
# --- Claims about ranks that the text states are verified here (the export stops if one fails) ---
ev_ind <- ev[at_indifference == TRUE]
meas <- list(sobol = setNames(sob$ST, sob$param),
             shap_ind = setNames(shi$shapley_effect, shi$param),
             shap_cop = setNames(shc$shapley_effect, shc$param),
             pawn = setNames(pawn$PAWN_median, pawn$param),
             evppi = setNames(ev_ind$EVPPI, ev_ind$param))
minor_p <- c("error_reduction", "cost_reduction", "tool_fixed_cost", "total_budget")
for (nm in names(meas)) {
  v <- meas[[nm]]; r <- rank(-v)
  stopifnot(r[["cycle_compression"]] == 1)                         # compression first on every measure
  stopifnot(r[["adoption_ceiling"]] == 4)                          # adoption ceiling fourth on every measure
  stopifnot(max(v[minor_p]) < min(v[setdiff(names(v), minor_p)]))  # four minor inputs below all others
  stopifnot(all(r[c("discount_rate", "bean_price")] %in% c(2, 3))) # discount rate and price share ranks 2 and 3
}
put("EVPPIzeroHurdleMax", money_M(max(ev[abs(hurdle_mult - 0) < 1e-9]$EVPPI, 0)), "%.1f", "evppi.csv (all inputs, hurdle 0)")
stopifnot(max(ev[abs(hurdle_mult - 0) < 1e-9]$EVPPI) < 1e-6)
put("EVPPIratioCycErr", ev_ind[param == "cycle_compression"]$EVPPI / ev_ind[param == "error_reduction"]$EVPPI, "%.0f", "evppi.csv")
put("NegNPVnonzeroComp", int(sum(lk$npv < 0 & lk$cycle_compression > 0)), "%d", "lookup_table.parquet")
put("NgridTech", int(nrow(unique(l1[, .(error_reduction, cost_reduction, tool_fixed_cost, total_budget)]))), src = l1_path)
put("NgridEcon", int(nrow(unique(lk[, .(bean_price, discount_rate, adoption_ceiling)]))), src = "lookup_table.parquet")

# --- Break-even running cost = equivalent annual value x area (level annual cost at which NPV is zero) ---
for (k in 0:3) {
  s_k <- base[cycle_compression == k]
  put(paste0("BECostCyc", num_word[k + 1]), median(s_k$eaa_per_ha_yr) * base$total_area_ha[1] / 1e6, "%.2f", "baseline_slice() EAA x area")
}

# --- Manual-program gain in the robustness block (uncalibrated, % of mean yield per year) ---
rb_gain <- rb[, .(g = 100 * (mean(mean_dg_trad) / breed$cycle_years) / econ$mean_yield_kgha), by = .(h2_yield, gxe_cor_onstation)]
put("RobManualGainPctMin", min(rb_gain$g), "%.1f", rb_path)
put("RobManualGainPctMax", max(rb_gain$g), "%.1f", rb_path)
put("RobManualGainFoldMin", min(rb_gain$g) / econ$gain_target_pct_per_yr, "%.0f", paste0(rb_path, " / gain target"))

# --- Break-even compression: compression at which the break-even running cost equals an annual cost level ---
be_c <- sapply(0:3, function(k) median(base[cycle_compression == k]$eaa_per_ha_yr) * base$total_area_ha[1] / 1e6)
for (nm in c("Low", "Mid", "High")) {
  cst <- econ$tool_running_cost_scenarios_usd_per_yr[match(nm, c("Low", "Mid", "High"))] / 1e6
  ci <- approx(be_c, 0:3, xout = cst)$y
  put(paste0("BECompCost", nm), if (is.na(ci)) "$>$3" else sprintf("%.1f", ci), src = "interpolation of the break-even running cost between integer compression levels, baseline state")
}

# --- Evidence routes for error reduction (meta_results.csv) ---
er_ids <- c(TOne = "MA2.ER.T1", TTwo = "MA2.ER.T2", TTwoV = "MA2.ER.T2v", Htwo = "MA2.ER.h2.median")
for (nm in names(er_ids)) {
  r <- m(er_ids[[nm]])
  put(paste0("ERint", nm, "Lo"), pct(r$ci_lb), "%.0f", paste0("meta_results.csv ", er_ids[[nm]]))
  put(paste0("ERint", nm, "Hi"), pct(r$ci_ub), "%.0f", paste0("meta_results.csv ", er_ids[[nm]]))
  if (nm != "Htwo")
    put(paste0("ERpr", nm, "Pos"), 100 * as.numeric(sub(".*P\\(ER>0\\)=([0-9.]+).*", "\\1", r$notes)), "%.0f", paste0("meta_results.csv ", er_ids[[nm]], " notes"))
}
put("ERtTwoV", pct(m("MA2.ER.T2v")$estimate), "%.0f", "meta_results.csv MA2.ER.T2v (median of Monte Carlo draws)")
put("MAtwoManualRelV",   m("MA2.manualrel.validity")$estimate, "%.2f", "meta_results.csv MA2.manualrel.validity")
put("MAtwoManualRelVLo", m("MA2.manualrel.validity")$ci_lb,    "%.2f", "meta_results.csv MA2.manualrel.validity")
put("MAtwoManualRelVHi", m("MA2.manualrel.validity")$ci_ub,    "%.2f", "meta_results.csv MA2.manualrel.validity")
put("MAtwoManualRelLo",  m("MA2.manualrel")$ci_lb,             "%.2f", "meta_results.csv MA2.manualrel")
put("MAtwoManualRelHi",  m("MA2.manualrel")$ci_ub,             "%.2f", "meta_results.csv MA2.manualrel")
put("MAtwoBreakEvenRel", sqrt(m("MA2.primary")$estimate), "%.2f", "sqrt of pooled correlation (break-even manual validity under T2)")
put("MAtwoManualTruthOnly", m("MA2.manualrel.truthonly")$estimate, "%.2f", "meta_results.csv MA2.manualrel.truthonly")

# --- Robustness scenarios (scenario_summary.csv; 18_revision_scenarios.R) ---
sc <- fread(file.path(out_dir, "scenario_summary.csv"))
sv <- function(id, met) { r <- sc[scenario == id & metric == met]$value; stopifnot(length(r) == 1); r }
scen_key <- c(Base = "Base", PriorZeroOne = "PriorZeroOne", PriorZeroMass = "PriorZeroMass", PriorZeroHeavy = "PriorZeroHeavy",
              PriorSeason = "PriorSeason", SeasonRef = "SeasonRef", RampStep = "RampStep", Recurrent = "Recurrent",
              HorizonForty = "HorizonForty", KTen = "KTen", KFifteen = "KFifteen", Conservative = "Conservative",
              ConservativeZeroOne = "ConservativeZeroOne", TargetQuarter = "TargetQuarter", TargetOne = "TargetOne",
              RunCost0.5 = "CostLow", RunCost1.0 = "CostMid", RunCost2.0 = "CostHigh")
for (id in names(scen_key)) {
  k <- scen_key[[id]]
  put(paste0("RvEV", k),        sv(id, "EV"),    "%.1f", paste0("scenario_summary.csv ", id))
  put(paste0("RvPrPos", k),     sv(id, "PrPos"), "%.1f", paste0("scenario_summary.csv ", id))
  put(paste0("RvSTCyc", k),     sv(id, "ST_cycle_compression"), "%.2f", paste0("scenario_summary.csv ", id))
  put(paste0("RvSTErr", k),     f3(sv(id, "ST_error_reduction")), src = paste0("scenario_summary.csv ", id))
  put(paste0("RvSTDisc", k),    sv(id, "ST_discount_rate"),     "%.2f", paste0("scenario_summary.csv ", id))
  put(paste0("RvEVPPICyc", k),  sv(id, "EVPPIind_cycle_compression"), "%.1f", paste0("scenario_summary.csv ", id))
  put(paste0("RvEVPPIErr", k),  sv(id, "EVPPIind_error_reduction"),   "%.1f", paste0("scenario_summary.csv ", id))
}
for (id in c("Base", "RampStep", "Recurrent", "HorizonForty", "KTen", "KFifteen", "Conservative", "TargetQuarter", "TargetOne")) {
  k <- scen_key[[id]]
  put(paste0("RvNPVCycZero", k), sv(id, "NPVbase_c0"), "%.1f", paste0("scenario_summary.csv ", id))
  put(paste0("RvNPVCycOne", k),  sv(id, "NPVbase_c1"), "%.1f", paste0("scenario_summary.csv ", id))
  put(paste0("RvNPVCycThree", k), sv(id, "NPVbase_c3"), "%.1f", paste0("scenario_summary.csv ", id))
}
# compression keeps the first rank on the total-order index in every scenario; range of the indices across scenarios
sc_ids <- names(scen_key)
st_cyc <- vapply(sc_ids, function(i) sv(i, "ST_cycle_compression"), 0)
st_err <- vapply(sc_ids, function(i) sv(i, "ST_error_reduction"), 0)
for (i in sc_ids) {
  stm <- sc[scenario == i & grepl("^ST_", metric)]
  stopifnot(stm[which.max(value)]$metric == "ST_cycle_compression")
}
put("RvNscenarios", length(sc_ids), "%d", "scenario_summary.csv")
put("RvSTCycMin", min(st_cyc), "%.2f", "scenario_summary.csv"); put("RvSTCycMax", max(st_cyc), "%.2f", "scenario_summary.csv")
put("RvSTErrMin", f3(min(st_err)), src = "scenario_summary.csv"); put("RvSTErrMax", max(st_err), "%.3f", "scenario_summary.csv")
put("RvRatioMin", min(st_cyc / st_err), "%.0f", "scenario_summary.csv (smallest ratio of the two total-order indices)")
put("RvGainScaleTargetQuarter", sv("TargetQuarter", "GainScale"), "%.3f", "scenario_summary.csv")
put("RvGainScaleTargetOne", sv("TargetOne", "GainScale"), "%.2f", "scenario_summary.csv")
for (id in c("RunCost0.5", "RunCost1.0", "RunCost2.0")) {
  k <- scen_key[[id]]
  put(paste0("CfgRun", k), sv(id, "RunCostM"), "%.1f", "economic_params.yaml tool_running_cost_scenarios_usd_per_yr")
  put(paste0("RvEVPPIzeroCyc", k),  sv(id, "EVPPIzero_cycle_compression"), "%.1f", paste0("scenario_summary.csv ", id))
  put(paste0("RvEVPPIzeroDisc", k), max(sv(id, "EVPPIzero_discount_rate"), 0),     "%.1f", paste0("scenario_summary.csv ", id))
  put(paste0("RvEVPPIzeroPrice", k), max(sv(id, "EVPPIzero_bean_price"), 0),       "%.1f", paste0("scenario_summary.csv ", id))
  put(paste0("RvEVPPIzeroErr", k),  max(sv(id, "EVPPIzero_error_reduction"), 0), "%.1f", paste0("scenario_summary.csv ", id))
  put(paste0("RvNPVCycOne", k),     sv(id, "NPVbase_c1"), "%.1f", paste0("scenario_summary.csv ", id))
  put(paste0("RvNPVCycThree", k),   sv(id, "NPVbase_c3"), "%.1f", paste0("scenario_summary.csv ", id))
}
for (e in c("0", "0.5")) {
  put(paste0("RvNPVErr", if (e == "0") "Zero" else "Fifty"), sv("ValueOfAccuracy", paste0("NPV_err_", e)), "%.2f", "scenario_summary.csv ValueOfAccuracy")
}
put("RvNPVCycOneErrZero", sv("ValueOfAccuracy", "NPV_c1_err0"), "%.1f", "scenario_summary.csv ValueOfAccuracy")

scnc <- load_config("revision_scenarios")
put("CfgHorizonLong", scnc$horizon_long_years, "%d", "revision_scenarios.yaml")
put("CfgKLow",  scnc$adoption_rate_k[[1]], "%.2f", "revision_scenarios.yaml")
put("CfgKMid",  scnc$adoption_rate_k[[2]], "%.2f", "revision_scenarios.yaml")
put("CfgTargetLow",  scnc$gain_targets_pct_per_yr[[1]], "%.2f", "revision_scenarios.yaml")
put("CfgTargetHigh", scnc$gain_targets_pct_per_yr[[2]], "%.1f", "revision_scenarios.yaml")
put("CfgZeroMassA", scnc$compression_zero_mass[[1]], "%.1f", "revision_scenarios.yaml")
put("CfgZeroMassB", scnc$compression_zero_mass[[2]], "%.1f", "revision_scenarios.yaml")
put("CfgSeasonStep", scnc$seasonal_time_step_years, "%.1f", "revision_scenarios.yaml")
put("CfgSeasonMid",  scnc$seasonal_compression_years[[2]], "%.1f", "revision_scenarios.yaml")
mi_used <- fread("paper/lit/meta_inputs_used.csv")
put("MAthreeWalterTotal", mi_used[unit == "Walter2019_handheld_total"]$pct_reduction[1], "%.1f", "meta_inputs_used.csv Walter2019_handheld_total (image-analysis time counted)")

# --- Stress block: a tool noisier than manual scoring (sweep_v2_stress) ---
stp <- file.path("layer1_genetic_sim/outputs/sweep_v2_stress/cell_summaries.csv")
st_cells <- fread(stp); st_sum <- fread(file.path(out_dir, "stress_block_summary.csv"))
put("CfgStressErr", pct(unique(st_cells$error_reduction)), "%.0f", stp)
put("NcellsStress", int(nrow(st_cells)), "%d", stp)
put("StressUpliftPct", st_sum$uplift_pct[1], "%.1f", "stress_block_summary.csv")
put("StressNPVCycZero", st_sum[compression == 0]$medNPV, "%.1f", "stress_block_summary.csv")
put("StressPrPosCycZero", st_sum[compression == 0]$PrPos, "%.0f", "stress_block_summary.csv")
put("StressNPVCycOne", st_sum[compression == 1]$medNPV, "%.1f", "stress_block_summary.csv")
put("StressPrPosCycOne", st_sum[compression == 1]$PrPos, "%.0f", "stress_block_summary.csv")
put("StressDiffMin", min(st_cells$mean_dg_diff / st_cells$mean_dg_trad) * 100, "%.1f", stp)
put("StressDiffMax", max(st_cells$mean_dg_diff / st_cells$mean_dg_trad) * 100, "%.1f", stp)

# --- Write ------------------------------------------------------------------
tab <- rbindlist(reg$rows)
dir.create(dirname(tex_path), showWarnings = FALSE, recursive = TRUE)
# TeX form: a leading hyphen becomes a minus sign, and \xspace keeps the space after the macro
# (the manuscript preamble loads xspace and registers \% as an exception)
val_tex <- sub("^-", "\\\\ensuremath{-}", tab$value)
writeLines(c("% AUTO-GENERATED by paper/R/export_numbers.R via `make numbers`. Do not edit.",
             sprintf("\\newcommand{\\Num%s}{%s\\xspace}", tab$name, val_tex)),
           tex_path)
fwrite(tab, csv_path)
cat(sprintf("Wrote %d numbers → %s, %s\n", nrow(tab), tex_path, csv_path))
