# ============================================================================
# Layer 1 — Breeder's-equation prediction versus the stochastic simulation
# ============================================================================
# For every (h2, kappa, error reduction, cost reduction) of the h2 x GxE robustness
# block, the per-cycle uplift of the tool over manual phenotyping is predicted
# with the multi-stage breeder's equation that drives the allocation search
# (approx_delta_g(), sum over stages of i_s * sqrt(h2_eff,s)), and compared with
# the simulated uplift of the same cells.
#
#   be_fixed  accuracy only: tool error reduction applied to the manual program's
#             optimal allocation (no reallocation, no cost channel)
#   be_opt    full tool profile: both programs re-optimized under the budget
#             (error reduction and per-plot cost reduction, tool fixed cost)
#   sim       simulated paired uplift, averaged over cycle compression, on two
#             bases: the manual gain over the founder mean (vsbase; the basis of
#             the breeder's equation) and over the check varieties (check; the
#             basis used in the economic layer)
#   sim_acc   simulated accuracy effect: sim at the given error reduction minus
#             sim at zero error reduction and the same cost reduction (vsbase)
#
# Usage (project root): Rscript layer1_genetic_sim/R/08_breeders_equation_check.R
# Input block: TOOL_ROBUST_SUBDIR (default sweep_v3_gxe). Output: be_vs_sim.csv
# in that block's directory.

source(file.path(getwd(), "layer1_genetic_sim/R/03_budget_optimizer.R"))
source(file.path(getwd(), "layer1_genetic_sim/R/05_monte_carlo.R"))   # apply_scenario_override()

sub_dir <- Sys.getenv("TOOL_ROBUST_SUBDIR", unset = "sweep_v3_gxe")
blk_dir <- file.path(getwd(), "layer1_genetic_sim/outputs", sub_dir)
rb <- fread(file.path(blk_dir, "cell_summaries.csv"))
stopifnot(uniqueN(rb$tool_fixed_cost) == 1, uniqueN(rb$total_budget) == 1)

master <- load_config("master_config")
gen0   <- load_config("genetic_params")
breed  <- load_config("breeding_params")

# Simulated uplift per (h2, GxE, error, cost), averaged over compression (which
# leaves the per-cycle selection unchanged)
sim <- rb[, .(sim_vsbase = 100 * mean(mean_dg_diff) / mean(mean_dg_vs_base_trad),
              sim_check  = 100 * mean(mean_dg_diff) / mean(mean_dg_trad),
              sim_se     = 100 * sqrt(sum(sd_dg_diff^2 / n_reps)) / .N / mean(mean_dg_vs_base_trad)),
          by = .(h2_yield, gxe_cor_onstation, error_reduction, cost_reduction)]
sim[, sim_acc := sim_vsbase - sim_vsbase[error_reduction == 0], by = .(h2_yield, gxe_cor_onstation, cost_reduction)]

be <- sim[, {
  gen <- gen0
  gen$h2_yield <- h2_yield
  gen$gxe_cor_onstation <- gxe_cor_onstation
  m <- apply_scenario_override(master, list(
    total_budget = rb$total_budget[1], tool_fixed_cost = rb$tool_fixed_cost[1],
    error_reduction = error_reduction, cost_reduction = cost_reduction, cycle_compression = 0))
  trad <- optimize_allocation(m, gen, breed, "traditional")
  tool <- optimize_allocation(m, gen, breed, "tool")
  sigma_a <- sqrt(gen$var_yield_additive)
  dg_fix <- approx_delta_g(trad$best_allocation, apply_error_reduction(h2_yield, error_reduction),
                           sigma_a, gxe_cor_onstation, cycle_years = 1)
  dg_man <- approx_delta_g(trad$best_allocation, h2_yield, sigma_a, gxe_cor_onstation, cycle_years = 1)
  .(be_fixed = 100 * (dg_fix / dg_man - 1),
    be_opt   = 100 * (tool$best_delta_g / trad$best_delta_g - 1),
    pyt_envs_trad = trad$best_allocation$PYT$envs, pyt_envs_tool = tool$best_allocation$PYT$envs)
}, by = .(h2_yield, gxe_cor_onstation, error_reduction, cost_reduction)]

out <- merge(sim, be, by = c("h2_yield", "gxe_cor_onstation", "error_reduction", "cost_reduction"))
out[, kappa := 1 - gxe_cor_onstation]
setcolorder(out, c("h2_yield", "kappa"))
fwrite(out, file.path(blk_dir, "be_vs_sim.csv"))
cat("Wrote", nrow(out), "rows ->", file.path(blk_dir, "be_vs_sim.csv"), "\n")
