# ============================================================================
# Layer 1 — Breeding Pipeline Simulation
# ============================================================================
# Simulates a complete common bean breeding cycle from crosses to release.

source(file.path(getwd(), "layer1_genetic_sim/R/utils.R"))

# --- Run a single trial stage (phenotype + select) ---

run_trial_stage <- function(pop, SP, stage_cfg, h2_base, gxe_cor, var_ref = NULL) {
  n_reps <- stage_cfg$reps
  n_envs <- stage_cfg$envs

  # Effective heritability on genotype-mean basis (already accounts for
  # reps, envs and GxE — phenotype_and_select must NOT average reps again).
  # var_ref anchors the error variance to the F5 pool entering ON (see
  # phenotype_and_select), so h2_eff is a heritability on the F5-line basis.
  h2_eff <- calc_h2_eff(h2_base, n_reps, n_envs, gxe_cor)

  # Determine how many to evaluate (cap at available pop)
  n_eval <- min(stage_cfg$n_genotypes, pop@nInd)
  if (pop@nInd > n_eval) {
    idx <- sample(pop@nInd, n_eval)
    pop <- pop[idx]
  }

  # Determine how many to select
  n_select <- if (!is.null(stage_cfg$n_select)) {
    min(stage_cfg$n_select, pop@nInd)
  } else {
    max(1, round(pop@nInd * stage_cfg$selection_intensity))
  }

  # Phenotype and select
  phenotype_and_select(pop, SP, trait = 1, n_select = n_select,
                        h2_eff = h2_eff, var_ref = var_ref)
}


# --- Full breeding pipeline ---

run_breeding_pipeline <- function(basePop, SP, breed_cfg, scenario_cfg,
                                   gen_cfg, budget_alloc = NULL) {
  h2_base <- gen_cfg$h2_yield
  gxe_cor <- gen_cfg$gxe_cor_onstation

  # Apply tool error reduction on the single-plot basis (sigma2_e scaled by
  # 1 - error_reduction). This happens BEFORE calc_h2_eff() in run_trial_stage,
  # so reps/envs averaging acts on the already-reduced plot error. h2_base is
  # the post-reduction PLOT heritability (reported as h2_plot), NOT h2_eff.
  h2_base <- apply_error_reduction(h2_base, scenario_cfg$error_reduction)

  cycle_years <- breed_cfg$cycle_years - scenario_cfg$cycle_compression

  # --- Step 1: Make crosses ---
  n_crosses <- breed_cfg$n_crosses_per_cycle
  n_per_cross <- breed_cfg$n_progeny_per_cross
  n_parents <- basePop@nInd
  if (n_parents < 2) stop("Need at least 2 parents for crossing")

  cross_plan <- matrix(NA, nrow = n_crosses, ncol = 2)
  for (i in seq_len(n_crosses)) {
    cross_plan[i, ] <- sample(n_parents, 2, replace = FALSE)
  }

  F1 <- makeCross(basePop, crossPlan = cross_plan, nProgeny = n_per_cross,
                   simParam = SP)

  # --- Step 2: Generation advance via selfing (SSD) ---
  pop <- F1
  for (gen in seq_len(breed_cfg$n_ssd_generations)) {
    pop <- self(pop, nProgeny = 1, simParam = SP)
  }

  # --- Step 3: Trial stages ---
  # Reference variance for plot error: genetic variance of the F5 pool entering
  # ON, computed once. Every stage uses varE = var_ref * (1/h2_eff - 1), so plot
  # error stays physically constant while selection depletes the variance. The
  # F5 pool is identical in both scenarios (common random numbers), and
  # error_reduction scales the plot error by exactly (1 - error_reduction)
  # through h2_base above.
  var_ref <- var(gv(pop)[, 1])

  # Use budget_alloc if provided (from optimizer), else fall back to breed_cfg
  stage_names <- c("ON", "PYT", "AYT", "NPT")
  for (st_name in stage_names) {
    if (!is.null(budget_alloc) && !is.null(budget_alloc[[st_name]])) {
      st_cfg <- budget_alloc[[st_name]]
    } else {
      st_cfg <- breed_cfg$stages[[st_name]]
    }
    pop <- run_trial_stage(pop, SP, st_cfg, h2_base, gxe_cor, var_ref = var_ref)
  }

  # --- Results ---
  # Gain is measured against the CHECK varieties farmers already grow: the
  # mean true genetic value of the top-k founders, k = number of released
  # lines (gain_baseline = check_topk, default). The legacy baseline, the
  # founder MEAN (gain_baseline = founder_mean), is always reported too as
  # delta_g_vs_base so v1 and v2 results stay comparable.
  released_gv <- mean(gv(pop))
  base_gv     <- mean(gv(basePop))
  n_released  <- pop@nInd
  founder_gv  <- sort(gv(basePop)[, 1], decreasing = TRUE)
  check_gv    <- mean(founder_gv[seq_len(min(n_released, length(founder_gv)))])

  gain_baseline <- if (is.null(gen_cfg$gain_baseline)) "check_topk" else gen_cfg$gain_baseline
  ref_gv <- switch(gain_baseline,
    check_topk   = check_gv,
    founder_mean = base_gv,
    stop("Unknown gain_baseline: ", gain_baseline,
         " (use 'check_topk' or 'founder_mean')")
  )

  list(
    delta_g         = released_gv - ref_gv,      # primary (vs configured baseline)
    delta_g_vs_base = released_gv - base_gv,     # legacy v1 definition
    released_gv     = released_gv,
    base_gv         = base_gv,
    check_gv        = check_gv,
    n_released      = n_released,
    cycle_years     = cycle_years,
    h2_plot         = h2_base                    # post-error-reduction plot h2
  )
}
