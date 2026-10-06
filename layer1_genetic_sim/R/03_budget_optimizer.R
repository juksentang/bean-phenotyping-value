# ============================================================================
# Layer 1 — Budget-Constrained Resource Allocation Optimizer
# ============================================================================
# Grid search over (n_genotypes, n_envs) at PYT stage to maximize ΔG
# under a fixed total budget constraint.
#
# Key insight: downstream stages have FIXED CAPACITY (AYT~50, NPT~15).
# So selection proportion at PYT = n_ayt_fixed / n_pyt, which decreases
# (= stronger selection) as n_pyt grows. This correctly captures the
# benefit of evaluating more genotypes.

source(file.path(getwd(), "layer1_genetic_sim/R/utils.R"))

# --- Calculate total trial cost for a scenario ---

calc_trial_cost <- function(stages, cost_per_plot, fixed_cost = 0) {
  total <- fixed_cost
  for (st_name in names(stages)) {
    st <- stages[[st_name]]
    total <- total + st$n_genotypes * st$envs * st$reps * cost_per_plot
  }
  total
}


# --- Approximate ΔG for a given allocation using breeder's equation ---

approx_delta_g <- function(stages, h2_base, sigma_a, gxe_cor, cycle_years) {
  total_dg <- 0
  for (st_name in names(stages)) {
    st <- stages[[st_name]]
    h2_eff <- calc_h2_eff(h2_base, st$reps, st$envs, gxe_cor)
    sel_prop <- st$n_select / st$n_genotypes
    dg_stage <- predict_delta_g(h2_eff, sigma_a, sel_prop, cycle_years = 1)
    total_dg <- total_dg + dg_stage
  }
  total_dg / cycle_years
}


# --- Grid search optimizer ---

optimize_allocation <- function(master_cfg, gen_cfg, breed_cfg,
                                 scenario_name = "traditional") {
  scenario <- master_cfg$scenarios[[scenario_name]]
  budget <- master_cfg$total_budget
  cost_per_plot <- scenario$cost_per_plot
  fixed_cost <- scenario$fixed_cost
  h2_base <- gen_cfg$h2_yield
  sigma_a <- sqrt(gen_cfg$var_yield_additive)
  gxe_cor <- gen_cfg$gxe_cor_onstation
  cycle_years <- breed_cfg$cycle_years - scenario$cycle_compression

  # Apply error reduction to h2 (same transform as run_breeding_pipeline)
  h2_base <- apply_error_reduction(h2_base, scenario$error_reduction)

  base_stages <- breed_cfg$stages

  # Fixed downstream capacities (number selected INTO the next stage)
  n_select_into_pyt <- base_stages$PYT$n_genotypes   # ~300
  n_select_into_ayt <- base_stages$AYT$n_genotypes   # ~50
  n_select_into_npt <- base_stages$NPT$n_genotypes   # ~15
  n_released        <- max(3, round(base_stages$NPT$n_genotypes *
                                     base_stages$NPT$selection_intensity))  # ~5

  # F5 lines the pipeline can actually evaluate at ON (crosses x progeny per
  # cross). run_trial_stage cannot phenotype more than this, so ON plots (and
  # cost) are capped here; otherwise the budget pays for phantom plots and the
  # ON selection proportion is overstated (n_pyt / n_on vs n_pyt / n_pool).
  n_pool <- breed_cfg$n_crosses_per_cycle * breed_cfg$n_progeny_per_cross

  # Grid: vary PYT n_genotypes and PYT envs
  n_grid <- seq(100, 3000, by = 50)
  e_grid <- seq(1, 20, by = 1)

  best_dg <- -Inf
  best_alloc <- NULL
  results <- data.table(
    n_pyt = integer(), e_pyt = integer(),
    total_cost = numeric(), delta_g = numeric()
  )

  for (n_pyt in n_grid) {
    for (e_pyt in e_grid) {
      # ON must supply n_pyt genotypes; ON evaluates more and selects down
      n_on <- min(n_pyt * 10, 5000, n_pool)  # ON is 10x PYT typically, <= F5 pool
      if (n_pyt > n_on) next                 # pool too small to supply PYT

      # Build stages with fixed downstream sizes and computed selection props
      stages <- list(
        ON  = list(n_genotypes = n_on,  n_select = n_pyt,
                   envs = 1, reps = 1),
        PYT = list(n_genotypes = n_pyt, n_select = n_select_into_ayt,
                   envs = e_pyt, reps = base_stages$PYT$reps),
        AYT = list(n_genotypes = n_select_into_ayt, n_select = n_select_into_npt,
                   envs = base_stages$AYT$envs, reps = base_stages$AYT$reps),
        NPT = list(n_genotypes = n_select_into_npt, n_select = n_released,
                   envs = base_stages$NPT$envs, reps = base_stages$NPT$reps)
      )

      # Check budget
      total_cost <- calc_trial_cost(stages, cost_per_plot, fixed_cost)
      if (total_cost > budget) next

      # Approximate ΔG
      dg <- approx_delta_g(stages, h2_base, sigma_a, gxe_cor, cycle_years)

      results <- rbindlist(list(results, data.table(
        n_pyt = n_pyt, e_pyt = e_pyt,
        total_cost = total_cost, delta_g = dg
      )))

      if (dg > best_dg) {
        best_dg <- dg
        best_alloc <- stages
      }
    }
  }

  list(
    best_allocation = best_alloc,
    best_delta_g    = best_dg,
    scenario        = scenario_name,
    budget          = budget,
    cost_used       = calc_trial_cost(best_alloc, cost_per_plot, fixed_cost),
    grid_results    = results
  )
}
