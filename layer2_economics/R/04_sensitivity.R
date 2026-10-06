# ============================================================================
# Layer 2 — Sensitivity Analysis
# ============================================================================
# Tornado diagram (one-at-a-time) and Monte Carlo sensitivity.

source(file.path(getwd(), "layer2_economics/R/03_npv_irr.R"))

# --- Tornado: evaluate NPV at low/high for each parameter ---

run_tornado <- function(mc_results, econ_cfg = NULL, master_cfg = NULL) {
  if (is.null(econ_cfg))   econ_cfg   <- load_config("economic_params")
  if (is.null(master_cfg)) master_cfg <- load_config("master_config")
  sens <- master_cfg$sensitivity

  # Base case NPV (using median MC rep)
  base_eval <- evaluate_all_reps(mc_results, econ_cfg, master_cfg)
  base_npv <- base_eval$summary$median_npv

  # Parameters that can be varied in Layer 2 without re-running Layer 1.
  # tool-specific parameters (accuracy, cost reduction, cycle compression)
  # require Layer 1 re-simulation; they are covered by the MC sensitivity
  # across Layer 1 reps instead.
  params <- list(
    list(name = "Discount Rate", field = "discount_rate",
         lo = sens$discount_rate[1], hi = sens$discount_rate[2],
         base = econ_cfg$discount_rate),
    list(name = "Bean Price (USD/kg)", field = "bean_price_usd_per_kg",
         lo = sens$bean_price[1], hi = sens$bean_price[2],
         base = econ_cfg$bean_price_usd_per_kg),
    list(name = "Adoption Ceiling", field = "adoption_ceiling",
         lo = sens$ceiling_adoption[1], hi = sens$ceiling_adoption[2],
         base = econ_cfg$adoption_ceiling),
    list(name = "Total Production Area (ha)", field = "total_area_ha",
         lo = 3000000, hi = 7000000,
         base = econ_cfg$total_area_ha),
    list(name = "Diffusion Rate (k)", field = "diffusion_rate_k",
         lo = 0.15, hi = 0.50,
         base = econ_cfg$diffusion_rate_k),
    list(name = "Diffusion Midpoint (yrs)", field = "diffusion_midpoint_years",
         lo = 8, hi = 20,
         base = econ_cfg$diffusion_midpoint_years)
  )

  tornado_dt <- data.table(
    param = character(), base_val = numeric(),
    lo_val = numeric(), hi_val = numeric(),
    npv_lo = numeric(), npv_hi = numeric(),
    npv_range = numeric()
  )

  for (p in params) {
    for (side in c("lo", "hi")) {
      cfg_mod <- copy(econ_cfg)
      val <- if (side == "lo") p$lo else p$hi
      cfg_mod[[p$field]] <- val

      # For parameters that affect Layer 1 (accuracy, cost, cycle),
      # we can't re-run Layer 1, so we approximate by adjusting the
      # economic model parameters
      ev <- evaluate_all_reps(mc_results, cfg_mod, master_cfg)
      if (side == "lo") npv_lo <- ev$summary$median_npv
      else npv_hi <- ev$summary$median_npv
    }

    tornado_dt <- rbindlist(list(tornado_dt, data.table(
      param = p$name, base_val = p$base,
      lo_val = p$lo, hi_val = p$hi,
      npv_lo = npv_lo, npv_hi = npv_hi,
      npv_range = abs(npv_hi - npv_lo)
    )))
  }

  # Sort by impact range
  tornado_dt <- tornado_dt[order(-npv_range)]

  list(tornado = tornado_dt, base_npv = base_npv)
}


# --- Monte Carlo sensitivity (joint parameter sampling) ---

run_mc_sensitivity <- function(mc_results, n_draws = 10000,
                                econ_cfg = NULL, master_cfg = NULL) {
  if (is.null(econ_cfg))   econ_cfg   <- load_config("economic_params")
  if (is.null(master_cfg)) master_cfg <- load_config("master_config")
  sens <- master_cfg$sensitivity

  set.seed(12345)

  # Use a single median MC rep for speed
  med_rep <- mc_results[rep == median(unique(mc_results$rep))]

  results <- numeric(n_draws)
  param_draws <- data.table(
    draw = integer(),
    discount_rate = numeric(),
    bean_price = numeric(),
    adoption_ceiling = numeric(),
    tool_fixed_cost = numeric()
  )

  for (d in seq_len(n_draws)) {
    cfg_d <- copy(econ_cfg)

    # Sample each parameter from triangular/uniform distribution
    cfg_d$discount_rate <- runif(1, sens$discount_rate[1], sens$discount_rate[2])
    cfg_d$bean_price_usd_per_kg <- runif(1, sens$bean_price[1], sens$bean_price[2])
    cfg_d$adoption_ceiling <- runif(1, sens$ceiling_adoption[1], sens$ceiling_adoption[2])
    cfg_d$tool_fixed_cost <- runif(1, sens$tool_fixed_cost[1], sens$tool_fixed_cost[2])

    ev <- evaluate_all_reps(med_rep, cfg_d, master_cfg)
    results[d] <- ev$summary$median_npv

    param_draws <- rbindlist(list(param_draws, data.table(
      draw = d,
      discount_rate = cfg_d$discount_rate,
      bean_price = cfg_d$bean_price_usd_per_kg,
      adoption_ceiling = cfg_d$adoption_ceiling,
      tool_fixed_cost = cfg_d$tool_fixed_cost
    )))
  }

  list(
    npv_distribution = results,
    param_draws      = param_draws,
    prob_positive    = mean(results > 0),
    mean_npv         = mean(results),
    sd_npv           = sd(results),
    q025             = quantile(results, 0.025),
    q975             = quantile(results, 0.975)
  )
}
