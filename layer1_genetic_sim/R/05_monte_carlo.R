# ============================================================================
# Layer 1 — Monte Carlo Simulation Framework
# ============================================================================
# Runs N replications of the breeding pipeline for both scenarios,
# collecting genetic gain distributions.

source(file.path(getwd(), "layer1_genetic_sim/R/01_founder_pop.R"))
source(file.path(getwd(), "layer1_genetic_sim/R/02_breeding_pipeline.R"))
source(file.path(getwd(), "layer1_genetic_sim/R/03_budget_optimizer.R"))

# --- Scenario override for parameter sweeps ---
# Patches master_cfg in-place for a single sweep cell. Fields in `override`:
#   total_budget      — USD, replaces master_cfg$total_budget
#   error_reduction   — fraction, tool measurement-error reduction
#   cost_reduction    — fraction, tool per-plot cost reduction (sets tool cost_per_plot
#                       = traditional cost_per_plot * (1 - cost_reduction))
#   cycle_compression — years saved by the tool
#   tool_fixed_cost    — USD fixed cost for the tool scenario
# Any NULL field leaves the corresponding master_cfg value untouched.

apply_scenario_override <- function(master_cfg, override) {
  if (is.null(override)) return(master_cfg)
  if (!is.null(override$total_budget))
    master_cfg$total_budget <- override$total_budget
  if (!is.null(override$error_reduction))
    master_cfg$scenarios$tool$error_reduction <- override$error_reduction
  if (!is.null(override$cycle_compression))
    master_cfg$scenarios$tool$cycle_compression <- override$cycle_compression
  if (!is.null(override$tool_fixed_cost))
    master_cfg$scenarios$tool$fixed_cost <- override$tool_fixed_cost
  if (!is.null(override$cost_reduction)) {
    trad_cost <- master_cfg$scenarios$traditional$cost_per_plot
    master_cfg$scenarios$tool$cost_per_plot <- trad_cost * (1 - override$cost_reduction)
  }
  master_cfg
}


# --- Single MC replication ---

run_single_rep <- function(rep_id, master_cfg, gen_cfg, breed_cfg,
                            precomputed_alloc = NULL) {
  set.seed(master_cfg$mc_seed + rep_id)

  # Create fresh founder population (built ONCE per rep and shared by both
  # scenarios, so they start from the identical basePop)
  founder <- create_founder_population(gen_cfg)
  basePop <- founder$pop
  SP      <- founder$SP

  # Common random numbers (CRN): derive one pipeline seed from this rep's
  # stream (cell seed + rep index) and re-seed with it immediately before EACH
  # scenario. Crosses, SSD and as much downstream sampling as possible are then
  # identical in traditional and tool, so noise cancels in the tool - traditional
  # difference. (Draws diverge once the scenarios' trial sizes / h2 differ.)
  pipeline_seed <- sample.int(.Machine$integer.max, 1L)

  results <- list()
  for (sc_name in names(master_cfg$scenarios)) {
    scenario <- master_cfg$scenarios[[sc_name]]
    set.seed(pipeline_seed)

    # Use precomputed optimal allocation if available
    alloc <- if (!is.null(precomputed_alloc)) precomputed_alloc[[sc_name]] else NULL

    res <- run_breeding_pipeline(
      basePop    = basePop,
      SP         = SP,
      breed_cfg  = breed_cfg,
      scenario_cfg = scenario,
      gen_cfg    = gen_cfg,
      budget_alloc = alloc
    )

    results[[sc_name]] <- data.table(
      rep             = rep_id,
      scenario        = sc_name,
      delta_g         = res$delta_g,           # vs configured baseline (check_topk)
      delta_g_vs_base = res$delta_g_vs_base,   # vs founder mean (v1 definition)
      released_gv     = res$released_gv,
      base_gv         = res$base_gv,
      check_gv        = res$check_gv,
      n_released      = res$n_released,
      cycle_years     = res$cycle_years,
      h2_plot         = res$h2_plot            # post-error-reduction PLOT h2 (not h2_eff)
    )
  }

  rbindlist(results)
}


# --- Full MC run ---

run_monte_carlo <- function(master_cfg = NULL, gen_cfg = NULL, breed_cfg = NULL,
                             n_reps = NULL, use_parallel = TRUE,
                             scenario_override = NULL, verbose = TRUE) {
  if (is.null(master_cfg)) master_cfg <- load_config("master_config")
  if (is.null(gen_cfg))    gen_cfg    <- load_config("genetic_params")
  if (is.null(breed_cfg))  breed_cfg  <- load_config("breeding_params")

  master_cfg <- apply_scenario_override(master_cfg, scenario_override)
  if (is.null(n_reps))     n_reps     <- master_cfg$mc_reps

  if (verbose) {
    cat("=== Layer 1: Monte Carlo Simulation ===\n")
    cat("Scenarios:", paste(names(master_cfg$scenarios), collapse = ", "), "\n")
    cat("MC reps:", n_reps, "\n")
    if (!is.null(scenario_override)) {
      cat("Override: budget=$", master_cfg$total_budget,
          " tool{err=", master_cfg$scenarios$tool$error_reduction,
          ", cost/plot=$", master_cfg$scenarios$tool$cost_per_plot,
          ", cycle-=", master_cfg$scenarios$tool$cycle_compression,
          ", fixed=$", master_cfg$scenarios$tool$fixed_cost, "}\n", sep = "")
    }
  }

  # --- Step 1: Optimize budget allocation (once, using breeder's eq approx) ---
  if (verbose) cat("\nOptimizing resource allocation...\n")
  precomputed_alloc <- list()
  for (sc_name in names(master_cfg$scenarios)) {
    opt <- optimize_allocation(master_cfg, gen_cfg, breed_cfg, sc_name)
    precomputed_alloc[[sc_name]] <- opt$best_allocation
    if (verbose) {
      cat(sprintf("  %s: approx ΔG = %.2f kg/ha/yr, budget used = $%.0f / $%.0f\n",
                  sc_name, opt$best_delta_g, opt$cost_used, master_cfg$total_budget))
      pyt <- opt$best_allocation$PYT
      cat(sprintf("    PYT: %d genotypes x %d envs x %d reps\n",
                  pyt$n_genotypes, pyt$envs, pyt$reps))
    }
  }

  # --- Step 2: Run MC replications ---
  if (verbose) cat("\nRunning", n_reps, "MC replications...\n")
  t_start <- Sys.time()

  if (use_parallel && n_reps > 1) {
    # Respect Slurm cgroup limit first (parallel::detectCores reports the
    # full host's cores inside a cgroup, leading to massive over-forking).
    slurm_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA_character_)))
    if (!is.na(slurm_cpus) && slurm_cpus >= 1) {
      n_cores <- slurm_cpus
    } else {
      n_cores <- max(1, parallel::detectCores() - 1)
    }
    if (verbose) cat("Using", n_cores, "cores (SLURM_CPUS_PER_TASK=",
                     Sys.getenv("SLURM_CPUS_PER_TASK"), ")\n", sep="")
    all_results <- parallel::mclapply(
      seq_len(n_reps),
      function(i) {
        run_single_rep(i, master_cfg, gen_cfg, breed_cfg, precomputed_alloc)
      },
      mc.cores = n_cores
    )
    all_results <- rbindlist(all_results)
  } else {
    all_results <- rbindlist(lapply(seq_len(n_reps), function(i) {
      if (verbose && (i %% 10 == 0 || i == 1)) cat("  rep", i, "/", n_reps, "\n")
      run_single_rep(i, master_cfg, gen_cfg, breed_cfg, precomputed_alloc)
    }))
  }

  elapsed <- difftime(Sys.time(), t_start, units = "mins")
  if (verbose) cat(sprintf("\nCompleted in %.1f minutes\n", as.numeric(elapsed)))

  # --- Summary ---
  summary_dt <- all_results[, .(
    mean_dg    = mean(delta_g),
    mean_dg_vs_base = mean(delta_g_vs_base),
    sd_dg      = sd(delta_g),
    median_dg  = median(delta_g),
    q025_dg    = quantile(delta_g, 0.025),
    q975_dg    = quantile(delta_g, 0.975),
    mean_cycle = mean(cycle_years)
  ), by = scenario]

  if (verbose) {
    cat("\n=== Summary ===\n")
    print(summary_dt)
    if ("traditional" %in% all_results$scenario && "tool" %in% all_results$scenario) {
      paired <- merge(
        all_results[scenario == "traditional", .(rep, dg_trad = delta_g)],
        all_results[scenario == "tool", .(rep, dg_tool = delta_g)],
        by = "rep"
      )
      paired[, dg_diff := dg_tool - dg_trad]
      cat(sprintf("\nTOOL advantage: %.2f +/- %.2f kg/ha (mean +/- SD)\n",
                  mean(paired$dg_diff), sd(paired$dg_diff)))
      cat(sprintf("P(tool > traditional): %.1f%%\n",
                  100 * mean(paired$dg_diff > 0)))
    }
  }

  list(
    results         = all_results,
    summary         = summary_dt,
    precomputed_alloc = precomputed_alloc,
    n_reps          = n_reps,
    elapsed_mins    = as.numeric(elapsed)
  )
}
