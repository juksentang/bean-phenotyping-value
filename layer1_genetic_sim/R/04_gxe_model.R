# ============================================================================
# Layer 1 — GxE Interaction Modeling
# ============================================================================
# Handles mixed on-station / on-farm environment structures.
# Tool allows more on-farm environments (lower cost), improving GxE estimation.

source(file.path(getwd(), "layer1_genetic_sim/R/utils.R"))

# --- Build a heterogeneous correlation matrix ---
# Combines on-station and on-farm environments with different correlations.

build_gxe_cor_matrix <- function(n_onstation, n_onfarm,
                                  cor_onstation = 0.60,
                                  cor_onfarm = 0.40,
                                  cor_cross = NULL) {
  # cor_cross: correlation between on-station and on-farm environments
  # default: geometric mean of the two
  if (is.null(cor_cross)) {
    cor_cross <- sqrt(cor_onstation * cor_onfarm)
  }

  n_total <- n_onstation + n_onfarm
  cormat <- matrix(cor_cross, nrow = n_total, ncol = n_total)

  # On-station block
  if (n_onstation > 1) {
    idx_os <- seq_len(n_onstation)
    cormat[idx_os, idx_os] <- cor_onstation
  }

  # On-farm block
  if (n_onfarm > 1) {
    idx_of <- (n_onstation + 1):n_total
    cormat[idx_of, idx_of] <- cor_onfarm
  }

  diag(cormat) <- 1.0
  cormat
}


# --- Effective heritability for multi-environment BLUP ---
# Uses the average genetic correlation across environments

calc_h2_eff_gxe <- function(h2_base, n_reps, cor_matrix) {
  n_envs <- nrow(cor_matrix)
  # Average genetic correlation (excluding diagonal)
  if (n_envs == 1) return(calc_h2_eff(h2_base, n_reps, 1, 1.0))

  avg_cor <- (sum(cor_matrix) - n_envs) / (n_envs * (n_envs - 1))
  calc_h2_eff(h2_base, n_reps, n_envs, avg_cor)
}


# --- Determine environment allocation per scenario ---

get_env_allocation <- function(scenario_name, gen_cfg, breed_cfg, master_cfg) {
  scenario <- master_cfg$scenarios[[scenario_name]]

  if (scenario_name == "traditional") {
    # Traditional: all on-station
    list(
      n_onstation = breed_cfg$stages$AYT$envs,
      n_onfarm    = 0
    )
  } else {
    # tool: cost savings allow adding on-farm environments
    # With 50% cost reduction, roughly double the environments
    # Split: keep same on-station count, add on-farm
    n_os <- breed_cfg$stages$AYT$envs
    # Additional on-farm envs funded by cost savings
    cost_saving_per_env <- breed_cfg$stages$AYT$n_genotypes *
      breed_cfg$stages$AYT$reps *
      (master_cfg$scenarios$traditional$cost_per_plot - scenario$cost_per_plot)
    available_budget <- cost_saving_per_env * n_os  # savings from cheaper on-station
    cost_per_onfarm_env <- breed_cfg$stages$AYT$n_genotypes *
      breed_cfg$stages$AYT$reps * scenario$cost_per_plot
    n_of <- max(0, floor(available_budget / cost_per_onfarm_env))

    list(
      n_onstation = n_os,
      n_onfarm    = n_of
    )
  }
}
