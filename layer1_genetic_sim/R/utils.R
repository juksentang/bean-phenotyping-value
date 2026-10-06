# ============================================================================
# Layer 1 — Shared Utilities
# ============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(AlphaSimR)
  library(data.table)
})

# --- Config loading ---

load_config <- function(name, config_dir = file.path(getwd(), "config")) {
  path <- file.path(config_dir, paste0(name, ".yaml"))
  if (!file.exists(path)) stop("Config not found: ", path)
  yaml::read_yaml(path)
}

# --- Effective heritability ---
# h2_eff adjusts base heritability for number of reps and environments:
#   h2_eff = sigma2_A / (sigma2_A + sigma2_GE/e + sigma2_e/(r*e))
# where sigma2_e = sigma2_A * (1 - h2) / h2
# and   sigma2_GE is derived from genetic correlation across environments

calc_h2_eff <- function(h2_base, n_reps, n_envs, gxe_cor = 1.0) {
  # Approximate effective heritability on a genotype-mean basis
  # gxe_cor = genetic correlation between environments
  # sigma2_GE / sigma2_A ~ (1 - gxe_cor) for uniformly correlated envs
  ratio_ge <- max(0, 1 - gxe_cor)           # sigma2_GE / sigma2_A
  ratio_e  <- (1 - h2_base) / h2_base       # sigma2_e / sigma2_A (single plot)
  h2_eff <- 1 / (1 + ratio_ge / n_envs + ratio_e / (n_reps * n_envs))
  h2_eff
}

# --- tool measurement-error reduction on the single-plot basis ---
# Normalise sigma2_e = 1 on the single-plot basis, so sigma2_A = h2/(1 - h2).
# Tool shrinks sigma2_e by (1 - error_reduction) BEFORE any reps/envs averaging;
# calc_h2_eff() then divides this reduced sigma2_e by r*e. Returns the
# post-reduction single-plot heritability (h2_plot). Shared by the pipeline and
# the budget optimizer so both apply the identical transform.

apply_error_reduction <- function(h2_base, error_reduction = 0) {
  # error_reduction < 0 is a tool noisier than manual scoring (stress block sweep_v2_stress);
  # values >= 0 behave exactly as in the main sweeps
  if (is.null(error_reduction) || error_reduction == 0) return(h2_base)
  sigma2_a <- h2_base / (1 - h2_base)
  sigma2_e <- 1.0 * (1 - error_reduction)
  sigma2_a / (sigma2_a + sigma2_e)
}

# --- Selection intensity from proportion selected ---

selection_intensity <- function(prop) {
  # i = phi(z_p) / (1 - Phi(z_p))  where z_p = qnorm(1 - prop)
  if (prop >= 1) return(0)
  if (prop <= 0) return(Inf)
  z <- qnorm(1 - prop)
  dnorm(z) / prop
}

# --- Breeder's equation prediction ---

predict_delta_g <- function(h2_eff, sigma_a, sel_prop, cycle_years) {
  i <- selection_intensity(sel_prop)
  h  <- sqrt(h2_eff)
  i * h * sigma_a / cycle_years
}

# --- BLUP-like selection from AlphaSimR ---
# Simulate phenotyping with given h2_eff, then select top fraction.
#
# h2_eff must be the heritability of the genotype MEAN over n_reps x n_envs
# (calc_h2_eff already folds in reps, envs and GxE), so ONE phenotype draw at
# that h2 is the correct selection criterion. Averaging n_reps independent
# draws on top of it would shrink the error a second time (v1 bug: inflated
# accuracy, tool error_reduction nearly irrelevant). `n_reps` is therefore
# unused here and only kept so existing callers do not break.
#
# Variance basis: SP$setVarE(h2 = h2_eff) / setPheno(h2 = ) are FOUNDER-
# anchored (varE = varA_founder/h2 - varG_founder, with varA_founder fixed when
# the trait was added), NOT relative to the pool being phenotyped. F5 inbred
# lines carry ~2x the founder variance and truncation selection shrinks it at
# every later stage, so realized h2 would drift away from h2_eff (inflated at
# ON, deflated at later stages). Instead the error variance is set explicitly
#   varE = var_ref * (1 / h2_eff - 1)
# where var_ref is the genetic variance of the F5 pool entering ON (passed by
# run_breeding_pipeline). Plot error is then physically constant across
# stages (realized h2 falls as selection depletes variance), and h2_eff keeps
# its genotype-mean meaning on the F5-line basis that literature h2 refers to.
# var_ref = NULL falls back to the variance of the pool being phenotyped
# (realized h2 == h2_eff at that stage). Single-trait SP only.

phenotype_and_select <- function(pop, SP, trait = 1, n_select,
                                  h2_eff, n_reps = 1, var_ref = NULL) {
  if (is.null(var_ref)) var_ref <- var(gv(pop)[, trait])
  # Genotype-mean error variance on the var_ref basis
  var_e <- var_ref * (1 / h2_eff - 1)
  # Single phenotype draw with explicit varE (no SP state is touched, no extra
  # averaging over reps)
  pop <- setPheno(pop, varE = var_e, simParam = SP)
  # Select top n_select based on target trait
  n_select <- min(n_select, pop@nInd)
  ord <- order(pheno(pop)[, trait], decreasing = TRUE)
  selected_idx <- ord[seq_len(n_select)]
  pop[selected_idx]
}
