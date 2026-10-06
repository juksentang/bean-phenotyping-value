# ============================================================================
# Layer 1 — Founder Population Initialization
# ============================================================================
# Creates the base population and SimParam for common bean simulation.

source(file.path(getwd(), "layer1_genetic_sim/R/utils.R"))

create_founder_population <- function(gen_cfg = NULL) {
  if (is.null(gen_cfg)) gen_cfg <- load_config("genetic_params")

  # --- Generate founder haplotypes via coalescent (MaCS) ---
  founderPop <- runMacs(
    nInd     = gen_cfg$n_founders,
    nChr     = gen_cfg$n_chr,
    segSites = gen_cfg$seg_sites_per_chr,
    species  = "GENERIC"
  )

  # --- Define simulation parameters ---
  SP <- SimParam$new(founderPop)

  # Single-threaded AlphaSimR: with nThreads > 1, makeCross()/self() are NOT
  # reproducible under set.seed(), which would defeat the common-random-numbers
  # re-seeding in run_single_rep(). Parallelism comes from mclapply across reps.
  SP$nThreads <- 1L

  # Add yield trait (single additive trait for base model)
  SP$addTraitA(
    nQtlPerChr = gen_cfg$n_qtl_per_chr,
    mean        = gen_cfg$mean_yield_kgha,
    var         = gen_cfg$var_yield_additive
  )

  # Set default heritability for phenotyping. Note setVarE(h2 = ) is
  # FOUNDER-anchored (varE = varA_founder/h2 - varG_founder); the breeding
  # pipeline does not rely on it and sets varE explicitly per stage (see
  # phenotype_and_select in utils.R).
  SP$setVarE(h2 = gen_cfg$h2_yield)

  # Create base population
  basePop <- newPop(founderPop, simParam = SP)

  list(pop = basePop, SP = SP, gen_cfg = gen_cfg)
}


# --- Multi-environment founder (with GxE) ---

create_founder_population_gxe <- function(gen_cfg = NULL, n_envs = NULL,
                                           env_type = "onstation") {
 if (is.null(gen_cfg)) gen_cfg <- load_config("genetic_params")
  if (is.null(n_envs)) n_envs <- gen_cfg$n_env_default

  # Choose genetic correlation based on environment type
  gxe_cor <- if (env_type == "onfarm") {
    gen_cfg$gxe_cor_onfarm
  } else {
    gen_cfg$gxe_cor_onstation
  }

  # Build genetic correlation matrix (compound symmetry)
  corA <- matrix(gxe_cor, nrow = n_envs, ncol = n_envs)
  diag(corA) <- 1.0

  # Mean and variance vectors (same trait in each environment)
  mean_vec <- rep(gen_cfg$mean_yield_kgha, n_envs)
  var_vec  <- rep(gen_cfg$var_yield_additive, n_envs)

  # Generate founder haplotypes
  founderPop <- runMacs(
    nInd     = gen_cfg$n_founders,
    nChr     = gen_cfg$n_chr,
    segSites = gen_cfg$seg_sites_per_chr,
    species  = "GENERIC"
  )

  SP <- SimParam$new(founderPop)
  SP$nThreads <- 1L   # reproducible under set.seed() (see create_founder_population)

  # Add correlated traits (one per environment)
  SP$addTraitA(
    nQtlPerChr = gen_cfg$n_qtl_per_chr,
    mean        = mean_vec,
    var         = var_vec,
    corA        = corA
  )

  # Set heritability (per-environment)
  h2_vec <- rep(gen_cfg$h2_yield, n_envs)
  SP$setVarE(h2 = h2_vec)

  basePop <- newPop(founderPop, simParam = SP)

  list(pop = basePop, SP = SP, gen_cfg = gen_cfg,
       n_envs = n_envs, gxe_cor = gxe_cor)
}
