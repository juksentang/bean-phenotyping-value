# ============================================================================
# Empirical calibration — shared variance-component helpers
# ============================================================================
# Definitions used throughout (plot level, single trait = grain yield):
#   Vg   genotype variance
#   Vge  genotype x environment interaction variance
#   Ve   plot residual variance (mean over environments if heterogeneous)
#   h2_w     = Vg / (Vg + Ve)            within-environment plot heritability
#   kappa    = Vge / Vg                  GxE-to-genetic variance ratio
#   rho_env  = Vg / (Vg + Vge)           between-environment genetic correlation
#                                        under compound symmetry
#   h2_total = Vg / (Vg + Vge + Ve)
#   H2_mean  = Vg / (Vg + Vge / e + Ve / (e r))   entry-mean heritability
#              with e, r the harmonic means of environments per genotype and
#              plots per genotype-environment cell of the data set analyzed.

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
})
HAS_GLMMTMB <- requireNamespace("glmmTMB", quietly = TRUE)

derived <- function(Vg, Vge, Ve) {
  c(h2_w     = Vg / (Vg + Ve),
    kappa    = if (Vg > 0) Vge / Vg else NA_real_,
    rho_env  = Vg / (Vg + Vge),
    h2_total = Vg / (Vg + Vge + Ve))
}

has_blocks <- function(d) d[, uniqueN(na.omit(block)) > 1]

# --- Single trial -------------------------------------------------------------
# y ~ rep + (1|genotype) [+ (1|rep:block)]
fit_trial <- function(d) {
  d <- copy(d)[, `:=`(rep = factor(rep), genotype = factor(genotype))]
  blk <- has_blocks(d) && d[, all(!is.na(block))]
  if (blk) d[, repblock := factor(paste(rep, block))]
  f <- if (blk) y ~ rep + (1 | genotype) + (1 | repblock) else y ~ rep + (1 | genotype)
  if (nlevels(d$rep) < 2) f <- update(f, . ~ . - rep)
  fit <- suppressMessages(lmer(f, data = d, REML = TRUE))
  vc  <- as.data.table(VarCorr(fit))
  Vg  <- vc[grp == "genotype", vcov]
  Ve  <- vc[grp == "Residual", vcov]
  Vb  <- if (blk) vc[grp == "repblock", vcov] else NA_real_
  reps_per_geno <- d[, .N, by = genotype]$N
  r_h <- 1 / mean(1 / reps_per_geno)
  data.table(n_plots = nrow(d), n_geno = nlevels(d$genotype), n_rep = nlevels(d$rep),
             n_block = if (blk) uniqueN(d$repblock) else 0L, blocks_fitted = blk,
             mean_kg_ha = mean(d$y), Vg = Vg, Ve = Ve, Vblock = Vb,
             h2_w = Vg / (Vg + Ve), H2_mean = Vg / (Vg + Ve / r_h),
             cv_e_pct = 100 * sqrt(Ve) / mean(d$y),
             vg_boundary = Vg < 1e-6 * Ve, singular = isSingular(fit))
}

# --- Connected sets -----------------------------------------------------------
# Trials are linked when they share at least `min_shared` genotypes and at least
# `min_frac` of the genotypes of the smaller trial; sets are connected components.
connected_sets <- function(d, min_shared = 10, min_frac = 0.5) {
  g  <- split(d$genotype, d$env); g <- lapply(g, unique)
  ev <- names(g); n <- length(ev)
  A  <- matrix(FALSE, n, n, dimnames = list(ev, ev))
  for (i in seq_len(n)) for (j in seq_len(n)) if (i < j) {
    s <- length(intersect(g[[i]], g[[j]]))
    A[i, j] <- A[j, i] <- s >= min_shared && s >= min_frac * min(length(g[[i]]), length(g[[j]]))
  }
  comp <- rep(NA_integer_, n); k <- 0L
  for (i in seq_len(n)) if (is.na(comp[i])) {
    k <- k + 1L; stack <- i
    while (length(stack)) {
      v <- stack[1]; stack <- stack[-1]
      if (!is.na(comp[v])) next
      comp[v] <- k
      stack <- c(stack, which(A[v, ] & is.na(comp)))
    }
  }
  data.table(env = ev, component = comp)
}

set_design <- function(d) {
  ge  <- d[, .N, by = .(genotype, env)]
  e_g <- ge[, .N, by = genotype]$N
  envs <- unique(d$env)
  common <- Reduce(intersect, lapply(split(d$genotype, d$env), unique))
  list(n_env = length(envs), n_geno = uniqueN(d$genotype), n_geno_common = length(common),
       e_h = 1 / mean(1 / e_g), r_h = 1 / mean(1 / ge$N), n_plots = nrow(d))
}

met_formula <- function(blk) {
  if (blk) y ~ env + (1 | genotype) + (1 | genotype:env) + (1 | env:rep) + (1 | env:rep:block)
  else     y ~ env + (1 | genotype) + (1 | genotype:env) + (1 | env:rep)
}

prep_met <- function(d) {
  d <- copy(d)
  d[, blk_ok := !is.na(block) & uniqueN(na.omit(block)) > 1, by = env]
  # Environments without incomplete blocks get a single block per replicate,
  # so the env:rep:block term collapses to env:rep for them.
  d[blk_ok == FALSE, block := 0L]
  d[, `:=`(env = factor(env), genotype = factor(genotype), rep = factor(rep), block = factor(block))]
  d
}

# --- Multi-environment, homogeneous residual (lme4) -----------------------------
fit_met_hom <- function(d) {
  blk <- d[, any(blk_ok)]
  fit <- suppressMessages(lmer(met_formula(blk), data = d, REML = TRUE))
  fit
}

vc_hom <- function(fit) {
  vc <- as.data.table(VarCorr(fit))
  c(Vg = vc[grp == "genotype", vcov], Vge = vc[grp == "genotype:env", vcov],
    Ve = vc[grp == "Residual", vcov])
}

# --- Multi-environment, heterogeneous residual by environment (glmmTMB) -----------
fit_met_het <- function(d, start_fit = NULL) {
  if (!HAS_GLMMTMB) return(NULL)
  blk <- d[, any(blk_ok)]
  # Starting values from the homogeneous lme4 fit (log standard deviations, in
  # the order of the random terms); without them the optimizer can stall at a
  # degenerate solution with all genetic variance at zero.
  start <- NULL
  if (!is.null(start_fit)) {
    vc <- as.data.table(VarCorr(start_fit))
    terms <- c("genotype", "genotype:env", "env:rep", if (blk) "env:rep:block")
    v0 <- vc[match(terms, grp), vcov]
    v0 <- pmax(v0, 1e-3 * vc[grp == "Residual", vcov])
    start <- list(theta = log(sqrt(v0)))
  }
  fit <- tryCatch(suppressWarnings(glmmTMB::glmmTMB(met_formula(blk), dispformula = ~ env,
                                                    data = d, REML = TRUE, start = start)),
                  error = function(e) NULL)
  # A non-positive-definite Hessian is accepted when the optimizer converged:
  # it arises when a variance component sits at zero, and only affects SEs.
  if (is.null(fit) || fit$fit$convergence != 0) return(NULL)
  fit
}

vc_het <- function(fit, d) {
  vc  <- glmmTMB::VarCorr(fit)$cond
  Vg  <- as.numeric(vc[["genotype"]])
  Vge <- as.numeric(vc[["genotype:env"]])
  # glmmTMB's gaussian dispersion is the residual standard deviation
  sd_env <- tapply(predict(fit, type = "disp"), d$env, mean)
  Ve_env <- as.numeric(sd_env)^2
  # Delta-method standard errors of Vg and Vge from the log-SD parameters
  th <- fit$fit$par[names(fit$fit$par) == "theta"]
  se_th <- if (isTRUE(fit$sdr$pdHess)) {
    V <- vcov(fit, full = TRUE); sqrt(diag(V)[grep("^theta", rownames(V))])
  } else rep(NA_real_, length(th))
  se_var <- 2 * exp(2 * th) * se_th
  names(se_var) <- names(vc)
  list(Vg = Vg, Vge = Vge, Ve_env = setNames(Ve_env, names(sd_env)),
       se_Vg = unname(se_var["genotype"]), se_Vge = unname(se_var["genotype:env"]))
}

# --- Parametric bootstrap of the homogeneous model ---------------------------------
boot_met <- function(fit, nsim = 200, seed = 20261018, ncpus = 8) {
  f <- function(m) {
    v <- vc_hom(m)
    c(v, derived(v[["Vg"]], v[["Vge"]], v[["Ve"]]))
  }
  set.seed(seed)
  b <- suppressMessages(suppressWarnings(
    bootMer(fit, f, nsim = nsim, use.u = FALSE, type = "parametric",
            parallel = if (ncpus > 1) "multicore" else "no", ncpus = ncpus)))
  t <- b$t
  t <- t[stats::complete.cases(t[, c("Vg", "Vge", "Ve")]), , drop = FALSE]
  q <- apply(t, 2, quantile, c(0.025, 0.975), na.rm = TRUE)
  s <- apply(t, 2, sd, na.rm = TRUE)
  # kappa = 1/rho_env - 1 is monotone in rho_env, so its percentile interval
  # follows from that of rho_env (and is unbounded above when rho_env can be 0).
  q[, "kappa"] <- rev(1 / q[, "rho_env"] - 1)
  list(lo = q[1, ], hi = q[2, ], sd = s, n_ok = nrow(t))
}

# --- Parametric bootstrap of the heterogeneous-residual model --------------------
boot_het <- function(fit, d, nsim = 200, seed = 20261018, ncpus = 8) {
  set.seed(seed)
  sims <- simulate(fit, nsim = nsim)
  blk  <- d[, any(blk_ok)]
  start <- list(theta = fit$fit$par[names(fit$fit$par) == "theta"],
                betadisp = fit$fit$par[names(fit$fit$par) == "betadisp"])
  one <- function(i) {
    di <- copy(d); di[, y := sims[[i]]]
    fi <- tryCatch(suppressWarnings(glmmTMB::glmmTMB(met_formula(blk), dispformula = ~ env,
                                                     data = di, REML = TRUE, start = start,
                                                     control = glmmTMB::glmmTMBControl(parallel = 1))),
                   error = function(e) NULL)
    if (is.null(fi) || fi$fit$convergence != 0) return(rep(NA_real_, 7))
    v  <- vc_het(fi, di)
    Ve <- mean(v$Ve_env)
    c(Vg = v$Vg, Vge = v$Vge, Ve = Ve, derived(v$Vg, v$Vge, Ve))
  }
  res <- parallel::mclapply(seq_len(nsim), one, mc.cores = ncpus)
  t <- do.call(rbind, res)
  colnames(t) <- c("Vg", "Vge", "Ve", "h2_w", "kappa", "rho_env", "h2_total")
  t <- t[stats::complete.cases(t[, c("Vg", "Vge", "Ve")]), , drop = FALSE]
  q <- apply(t, 2, quantile, c(0.025, 0.975), na.rm = TRUE)
  q[, "kappa"] <- rev(1 / q[, "rho_env"] - 1)
  list(lo = q[1, ], hi = q[2, ], sd = apply(t, 2, sd, na.rm = TRUE), n_ok = nrow(t))
}

# --- Analyze one connected set --------------------------------------------------
analyze_set <- function(d, set_name, nsim = 200, ncpus = 8, boot = TRUE) {
  d   <- prep_met(d)
  des <- set_design(d)
  fh  <- fit_met_hom(d)
  vh  <- vc_hom(fh)
  ft  <- fit_met_het(d, start_fit = fh)
  # The heterogeneous model nests the homogeneous one; a lower restricted
  # likelihood indicates a failed optimization, in which case it is not used.
  ll_gain <- if (!is.null(ft)) as.numeric(logLik(ft)) - as.numeric(logLik(fh)) else NA_real_
  if (!is.null(ft) && ll_gain < -1) ft <- NULL
  if (!is.null(ft)) {
    vt <- vc_het(ft, d)
    Vg <- vt$Vg; Vge <- vt$Vge; Ve <- mean(vt$Ve_env)
    resid_model <- "heterogeneous_by_env"
    pd_hess <- isTRUE(ft$sdr$pdHess)
    se_Vg <- vt$se_Vg; se_Vge <- vt$se_Vge
    ve_rng <- range(vt$Ve_env)
  } else {
    Vg <- vh[["Vg"]]; Vge <- vh[["Vge"]]; Ve <- vh[["Ve"]]
    resid_model <- "homogeneous"; se_Vg <- se_Vge <- NA_real_; pd_hess <- NA
    ve_rng <- c(NA_real_, NA_real_)
  }
  dv <- derived(Vg, Vge, Ve)
  H2 <- Vg / (Vg + Vge / des$e_h + Ve / (des$e_h * des$r_h))
  out <- data.table(set = set_name, n_env = des$n_env, n_geno = des$n_geno,
                    n_geno_common = des$n_geno_common, n_plots = des$n_plots,
                    e_harm = des$e_h, r_harm = des$r_h,
                    blocks_fitted = d[, any(blk_ok)], resid_model = resid_model, pd_hessian = pd_hess,
                    reml_ll_gain_het = ll_gain,
                    Vg = Vg, Vge = Vge, Ve = Ve, Ve_env_min = ve_rng[1], Ve_env_max = ve_rng[2],
                    se_Vg = se_Vg, se_Vge = se_Vge,
                    h2_w = dv[["h2_w"]], kappa = dv[["kappa"]], rho_env = dv[["rho_env"]],
                    h2_total = dv[["h2_total"]], H2_mean = H2,
                    Vg_hom = vh[["Vg"]], Vge_hom = vh[["Vge"]], Ve_hom = vh[["Ve"]],
                    h2_w_hom = derived(vh[["Vg"]], vh[["Vge"]], vh[["Ve"]])[["h2_w"]],
                    kappa_hom = derived(vh[["Vg"]], vh[["Vge"]], vh[["Ve"]])[["kappa"]],
                    singular_hom = isSingular(fh))
  if (boot) {
    b <- if (!is.null(ft)) boot_het(ft, d, nsim = nsim, ncpus = ncpus)
         else boot_met(fh, nsim = nsim, ncpus = ncpus)
    out[, boot_model := resid_model]
    for (k in c("Vg", "Vge", "Ve", "h2_w", "kappa", "rho_env", "h2_total")) {
      out[, paste0(k, "_lo") := b$lo[[k]]]
      out[, paste0(k, "_hi") := b$hi[[k]]]
    }
    out[, n_boot := b$n_ok]
  }
  out
}

# --- Distribution summary ------------------------------------------------------
dist_summary <- function(x) {
  x <- x[is.finite(x)]
  data.table(n = length(x), median = median(x), q25 = quantile(x, 0.25), q75 = quantile(x, 0.75),
             min = min(x), max = max(x), mean = mean(x))
}

fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
fmt_k <- function(x) ifelse(!is.finite(x) | x > 99, "inf", fmt(x))
fmt0 <- function(x) formatC(x, format = "f", digits = 0, big.mark = ",")
