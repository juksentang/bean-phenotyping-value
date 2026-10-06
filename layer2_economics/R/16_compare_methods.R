# ============================================================================
# Layer 2 — Combined comparison of sensitivity methods
# ============================================================================
# One-stop chart: Sobol ST, Shapley (indep + copula), PAWN (median + max),
# EVPPI (at tau = E[NPV]) — for all 8 parameters.
# Output: report/figures/method_comparison.png

source(file.path(getwd(), "layer2_economics/R/00_paths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

out_dir <- l2_out_dir()
fig_dir <- l2_fig_dir()

sobol   <- fread(file.path(out_dir, "sobol_indices_full.csv"))
sh_i    <- fread(file.path(out_dir, "shapley_indices_indep.csv"))
sh_c    <- fread(file.path(out_dir, "shapley_indices_copula.csv"))
pawn    <- fread(file.path(out_dir, "pawn_indices.csv"))
evppi   <- fread(file.path(out_dir, "evppi.csv"))

# EVPPI at the peak-information hurdle
ev_peak <- evppi[at_indifference == TRUE, .(param, EVPPI = EVPPI / 1e6)]  # tau = E[NPV]

sobol_tidy <- sobol[, .(param, value = ST,             method = "Sobol ST")]
shi_tidy   <- sh_i[, .(param, value = shapley_effect, method = "Shapley (indep)")]
shc_tidy   <- sh_c[, .(param, value = shapley_effect, method = "Shapley (copula)")]
pawn_med   <- pawn[, .(param, value = PAWN_median,    method = "PAWN (median)")]
pawn_max_  <- pawn[, .(param, value = PAWN_max,       method = "PAWN (max)")]

all_dt <- rbindlist(list(sobol_tidy, shi_tidy, shc_tidy, pawn_med, pawn_max_))

# Order params by Sobol ST descending
param_order <- sobol$param[order(-sobol$ST)]
all_dt[, param := factor(param, levels = rev(param_order))]
all_dt[, method := factor(method, levels = c(
  "Sobol ST", "Shapley (indep)", "Shapley (copula)",
  "PAWN (median)", "PAWN (max)"
))]

p_comp <- ggplot(all_dt, aes(x = param, y = value, fill = method)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Sensitivity indices — method comparison",
    subtitle = "All normalised so values are comparable within their own scale",
    x = NULL, y = "sensitivity index"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "method_comparison.png"),
       p_comp, width = 10, height = 6, dpi = 150)

# Separate panel for EVPPI because units are $M, not [0,1]
ev_peak[, param := factor(param, levels = rev(param_order))]
p_ev <- ggplot(ev_peak, aes(x = param, y = EVPPI)) +
  geom_col(fill = "#d95f02") +
  geom_text(aes(label = sprintf("$%.1fM", EVPPI)),
             hjust = -0.05, size = 3) +
  coord_flip() +
  labs(
    title = "EVPPI at the indifference hurdle (tau = E[NPV])",
    subtitle = "Dollar value of perfectly measuring each parameter before deciding",
    x = NULL, y = "EVPPI ($M)"
  ) +
  expand_limits(y = max(ev_peak$EVPPI) * 1.15) +
  theme_minimal(base_size = 12)

ggsave(file.path(fig_dir, "evppi_comparison.png"),
       p_ev, width = 8, height = 5, dpi = 150)

cat("Wrote:\n  ",
    file.path(fig_dir, "method_comparison.png"), "\n  ",
    file.path(fig_dir, "evppi_comparison.png"), "\n")
