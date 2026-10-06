# ============================================================================
# Layer 2 — Visualization
# ============================================================================

suppressPackageStartupMessages(library(ggplot2))

# --- Tornado diagram ---

plot_tornado <- function(tornado_result, output_path = NULL) {
  dt <- tornado_result$tornado
  base_npv <- tornado_result$base_npv

  # Reorder by range
  dt$param <- factor(dt$param, levels = rev(dt$param))

  p <- ggplot(dt) +
    geom_segment(aes(x = npv_lo, xend = npv_hi, y = param, yend = param),
                  linewidth = 6, color = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = base_npv, linetype = "dashed", color = "red") +
    labs(title = "Tornado Diagram: NPV Sensitivity",
         x = "NPV (USD)", y = NULL,
         caption = paste0("Base NPV = $", formatC(base_npv, format = "f",
                                                     big.mark = ",", digits = 0))) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  if (!is.null(output_path)) {
    ggsave(output_path, p, width = 10, height = 6, dpi = 150)
    cat("Tornado plot saved:", output_path, "\n")
  }
  p
}


# --- Genetic gain trajectory comparison ---

plot_gain_trajectories <- function(mc_results, output_path = NULL) {
  dt <- copy(mc_results)
  dt[, annual_dg := delta_g / cycle_years]

  summary_dt <- dt[, .(
    mean_dg = mean(annual_dg),
    lo = quantile(annual_dg, 0.025),
    hi = quantile(annual_dg, 0.975)
  ), by = scenario]

  p <- ggplot(summary_dt, aes(x = scenario, y = mean_dg, fill = scenario)) +
    geom_col(alpha = 0.7, width = 0.5) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
    labs(title = "Annual Genetic Gain: Traditional vs tool",
         x = NULL, y = "Annual Genetic Gain (kg/ha/yr)") +
    scale_fill_manual(values = c("traditional" = "grey60", "tool" = "steelblue")) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))

  if (!is.null(output_path)) {
    ggsave(output_path, p, width = 7, height = 5, dpi = 150)
    cat("Gain trajectory plot saved:", output_path, "\n")
  }
  p
}


# --- NPV distribution histogram ---

plot_npv_distribution <- function(npv_values, output_path = NULL) {
  dt <- data.table(npv = npv_values)
  prob_pos <- mean(npv_values > 0)

  p <- ggplot(dt, aes(x = npv)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    labs(title = "NPV Distribution (Monte Carlo Sensitivity)",
         x = "NPV (USD)", y = "Count",
         caption = sprintf("P(NPV > 0) = %.1f%%", prob_pos * 100)) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  if (!is.null(output_path)) {
    ggsave(output_path, p, width = 8, height = 5, dpi = 150)
    cat("NPV distribution plot saved:", output_path, "\n")
  }
  p
}


# --- 2D response heatmap (from lookup table or break-even marginal) ---

plot_response_heatmap <- function(marginal_dt,
                                    x_name = "x_val",
                                    y_name = "y_val",
                                    fill_name = "npv",
                                    title = "NPV response surface",
                                    x_label = NULL,
                                    y_label = NULL,
                                    contour_dt = NULL,
                                    output_path = NULL,
                                    fill_label = "NPV (M USD)",
                                    fmt = "%.0f",
                                    rescale = NULL) {
  if (is.null(x_label)) x_label <- x_name
  if (is.null(y_label)) y_label <- y_name
  dt <- copy(marginal_dt)
  # Default rescale: NPV → millions. If fill is already per-ha/yr, pass rescale=1.
  if (is.null(rescale)) {
    rescale <- if (grepl("^npv", fill_name)) 1e6 else 1
  }
  dt$npv_M <- dt[[fill_name]] / rescale
  # Discretize: gridded data plots much better with factor axes
  x_levels <- sort(unique(dt[[x_name]]))
  y_levels <- sort(unique(dt[[y_name]]))
  dt$x_fct <- factor(dt[[x_name]], levels = x_levels)
  dt$y_fct <- factor(dt[[y_name]], levels = y_levels)

  p <- ggplot(dt, aes(x = x_fct, y = y_fct, fill = npv_M)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf(fmt, npv_M)), size = 3.2, color = "black") +
    scale_fill_gradient2(low = "firebrick", mid = "white", high = "steelblue",
                          midpoint = 0, name = fill_label) +
    labs(title = title, x = x_label, y = y_label) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid = element_blank())

  if (!is.null(contour_dt) && nrow(contour_dt) > 0) {
    # Contour lines are in the original continuous coord space; map them
    # to the factor positions via the rank transform.
    cd <- copy(contour_dt)
    cd[, x_fct := factor(findInterval(x, c(-Inf, (x_levels[-1] + x_levels[-length(x_levels)]) / 2, Inf)),
                          levels = seq_along(x_levels), labels = as.character(x_levels))]
    cd[, y_fct := factor(findInterval(y, c(-Inf, (y_levels[-1] + y_levels[-length(y_levels)]) / 2, Inf)),
                          levels = seq_along(y_levels), labels = as.character(y_levels))]
    p <- p + geom_path(data = cd,
                        aes(x = x_fct, y = y_fct, group = segment),
                        inherit.aes = FALSE,
                        color = "black", linewidth = 0.8, linetype = "dashed")
  }

  if (!is.null(output_path)) {
    ggsave(output_path, p, width = 8, height = 6, dpi = 150)
    cat("Response heatmap saved:", output_path, "\n")
  }
  p
}


# --- Break-even contour overlay ---

plot_break_even_contour <- function(marginal_dt, contour_dt = NULL,
                                     title = "Break-even surface",
                                     output_path = NULL) {
  plot_response_heatmap(
    marginal_dt,
    title       = title,
    contour_dt  = contour_dt,
    output_path = output_path
  )
}


# --- Sobol indices bar chart ---

plot_sobol_bars <- function(sobol_result, output_path = NULL) {
  dt <- if (is.list(sobol_result) && !is.null(sobol_result$indices)) sobol_result$indices else sobol_result
  if (!is.data.table(dt)) dt <- as.data.table(dt)

  long <- rbind(
    dt[, .(param, group, idx = "S1 (first-order)", val = pmax(S1, 0))],
    dt[, .(param, group, idx = "ST (total-order)", val = pmax(ST, 0))]
  )
  long[, param := factor(param, levels = dt[order(ST), param])]

  p <- ggplot(long, aes(x = val, y = param, fill = idx)) +
    geom_col(position = "dodge", alpha = 0.85) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = c("S1 (first-order)" = "steelblue",
                                   "ST (total-order)" = "firebrick"),
                       name = NULL) +
    labs(title = "Global Sobol sensitivity indices",
         x = "Variance contribution", y = NULL,
         caption = "Values clipped at 0 (small-sample noise below zero not shown)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "top",
          strip.text.y = element_text(angle = 0, face = "bold"))

  if (!is.null(output_path)) {
    ggsave(output_path, p, width = 9, height = 6, dpi = 150)
    cat("Sobol bar chart saved:", output_path, "\n")
  }
  p
}


# --- Profitability surface over (price, discount_rate) ---
# Shows median NPV (in millions) on a log scale across the economic grid,
# with text annotations for the fraction of L1 cells with NPV>0.
# At tool baseline parameters the sign of NPV is dominated by the L1 grid
# so frac_positive is often nearly constant; the log-NPV magnitude is the
# informative signal in this slice.

plot_profit_surface <- function(profit_frac, output_path = NULL) {
  area_base <- median(unique(profit_frac$total_area_ha))
  A_base    <- median(unique(profit_frac$adoption_ceiling))
  dt <- profit_frac[total_area_ha == area_base & adoption_ceiling == A_base]
  dt[, log_npv := log10(pmax(median_npv, 1))]
  dt[, npv_M   := median_npv / 1e6]

  dt[, price_fct := factor(bean_price, levels = sort(unique(bean_price)))]
  dt[, r_fct     := factor(discount_rate, levels = sort(unique(discount_rate)))]

  p <- ggplot(dt, aes(x = price_fct, y = r_fct, fill = log_npv)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("$%.0fM\n%.0f%%",
                                     npv_M, frac_positive * 100)),
               size = 3) +
    scale_fill_gradient(low = "khaki", high = "darkgreen",
                         name = "log10(median NPV)") +
    labs(title = "Tool value surface over market conditions",
         subtitle = sprintf("Adoption = %.2f, Area = %.0f ha  |  numbers: median NPV (M USD) / cells with NPV>0",
                             A_base, area_base),
         x = "Bean price (USD/kg)", y = "Discount rate") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid = element_blank())

  if (!is.null(output_path)) {
    ggsave(output_path, p, width = 9, height = 6, dpi = 150)
    cat("Profit surface saved:", output_path, "\n")
  }
  p
}


# --- Adoption curve comparison ---

plot_adoption_curves <- function(horizon = 25, econ_cfg = NULL, output_path = NULL) {
  if (is.null(econ_cfg)) econ_cfg <- load_config("economic_params")

  A_max <- econ_cfg$adoption_ceiling
  k     <- econ_cfg$diffusion_rate_k
  t_mid <- econ_cfg$diffusion_midpoint_years

  # Traditional releases at cycle_years (10); tool at cycle_years - compression (9)
  cycle_trad <- 10
  cycle_tool  <- cycle_trad - econ_cfg$tool_cycle_compression

  dt <- compare_adoption(horizon, A_max, k, t_mid, cycle_trad, cycle_tool)

  p <- ggplot(dt, aes(x = year, y = adoption, color = scenario)) +
    geom_line(linewidth = 1.2) +
    labs(title = "Variety Adoption Trajectories",
         x = "Year", y = "Adoption Rate (fraction of area)",
         color = "Scenario") +
    scale_color_manual(values = c("traditional" = "grey60", "tool" = "steelblue")) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  if (!is.null(output_path)) {
    ggsave(output_path, p, width = 8, height = 5, dpi = 150)
    cat("Adoption curve plot saved:", output_path, "\n")
  }
  p
}
