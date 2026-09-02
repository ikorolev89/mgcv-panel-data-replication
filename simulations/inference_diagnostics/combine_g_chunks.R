#!/usr/bin/env Rscript

# Combine four M=250 chunks into the final M=1000 n=500,T=4 g design.
# Usage:
#   Rscript simulations/inference_diagnostics/combine_g_chunks.R hetero
#   Rscript simulations/inference_diagnostics/combine_g_chunks.R hetero_ar1

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
error_design <- if (length(args)) args[[1]] else stop("Supply an error design.")
if (!error_design %in% c("hetero", "hetero_ar1")) {
  stop("Error design must be hetero or hetero_ar1.")
}

out_dir <- file.path(getwd(), "simulations", "inference_diagnostics", "outputs")
chunk_prefix <- if (error_design == "hetero") "h" else "a"
error_suffix <- if (error_design == "hetero") {
  "_hetero"
} else {
  "_hetero_ar1_rho0p5"
}

chunk_paths <- vapply(seq_len(4L), function(i) {
  file.path(
    out_dir,
    paste0(
      "panel_fe_g_diagnostics_grid_M250_n500_T4_sin_2pi_np_pl_J50",
      error_suffix,
      "_dfsubstantive_chunk", chunk_prefix, i, ".csv.gz"
    )
  )
}, character(1))

missing_paths <- chunk_paths[!file.exists(chunk_paths)]
if (length(missing_paths)) {
  stop("Missing chunk files: ", paste(basename(missing_paths), collapse = ", "))
}

chunk_results <- lapply(seq_along(chunk_paths), function(i) {
  d <- read.csv(gzfile(chunk_paths[i]), check.names = FALSE)
  if (length(unique(d$sim)) != 250L) {
    stop("Chunk does not contain 250 simulations: ", basename(chunk_paths[i]))
  }
  d$sim <- d$sim + 250L * (i - 1L)
  d
})
results <- do.call(rbind, chunk_results)
row.names(results) <- NULL

expected_rows <- 1000L * 2L * 3L * 5L * 50L
if (nrow(results) != expected_rows) {
  stop("Expected ", expected_rows, " grid rows but found ", nrow(results), ".")
}
if (length(unique(results$sim)) != 1000L) stop("Combined simulation IDs are not unique.")
if (!all(is.finite(results$estimate)) || !all(is.finite(results$se))) {
  stop("Combined results contain nonfinite estimates or standard errors.")
}

results$M <- 1000L
results$n <- 500L
results$T <- 4L
results$K <- 20L
results$dgp <- "sin_2pi"
results$error_design <- error_design
results$rho <- if (error_design == "hetero_ar1") 0.5 else NA_real_
results$df_correction <- "substantive"

vcov_types <- c(
  "classic_cluster",
  "penalty_cluster",
  "penalty_cluster_delta",
  "mgcv_sandwich_freq",
  "mgcv_sandwich_default"
)
x_grid <- sort(unique(results$x))
selected_points <- c(0.25, 0.50, 0.75)
selected_grid_idx <- vapply(
  selected_points,
  function(x0) which.min(abs(x_grid - x0)),
  integer(1)
)

by_x_mean <- aggregate(
  cbind(point_covered, estimate, se) ~ model + method + vcov_type + x,
  data = results,
  FUN = mean
)
names(by_x_mean)[names(by_x_mean) == "point_covered"] <- "point_coverage"
names(by_x_mean)[names(by_x_mean) == "estimate"] <- "mean_estimate"
names(by_x_mean)[names(by_x_mean) == "se"] <- "mean_reported_se"

by_x_sd <- aggregate(
  estimate ~ model + method + vcov_type + x,
  data = results,
  FUN = sd
)
names(by_x_sd)[names(by_x_sd) == "estimate"] <- "empirical_sd"

coverage_by_x <- merge(
  by_x_mean,
  by_x_sd,
  by = c("model", "method", "vcov_type", "x")
)
truth_by_x <- unique(results[, c("x", "truth")])
coverage_by_x$truth <- truth_by_x$truth[match(coverage_by_x$x, truth_by_x$x)]
coverage_by_x$bias <- coverage_by_x$mean_estimate - coverage_by_x$truth
coverage_by_x$se_to_sd <- coverage_by_x$mean_reported_se / coverage_by_x$empirical_sd
coverage_by_x$abs_bias <- abs(coverage_by_x$bias)
coverage_by_x$error_design <- error_design
coverage_by_x$rho <- if (error_design == "hetero_ar1") 0.5 else NA_real_
coverage_by_x$df_correction <- "substantive"
coverage_by_x$M <- 1000L
coverage_by_x$n <- 500L
coverage_by_x$T <- 4L
coverage_by_x$K <- 20L
coverage_by_x$dgp <- "sin_2pi"

selected_coverage <- do.call(
  rbind,
  lapply(seq_along(selected_points), function(j) {
    d <- coverage_by_x[abs(coverage_by_x$x - x_grid[selected_grid_idx[j]]) < 1e-12, ]
    d$selected_x <- selected_points[j]
    d[, c("model", "method", "vcov_type", "selected_x", "point_coverage")]
  })
)
selected_wide <- reshape(
  selected_coverage,
  idvar = c("model", "method", "vcov_type"),
  timevar = "selected_x",
  direction = "wide"
)

summary_base <- aggregate(
  cbind(point_covered, point_width, uniform_covered, uniform_width, sup_crit) ~
    model + method + vcov_type,
  data = results,
  FUN = mean
)
names(summary_base)[names(summary_base) == "point_covered"] <- "avg_pointwise_coverage"
names(summary_base)[names(summary_base) == "point_width"] <- "avg_pointwise_width"
names(summary_base)[names(summary_base) == "uniform_covered"] <- "uniform_coverage"
names(summary_base)[names(summary_base) == "uniform_width"] <- "avg_uniform_width"
names(summary_base)[names(summary_base) == "sup_crit"] <- "avg_sup_crit"

min_coverage <- aggregate(
  point_coverage ~ model + method + vcov_type,
  data = coverage_by_x,
  FUN = min
)
names(min_coverage)[names(min_coverage) == "point_coverage"] <- "min_pointwise_coverage"

diagnostic_means <- aggregate(
  cbind(abs_bias, empirical_sd, mean_reported_se, se_to_sd) ~
    model + method + vcov_type,
  data = coverage_by_x,
  FUN = mean
)
max_abs_bias <- aggregate(
  abs_bias ~ model + method + vcov_type,
  data = coverage_by_x,
  FUN = max
)
names(max_abs_bias)[names(max_abs_bias) == "abs_bias"] <- "max_abs_bias"

n_success <- aggregate(
  sim ~ model + method + vcov_type,
  data = unique(results[, c("sim", "model", "method", "vcov_type")]),
  FUN = length
)
names(n_success)[names(n_success) == "sim"] <- "n_success"

summary_results <- Reduce(
  function(x, y) merge(x, y, by = c("model", "method", "vcov_type"), all.x = TRUE),
  list(summary_base, min_coverage, selected_wide, diagnostic_means, max_abs_bias, n_success)
)
names(summary_results) <- sub("point_coverage\\.", "coverage_at_", names(summary_results))
summary_results <- summary_results[
  order(summary_results$model, summary_results$method, match(summary_results$vcov_type, vcov_types)),
]
row.names(summary_results) <- NULL
summary_results$error_design <- error_design
summary_results$rho <- if (error_design == "hetero_ar1") 0.5 else NA_real_
summary_results$df_correction <- "substantive"
summary_results$M <- 1000L
summary_results$n <- 500L
summary_results$T <- 4L
summary_results$K <- 20L
summary_results$dgp <- "sin_2pi"

file_suffix <- paste0(
  "M1000_n500_T4_sin_2pi_np_pl_J50",
  error_suffix,
  "_dfsubstantive"
)
results_path <- file.path(out_dir, paste0("panel_fe_g_diagnostics_grid_", file_suffix, ".csv.gz"))
coverage_x_path <- file.path(out_dir, paste0("panel_fe_g_diagnostics_by_x_", file_suffix, ".csv"))
summary_path <- file.path(out_dir, paste0("panel_fe_g_diagnostics_summary_", file_suffix, ".csv"))
plot_path <- file.path(out_dir, paste0("panel_fe_g_diagnostics_coverage_curve_", file_suffix, ".png"))
se_ratio_plot_path <- file.path(out_dir, paste0("panel_fe_g_diagnostics_se_ratio_", file_suffix, ".png"))

write.csv(results, gzfile(results_path), row.names = FALSE)
write.csv(coverage_by_x, coverage_x_path, row.names = FALSE)
write.csv(summary_results, summary_path, row.names = FALSE)

plot_theme <- theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title.position = "plot",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

coverage_plot <- ggplot(
  coverage_by_x,
  aes(x = x, y = point_coverage, color = vcov_type)
) +
  geom_hline(yintercept = 0.95, color = "grey45", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  facet_grid(model ~ method) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Pointwise coverage for centered g(x)",
    subtitle = paste0("n = 500; T = 4; errors: ", error_design, "; M = 1000"),
    x = "x", y = "pointwise coverage", color = "covariance"
  ) + plot_theme

se_ratio_plot <- ggplot(
  coverage_by_x,
  aes(x = x, y = se_to_sd, color = vcov_type)
) +
  geom_hline(yintercept = 1, color = "grey45", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  facet_grid(model ~ method) +
  labs(
    title = "Reported standard error relative to Monte Carlo standard deviation",
    subtitle = paste0("n = 500; T = 4; errors: ", error_design, "; M = 1000"),
    x = "x", y = "mean reported SE / empirical SD", color = "covariance"
  ) + plot_theme

ggsave(plot_path, coverage_plot, width = 9, height = 6.5, dpi = 180, bg = "white")
ggsave(se_ratio_plot_path, se_ratio_plot, width = 9, height = 6.5, dpi = 180, bg = "white")

cat("Combined", error_design, "chunks into", nrow(results), "grid rows.\n")
cat("Wrote:", summary_path, "\n")
