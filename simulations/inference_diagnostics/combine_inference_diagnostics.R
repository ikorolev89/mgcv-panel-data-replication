#!/usr/bin/env Rscript

# Combine completed inference-diagnostic jobs and make comparison plots.
# Usage:
#   Rscript simulations/inference_diagnostics/combine_inference_diagnostics.R 1000

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
M_requested <- if (length(args) >= 1L) as.integer(args[[1]]) else 1000L
if (is.na(M_requested) || M_requested < 1L) stop("M must be positive.")

root <- getwd()
out_dir <- file.path(root, "simulations", "inference_diagnostics", "outputs")
combined_dir <- file.path(out_dir, "combined")
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)

read_matching <- function(pattern) {
  paths <- list.files(out_dir, pattern = pattern, full.names = TRUE)
  paths <- paths[grepl(paste0("_M", M_requested, "_"), basename(paths), fixed = TRUE)]
  if (!length(paths)) stop("No files matched: ", pattern)
  out <- do.call(rbind, lapply(paths, read.csv, check.names = FALSE))
  row.names(out) <- NULL
  out[out$M == M_requested, , drop = FALSE]
}

beta <- read_matching("^panel_fe_beta_diagnostics_summary_.*\\.csv$")
g_summary <- read_matching("^panel_fe_g_diagnostics_summary_.*\\.csv$")
g_by_x <- read_matching("^panel_fe_g_diagnostics_by_x_.*\\.csv$")

expected_beta_rows <- 3L * 2L * 2L * 3L * 5L
expected_g_summary_rows <- 3L * 2L * 2L * 3L * 5L
expected_g_by_x_rows <- expected_g_summary_rows * 50L

if (nrow(beta) != expected_beta_rows) {
  warning("Expected ", expected_beta_rows, " beta-summary rows but found ", nrow(beta), ".")
}
if (nrow(g_summary) != expected_g_summary_rows) {
  warning(
    "Expected ", expected_g_summary_rows,
    " g-summary rows but found ", nrow(g_summary), "."
  )
}
if (nrow(g_by_x) != expected_g_by_x_rows) {
  warning("Expected ", expected_g_by_x_rows, " g-by-x rows but found ", nrow(g_by_x), ".")
}

sample_label <- function(n, T) paste0("n=", n, ", T=", T)
error_label <- function(x) {
  ifelse(x == "hetero", "Heteroskedastic, serially independent", "Heteroskedastic AR(1), rho=0.5")
}

beta$sample <- factor(
  sample_label(beta$n, beta$T),
  levels = c("n=200, T=4", "n=200, T=8", "n=500, T=4")
)
beta$error_label <- error_label(beta$error_design)
g_summary$sample <- factor(
  sample_label(g_summary$n, g_summary$T),
  levels = levels(beta$sample)
)
g_summary$error_label <- error_label(g_summary$error_design)
g_by_x$sample <- factor(
  sample_label(g_by_x$n, g_by_x$T),
  levels = levels(beta$sample)
)
g_by_x$error_label <- error_label(g_by_x$error_design)

vcov_levels <- c(
  "classic_cluster",
  "penalty_cluster",
  "penalty_cluster_delta",
  "mgcv_sandwich_freq",
  "mgcv_sandwich_default"
)
method_levels <- c("bam_fe", "bam_re", "bam_fd")

for (object_name in c("beta", "g_summary", "g_by_x")) {
  object <- get(object_name)
  object$vcov_type <- factor(object$vcov_type, levels = vcov_levels)
  object$method <- factor(object$method, levels = method_levels)
  assign(object_name, object)
}

write.csv(
  beta,
  file.path(combined_dir, paste0("beta_diagnostics_all_M", M_requested, ".csv")),
  row.names = FALSE
)
write.csv(
  g_summary,
  file.path(combined_dir, paste0("g_diagnostics_all_M", M_requested, ".csv")),
  row.names = FALSE
)
write.csv(
  g_by_x,
  file.path(combined_dir, paste0("g_diagnostics_by_x_all_M", M_requested, ".csv")),
  row.names = FALSE
)

theme_diagnostics <- function() {
  theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      plot.title.position = "plot",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.text.x = element_text(angle = 25, hjust = 1)
    )
}

beta_size <- beta[abs(beta$beta_true - beta$beta_null) < 1e-12, ]
beta_power <- beta[abs(beta$beta_true - beta$beta_null) >= 1e-12, ]

size_plot <- ggplot(
  beta_size,
  aes(x = sample, y = rejection_frequency, color = vcov_type, group = vcov_type)
) +
  geom_hline(yintercept = 0.05, color = "grey45", linewidth = 0.4) +
  geom_point(size = 1.7) +
  geom_line(linewidth = 0.7) +
  facet_grid(error_label ~ method) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    title = "Empirical size of the beta test",
    subtitle = paste0("Nominal size 0.05; M = ", M_requested, "; K = 20"),
    x = NULL,
    y = "rejection frequency",
    color = "covariance"
  ) +
  theme_diagnostics()

power_plot <- ggplot(
  beta_power,
  aes(x = sample, y = rejection_frequency, color = vcov_type, group = vcov_type)
) +
  geom_point(size = 1.7) +
  geom_line(linewidth = 0.7) +
  facet_grid(error_label ~ method) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Power of the beta test against beta = 1.075",
    subtitle = paste0("Null beta = 1; M = ", M_requested, "; K = 20"),
    x = NULL,
    y = "rejection frequency",
    color = "covariance"
  ) +
  theme_diagnostics()

beta_ratio_plot <- ggplot(
  beta_size,
  aes(x = sample, y = se_to_sd, color = vcov_type, group = vcov_type)
) +
  geom_hline(yintercept = 1, color = "grey45", linewidth = 0.4) +
  geom_point(size = 1.7) +
  geom_line(linewidth = 0.7) +
  facet_grid(error_label ~ method) +
  labs(
    title = "Beta standard-error diagnostic",
    subtitle = paste0("Mean reported SE / Monte Carlo SD; M = ", M_requested),
    x = NULL,
    y = "SE / SD",
    color = "covariance"
  ) +
  theme_diagnostics()

ggsave(
  file.path(combined_dir, paste0("beta_size_M", M_requested, ".png")),
  size_plot, width = 12, height = 7, dpi = 180, bg = "white"
)
ggsave(
  file.path(combined_dir, paste0("beta_power_M", M_requested, ".png")),
  power_plot, width = 12, height = 7, dpi = 180, bg = "white"
)
ggsave(
  file.path(combined_dir, paste0("beta_se_ratio_M", M_requested, ".png")),
  beta_ratio_plot, width = 12, height = 7, dpi = 180, bg = "white"
)

for (err in unique(g_by_x$error_design)) {
  for (sample_i in levels(g_by_x$sample)) {
    d <- g_by_x[g_by_x$error_design == err & g_by_x$sample == sample_i, ]
    if (!nrow(d)) next

    coverage_plot <- ggplot(
      d,
      aes(x = x, y = point_coverage, color = vcov_type)
    ) +
      geom_hline(yintercept = 0.95, color = "grey45", linewidth = 0.4) +
      geom_line(linewidth = 0.75) +
      facet_grid(model ~ method) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(
        title = "Pointwise coverage for centered g(x)",
        subtitle = paste0(sample_i, "; ", error_label(err), "; M = ", M_requested),
        x = "x",
        y = "coverage",
        color = "covariance"
      ) +
      theme_diagnostics() +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

    ratio_plot <- ggplot(
      d,
      aes(x = x, y = se_to_sd, color = vcov_type)
    ) +
      geom_hline(yintercept = 1, color = "grey45", linewidth = 0.4) +
      geom_line(linewidth = 0.75) +
      facet_grid(model ~ method) +
      labs(
        title = "Standard-error diagnostic for centered g(x)",
        subtitle = paste0(sample_i, "; ", error_label(err), "; M = ", M_requested),
        x = "x",
        y = "mean reported SE / empirical SD",
        color = "covariance"
      ) +
      theme_diagnostics() +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

    file_tag <- paste0(
      gsub("[^A-Za-z0-9]+", "_", sample_i), "_", err, "_M", M_requested
    )
    ggsave(
      file.path(combined_dir, paste0("g_coverage_", file_tag, ".png")),
      coverage_plot, width = 11, height = 6.5, dpi = 180, bg = "white"
    )
    ggsave(
      file.path(combined_dir, paste0("g_se_ratio_", file_tag, ".png")),
      ratio_plot, width = 11, height = 6.5, dpi = 180, bg = "white"
    )
  }
}

cat("Combined beta rows:", nrow(beta), "\n")
cat("Combined g-summary rows:", nrow(g_summary), "\n")
cat("Combined g-by-x rows:", nrow(g_by_x), "\n")
cat("Wrote combined results to:", combined_dir, "\n")
