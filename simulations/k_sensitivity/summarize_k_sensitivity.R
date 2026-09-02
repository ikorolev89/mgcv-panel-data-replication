#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  getwd()
}
setwd(script_dir)

input_dir <- "outputs"
combined_dir <- "combined"
dir.create(combined_dir, showWarnings = FALSE, recursive = TRUE)

extract_integer <- function(x, pattern) {
  as.integer(sub(paste0(".*", pattern, "([0-9]+).*"), "\\1", x))
}

read_with_source <- function(files) {
  do.call(rbind, lapply(files, function(f) {
    d <- read.csv(f, stringsAsFactors = FALSE)
    d$source_file <- basename(f)
    d
  }))
}

mse_files <- list.files(
  input_dir,
  pattern = "^panel_fe_test_mse_summary_M1000_.*\\.csv$",
  full.names = TRUE
)
beta_files <- list.files(
  input_dir,
  pattern = "^panel_fe_beta_inference_summary_M1000_.*\\.csv$",
  full.names = TRUE
)
coverage_files <- list.files(
  input_dir,
  pattern = "^panel_fe_g_inference_summary_M1000_.*\\.csv$",
  full.names = TRUE
)

if (length(mse_files) != 9L) stop("Expected 9 MSE summaries; found ", length(mse_files), ".")
if (length(beta_files) != 18L) stop("Expected 18 beta summaries; found ", length(beta_files), ".")
if (length(coverage_files) != 9L) stop("Expected 9 coverage summaries; found ", length(coverage_files), ".")

mse <- read_with_source(mse_files)
mse$n <- extract_integer(mse$source_file, "_n")
mse$T <- extract_integer(mse$source_file, "_T")
mse$K_smooth_file <- extract_integer(mse$source_file, "_Kbam")
mse$configuration <- paste0("n=", mse$n, ", T=", mse$T)
mse$rmse_times_100 <- 100 * mse$mean_rmse

beta <- read_with_source(beta_files)
beta$n <- extract_integer(beta$source_file, "_n")
beta$T <- extract_integer(beta$source_file, "_T")
beta$beta_true <- ifelse(grepl("_btrue1p075_", beta$source_file), 1.075, 1)
beta$configuration <- paste0("n=", beta$n, ", T=", beta$T)

coverage <- read_with_source(coverage_files)
coverage$n <- extract_integer(coverage$source_file, "_n")
coverage$T <- extract_integer(coverage$source_file, "_T")
coverage$configuration <- paste0("n=", coverage$n, ", T=", coverage$T)

write.csv(mse, file.path(combined_dir, "mse_summary_all_K.csv"), row.names = FALSE)
write.csv(beta, file.path(combined_dir, "beta_summary_all_K.csv"), row.names = FALSE)
write.csv(coverage, file.path(combined_dir, "coverage_summary_all_K.csv"), row.names = FALSE)

mse_bam <- mse[grepl("mgcv::bam", mse$method), ]
mse_stability <- do.call(rbind, lapply(
  split(mse_bam, list(mse_bam$configuration, mse_bam$dgp, mse_bam$method), drop = TRUE),
  function(d) {
    baseline <- d$mean_rmse[d$K_smooth == 20]
    data.frame(
      configuration = d$configuration[1],
      dgp = d$dgp[1],
      method = d$method[1],
      min_rmse_times_100 = 100 * min(d$mean_rmse),
      max_rmse_times_100 = 100 * max(d$mean_rmse),
      range_rmse_times_100 = 100 * diff(range(d$mean_rmse)),
      relative_range_vs_K20 = diff(range(d$mean_rmse)) / baseline
    )
  }
))
row.names(mse_stability) <- NULL

beta_stability <- do.call(rbind, lapply(
  split(beta, list(beta$configuration, beta$beta_true, beta$method), drop = TRUE),
  function(d) data.frame(
    configuration = d$configuration[1],
    beta_true = d$beta_true[1],
    method = d$method[1],
    min_rejection_penalty = min(d$rejection_penalty_cluster),
    max_rejection_penalty = max(d$rejection_penalty_cluster),
    range_rejection_penalty = diff(range(d$rejection_penalty_cluster)),
    range_mean_beta_hat = diff(range(d$mean_beta_hat))
  )
))
row.names(beta_stability) <- NULL

coverage_penalty <- coverage[coverage$se_type == "penalty", ]
coverage_stability <- do.call(rbind, lapply(
  split(
    coverage_penalty,
    list(coverage_penalty$configuration, coverage_penalty$model, coverage_penalty$method),
    drop = TRUE
  ),
  function(d) data.frame(
    configuration = d$configuration[1],
    model = d$model[1],
    method = d$method[1],
    range_avg_pointwise_coverage = diff(range(d$avg_pointwise_coverage)),
    range_min_pointwise_coverage = diff(range(d$min_pointwise_coverage)),
    range_uniform_coverage = diff(range(d$uniform_coverage)),
    range_avg_pointwise_width = diff(range(d$avg_pointwise_width)),
    range_avg_uniform_width = diff(range(d$avg_uniform_width))
  )
))
row.names(coverage_stability) <- NULL

write.csv(mse_stability, file.path(combined_dir, "mse_stability_ranges.csv"), row.names = FALSE)
write.csv(beta_stability, file.path(combined_dir, "beta_stability_ranges.csv"), row.names = FALSE)
write.csv(coverage_stability, file.path(combined_dir, "coverage_stability_ranges.csv"), row.names = FALSE)

theme_set(theme_minimal(base_size = 11))

p_mse <- ggplot(mse_bam, aes(K_smooth, rmse_times_100, color = method)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.7) +
  facet_grid(dgp ~ configuration, scales = "free_y") +
  scale_x_continuous(breaks = c(10, 20, 40)) +
  labs(x = "Basis dimension K", y = "100 x mean out-of-sample RMSE", color = NULL)
ggsave(file.path(combined_dir, "mse_by_K.png"), p_mse, width = 10, height = 8, dpi = 180, bg = "white")

p_beta <- ggplot(
  beta,
  aes(K_smooth, rejection_penalty_cluster, color = method)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.7) +
  facet_grid(beta_true ~ configuration, labeller = label_both) +
  scale_x_continuous(breaks = c(10, 20, 40)) +
  labs(x = "Basis dimension K", y = "Rejection frequency", color = NULL)
ggsave(file.path(combined_dir, "beta_rejection_by_K.png"), p_beta, width = 10, height = 6, dpi = 180, bg = "white")

coverage_long <- rbind(
  transform(coverage_penalty, estimand = "Average pointwise", coverage_value = avg_pointwise_coverage),
  transform(coverage_penalty, estimand = "Uniform", coverage_value = uniform_coverage)
)
p_coverage <- ggplot(
  coverage_long,
  aes(K_smooth, coverage_value, color = method)
) +
  geom_hline(yintercept = 0.95, color = "grey55", linewidth = 0.35) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  facet_grid(interaction(model, estimand) ~ configuration, scales = "free_y") +
  scale_x_continuous(breaks = c(10, 20, 40)) +
  labs(x = "Basis dimension K", y = "Coverage", color = NULL)
ggsave(file.path(combined_dir, "coverage_by_K.png"), p_coverage, width = 10, height = 9, dpi = 180, bg = "white")

edf_mse <- mse_bam[, c(
  "configuration", "dgp", "method", "K_smooth", "mean_edf",
  "mean_edf_ceiling", "mean_edf_ratio", "fraction_binding_90", "fraction_binding_95"
)]
edf_mse$block <- "mse"
edf_mse$model <- edf_mse$dgp
edf_mse$dgp <- NULL

edf_beta <- beta[, c(
  "configuration", "method", "K_smooth", "mean_edf", "mean_edf_ceiling",
  "mean_edf_ratio", "fraction_binding_90", "fraction_binding_95"
)]
edf_beta$block <- "beta"
edf_beta$model <- "pl"
edf_beta <- unique(edf_beta)

edf_coverage <- coverage_penalty[, c(
  "configuration", "model", "method", "K_smooth", "mean_edf",
  "mean_edf_ceiling", "mean_edf_ratio", "fraction_binding_90", "fraction_binding_95"
)]
edf_coverage$block <- "coverage"

edf <- rbind(
  edf_mse[, names(edf_coverage)],
  edf_beta[, names(edf_coverage)],
  edf_coverage
)
edf <- edf[order(edf$block, edf$configuration, edf$model, edf$method, edf$K_smooth), ]
write.csv(edf, file.path(combined_dir, "edf_binding_summary.csv"), row.names = FALSE)

cat("Wrote combined summaries and plots to", normalizePath(combined_dir), "\n")

