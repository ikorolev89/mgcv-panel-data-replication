#!/usr/bin/env Rscript
# Fast numerical verification from saved replication records and grid summaries.
close_enough <- function(actual, expected, label) {
  if (!isTRUE(all.equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-10))) {
    stop("Mismatch: ", label)
  }
}
get_files <- function(directory, pattern, expected) {
  paths <- list.files(directory, pattern, full.names = TRUE)
  if (length(paths) != expected) stop("Expected ", expected, " files matching ", directory, "/", pattern, "; found ", length(paths))
  paths
}
read <- function(path) read.csv(path, check.names = FALSE)
group_key <- function(d, keys) do.call(paste, c(d[keys], sep = "|"))
check_reps <- function(d, keys) {
  stopifnot(!anyDuplicated(d[c(keys, "sim")]))
  groups <- split(d$sim, group_key(d, keys))
  stopifnot(all(vapply(groups, function(x) identical(sort(as.integer(x)), seq_len(1000L)), logical(1))))
}
compare_groups <- function(d, s, keys, mapping) {
  sk <- group_key(s, keys)
  dk <- group_key(d, keys)
  stopifnot(!anyDuplicated(sk), setequal(unique(dk), sk))
  for (i in seq_len(nrow(s))) {
    z <- d[dk == sk[i], , drop = FALSE]
    for (column in names(mapping)) close_enough(mapping[[column]](z), s[[column]][i], paste(sk[i], column))
  }
}
edf_mapping <- list(
  mean_edf = function(d) mean(d$smooth_edf),
  mean_edf_ceiling = function(d) mean(d$edf_ceiling),
  mean_edf_ratio = function(d) mean(d$edf_ratio),
  fraction_binding_90 = function(d) mean(d$binding_90),
  fraction_binding_95 = function(d) mean(d$binding_95)
)
mse_paths <- c(
  get_files("outputs", "^panel_fe_test_mse_losses_M1000_.*_no_forests_hetero_ar1_rho0p5.csv$", 3L),
  get_files("simulations/k_sensitivity/outputs", "^panel_fe_test_mse_losses_M1000_.*_no_forests_hetero_ar1_rho0p5.csv$", 9L)
)
for (path in mse_paths) {
  d <- read(path); s <- read(sub("_losses_", "_summary_", path))
  check_reps(d, c("dgp", "method"))
  stopifnot(nrow(s) == 21L, all(is.finite(d$mse)), all(is.finite(d$rmse)))
  mapping <- list(mean_mse = function(d) mean(d$mse), sd_mse = function(d) sd(d$mse),
    se_mse = function(d) sd(d$mse) / sqrt(nrow(d)), median_mse = function(d) median(d$mse), mean_rmse = function(d) mean(d$rmse))
  compare_groups(d, s, c("dgp", "method"), mapping)
  if ("smooth_edf" %in% names(d)) {
    compare_groups(d[!is.na(d$smooth_edf), ], s[!is.na(s$mean_edf), ], c("dgp", "method"), edf_mapping)
  }
}
cat("Verified all 12 estimation-accuracy runs.\n")
beta_paths <- c(
  get_files("outputs", "^panel_fe_beta_inference_results_M1000_.*_dfsubstantive.csv$", 6L),
  get_files("simulations/k_sensitivity/outputs", "^panel_fe_beta_inference_results_M1000_.*_dfsubstantive.csv$", 18L)
)
for (path in beta_paths) {
  d <- read(path); s <- read(sub("_results_", "_summary_", path))
  check_reps(d, "method")
  stopifnot(nrow(s) == 3L, all(s$n_success == 1000L), all(is.finite(d$beta_hat)))
  mapping <- list(rejection_classic_cluster = function(d) mean(d$reject_classic_cluster),
    rejection_penalty_cluster = function(d) mean(d$reject_penalty_cluster),
    mean_beta_hat = function(d) mean(d$beta_hat), sd_beta_hat = function(d) sd(d$beta_hat),
    mean_se_classic_cluster = function(d) mean(d$se_classic_cluster), mean_se_penalty_cluster = function(d) mean(d$se_penalty_cluster))
  if ("smooth_edf" %in% names(d)) mapping <- c(mapping, edf_mapping)
  compare_groups(d, s, "method", mapping)
}
diag_beta <- get_files("simulations/inference_diagnostics/outputs", "^panel_fe_beta_diagnostics_results_M1000_.*_dfsubstantive.csv$", 12L)
for (path in diag_beta) {
  d <- read(path); s <- read(sub("_results_", "_summary_", path))
  check_reps(d, c("method", "vcov_type"))
  stopifnot(nrow(s) == 15L, all(s$n_success == 1000L))
  compare_groups(d, s, c("method", "vcov_type"), list(
    rejection_frequency = function(d) mean(d$reject), beta_coverage = function(d) mean(d$beta_covered),
    mean_beta_hat = function(d) mean(d$beta_hat), sd_beta_hat = function(d) sd(d$beta_hat),
    mean_se = function(d) mean(d$se), se_to_sd = function(d) mean(d$se) / sd(d$beta_hat),
    mean_delta_beta_variance = function(d) mean(d$delta_beta_variance)))
}
cat("Verified all 36 beta-inference runs.\n")
g_paths <- get_files("replication/compact", "^panel_fe_g_.*_replications_M1000_.*.csv.gz$", 18L)
for (path in g_paths) {
  d <- read(path)
  original <- unique(d$source_file); stopifnot(length(original) == 1L)
  s <- read(sub(".csv.gz", ".csv", sub("_grid_", "_summary_", original), fixed = TRUE))
  diagnostic <- "vcov_type" %in% names(d)
  xp <- sub("_grid_", if (diagnostic) "_by_x_" else "_coverage_by_x_", original)
  x <- read(sub(".csv.gz", ".csv", xp, fixed = TRUE))
  keys <- c("model", "method", if (diagnostic) "vcov_type" else "se_type")
  check_reps(d, keys)
  stopifnot(nrow(s) == if (diagnostic) 30L else 12L, all(d$grid_points == 50L),
            all(d$uniform_covered %in% c(0, 1)), !anyDuplicated(x[c(keys, "x")]))
  mapping <- list(avg_pointwise_coverage = function(d) mean(d$point_covered),
    avg_pointwise_width = function(d) mean(d$point_width), uniform_coverage = function(d) mean(d$uniform_covered),
    avg_uniform_width = function(d) mean(d$uniform_width), avg_sup_crit = function(d) mean(d$sup_crit))
  if ("smooth_edf" %in% names(d)) mapping <- c(mapping, edf_mapping)
  compare_groups(d, s, keys, mapping)
  compare_groups(x, s, keys, list(avg_pointwise_coverage = function(d) mean(d$point_coverage),
    min_pointwise_coverage = function(d) min(d$point_coverage)))
  if (diagnostic) {
    stopifnot(all(s$n_success == 1000L))
    compare_groups(x, s, keys, list(abs_bias = function(d) mean(d$abs_bias),
      empirical_sd = function(d) mean(d$empirical_sd), mean_reported_se = function(d) mean(d$mean_reported_se),
      se_to_sd = function(d) mean(d$se_to_sd), max_abs_bias = function(d) max(d$abs_bias)))
  }
}
cat("Verified all 18 coverage runs from compact records and point-specific summaries.\n")
seed_dir <- "simulations/inference_diagnostics"
seed_plan <- read(file.path(seed_dir, "seeded_g_schedule.csv"))
stopifnot(nrow(seed_plan) == 8L, !anyDuplicated(seed_plan$mc_seed),
          !anyDuplicated(seed_plan$chunk), all(seed_plan$M == 250L),
          all(seed_plan$gaussian_seed == 20260506L))
for (err in c("hetero", "hetero_ar1")) {
  plan <- seed_plan[seed_plan$error_design == err, ]
  records <- read(file.path(seed_dir, "outputs", paste0("seeded_g_manifest_", err, ".csv")))
  stopifnot(nrow(records) == 4L, all(records$status == "complete"), all(records$exit_status == 0L))
  for (column in names(plan)) stopifnot(identical(records[[column]], plan[[column]]))
  for (spec in list(c("chunk_script_md5", "panel_fe_g_inference_diagnostics_chunk.R"),
                    c("combine_script_md5", "combine_g_chunks.R"),
                    c("schedule_md5", "seeded_g_schedule.csv"))) {
    stopifnot(all(records[[spec[1]]] == unname(tools::md5sum(file.path(seed_dir, spec[2])))))
  }
  suffix <- paste0("M1000_n500_T4_sin_2pi_np_pl_J50_", if (err == "hetero") "hetero" else "hetero_ar1_rho0p5", "_dfsubstantive")
  compact <- read(file.path("replication/compact", paste0("panel_fe_g_diagnostics_replications_", suffix, ".csv.gz")))
  for (i in seq_len(nrow(plan))) {
    chunk <- plan[i, ]
    summary_path <- file.path(seed_dir, "outputs", paste0("panel_fe_g_diagnostics_summary_",
      sub("M1000", "M250", suffix), "_chunk", chunk$chunk, ".csv"))
    stopifnot(unname(tools::md5sum(summary_path)) == records$summary_md5[i])
    d <- compact[compact$sim >= chunk$sim_first & compact$sim <= chunk$sim_last, ]
    d$sim <- d$sim - chunk$sim_first + 1L
    s <- read(summary_path)
    stopifnot(nrow(d) == 250L * 30L, nrow(s) == 30L, all(s$n_success == 250L))
    compare_groups(d, s, c("model", "method", "vcov_type"), list(
      avg_pointwise_coverage = function(d) mean(d$point_covered),
      avg_pointwise_width = function(d) mean(d$point_width),
      uniform_coverage = function(d) mean(d$uniform_covered),
      avg_uniform_width = function(d) mean(d$uniform_width),
      avg_sup_crit = function(d) mean(d$sup_crit)))
  }
}
cat("Verified the eight recorded chunk seeds, code fingerprints, and chunk-to-combined results.\n")
manifest <- read("simulations/inference_diagnostics/outputs/combined/final_manifest_M1000.csv")
stopifnot(nrow(manifest) == 18L, all(manifest$min_success == 1000L), all(manifest$max_success == 1000L))
status <- system2(file.path(R.home("bin"), "Rscript"), c("replication/generate_main_tables.R", "--check"))
stopifnot(status == 0L)
cat("PASS: 66 saved simulation runs; 1,000 replications per design-method cell.\n")
