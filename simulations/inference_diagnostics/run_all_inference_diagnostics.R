#!/usr/bin/env Rscript

# Resumable runner for the inference diagnostic exercise.
#
# Usage:
#   Rscript simulations/inference_diagnostics/run_all_inference_diagnostics.R 1000 6
#
# Arguments:
#   1. Monte Carlo repetitions per design (default 1000)
#   2. maximum number of concurrent R jobs (default 6)

args <- commandArgs(trailingOnly = TRUE)
M <- if (length(args) >= 1L) as.integer(args[[1]]) else 1000L
n_workers <- if (length(args) >= 2L) as.integer(args[[2]]) else 6L

if (is.na(M) || M < 1L) stop("M must be a positive integer.")
if (is.na(n_workers) || n_workers < 1L) stop("n_workers must be positive.")

root <- getwd()
script_dir <- file.path(root, "simulations", "inference_diagnostics")
out_dir <- file.path(script_dir, "outputs")
log_dir <- file.path(script_dir, "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

configs <- data.frame(
  n = c(200L, 200L, 500L),
  T = c(4L, 8L, 4L)
)
error_designs <- c("hetero", "hetero_ar1")
beta_values <- c(1, 1.075)

format_suffix_num <- function(x) {
  out <- gsub("\\.", "p", format(x, scientific = FALSE, trim = TRUE))
  gsub("-", "m", out, fixed = TRUE)
}

jobs <- list()
job_index <- 0L

add_job <- function(kind, label, script, script_args, expected) {
  job_index <<- job_index + 1L
  jobs[[job_index]] <<- list(
    kind = kind,
    label = label,
    script = script,
    args = script_args,
    expected = expected,
    log = file.path(log_dir, paste0(label, ".log"))
  )
}

for (i in seq_len(nrow(configs))) {
  n_i <- configs$n[i]
  T_i <- configs$T[i]

  for (err in error_designs) {
    err_suffix <- if (err == "hetero") "_hetero" else "_hetero_ar1_rho0p5"

    g_suffix <- paste0(
      "M", M, "_n", n_i, "_T", T_i,
      "_sin_2pi_np_pl_J50", err_suffix, "_dfsubstantive"
    )
    g_label <- paste0("g_n", n_i, "_T", T_i, "_", err)
    add_job(
      kind = "g",
      label = g_label,
      script = file.path(script_dir, "panel_fe_g_inference_diagnostics.R"),
      script_args = c(M, n_i, T_i, "sin_2pi", "both", 50, 499, err, 0.5, "substantive"),
      expected = file.path(out_dir, paste0("panel_fe_g_diagnostics_summary_", g_suffix, ".csv"))
    )

    for (beta_true in beta_values) {
      beta_suffix <- paste0(
        "M", M, "_B0_n", n_i, "_T", T_i,
        "_baseline_btrue", format_suffix_num(beta_true),
        "_bnull1", err_suffix, "_dfsubstantive"
      )
      beta_label <- paste0(
        "beta_n", n_i, "_T", T_i, "_b", format_suffix_num(beta_true), "_", err
      )
      add_job(
        kind = "beta",
        label = beta_label,
        script = file.path(script_dir, "panel_fe_beta_inference_diagnostics.R"),
        script_args = c(M, 0, n_i, T_i, "baseline", beta_true, 1, err, 0.5, "substantive"),
        expected = file.path(out_dir, paste0("panel_fe_beta_diagnostics_summary_", beta_suffix, ".csv"))
      )
    }
  }
}

run_job <- function(job) {
  if (file.exists(job$expected) && file.info(job$expected)$size > 0) {
    cat(format(Sys.time()), "SKIP", job$label, "(completed output exists)\n")
    return(data.frame(label = job$label, status = 0L, skipped = TRUE))
  }

  cat(format(Sys.time()), "START", job$label, "\n")
  status <- system2(
    "Rscript",
    args = c(shQuote(job$script), shQuote(as.character(job$args))),
    stdout = job$log,
    stderr = job$log
  )
  complete <- identical(status, 0L) &&
    file.exists(job$expected) && file.info(job$expected)$size > 0
  cat(
    format(Sys.time()),
    if (complete) "DONE" else "FAILED",
    job$label,
    "status", status,
    "\n"
  )
  data.frame(label = job$label, status = if (complete) 0L else if (status == 0L) 1L else status, skipped = FALSE)
}

cat(
  "Running", length(jobs), "jobs with M =", M,
  "and up to", n_workers, "concurrent workers.\n"
)

if (.Platform$OS.type == "unix" && n_workers > 1L) {
  result_list <- parallel::mclapply(
    jobs,
    run_job,
    mc.cores = min(n_workers, length(jobs)),
    mc.preschedule = FALSE
  )
} else {
  result_list <- lapply(jobs, run_job)
}

run_results <- do.call(rbind, result_list)
manifest_path <- file.path(
  out_dir,
  paste0("inference_diagnostics_run_manifest_M", M, ".csv")
)
write.csv(run_results, manifest_path, row.names = FALSE)

if (any(run_results$status != 0L)) {
  print(run_results[run_results$status != 0L, ], row.names = FALSE)
  stop("At least one simulation job failed; inspect the corresponding log files.")
}

cat("All jobs completed successfully.\n")
cat("Manifest:", manifest_path, "\n")
