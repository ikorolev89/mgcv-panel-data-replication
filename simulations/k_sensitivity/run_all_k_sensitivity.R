#!/usr/bin/env Rscript

# Run the paper's simulation designs at several penalized-smooth basis dimensions.
# Completed jobs are marked and skipped on a later invocation, so this runner is
# safe to resume after an interruption.

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  getwd()
}
setwd(script_dir)

trailing <- commandArgs(trailingOnly = TRUE)
n_workers <- if (length(trailing)) as.integer(trailing[[1]]) else 9L
if (is.na(n_workers) || n_workers < 1L) stop("n_workers must be positive.")

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE, recursive = TRUE)
dir.create("status", showWarnings = FALSE, recursive = TRUE)

M <- 1000L
Ks <- c(10L, 20L, 40L)
configs <- data.frame(
  n = c(200L, 200L, 500L),
  T = c(4L, 8L, 4L)
)

jobs <- list()
add_job <- function(id, block, script, args) {
  jobs[[length(jobs) + 1L]] <<- list(
    id = id,
    block = block,
    script = script,
    args = as.character(args)
  )
}

for (ii in seq_len(nrow(configs))) {
  n <- configs$n[ii]
  TT <- configs$T[ii]
  for (K in Ks) {
    add_job(
      sprintf("mse_n%d_T%d_K%d", n, TT, K),
      "mse",
      "panel_fe_test_mse_k_sensitivity.R",
      c(M, 7, 500, TT, "paper_main", n, "FALSE", "hetero_ar1", 0.5, K)
    )

    for (beta_true in c(1, 1.075)) {
      beta_tag <- gsub("\\.", "p", format(beta_true, trim = TRUE))
      add_job(
        sprintf("beta_n%d_T%d_b%s_K%d", n, TT, beta_tag, K),
        "beta",
        "panel_fe_beta_inference_k_sensitivity.R",
        c(M, 0, n, TT, "sine_g", beta_true, 1, "hetero_ar1", 0.5, "substantive", K)
      )
    }

    add_job(
      sprintf("coverage_n%d_T%d_K%d", n, TT, K),
      "coverage",
      "panel_fe_g_inference_k_sensitivity.R",
      c(M, n, TT, "sin_2pi", "both", 50, 499, "hetero_ar1", 0.5, "substantive", K)
    )
  }
}

run_job <- function(job) {
  done_path <- file.path("status", paste0(job$id, ".done"))
  failed_path <- file.path("status", paste0(job$id, ".failed"))
  log_path <- file.path("logs", paste0(job$id, ".log"))

  if (file.exists(done_path) && job$block != "coverage") {
    cat(sprintf("SKIP  %s (already complete)\n", job$id))
    return(data.frame(
      id = job$id, block = job$block, status = "skipped_complete",
      elapsed_seconds = NA_real_, log = log_path
    ))
  }

  if (file.exists(failed_path)) unlink(failed_path)
  cat(sprintf("START %s\n", job$id))
  start <- proc.time()[["elapsed"]]
  status <- system2(
    "Rscript",
    c(job$script, job$args),
    stdout = log_path,
    stderr = log_path
  )
  elapsed <- proc.time()[["elapsed"]] - start

  if (identical(status, 0L)) {
    writeLines(
      c(
        paste0("completed_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
        paste0("elapsed_seconds=", elapsed),
        paste0("log=", log_path)
      ),
      done_path
    )
    final_status <- "complete"
  } else {
    writeLines(
      c(
        paste0("failed_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
        paste0("exit_status=", status),
        paste0("elapsed_seconds=", elapsed),
        paste0("log=", log_path)
      ),
      failed_path
    )
    final_status <- paste0("failed_", status)
  }

  cat(sprintf("END   %s [%s, %.1f sec]\n", job$id, final_status, elapsed))
  data.frame(
    id = job$id, block = job$block, status = final_status,
    elapsed_seconds = elapsed, log = log_path
  )
}

cat(sprintf(
  "Running %d jobs with %d workers (M=%d; K=%s)\n",
  length(jobs), n_workers, M, paste(Ks, collapse = ",")
))

if (.Platform$OS.type == "unix" && n_workers > 1L) {
  run_results <- parallel::mclapply(
    jobs,
    run_job,
    mc.cores = n_workers,
    mc.preschedule = FALSE
  )
} else {
  run_results <- lapply(jobs, run_job)
}

run_results <- do.call(rbind, run_results)
write.csv(run_results, "run_manifest.csv", row.names = FALSE)
print(run_results, row.names = FALSE)

if (any(startsWith(run_results$status, "failed"))) {
  quit(status = 1L)
}
