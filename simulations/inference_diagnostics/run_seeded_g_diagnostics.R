#!/usr/bin/env Rscript
# Reproduce the n=500,T=4 appendix function-coverage results using a fixed,
# versioned seed schedule. Every requested chunk is rerun; no old output is reused.
# Usage, from the repository root:
#   Rscript simulations/inference_diagnostics/run_seeded_g_diagnostics.R both 8
# Arguments: error design (both, hetero, hetero_ar1), concurrent workers,
# optional output root (default current working directory).
args <- commandArgs(trailingOnly = TRUE)
design <- if (length(args)) args[1] else "both"
workers <- if (length(args) >= 2L) as.integer(args[2]) else 2L
if (!design %in% c("both", "hetero", "hetero_ar1")) stop("Unknown error design.")
if (is.na(workers) || workers < 1L) stop("workers must be positive.")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Run with Rscript.")
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
output_root <- if (length(args) >= 3L) args[3] else getwd()
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
setwd(normalizePath(output_root))
out_dir <- file.path("simulations", "inference_diagnostics", "outputs")
log_dir <- file.path("simulations", "inference_diagnostics", "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
schedule_path <- file.path(script_dir, "seeded_g_schedule.csv")
plan <- read.csv(schedule_path, stringsAsFactors = FALSE)
stopifnot(nrow(plan) == 8L, !anyDuplicated(plan$chunk), !anyDuplicated(plan$mc_seed),
          all(plan$M == 250L), all(plan$n == 500L), all(plan$T == 4L),
          all(plan$K == 20L), all(plan$grid_points == 50L),
          all(plan$gaussian_draws == 499L), all(plan$gaussian_seed == 20260506L))
for (err in c("hetero", "hetero_ar1")) {
  p <- plan[plan$error_design == err, ]
  stopifnot(nrow(p) == 4L, identical(p$sim_first, c(1L, 251L, 501L, 751L)),
            identical(p$sim_last, c(250L, 500L, 750L, 1000L)))
}
if (design != "both") plan <- plan[plan$error_design == design, ]
chunk_script <- file.path(script_dir, "panel_fe_g_inference_diagnostics_chunk.R")
combine_script <- file.path(script_dir, "combine_g_chunks.R")
script_md5 <- unname(tools::md5sum(chunk_script))
combine_md5 <- unname(tools::md5sum(combine_script))
schedule_md5 <- unname(tools::md5sum(schedule_path))
rscript <- file.path(R.home("bin"), "Rscript")
run_chunk <- function(i) {
  job <- plan[i, , drop = FALSE]
  err_suffix <- if (job$error_design == "hetero") "hetero" else "hetero_ar1_rho0p5"
  suffix <- paste0("M250_n500_T4_sin_2pi_np_pl_J50_", err_suffix, "_dfsubstantive_chunk", job$chunk)
  summary_path <- file.path(out_dir, paste0("panel_fe_g_diagnostics_summary_", suffix, ".csv"))
  grid_path <- file.path(out_dir, paste0("panel_fe_g_diagnostics_grid_", suffix, ".csv.gz"))
  logfile <- file.path(log_dir, paste0("seeded_g_", job$chunk, ".log"))
  command_args <- c(250, 500, 4, "sin_2pi", "both", 50, 499, job$error_design,
                    0.5, "substantive", job$mc_seed, job$chunk)
  started <- Sys.time()
  cat(sprintf("%s START %s seed %d\n", format(started, tz = "UTC", usetz = TRUE), job$chunk, job$mc_seed))
  status <- system2(rscript, c(shQuote(chunk_script), shQuote(as.character(command_args))),
                    stdout = logfile, stderr = logfile)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  summary <- if (status == 0L && file.exists(summary_path)) read.csv(summary_path) else NULL
  complete <- !is.null(summary) && nrow(summary) == 30L && all(summary$n_success == 250L) &&
    file.exists(grid_path) && file.info(grid_path)$size > 0
  job$status <- if (complete) "complete" else "failed"
  job$exit_status <- if (complete) 0L else if (status == 0L) 1L else status
  job$started_utc <- format(started, tz = "UTC", usetz = TRUE)
  job$finished_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  job$elapsed_seconds <- elapsed
  job$R_version <- R.version.string
  job$mgcv_version <- as.character(packageVersion("mgcv"))
  job$ggplot2_version <- as.character(packageVersion("ggplot2"))
  job$rng_kind <- "Mersenne-Twister / Inversion / Rejection"
  job$chunk_script_md5 <- script_md5
  job$combine_script_md5 <- combine_md5
  job$schedule_md5 <- schedule_md5
  job$grid_md5 <- if (complete) unname(tools::md5sum(grid_path)) else NA_character_
  job$summary_md5 <- if (complete) unname(tools::md5sum(summary_path)) else NA_character_
  job$command <- paste("Rscript simulations/inference_diagnostics/panel_fe_g_inference_diagnostics_chunk.R",
                       paste(command_args, collapse = " "))
  write.csv(job, file.path(out_dir, paste0("seeded_g_chunk_", job$chunk, ".csv")), row.names = FALSE)
  cat(sprintf("%s %s %s (%.1f seconds)\n", format(Sys.time(), tz = "UTC", usetz = TRUE),
              toupper(job$status), job$chunk, elapsed))
  job
}
if (.Platform$OS.type == "unix" && workers > 1L) {
  records <- parallel::mclapply(seq_len(nrow(plan)), run_chunk,
    mc.cores = min(workers, nrow(plan)), mc.preschedule = FALSE, mc.set.seed = FALSE)
} else {
  records <- lapply(seq_len(nrow(plan)), run_chunk)
}
if (any(vapply(records, inherits, logical(1), "try-error"))) stop("A worker failed; inspect logs.")
records <- do.call(rbind, records)
for (err in unique(records$error_design)) {
  write.csv(records[records$error_design == err, ],
    file.path(out_dir, paste0("seeded_g_manifest_", err, ".csv")), row.names = FALSE)
}
if (any(records$status != "complete")) stop("At least one chunk failed; no combined output was produced.")
for (err in unique(plan$error_design)) {
  status <- system2(rscript, c(shQuote(combine_script), shQuote(err)))
  if (status != 0L) stop("Combination failed for ", err)
}
cat("All requested seeded coverage designs completed and combined successfully.\n")
