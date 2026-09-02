#!/usr/bin/env Rscript
# Exact arguments for Tables 1--5. Run from the repository root.
# Optional first argument: Monte Carlo repetitions (default 1000).
args <- commandArgs(trailingOnly = TRUE)
M <- if (length(args)) as.integer(args[1]) else 1000L
if (is.na(M) || M < 1L) stop("M must be a positive integer.")
run <- function(script, args) {
  message("Running ", script, " ", paste(args, collapse = " "))
  status <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(script), shQuote(as.character(args))))
  if (status != 0L) stop("Simulation failed: ", script)
}
for (cfg in list(c(200, 4), c(200, 8), c(500, 4))) {
  n <- cfg[1]; TT <- cfg[2]
  run("simulations/panel_fe_test_mse_first_pass.R", c(M, 7, 500, TT, "paper_main", n, "FALSE", "hetero_ar1", 0.5))
  for (beta in c(1, 1.075)) {
    run("simulations/panel_fe_beta_inference_bootstrap.R", c(M, 0, n, TT, "baseline", beta, 1, "hetero_ar1", 0.5, "substantive"))
  }
  run("simulations/panel_fe_g_inference_coverage.R", c(M, n, TT, "sin_2pi", "both", 50, 499, "hetero_ar1", 0.5, "substantive"))
}
