#!/usr/bin/env Rscript
# Run small fits from every simulation family in a temporary directory.
# Saved paper results are not overwritten. This checks execution, not accuracy.
main <- function() {
  root <- normalizePath(getwd())
  scratch <- tempfile("replication smoke-")
  dir.create(scratch)
  on.exit(setwd(root), add = TRUE)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  cases <- list(
    list("simulations/panel_fe_test_mse_first_pass.R", c(2, 7, 50, 4, "paper_main", 40, "FALSE", "hetero_ar1", 0.5)),
    list("simulations/panel_fe_beta_inference_bootstrap.R", c(2, 0, 40, 4, "baseline", 1, 1, "hetero_ar1", 0.5, "substantive")),
    list("simulations/panel_fe_g_inference_coverage.R", c(2, 40, 4, "sin_2pi", "both", 10, 99, "hetero_ar1", 0.5, "substantive")),
    list("simulations/k_sensitivity/panel_fe_test_mse_k_sensitivity.R", c(2, 7, 50, 4, "paper_main", 40, "FALSE", "hetero_ar1", 0.5, 10)),
    list("simulations/k_sensitivity/panel_fe_beta_inference_k_sensitivity.R", c(2, 0, 40, 4, "sine_g", 1, 1, "hetero_ar1", 0.5, "substantive", 10)),
    list("simulations/k_sensitivity/panel_fe_g_inference_k_sensitivity.R", c(2, 40, 4, "sin_2pi", "both", 10, 99, "hetero_ar1", 0.5, "substantive", 10)),
    list("simulations/inference_diagnostics/panel_fe_beta_inference_diagnostics.R", c(2, 0, 40, 4, "baseline", 1, 1, "hetero_ar1", 0.5, "substantive")),
    list("simulations/inference_diagnostics/panel_fe_g_inference_diagnostics.R", c(2, 40, 4, "sin_2pi", "both", 10, 99, "hetero_ar1", 0.5, "substantive"))
  )
  for (i in seq_along(cases)) {
    case <- cases[[i]]
    script <- file.path(root, case[[1]])
    # Each script writes relative to its working directory.
    setwd(scratch)
    logfile <- file.path(scratch, paste0("case-", i, ".log"))
    status <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(script), shQuote(case[[2]])), stdout = logfile, stderr = logfile)
    if (status != 0L) stop(paste(readLines(logfile, warn = FALSE), collapse = "\n"))
    cat("PASS", case[[1]], "\n")
  }
  cat("All eight simulation smoke cases passed.\n")
}
main()
