#!/usr/bin/env Rscript
# Retain the replication-level quantities needed for the reported coverage
# summaries without distributing redundant 50-point interval endpoints.
# Run after the simulation scripts when refreshing the replication package.
dir.create("replication/compact", recursive = TRUE, showWarnings = FALSE)
paths <- c(
  list.files("outputs", "^panel_fe_g_inference_grid_M1000_.*_dfsubstantive.csv$", full.names = TRUE),
  list.files("simulations/k_sensitivity/outputs", "^panel_fe_g_inference_grid_M1000_.*_dfsubstantive.csv$", full.names = TRUE),
  list.files("simulations/inference_diagnostics/outputs", "^panel_fe_g_diagnostics_grid_M1000_.*_dfsubstantive.csv.gz$", full.names = TRUE)
)
if (length(paths) != 18L) stop("Expected 18 full coverage traces; found ", length(paths))
for (path in paths) {
  message("Compacting ", basename(path))
  input <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  d <- read.csv(input, check.names = FALSE)
  close(input)
  covariance <- if ("se_type" %in% names(d)) "se_type" else "vcov_type"
  keys <- c("sim", "model", "method", covariance)
  values <- c("point_covered", "point_width", "uniform_covered", "uniform_width", "sup_crit")
  values <- c(values, intersect(c("smooth_edf", "edf_ceiling", "edf_ratio", "binding_90", "binding_95"), names(d)))
  if (any(!is.finite(as.matrix(d[values])))) stop("Nonfinite coverage result: ", path)
  out <- aggregate(d[values], d[keys], mean)
  counts <- aggregate(d["x"], d[keys], length)
  stopifnot(identical(out[keys], counts[keys]), all(counts$x == 50L),
            all(out$uniform_covered %in% c(0, 1)),
            !anyDuplicated(d[c(keys, "x")]))
  out$grid_points <- counts$x
  out$source_file <- path
  dest <- file.path("replication/compact", sub("\\.csv(\\.gz)?$", ".csv.gz", sub("_grid_", "_replications_", basename(path))))
  connection <- gzfile(dest, "wt")
  write.csv(out, connection, row.names = FALSE)
  close(connection)
  rm(d, out)
  gc(verbose = FALSE)
}
message("Wrote compact records for all 18 coverage runs.")
