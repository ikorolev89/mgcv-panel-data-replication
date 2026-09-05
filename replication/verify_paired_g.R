# Verify provenance and equality of overlapping g-coverage results.
# --record-output-hashes is a deliberate maintenance action after a checked rerun.
manifest <- read.csv("replication/paired_g_manifest.csv")
fingerprints <- read.csv("replication/paired_g_code_md5.csv")
stopifnot(nrow(manifest) == 96L, length(unique(manifest$id)) == 12L,
  all(manifest$status == "complete"), all(manifest$M == 1000L),
  all(manifest$J == 50L), all(manifest$B == 499L),
  all(manifest$master_seed == 20260903L),
  all(manifest$rng_scheme == "paired-replication-v1"),
  !anyDuplicated(manifest[c("id", "chunk")]),
  identical(unname(tools::md5sum(fingerprints$path)), fingerprints$md5))
for (id in unique(manifest$id)) {
  z <- manifest[manifest$id == id, ]; z <- z[order(z$chunk), ]
  stopifnot(identical(z$sim_first, seq.int(1L, 1000L, 125L)),
    identical(z$sim_last, seq.int(125L, 1000L, 125L)))
}
key_seed <- unique(manifest[c("n", "T", "error", "design_seed")])
stopifnot(nrow(key_seed) == 6L, !anyDuplicated(key_seed$design_seed))
expected_seed <- with(key_seed, as.integer((20260903 + 1000003 * n +
  10007 * T + 1009 * match(error, c("iid", "hetero", "hetero_ar1"))) %% 2147483646 + 1))
stopifnot(identical(expected_seed, key_seed$design_seed))

configs_paired <- data.frame(n = c(200L, 200L, 500L), T = c(4L, 8L, 4L))
common_columns <- c("sim", "model", "method", "se_type", "point_covered", "point_width",
  "uniform_covered", "uniform_width", "sup_crit", "smooth_edf", "edf_ceiling",
  "edf_ratio", "binding_90", "binding_95", "grid_points")
sort_common <- function(d) {
  d <- d[do.call(order, unname(d[c("sim", "model", "method", "se_type")])), common_columns]
  row.names(d) <- NULL; d
}
for (i in seq_len(nrow(configs_paired))) {
  suffix <- sprintf("M1000_n%d_T%d_sin_2pi_np_pl_J50_hetero_ar1_rho0p5_dfsubstantive.csv.gz",
    configs_paired$n[i], configs_paired$T[i])
  main <- read.csv(paste0("replication/compact/panel_fe_g_inference_replications_", suffix))
  kval <- read.csv(paste0("replication/compact/panel_fe_g_inference_replications_", sub("_J50_", "_J50_K20_", suffix)))
  diag <- read.csv(paste0("replication/compact/panel_fe_g_diagnostics_replications_", suffix))
  diag <- diag[diag$vcov_type %in% c("classic_cluster", "penalty_cluster"), ]
  diag$se_type <- sub("_cluster", "", diag$vcov_type)
  stopifnot(identical(sort_common(main), sort_common(kval)),
    identical(sort_common(main), sort_common(diag)))
}
cat("Verified exact replication-level equality: main, K20, and diagnostic Classic/Penalty results.\n")

hash_file <- "replication/paired_g_output_md5.csv"
if ("--record-output-hashes" %in% commandArgs(trailingOnly = TRUE)) {
  files <- c(list.files("replication/compact", "^panel_fe_g_.*M1000.*csv.gz$", full.names = TRUE),
    list.files("outputs", "^panel_fe_g_inference_(summary|coverage_by_x)_M1000.*_dfsubstantive.csv$", full.names = TRUE),
    list.files("simulations/k_sensitivity/outputs", "^panel_fe_g_inference_(summary|coverage_by_x)_M1000.*_dfsubstantive.csv$", full.names = TRUE),
    list.files("simulations/inference_diagnostics/outputs", "^panel_fe_g_diagnostics_(summary|by_x)_M1000.*dfsubstantive.csv$", full.names = TRUE))
  stopifnot(length(files) == 54L)
  write.csv(data.frame(path = files, md5 = unname(tools::md5sum(files))), hash_file, row.names = FALSE)
}
recorded <- read.csv(hash_file)
stopifnot(nrow(recorded) == 54L, !anyDuplicated(recorded$path),
  identical(unname(tools::md5sum(recorded$path)), recorded$md5))
cat("Verified paired RNG manifest, code fingerprints, and 54 saved-result checksums.\n")
