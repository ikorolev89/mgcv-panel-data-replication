#!/usr/bin/env Rscript
# Paired g-coverage rerun: --test, --run [workers], --summarize, or --promote.
# Run from the repository root. Chunk files are resumable only with matching code.
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "--test"
workers <- if (length(args) >= 2) as.integer(args[2]) else 6L
root <- normalizePath(getwd())
stage <- file.path(root, "tmp", "paired_g_20260903")
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
paths <- c(main = "simulations/panel_fe_g_inference_coverage.R",
  k = "simulations/k_sensitivity/panel_fe_g_inference_k_sensitivity.R",
  diag = "simulations/inference_diagnostics/panel_fe_g_inference_diagnostics.R")
configs <- data.frame(n = c(200L, 200L, 500L), T = c(4L, 8L, 4L))
cases <- do.call(rbind, lapply(seq_len(nrow(configs)), function(i) {
  data.frame(n = configs$n[i], T = configs$T[i],
    kind = c("diag", "diag", "k", "k"),
    error = c("hetero", "hetero_ar1", "hetero_ar1", "hetero_ar1"),
    K = c(20L, 20L, 10L, 40L))
}))
cases$id <- with(cases, paste(kind, n, T, error, K, sep = "_"))
code_files <- c(unname(paths), "simulations/g_coverage_shared.R",
  "replication/run_paired_g_simulations.R")
hashes <- tools::md5sum(file.path(root, code_files))
names(hashes) <- code_files

load_case <- function(case, M = 1000L, B = 499L) {
  e <- new.env(parent = .GlobalEnv)
  e$g_source_file <- file.path(root, paths[[case$kind]])
  cli <- c(M, case$n, case$T, "sin_2pi", "both", 50, B,
    case$error, 0.5, "substantive", case$K)
  e$commandArgs <- function(trailingOnly = TRUE) as.character(cli)
  lines <- readLines(e$g_source_file)
  marker <- which(lines == "# G_COVERAGE_EXECUTION_START")
  stopifnot(length(marker) == 1L)
  suppressMessages(eval(parse(text = lines[seq_len(marker - 1L)]), e))
  e$summary_code <- lines[seq(which(lines == "results <- do.call(rbind, results)"), length(lines))]
  e
}

if (mode == "--test") {
  case <- cases[2, ]; case$n <- 30L
  main_case <- case; main_case$kind <- "main"
  k_case <- case; k_case$kind <- "k"
  em <- load_case(main_case, 2L); ek <- load_case(k_case, 2L); ed <- load_case(case, 2L)
  main <- em$one_draw(1L); kval <- ek$one_draw(1L); diag <- ed$one_draw(1L)
  diag <- diag[diag$vcov_type %in% c("classic_cluster", "penalty_cluster"), ]
  diag$se_type <- sub("_cluster", "", diag$vcov_type)
  normalize <- function(d, cols) {
    d <- d[do.call(order, unname(d[c("sim", "model", "method", "se_type", "x")])), cols]
    row.names(d) <- NULL; d
  }
  cols <- names(main)
  stopifnot(isTRUE(all.equal(normalize(main, cols), normalize(kval, cols), tolerance = 1e-12)),
    isTRUE(all.equal(normalize(main, cols), normalize(diag, cols), tolerance = 1e-12)))
  first <- ed$one_draw(1L); second <- ed$one_draw(2L)
  stopifnot(identical(first, ed$one_draw(1L)), identical(second, ed$one_draw(2L)))
  eb <- load_case(case, 2L, 999L)
  ed$begin_g_replication(2L); data1 <- ed$simulate_panel(); z1 <- ed$sup_z
  eb$begin_g_replication(2L); data2 <- eb$simulate_panel()
  ed$begin_g_replication(1L); z2 <- ed$sup_z
  stopifnot(identical(data1, data2), !identical(z1, z2))
  point_cols <- setdiff(names(first), c("sup_crit", "uniform_lower", "uniform_upper", "uniform_covered", "uniform_width"))
  stopifnot(identical(first[point_cols], eb$one_draw(1L)[point_cols]))
  ed$vcov_types <- rev(ed$vcov_types)
  revfit <- ed$one_draw(1L)
  ord <- function(d) {d <- d[do.call(order, unname(d[c("model", "method", "vcov_type", "x")])), ]; row.names(d) <- NULL; d}
  stopifnot(identical(ord(first), ord(revfit)))
  message("PASS: main/K20/diagnostic overlap; repeat/order invariance; separate data and band RNG; fresh draws across replications.")
  quit(status = 0L)
}

jobs <- merge(data.frame(case = seq_len(nrow(cases))), data.frame(chunk = seq_len(8L)))
jobs <- jobs[order(jobs$chunk, jobs$case), ]
chunk_path <- function(case, chunk) file.path(stage, paste0(case$id, "_chunk", chunk, ".rds"))
expected_meta <- function(case, chunk) list(case = case, chunk = chunk,
  ids = seq.int((chunk - 1L) * 125L + 1L, chunk * 125L),
  M = 1000L, J = 50L, B = 499L, code_md5 = hashes,
  rng_scheme = "paired-replication-v1", master_seed = 20260903L,
  R = R.version.string, mgcv = as.character(packageVersion("mgcv")))
read_chunk <- function(case, chunk) {
  z <- readRDS(chunk_path(case, chunk))
  stopifnot(identical(z$meta, expected_meta(case, chunk)),
    identical(sort(unique(z$data$sim)), z$meta$ids),
    nrow(z$data) == 125L * 50L * 6L * if (case$kind == "diag") 5L else 2L,
    all(is.finite(z$data$estimate)), all(is.finite(z$data$se)))
  z
}
if (mode == "--run") {
  stopifnot(is.finite(workers), workers >= 1L)
  run_job <- function(j) {
    case <- cases[jobs$case[j], ]; chunk <- jobs$chunk[j]
    path <- chunk_path(case, chunk)
    if (file.exists(path)) {
      read_chunk(case, chunk)
      message("Verified cached ", basename(path)); return(TRUE)
    }
    e <- load_case(case)
    meta <- expected_meta(case, chunk)
    started <- Sys.time()
    data <- do.call(rbind, lapply(meta$ids, e$one_draw))
    payload <- list(meta = meta, data = data, design_seed = e$g_design_seed,
      started = format(started, tz = "UTC", usetz = TRUE),
      finished = format(Sys.time(), tz = "UTC", usetz = TRUE),
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
    temporary <- paste0(path, ".partial")
    saveRDS(payload, temporary, compress = "gzip")
    stopifnot(file.rename(temporary, path))
    read_chunk(case, chunk)
    message("Completed ", basename(path), " in ", round(payload$elapsed_seconds), " s")
    TRUE
  }
  ans <- parallel::mclapply(seq_len(nrow(jobs)), run_job,
    mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE)
  if (!all(vapply(ans, identical, logical(1), TRUE))) stop("At least one simulation chunk failed.")
  message("Completed all 12 unique designs: 12,000 replications and 72,000 fitted models.")
} else if (mode == "--summarize") {
  artifacts <- file.path(stage, "artifacts")
  dir.create(artifacts, recursive = TRUE, showWarnings = FALSE)
  manifests <- list()
  summarize_case <- function(case, data) {
    e <- load_case(case)
    relative <- switch(case$kind, main = "outputs", k = "simulations/k_sensitivity/outputs",
      diag = "simulations/inference_diagnostics/outputs")
    e$out_dir <- file.path(artifacts, relative)
    dir.create(e$out_dir, recursive = TRUE, showWarnings = FALSE)
    e$results <- list(data)
    eval(parse(text = e$summary_code), e)
    key <- c("sim", "model", "method", if (case$kind == "diag") "vcov_type" else "se_type", "x")
    stopifnot(!anyDuplicated(e$results[key]), identical(sort(unique(e$results$sim)), seq_len(1000L)))
    invisible(NULL)
  }
  for (i in seq_len(nrow(cases))) {
    case <- cases[i, ]
    chunks <- lapply(seq_len(8L), function(j) read_chunk(case, j))
    data <- do.call(rbind, lapply(chunks, `[[`, "data"))
    summarize_case(case, data)
    manifests[[i]] <- do.call(rbind, lapply(seq_len(8L), function(j) {
      z <- chunks[[j]]
      data.frame(case, chunk = j, sim_first = min(z$meta$ids), sim_last = max(z$meta$ids),
        M = 1000L, J = 50L, B = 499L, rng_scheme = z$meta$rng_scheme,
        master_seed = z$meta$master_seed, design_seed = z$design_seed,
        R = z$meta$R, mgcv = z$meta$mgcv, started = z$started, finished = z$finished,
        elapsed_seconds = z$elapsed_seconds,
        chunk_md5 = unname(tools::md5sum(chunk_path(case, j))), status = "complete")
    }))
    if (case$kind == "diag" && case$error == "hetero_ar1") {
      shared <- data[data$vcov_type %in% c("classic_cluster", "penalty_cluster"), ]
      shared$se_type <- sub("_cluster", "", shared$vcov_type)
      shared$vcov_type <- NULL
      for (kind in c("main", "k")) {derived <- case; derived$kind <- kind; summarize_case(derived, shared)}
    }
    rm(data, chunks); gc(verbose = FALSE)
  }
  dir.create(file.path(artifacts, "replication"), showWarnings = FALSE)
  write.csv(do.call(rbind, manifests), file.path(artifacts, "replication/paired_g_manifest.csv"), row.names = FALSE)
  write.csv(data.frame(path = names(hashes), md5 = unname(hashes)),
    file.path(artifacts, "replication/paired_g_code_md5.csv"), row.names = FALSE)
  files <- list.files(artifacts, recursive = TRUE, full.names = FALSE)
  write.csv(data.frame(path = files, md5 = unname(tools::md5sum(file.path(artifacts, files)))),
    file.path(stage, "artifact_md5.csv"), row.names = FALSE)
  message("Staged all 18 runs without changing production outputs.")
} else if (mode == "--promote") {
  artifacts <- file.path(stage, "artifacts")
  manifest <- read.csv(file.path(stage, "artifact_md5.csv"))
  stopifnot(identical(unname(tools::md5sum(file.path(artifacts, manifest$path))), manifest$md5))
  backup <- file.path(root, "tmp", "g_coverage_before_paired_20260903")
  for (p in manifest$path) {
    dest <- file.path(root, p); prior <- file.path(backup, p)
    if (file.exists(dest) && !file.exists(prior)) {
      dir.create(dirname(prior), recursive = TRUE, showWarnings = FALSE)
      stopifnot(file.copy(dest, prior))
    }
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(file.path(artifacts, p), dest, overwrite = TRUE))
  }
  message("Promoted validated staged outputs; prior files retained at ", backup)
} else stop("Unknown mode: ", mode)
