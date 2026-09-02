#!/usr/bin/env Rscript

# Reproducible computation benchmarks for the three bam estimators reported in
# the manuscript. Run from the project root with
#
#   Rscript simulations/computational_benchmarks/run_benchmarks.R
#
# The driver launches each method/sample-size combination in a fresh R process
# under /usr/bin/time -l. The worker records model-fitting and penalty-adjusted
# cluster-covariance times; the driver extracts peak resident memory.

suppressPackageStartupMessages(library(mgcv))

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
benchmark_dir <- dirname(script_path)
project_dir <- dirname(dirname(benchmark_dir))
paper_dir <- file.path(project_dir, "paper")
output_dir <- file.path(benchmark_dir, "outputs")
log_dir <- file.path(benchmark_dir, "logs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

arg_value <- function(flag, default = NULL) {
  pos <- match(flag, args)
  if (is.na(pos)) return(default)
  if (pos == length(args)) stop("Missing value after ", flag)
  args[[pos + 1L]]
}

normalized_error_variance <- function(index, sigma_e = 0.35) {
  raw_var <- 1 + 2 * exp(0.75 * index)
  sigma_e^2 * raw_var / mean(raw_var)
}

draw_panel_errors <- function(unit_id, t_periods, x, rho = 0.5) {
  error_var <- normalized_error_variance(x)
  epsilon <- numeric(length(unit_id))
  for (i in unique(unit_id)) {
    idx <- which(unit_id == i)
    ti <- length(idx)
    if (ti != t_periods) stop("Each unit must have T observations.")
    ar_corr <- outer(seq_len(ti), seq_len(ti), function(s, t) rho^abs(s - t))
    sd_i <- sqrt(error_var[idx])
    omega <- ar_corr * outer(sd_i, sd_i)
    epsilon[idx] <- as.numeric(t(chol(omega)) %*% rnorm(ti))
  }
  epsilon
}

simulate_panel <- function(n_units, t_periods) {
  unit_id <- rep(seq_len(n_units), each = t_periods)
  time <- rep(seq_len(t_periods), times = n_units)
  alpha_i <- rbeta(n_units, shape1 = 2.2, shape2 = 2.2)
  mu_i <- 1.1 * (alpha_i - mean(alpha_i)) + rnorm(n_units, sd = 1.25)
  x <- numeric(n_units * t_periods)
  for (i in seq_len(n_units)) {
    idx <- unit_id == i
    x[idx] <- pmin(pmax(alpha_i[i] + rnorm(t_periods, sd = 0.18), 0), 1)
  }
  epsilon <- draw_panel_errors(unit_id, t_periods, x)
  data.frame(
    unit_id = unit_id,
    unit = factor(unit_id),
    time = time,
    x = x,
    y = sin(2 * pi * x) + mu_i[unit_id] + epsilon
  )
}

make_first_differences <- function(df) {
  df <- df[order(df$unit_id, df$time), ]
  pieces <- lapply(split(df, df$unit_id), function(d) {
    data.frame(
      unit_id = d$unit_id[-1L],
      dy = diff(d$y),
      x_now = d$x[-1L],
      x_lag = d$x[-nrow(d)]
    )
  })
  do.call(rbind, pieces)
}

fit_method <- function(df, method, K = 20L) {
  if (method == "FE") {
    fit <- bam(
      y ~ s(x, k = K) + unit,
      data = df,
      method = "fREML",
      select = TRUE,
      discrete = TRUE
    )
    return(list(fit = fit, cluster = df$unit_id))
  }

  if (method == "RE") {
    fit <- bam(
      y ~ s(x, k = K) + s(unit, bs = "re"),
      data = df,
      method = "fREML",
      select = TRUE,
      discrete = TRUE
    )
    return(list(fit = fit, cluster = df$unit_id))
  }

  if (method == "FD") {
    d <- make_first_differences(df)
    X_pair <- cbind(d$x_now, d$x_lag)
    L_pair <- cbind(rep(1, nrow(d)), rep(-1, nrow(d)))
    fit <- bam(
      dy ~ s(X_pair, by = L_pair, k = K),
      data = list(dy = d$dy, X_pair = X_pair, L_pair = L_pair),
      method = "fREML",
      select = TRUE,
      discrete = TRUE
    )
    return(list(fit = fit, cluster = d$unit_id))
  }

  stop("Unknown method: ", method)
}

is_unit_effect_coef <- function(x) {
  startsWith(x, "unit") | startsWith(x, "s(unit)")
}

cluster_vcov_penalty <- function(fit, cluster) {
  R <- predict(fit, type = "lpmatrix")
  u <- residuals(fit, type = "response")
  scale <- summary(fit)$scale
  bread_inv <- vcov(fit, freq = FALSE) / scale
  score_g <- rowsum(
    R * as.numeric(u),
    group = as.factor(cluster),
    reorder = FALSE
  )
  meat <- crossprod(score_g)
  coef_names <- names(fit$edf)
  keep <- !is_unit_effect_coef(coef_names)
  df_used <- sum(fit$edf[keep])
  G <- nrow(score_g)
  N <- nrow(R)
  correction <- (G / (G - 1)) * ((N - 1) / max(N - df_used, 1))
  correction * bread_inv %*% meat %*% bread_inv
}

run_worker <- function() {
  method <- arg_value("--method")
  n_units <- as.integer(arg_value("--n"))
  t_periods <- as.integer(arg_value("--T"))
  repetitions <- as.integer(arg_value("--repetitions", "10"))
  seed <- as.integer(arg_value("--seed", "20260825"))
  output_file <- arg_value("--output")
  if (anyNA(c(n_units, t_periods, repetitions, seed)) || is.null(method)) {
    stop("Invalid worker arguments.")
  }

  set.seed(seed)
  fit_seconds <- numeric(repetitions)
  covariance_seconds <- numeric(repetitions)

  for (b in seq_len(repetitions)) {
    df <- simulate_panel(n_units, t_periods)
    gc(FALSE)
    start_fit <- proc.time()[["elapsed"]]
    fit_info <- fit_method(df, method)
    fit_seconds[b] <- proc.time()[["elapsed"]] - start_fit

    start_covariance <- proc.time()[["elapsed"]]
    V <- cluster_vcov_penalty(fit_info$fit, fit_info$cluster)
    covariance_seconds[b] <- proc.time()[["elapsed"]] - start_covariance
    if (any(!is.finite(diag(V)))) stop("Non-finite covariance benchmark result.")
    rm(V, fit_info, df)
  }

  result <- data.frame(
    method = method,
    n = n_units,
    T = t_periods,
    observations = n_units * t_periods,
    repetitions = repetitions,
    median_fit_seconds = median(fit_seconds),
    q25_fit_seconds = unname(quantile(fit_seconds, 0.25)),
    q75_fit_seconds = unname(quantile(fit_seconds, 0.75)),
    median_covariance_seconds = median(covariance_seconds),
    q25_covariance_seconds = unname(quantile(covariance_seconds, 0.25)),
    q75_covariance_seconds = unname(quantile(covariance_seconds, 0.75))
  )
  write.csv(result, output_file, row.names = FALSE)
}

extract_peak_rss <- function(time_file) {
  lines <- readLines(time_file, warn = FALSE)
  rss_line <- grep("maximum resident set size", lines, value = TRUE)
  if (length(rss_line) != 1L) {
    stop(
      "Could not extract peak resident memory from ", time_file,
      ". Run the driver in an environment where /usr/bin/time -l can read ",
      "process statistics."
    )
  }
  rss_bytes <- as.numeric(sub("maximum resident set size.*$", "", trimws(rss_line)))
  rss_bytes / 1024^2
}

latex_method <- function(x) paste0("\\code{bam} ", x)

write_latex_table <- function(results, metadata) {
  rows <- character()
  configurations <- unique(results[c("n", "T")])
  for (i in seq_len(nrow(configurations))) {
    n_i <- configurations$n[i]
    T_i <- configurations$T[i]
    d <- results[results$n == n_i & results$T == T_i, ]
    rows <- c(
      rows,
      sprintf(
        "\\multicolumn{4}{l}{\\textit{Panel: $n=%d$, $T=%d$ (%s observations)}} \\\\",
        n_i, T_i, format(n_i * T_i, big.mark = ",", scientific = FALSE)
      ),
      "\\midrule"
    )
    for (j in seq_len(nrow(d))) {
      rows <- c(
        rows,
        sprintf(
          "%s & %.3f & %.3f & %.1f \\\\",
          latex_method(d$method[j]),
          d$median_fit_seconds[j],
          d$median_covariance_seconds[j],
          d$peak_rss_mb[j]
        )
      )
    }
    if (i < nrow(configurations)) rows <- c(rows, "\\midrule")
  }

  tex <- c(
    "% Generated by simulations/computational_benchmarks/run_benchmarks.R.",
    "% Do not edit numerical entries by hand.",
    "",
    "\\begin{table}[H]",
    "\\centering",
    "\\caption{Computation time and memory use}",
    "\\label{tab:computation_benchmarks}",
    "\\begin{tabular}{lrrr}",
    "\\toprule",
    "Method & Median fit time (s) & Median covariance time (s) & Peak RSS (MB) \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{minipage}{0.94\\textwidth}",
    paste0(
      "\\footnotesize \\tablenote\\ Each entry is based on 10 independently simulated ",
      "panels from the $g(x)=\\sin(2\\pi x)$ heteroskedastic AR(1) design with ",
      "$\\rho=0.5$ and $K=20$. Fit and covariance columns report median elapsed ",
      "(wall-clock) time across the 10 draws. Covariance time covers construction of the penalty-adjusted ",
      "unit-cluster covariance matrix. Peak RSS is the maximum resident memory of ",
      "the fresh R process that performs the 10 repetitions and therefore includes ",
      "the R session and loaded packages. Benchmarks use a single R process with ",
      "default threading. The factor-FE column uses direct unit indicators in ",
      "\\code{bam} and does not benchmark specialized high-dimensional fixed-effect ",
      "absorption routines. Environment: ", metadata, "."
    ),
    "\\end{minipage}",
    "\\end{table}"
  )
  writeLines(tex, file.path(paper_dir, "computation_benchmarks.tex"))
}

run_driver <- function() {
  tasks <- expand.grid(
    method = c("FD", "FE", "RE"),
    configuration = c("n200_T4", "n200_T8", "n500_T4"),
    stringsAsFactors = FALSE
  )
  config <- list(
    n200_T4 = c(n = 200L, T = 4L),
    n200_T8 = c(n = 200L, T = 8L),
    n500_T4 = c(n = 500L, T = 4L)
  )

  result_list <- vector("list", nrow(tasks))
  rscript <- Sys.which("Rscript")
  if (!nzchar(rscript)) stop("Rscript not found.")

  for (i in seq_len(nrow(tasks))) {
    method <- tasks$method[i]
    cfg <- config[[tasks$configuration[i]]]
    stem <- sprintf("%s_n%d_T%d", method, cfg[["n"]], cfg[["T"]])
    result_file <- file.path(output_dir, paste0(stem, ".csv"))
    stdout_file <- file.path(log_dir, paste0(stem, ".out"))
    time_file <- file.path(log_dir, paste0(stem, ".time"))

    status <- system2(
      "/usr/bin/time",
      args = c(
        "-l", rscript, script_path, "--worker",
        "--method", method,
        "--n", cfg[["n"]],
        "--T", cfg[["T"]],
        "--repetitions", "10",
        "--seed", 20260825L + i,
        "--output", result_file
      ),
      stdout = stdout_file,
      stderr = time_file
    )
    if (status != 0L) {
      stop("Benchmark worker failed; see ", stdout_file, " and ", time_file)
    }
    d <- read.csv(result_file)
    d$peak_rss_mb <- extract_peak_rss(time_file)
    result_list[[i]] <- d
    message("Completed ", stem)
  }

  results <- do.call(rbind, result_list)
  method_order <- c("FD", "FE", "RE")
  results$method <- factor(results$method, levels = method_order)
  results <- results[order(results$n, results$T, results$method), ]
  results$method <- as.character(results$method)
  write.csv(
    results,
    file.path(output_dir, "computation_benchmarks.csv"),
    row.names = FALSE
  )

  hardware <- system2(
    "system_profiler",
    c("SPHardwareDataType", "-detailLevel", "mini"),
    stdout = TRUE,
    stderr = FALSE
  )
  chip <- trimws(sub("^.*Chip: ", "", grep("Chip:", hardware, value = TRUE)[1]))
  cores <- trimws(sub(
    "^.*Total Number of Cores: ", "",
    grep("Total Number of Cores:", hardware, value = TRUE)[1]
  ))
  cores <- sub("^([0-9]+) ", "\\1 cores ", cores)
  memory <- trimws(sub("^.*Memory: ", "", grep("Memory:", hardware, value = TRUE)[1]))
  macos <- system2("sw_vers", "-productVersion", stdout = TRUE)
  metadata <- paste0(
    chip, ", ", cores, ", ", memory, ", macOS ", macos,
    ", ", R.version.string, ", mgcv ", as.character(packageVersion("mgcv"))
  )
  writeLines(metadata, file.path(output_dir, "benchmark_environment.txt"))
  write_latex_table(results, metadata)

  message("Combined results: ", file.path(output_dir, "computation_benchmarks.csv"))
  message("LaTeX table: ", file.path(paper_dir, "computation_benchmarks.tex"))
}

if ("--worker" %in% args) {
  run_worker()
} else if ("--tables-only" %in% args) {
  results <- read.csv(file.path(output_dir, "computation_benchmarks.csv"))
  metadata <- readLines(file.path(output_dir, "benchmark_environment.txt"), warn = FALSE)
  write_latex_table(results, paste(metadata, collapse = " "))
  message("Regenerated the benchmark table from saved results.")
} else {
  run_driver()
}
