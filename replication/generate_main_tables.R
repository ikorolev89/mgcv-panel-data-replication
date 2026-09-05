#!/usr/bin/env Rscript
# Regenerate Tables 1--5 from the saved main-run summaries.
# --check compares numerical cells with the revision snapshot without writing.
check <- "--check" %in% commandArgs(trailingOnly = TRUE)
refresh <- "--refresh-coverage-snapshot" %in% commandArgs(trailingOnly = TRUE)
if (check && refresh) stop("Choose checking or deliberate snapshot refresh, not both.")
configs <- data.frame(n = c(200, 200, 500), T = c(4, 8, 4))
methods <- c("bam_fd", "bam_fe", "bam_re")
labels <- c(bam_fd = "\\code{bam} FD", bam_fe = "\\code{bam} FE", bam_re = "\\code{bam} RE")
read_one <- function(pattern) {
  p <- list.files("outputs", pattern, full.names = TRUE)
  if (length(p) != 1L) stop("Expected one source for ", pattern)
  read.csv(p, check.names = FALSE)
}
tables <- list()
values <- list()
add_table <- function(id, caption, header, body, numeric_rows) {
  tables[[id]] <<- c("\\begin{table}[htbp]", "\\centering", paste0("\\caption{", caption, "}"),
    paste0("\\label{", id, "}"), paste0("\\begin{tabular}{l", paste(rep("r", length(strsplit(header, " & ", fixed = TRUE)[[1]]) - 1L), collapse = ""), "}"),
    "\\toprule", paste0(header, " \\\\"), "\\midrule", body, "\\bottomrule", "\\end{tabular}", "\\end{table}", "")
  values[[id]] <<- data.frame(table = id, row = seq_along(numeric_rows), values = numeric_rows)
}
mse <- lapply(seq_len(nrow(configs)), function(i) read_one(sprintf(
  "^panel_fe_test_mse_summary_M1000_n%d_T%d_K7_.*_hetero_ar1_rho0p5.csv$", configs$n[i], configs$T[i])))
mse_methods <- c("linear FE", "linear FD", "spline FE, fixed K", "FD spline, fixed K", "mgcv::bam FD", "mgcv::bam FE", "mgcv::bam RE")
rows <- cells <- character()
for (dgp in c("linear", "quadratic", "sin_2pi")) for (method in mse_methods) {
  v <- vapply(mse, function(d) {
    z <- d$mean_rmse[d$dgp == dgp & d$method == method]
    if (length(z) != 1) stop("Missing RMSE cell")
    100 * z
  }, numeric(1))
  formatted <- paste(sprintf("%.3f", v), collapse = " & ")
  cells <- c(cells, formatted)
  label <- paste(dgp, method, sep = ": ")
  label <- gsub("_", "\\_", label, fixed = TRUE)
  rows <- c(rows, paste0(label, " & ", formatted, " \\\\"))
}
add_table("tab:estimation_rmse", "Out-of-sample RMSE for estimating $g$ ($\\times 100$)",
  "DGP / Method & $n=200,T=4$ & $n=200,T=8$ & $n=500,T=4$", rows, cells)
for (truth in c("1", "1p075")) {
  rows <- cells <- character()
  for (i in seq_len(nrow(configs))) {
    d <- read_one(sprintf("^panel_fe_beta_inference_summary_M1000_B0_n%d_T%d_baseline_btrue%s_bnull1_hetero_ar1_rho0p5_dfsubstantive.csv$", configs$n[i], configs$T[i], truth))
    rows <- c(rows, sprintf("\\multicolumn{6}{l}{$n=%d$, $T=%d$} \\\\", configs$n[i], configs$T[i]))
    for (method in methods) {
      z <- d[d$method == method, ]; stopifnot(nrow(z) == 1L)
      formatted <- paste(c(sprintf("%.3f", unlist(z[c("rejection_classic_cluster", "rejection_penalty_cluster")])),
        sprintf("%.4f", unlist(z[c("mean_beta_hat", "sd_beta_hat", "mean_se_penalty_cluster")]))), collapse = " & ")
      cells <- c(cells, formatted)
      rows <- c(rows, paste0(labels[method], " & ", formatted, " \\\\"))
    }
  }
  id <- if (truth == "1") "tab:beta_size" else "tab:beta_power"
  add_table(id, if (truth == "1") "Size of tests of $H_0: \\beta=1$" else "Power of tests of $H_0: \\beta=1$ at $\\beta=1.075$",
    "Method & Classic & Penalty & Mean $\\hat\\beta$ & SD $\\hat\\beta$ & Mean SE", rows, cells)
}
for (model in c("pl", "np")) {
  rows <- cells <- character()
  for (i in seq_len(nrow(configs))) {
    d <- read_one(sprintf("^panel_fe_g_inference_summary_M1000_n%d_T%d_sin_2pi_np_pl_J50_hetero_ar1_rho0p5_dfsubstantive.csv$", configs$n[i], configs$T[i]))
    rows <- c(rows, sprintf("\\multicolumn{5}{l}{$n=%d$, $T=%d$} \\\\", configs$n[i], configs$T[i]))
    for (method in methods) {
      z <- d[d$model == model & d$method == method & d$se_type == "penalty", ]; stopifnot(nrow(z) == 1L)
      formatted <- paste(c(sprintf("%.3f", unlist(z[c("avg_pointwise_coverage", "min_pointwise_coverage", "uniform_coverage")])), sprintf("%.4f", z$avg_pointwise_width)), collapse = " & ")
      cells <- c(cells, formatted)
      rows <- c(rows, paste0(labels[method], " & ", formatted, " \\\\"))
    }
  }
  add_table(paste0("tab:g_coverage_", model, "_penalty"), paste0("Penalty-adjusted coverage for centered $g$: ", model, " model"),
    "Method & Avg. pointwise & Min. pointwise & Uniform & Avg. width", rows, cells)
}
actual <- do.call(rbind, values); row.names(actual) <- NULL
reference <- read.csv("replication/paper_values.csv", check.names = FALSE)
if (refresh) {
  unaffected <- !grepl("^tab:g_coverage_", actual$table)
  stopifnot(identical(actual[unaffected, ], reference[unaffected, ]))
  write.csv(actual, "replication/paper_values.csv", row.names = FALSE)
  reference <- actual
}
if (!identical(actual, reference)) stop("Generated numerical cells differ from the September 2026 revision snapshot.")
if (!check) {
  dir.create("paper", showWarnings = FALSE)
  writeLines(c("% Generated by replication/generate_main_tables.R; requires booktabs.", unlist(tables)), "paper/main_simulation_tables.tex")
}
cat("All 225 numerical cells in Tables 1--5 match the revision.\n")
