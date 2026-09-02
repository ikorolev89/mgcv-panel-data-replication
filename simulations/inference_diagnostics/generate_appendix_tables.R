#!/usr/bin/env Rscript

# Generate the four appendix tables that compare the five inference methods.
# Run from the project root after combine_inference_diagnostics.R.

beta_path <- file.path(
  "simulations", "inference_diagnostics", "outputs", "combined",
  "beta_diagnostics_all_M1000.csv"
)
g_path <- file.path(
  "simulations", "inference_diagnostics", "outputs", "combined",
  "g_diagnostics_all_M1000.csv"
)
output_path <- file.path("paper", "inference_method_appendix.tex")

beta <- read.csv(beta_path, stringsAsFactors = FALSE)
g <- read.csv(g_path, stringsAsFactors = FALSE)

method_order <- c("bam_fd", "bam_fe", "bam_re")
method_label <- c(
  bam_fd = "\\code{bam} FD",
  bam_fe = "\\code{bam} FE",
  bam_re = "\\code{bam} RE"
)

vcov_order <- c(
  "classic_cluster",
  "penalty_cluster",
  "penalty_cluster_delta",
  "mgcv_sandwich_freq",
  "mgcv_sandwich_default"
)

designs <- data.frame(
  panel = LETTERS[1:6],
  error_design = rep(c("hetero_ar1", "hetero"), each = 3),
  n = rep(c(200L, 200L, 500L), 2),
  T = rep(c(4L, 8L, 4L), 2),
  stringsAsFactors = FALSE
)

design_label <- function(error_design, n, T) {
  error_text <- if (error_design == "hetero_ar1") {
    "heteroskedastic AR(1), $\\rho=0.5$"
  } else {
    "heteroskedastic serially independent errors"
  }
  sprintf("%s; $n=%d$, $T=%d$", error_text, n, T)
}

fmt <- function(x) sprintf("%.3f", x)

extract_one <- function(data, filters, variable) {
  keep <- rep(TRUE, nrow(data))
  for (name in names(filters)) {
    keep <- keep & data[[name]] == filters[[name]]
  }
  out <- data[keep, variable]
  if (length(out) != 1L) {
    stop("Expected one cell but found ", length(out), ": ",
         paste(names(filters), filters, sep = "=", collapse = ", "))
  }
  out
}

beta_table <- function(beta_true, caption, label, mc_note) {
  lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\renewcommand{\\arraystretch}{0.88}",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\begin{tabular}{lrrrrr}",
    "\\toprule",
    "Method & Classic & Penalty & Penalty $+\\widehat{\\Delta}$ & \\code{mgcv}-F & \\code{mgcv}-B \\\\",
    "\\midrule"
  )

  for (i in seq_len(nrow(designs))) {
    d <- designs[i, ]
    lines <- c(
      lines,
      sprintf(
        "\\multicolumn{6}{l}{\\textit{Panel %s: %s}} \\\\",
        d$panel, design_label(d$error_design, d$n, d$T)
      ),
      "\\midrule"
    )
    for (method in method_order) {
      values <- vapply(vcov_order, function(vcov_type) {
        extract_one(
          beta,
          list(
            beta_true = beta_true,
            error_design = d$error_design,
            n = d$n,
            T = d$T,
            method = method,
            vcov_type = vcov_type
          ),
          "rejection_frequency"
        )
      }, numeric(1))
      lines <- c(
        lines,
        sprintf("%s & %s \\\\", method_label[[method]], paste(fmt(values), collapse = " & "))
      )
    }
    if (i < nrow(designs)) lines <- c(lines, "\\midrule")
  }

  c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{minipage}{0.94\\textwidth}",
    "\\singlespacing",
    sprintf("\\scriptsize \\tablenote\\ Entries report empirical rejection frequencies based on 1,000 Monte Carlo replications. Classic is the unit-cluster sandwich that ignores the penalty; Penalty is the penalty-adjusted unit-cluster sandwich; and Penalty $+\\widehat{\\Delta}$ adds $\\widehat{\\Delta}=\\code{vcov}(\\widehat m,\\code{freq=FALSE})-\\code{vcov}(\\widehat m,\\code{freq=TRUE})$. The columns \\code{mgcv}-F and \\code{mgcv}-B use \\code{vcov}$(\\widehat m,\\code{sandwich=TRUE})$ with \\code{freq=TRUE} and \\code{freq=FALSE}, respectively; these are observation-level rather than unit-clustered sandwich estimators. %s", mc_note),
    "\\end{minipage}",
    "\\end{table}"
  )
}

g_table <- function(model, caption, label) {
  lines <- c(
    "\\begin{landscape}",
    "\\begin{table}[H]",
    "\\centering",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{3.5pt}",
    "\\renewcommand{\\arraystretch}{0.75}",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\begin{tabular}{lrrrrrrrrrr}",
    "\\toprule",
    "& \\multicolumn{2}{c}{Classic} & \\multicolumn{2}{c}{Penalty} & \\multicolumn{2}{c}{Penalty $+\\widehat{\\Delta}$} & \\multicolumn{2}{c}{\\code{mgcv}-F} & \\multicolumn{2}{c}{\\code{mgcv}-B} \\\\",
    "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9} \\cmidrule(lr){10-11}",
    "Method & Avg. point. & Uniform & Avg. point. & Uniform & Avg. point. & Uniform & Avg. point. & Uniform & Avg. point. & Uniform \\\\",
    "\\midrule"
  )

  for (i in seq_len(nrow(designs))) {
    d <- designs[i, ]
    lines <- c(
      lines,
      sprintf(
        "\\multicolumn{11}{l}{\\textit{Panel %s: %s}} \\\\",
        d$panel, design_label(d$error_design, d$n, d$T)
      ),
      "\\midrule"
    )
    for (method in method_order) {
      values <- unlist(lapply(vcov_order, function(vcov_type) {
        filters <- list(
          model = model,
          error_design = d$error_design,
          n = d$n,
          T = d$T,
          method = method,
          vcov_type = vcov_type
        )
        c(
          extract_one(g, filters, "avg_pointwise_coverage"),
          extract_one(g, filters, "uniform_coverage")
        )
      }), use.names = FALSE)
      lines <- c(
        lines,
        sprintf("%s & %s \\\\", method_label[[method]], paste(fmt(values), collapse = " & "))
      )
    }
    if (i < nrow(designs)) lines <- c(lines, "\\midrule")
  }

  c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{minipage}{0.92\\linewidth}",
    "\\singlespacing",
    "\\tiny \\tablenote\\ Entries report coverage probabilities based on 1,000 Monte Carlo replications. ``Avg. point.'' averages coverage over the 50 evaluation points; ``Uniform'' is coverage of the entire curve by a sup-$t$ band based on 499 Gaussian draws. The five inference methods are defined in the notes to Appendix Table~\\ref{tab:beta_size_methods_appendix}. The Monte Carlo standard error near 0.95 is approximately 0.007.",
    "\\end{minipage}",
    "\\end{table}",
    "\\end{landscape}"
  )
}

output <- c(
  "% This file is generated by simulations/inference_diagnostics/generate_appendix_tables.R.",
  "% Do not edit the table entries by hand.",
  "",
  beta_table(
    beta_true = 1,
    caption = "Size results across inference methods for tests of $H_0:\\beta=1$",
    label = "tab:beta_size_methods_appendix",
    mc_note = "The Monte Carlo standard error of a rejection frequency near 0.05 is approximately 0.007."
  ),
  "",
  beta_table(
    beta_true = 1.075,
    caption = "Power results across inference methods for tests of $H_0:\\beta=1$",
    label = "tab:beta_power_methods_appendix",
    mc_note = "The Monte Carlo standard error is at most approximately 0.016."
  ),
  "",
  g_table(
    model = "pl",
    caption = "Coverage results across inference methods for centered $g$: partially linear model",
    label = "tab:g_coverage_pl_methods_appendix"
  ),
  "",
  g_table(
    model = "np",
    caption = "Coverage results across inference methods for centered $g$: nonparametric model",
    label = "tab:g_coverage_np_methods_appendix"
  )
)

writeLines(output, output_path)
message("Wrote ", output_path)
