#!/usr/bin/env Rscript

# Manuscript-ready EmplUK application.
#
# This script estimates a linear two-way fixed-effects model, partially linear
# factor-FE and first-difference models with a smooth output effect, and an
# additive factor-FE GAM. It constructs penalty-adjusted firm-clustered
# inference and writes the tables and figures included by the manuscript.
#
# Run from the project root with:
#   Rscript empirical_application_empluk/run_manuscript_application.R

suppressPackageStartupMessages({
  library(plm)
  library(mgcv)
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
app_dir <- dirname(script_path)
project_dir <- dirname(app_dir)
paper_dir <- file.path(project_dir, "paper")
output_dir <- file.path(app_dir, "manuscript_outputs")
figure_dir <- file.path(app_dir, "manuscript_figures")
paper_figure_dir <- file.path(paper_dir, "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_figure_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(1991277)

data("EmplUK", package = "plm", envir = environment())
raw_df <- as.data.frame(EmplUK)
raw_df <- raw_df[order(raw_df$firm, raw_df$year), ]

df <- transform(
  raw_df,
  log_emp = log(emp),
  log_wage = log(wage),
  log_capital = log(capital),
  log_output = log(output),
  firm_fe = factor(firm),
  year_fe = factor(year)
)

if (anyNA(df)) stop("EmplUK unexpectedly contains missing values.")
if (nrow(df) != 1031L || nlevels(df$firm_fe) != 140L) {
  stop("Unexpected EmplUK sample dimensions.")
}

lag_within_firm <- function(x, firm) {
  ave(x, firm, FUN = function(z) c(NA, head(z, -1L)))
}

df$year_lag <- lag_within_firm(df$year, df$firm)
for (variable in c("log_emp", "log_wage", "log_capital", "log_output")) {
  df[[paste0(variable, "_lag")]] <- lag_within_firm(df[[variable]], df$firm)
}

df_fd <- subset(df, !is.na(year_lag) & year - year_lag == 1)
df_fd$d_log_emp <- df_fd$log_emp - df_fd$log_emp_lag
df_fd$d_log_wage <- df_fd$log_wage - df_fd$log_wage_lag
df_fd$d_log_capital <- df_fd$log_capital - df_fd$log_capital_lag
df_fd$year_fd <- factor(df_fd$year)

if (nrow(df_fd) != 891L || nlevels(droplevels(df_fd$firm_fe)) != 140L) {
  stop("Unexpected first-difference sample dimensions.")
}

K_main <- 20L
covariates <- c("log_wage", "log_capital", "log_output")
variable_labels <- c(
  log_wage = "Log wages",
  log_capital = "Log capital",
  log_output = "Log output"
)

linear_formula <-
  log_emp ~ log_wage + log_capital + log_output + firm_fe + year_fe
pl_formula <-
  log_emp ~ log_wage + log_capital +
    s(log_output, bs = "cr", k = K_main) + firm_fe + year_fe
gam_formula <-
  log_emp ~ s(log_wage, bs = "cr", k = K_main) +
    s(log_capital, bs = "cr", k = K_main) +
    s(log_output, bs = "cr", k = K_main) + firm_fe + year_fe

output_pair <- cbind(df_fd$log_output, df_fd$log_output_lag)
output_weights <- cbind(
  rep(1, nrow(df_fd)),
  rep(-1, nrow(df_fd))
)
pl_fd_formula <-
  d_log_emp ~ d_log_wage + d_log_capital + year_fd +
    s(output_pair, by = output_weights, bs = "cr", k = K_main)

linear_fe <- lm(linear_formula, data = df)
pl_fe <- mgcv::bam(
  pl_formula,
  data = df,
  method = "fREML",
  discrete = TRUE,
  select = TRUE
)
gam_fe <- mgcv::bam(
  gam_formula,
  data = df,
  method = "fREML",
  discrete = TRUE,
  select = TRUE
)
pl_fd <- mgcv::bam(
  pl_fd_formula,
  data = list(
    d_log_emp = df_fd$d_log_emp,
    d_log_wage = df_fd$d_log_wage,
    d_log_capital = df_fd$d_log_capital,
    year_fd = df_fd$year_fd,
    output_pair = output_pair,
    output_weights = output_weights
  ),
  method = "fREML",
  discrete = TRUE,
  select = TRUE
)

fit_sizes <- c(
  linear_fe = stats::nobs(linear_fe),
  pl_fe = stats::nobs(pl_fe),
  gam_fe = stats::nobs(gam_fe),
  pl_fd = stats::nobs(pl_fd)
)
if (any(fit_sizes[c("linear_fe", "pl_fe", "gam_fe")] != nrow(df)) ||
    fit_sizes[["pl_fd"]] != nrow(df_fd)) {
  stop("Unexpected empirical estimation samples: ",
       paste(names(fit_sizes), fit_sizes, sep = "=", collapse = ", "))
}

symmetrize <- function(A) (A + t(A)) / 2
rmse <- function(actual, fitted) sqrt(mean((actual - fitted)^2))

is_unit_effect_coef <- function(coef_names, unit_var) {
  startsWith(coef_names, unit_var) |
    startsWith(coef_names, paste0("s(", unit_var, ")."))
}

substantive_edf <- function(fit, unit_var) {
  coef_names <- names(fit$edf)
  sum(fit$edf[!is_unit_effect_coef(coef_names, unit_var)])
}

cluster_vcov_penalty <- function(fit, cluster, unit_var = "firm_fe") {
  R <- predict(fit, type = "lpmatrix")
  u <- residuals(fit, type = "response")
  bread_inv <- vcov(fit, freq = FALSE) / summary(fit)$scale
  score_g <- rowsum(
    R * as.numeric(u),
    group = as.factor(cluster),
    reorder = FALSE
  )
  meat <- crossprod(score_g)
  G <- nrow(score_g)
  N <- nrow(R)
  d <- substantive_edf(fit, unit_var)
  correction <- (G / (G - 1)) * ((N - 1) / max(N - d, 1))
  symmetrize(correction * bread_inv %*% meat %*% bread_inv)
}

cluster_vcov_linear <- function(fit, cluster, unit_var = "firm_fe") {
  R <- model.matrix(fit)
  u <- residuals(fit)
  bread_inv <- solve(crossprod(R))
  score_g <- rowsum(
    R * as.numeric(u),
    group = as.factor(cluster),
    reorder = FALSE
  )
  meat <- crossprod(score_g)
  G <- nrow(score_g)
  N <- nrow(R)
  keep <- !is_unit_effect_coef(colnames(R), unit_var)
  d <- qr(R[, keep, drop = FALSE])$rank
  correction <- (G / (G - 1)) * ((N - 1) / max(N - d, 1))
  symmetrize(correction * bread_inv %*% meat %*% bread_inv)
}

linear_vcov <- cluster_vcov_linear(linear_fe, df$firm_fe)
pl_vcov <- cluster_vcov_penalty(pl_fe, df$firm_fe)
gam_vcov <- cluster_vcov_penalty(gam_fe, df$firm_fe)
pl_fd_vcov <- cluster_vcov_penalty(pl_fd, df_fd$firm_fe)

coef_rows <- c("log_wage", "log_capital", "log_output")
pl_available <- coef_rows %in% names(coef(pl_fe))
coefficient_results <- data.frame(
  variable = coef_rows,
  label = unname(variable_labels[coef_rows]),
  linear_estimate = coef(linear_fe)[coef_rows],
  linear_cluster_se = sqrt(diag(linear_vcov))[coef_rows],
  pl_estimate = NA_real_,
  pl_cluster_se = NA_real_,
  fd_estimate = NA_real_,
  fd_cluster_se = NA_real_,
  stringsAsFactors = FALSE
)
coefficient_results$pl_estimate[pl_available] <-
  coef(pl_fe)[coef_rows[pl_available]]
coefficient_results$pl_cluster_se[pl_available] <-
  sqrt(diag(pl_vcov))[coef_rows[pl_available]]
fd_coef_rows <- c(
  log_wage = "d_log_wage",
  log_capital = "d_log_capital"
)
fd_available <- names(fd_coef_rows) %in% coefficient_results$variable
fd_positions <- match(names(fd_coef_rows)[fd_available], coefficient_results$variable)
coefficient_results$fd_estimate[fd_positions] <-
  coef(pl_fd)[unname(fd_coef_rows[fd_available])]
coefficient_results$fd_cluster_se[fd_positions] <-
  sqrt(diag(pl_fd_vcov))[unname(fd_coef_rows[fd_available])]
write.csv(
  coefficient_results,
  file.path(output_dir, "coefficient_results.csv"),
  row.names = FALSE
)

smooth_edf <- function(fit, variable) {
  smooth_table <- as.data.frame(summary(fit)$s.table)
  term <- paste0("s(", variable, ")")
  if (!term %in% rownames(smooth_table)) return(1)
  unname(smooth_table[term, "edf"])
}

fd_output_edf <- unname(as.data.frame(summary(pl_fd)$s.table)[1, "edf"])

model_comparison <- data.frame(
  model = c("Linear FE", "Partially linear FE", "Additive GAM FE", "Partially linear FD"),
  wage_edf = c(1, 1, smooth_edf(gam_fe, "log_wage"), 1),
  capital_edf = c(1, 1, smooth_edf(gam_fe, "log_capital"), 1),
  output_edf = c(
    1,
    smooth_edf(pl_fe, "log_output"),
    smooth_edf(gam_fe, "log_output"),
    fd_output_edf
  ),
  model_df = c(
    attr(logLik(linear_fe), "df"),
    attr(logLik(pl_fe), "df"),
    attr(logLik(gam_fe), "df"),
    attr(logLik(pl_fd), "df")
  ),
  RMSE = c(
    rmse(df$log_emp, fitted(linear_fe)),
    rmse(df$log_emp, fitted(pl_fe)),
    rmse(df$log_emp, fitted(gam_fe)),
    rmse(df_fd$d_log_emp, fitted(pl_fd))
  ),
  AIC = c(AIC(linear_fe), AIC(pl_fe), AIC(gam_fe), AIC(pl_fd)),
  BIC = c(BIC(linear_fe), BIC(pl_fe), BIC(gam_fe), BIC(pl_fd)),
  stringsAsFactors = FALSE
)
write.csv(
  model_comparison,
  file.path(output_dir, "model_comparison.csv"),
  row.names = FALSE
)

fit_basis_model <- function(model, K) {
  if (model == "Partially linear FE") {
    formula <- as.formula(sprintf(
      paste0(
        "log_emp ~ log_wage + log_capital + ",
        "s(log_output, bs = 'cr', k = %d) + firm_fe + year_fe"
      ),
      K
    ))
  } else {
    formula <- as.formula(sprintf(
      paste0(
        "log_emp ~ s(log_wage, bs = 'cr', k = %d) + ",
        "s(log_capital, bs = 'cr', k = %d) + ",
        "s(log_output, bs = 'cr', k = %d) + firm_fe + year_fe"
      ),
      K, K, K
    ))
  }
  fit <- mgcv::bam(
    formula,
    data = df,
    method = "fREML",
    discrete = TRUE,
    select = TRUE
  )
  variables <- if (model == "Partially linear FE") {
    "log_output"
  } else {
    covariates
  }
  data.frame(
    model = model,
    K = K,
    variable = variables,
    edf = vapply(variables, function(x) smooth_edf(fit, x), numeric(1)),
    RMSE = rmse(df$log_emp, fitted(fit)),
    AIC = AIC(fit),
    BIC = BIC(fit),
    stringsAsFactors = FALSE
  )
}

basis_sensitivity <- do.call(
  rbind,
  lapply(c("Partially linear FE", "Additive GAM FE"), function(model) {
    do.call(rbind, lapply(c(10L, 20L, 40L), function(K) {
      if (K == K_main) {
        fit <- if (model == "Partially linear FE") pl_fe else gam_fe
        variables <- if (model == "Partially linear FE") {
          "log_output"
        } else {
          covariates
        }
        data.frame(
          model = model,
          K = K,
          variable = variables,
          edf = vapply(
            variables,
            function(x) smooth_edf(fit, x),
            numeric(1)
          ),
          RMSE = rmse(df$log_emp, fitted(fit)),
          AIC = AIC(fit),
          BIC = BIC(fit),
          stringsAsFactors = FALSE
        )
      } else {
        fit_basis_model(model, K)
      }
    }))
  })
)
write.csv(
  basis_sensitivity,
  file.path(output_dir, "basis_sensitivity.csv"),
  row.names = FALSE
)

concurvity_estimate <- as.data.frame(
  mgcv::concurvity(gam_fe, full = FALSE)$estimate
)
concurvity_estimate$term_explained <- rownames(concurvity_estimate)
rownames(concurvity_estimate) <- NULL
write.csv(
  concurvity_estimate,
  file.path(output_dir, "concurvity.csv"),
  row.names = FALSE
)

reference_data <- data.frame(
  log_wage = median(df$log_wage),
  log_capital = median(df$log_capital),
  log_output = median(df$log_output),
  firm_fe = factor(
    levels(df$firm_fe)[1],
    levels = levels(df$firm_fe)
  ),
  year_fe = factor(
    levels(df$year_fe)[1],
    levels = levels(df$year_fe)
  )
)

make_grid <- function(variable, n = 180L) {
  limits <- as.numeric(
    quantile(df[[variable]], c(0.01, 0.99), names = FALSE)
  )
  seq(limits[1], limits[2], length.out = n)
}

smooth_columns <- function(fit, variable) {
  label <- paste0("s(", variable, ")")
  index <- which(vapply(fit$smooth, `[[`, "", "label") == label)
  if (length(index) != 1L) stop("Could not identify smooth ", label)
  sm <- fit$smooth[[index]]
  sm$first.para:sm$last.para
}

smooth_contrast <- function(fit, variable, x_grid) {
  newdata <- reference_data[rep(1L, length(x_grid)), ]
  newdata[[variable]] <- x_grid
  Rg <- predict(fit, newdata = newdata, type = "lpmatrix")
  cols <- smooth_columns(fit, variable)
  L <- Rg
  L[, -cols] <- 0
  sweep(L, 2, colMeans(L), "-")
}

uniform_critical_value <- function(C, draws = 9999L) {
  se <- sqrt(pmax(diag(C), 0))
  valid <- is.finite(se) & se > 1e-10
  Corr <- C[valid, valid, drop = FALSE] / outer(se[valid], se[valid])
  Corr <- symmetrize(Corr)
  eig <- eigen(Corr, symmetric = TRUE)
  A <- eig$vectors %*%
    diag(sqrt(pmax(eig$values, 0)), nrow = length(eig$values))
  z <- matrix(rnorm(draws * nrow(A)), nrow = draws)
  sim <- z %*% t(A)
  simulated_cutoff <- as.numeric(
    quantile(apply(abs(sim), 1, max), 0.95, names = FALSE)
  )
  # The Gaussian maximum's cutoff cannot be below a marginal normal cutoff.
  # Enforce that bound when simulation error matters for a nearly linear smooth.
  max(qnorm(0.975), simulated_cutoff)
}

output_grid <- make_grid("log_output")
L_output <- smooth_contrast(pl_fe, "log_output", output_grid)
output_estimate <- as.numeric(L_output %*% coef(pl_fe))
output_cov <- symmetrize(L_output %*% pl_vcov %*% t(L_output))
output_se <- sqrt(pmax(diag(output_cov), 0))
output_uniform_critical <- uniform_critical_value(output_cov)

output_inference <- data.frame(
  estimator = "Factor fixed effects",
  log_output = output_grid,
  estimate = output_estimate,
  point_low = output_estimate - qnorm(0.975) * output_se,
  point_high = output_estimate + qnorm(0.975) * output_se,
  uniform_low = output_estimate - output_uniform_critical * output_se,
  uniform_high = output_estimate + output_uniform_critical * output_se
)
write.csv(
  output_inference,
  file.path(output_dir, "output_smooth_inference.csv"),
  row.names = FALSE
)

fd_smooth_contrast <- function(x_grid) {
  newdata <- data.frame(
    d_log_wage = rep(0, length(x_grid)),
    d_log_capital = rep(0, length(x_grid)),
    year_fd = factor(
      rep(levels(df_fd$year_fd)[1], length(x_grid)),
      levels = levels(df_fd$year_fd)
    ),
    output_pair = x_grid,
    output_weights = rep(1, length(x_grid))
  )
  Rg <- predict(pl_fd, newdata = newdata, type = "lpmatrix")
  sm <- pl_fd$smooth[[1]]
  cols <- sm$first.para:sm$last.para
  L <- Rg
  L[, -cols] <- 0
  sweep(L, 2, colMeans(L), "-")
}

L_output_fd <- fd_smooth_contrast(output_grid)
output_estimate_fd <- as.numeric(L_output_fd %*% coef(pl_fd))
output_cov_fd <- symmetrize(L_output_fd %*% pl_fd_vcov %*% t(L_output_fd))
output_se_fd <- sqrt(pmax(diag(output_cov_fd), 0))
output_uniform_critical_fd <- uniform_critical_value(output_cov_fd)

output_inference_fd <- data.frame(
  estimator = "First differences",
  log_output = output_grid,
  estimate = output_estimate_fd,
  point_low = output_estimate_fd - qnorm(0.975) * output_se_fd,
  point_high = output_estimate_fd + qnorm(0.975) * output_se_fd,
  uniform_low = output_estimate_fd - output_uniform_critical_fd * output_se_fd,
  uniform_high = output_estimate_fd + output_uniform_critical_fd * output_se_fd
)
write.csv(
  output_inference_fd,
  file.path(output_dir, "output_smooth_inference_fd.csv"),
  row.names = FALSE
)

output_inference_both <- rbind(output_inference, output_inference_fd)
output_inference_both$estimator <- factor(
  output_inference_both$estimator,
  levels = c("Factor fixed effects", "First differences")
)

model_effect <- function(model, variable, x_grid) {
  if (model == "Linear FE") {
    estimate <- coef(linear_fe)[[variable]] * (x_grid - mean(x_grid))
  } else if (model == "Partially linear FE" && variable != "log_output") {
    estimate <- coef(pl_fe)[[variable]] * (x_grid - mean(x_grid))
  } else {
    fit <- if (model == "Partially linear FE") pl_fe else gam_fe
    L <- smooth_contrast(fit, variable, x_grid)
    estimate <- as.numeric(L %*% coef(fit))
  }
  data.frame(
    model = model,
    variable = variable,
    label = unname(variable_labels[[variable]]),
    x = x_grid,
    estimate = estimate,
    stringsAsFactors = FALSE
  )
}

models <- c("Linear FE", "Partially linear FE", "Additive GAM FE")
effect_comparison <- do.call(
  rbind,
  lapply(covariates, function(variable) {
    x_grid <- make_grid(variable)
    do.call(rbind, lapply(models, model_effect, variable = variable, x_grid = x_grid))
  })
)
effect_comparison$model <- factor(effect_comparison$model, levels = models)
effect_comparison$label <- factor(
  effect_comparison$label,
  levels = unname(variable_labels[covariates])
)
write.csv(
  effect_comparison,
  file.path(output_dir, "effect_comparison.csv"),
  row.names = FALSE
)

rug_data <- do.call(rbind, lapply(covariates, function(variable) {
  values <- df[[variable]]
  limits <- range(make_grid(variable))
  values <- values[values >= limits[1] & values <= limits[2]]
  if (length(values) > 350L) values <- sample(values, 350L)
  data.frame(
    variable = variable,
    label = unname(variable_labels[[variable]]),
    x = values
  )
}))
rug_data$label <- factor(
  rug_data$label,
  levels = unname(variable_labels[covariates])
)

paper_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title.position = "plot",
    plot.caption.position = "plot"
  )

comparison_plot <- ggplot(
  effect_comparison,
  aes(x = x, y = estimate, color = model, linetype = model)
) +
  geom_line(linewidth = 0.85) +
  geom_rug(
    data = rug_data,
    aes(x = x),
    inherit.aes = FALSE,
    alpha = 0.10,
    sides = "b"
  ) +
  facet_wrap(~label, scales = "free", ncol = 2) +
  scale_color_manual(values = c(
    "Linear FE" = "#6C757D",
    "Partially linear FE" = "#D55E00",
    "Additive GAM FE" = "#0072B2"
  )) +
  scale_linetype_manual(values = c(
    "Linear FE" = "dotted",
    "Partially linear FE" = "longdash",
    "Additive GAM FE" = "solid"
  )) +
  labs(
    title = "EmplUK: partial effects across specifications",
    subtitle = "Effects are centered over each displayed grid; axes cover the 1st--99th percentiles.",
    x = NULL,
    y = "Partial effect on log employment",
    color = NULL,
    linetype = NULL
  ) +
  paper_theme

band_data <- rbind(
  data.frame(
    estimator = output_inference_both$estimator,
    log_output = output_inference_both$log_output,
    low = output_inference_both$uniform_low,
    high = output_inference_both$uniform_high,
    band = "Uniform 95% band"
  ),
  data.frame(
    estimator = output_inference_both$estimator,
    log_output = output_inference_both$log_output,
    low = output_inference_both$point_low,
    high = output_inference_both$point_high,
    band = "Pointwise 95% intervals"
  )
)
band_data$band <- factor(
  band_data$band,
  levels = c("Uniform 95% band", "Pointwise 95% intervals")
)

output_plot <- ggplot() +
  geom_ribbon(
    data = band_data,
    aes(x = log_output, ymin = low, ymax = high, fill = band),
    alpha = 0.30
  ) +
  geom_line(
    data = output_inference_both,
    aes(x = log_output, y = estimate),
    color = "#163A5F",
    linewidth = 1
  ) +
  geom_hline(yintercept = 0, color = "grey65", linewidth = 0.45) +
  geom_rug(
    data = rbind(
      data.frame(
        estimator = factor(
          "Factor fixed effects",
          levels = levels(output_inference_both$estimator)
        ),
        log_output = df$log_output[
          df$log_output >= min(output_grid) &
            df$log_output <= max(output_grid)
        ]
      ),
      data.frame(
        estimator = factor(
          "First differences",
          levels = levels(output_inference_both$estimator)
        ),
        log_output = df_fd$log_output[
          df_fd$log_output >= min(output_grid) &
            df_fd$log_output <= max(output_grid)
        ]
      )
    ),
    aes(x = log_output),
    alpha = 0.10,
    sides = "b"
  ) +
  scale_fill_manual(values = c(
    "Uniform 95% band" = "#9ECAE1",
    "Pointwise 95% intervals" = "#3182BD"
  )) +
  facet_wrap(~estimator, nrow = 1) +
  labs(
    title = "Partially linear model: estimated output effect",
    subtitle = "Penalty-adjusted firm-clustered inference; each smooth is centered over the displayed grid.",
    x = "Log output",
    y = "Partial effect on log employment",
    fill = NULL
  ) +
  paper_theme

save_plot <- function(plot, stem, width, height) {
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(
    pdf_path,
    plot,
    width = width,
    height = height,
    device = grDevices::pdf,
    bg = "white"
  )
  file.copy(
    pdf_path,
    file.path(paper_figure_dir, paste0(stem, ".pdf")),
    overwrite = TRUE
  )
}

save_plot(comparison_plot, "empluk_effect_comparison", 9, 6.6)
save_plot(output_plot, "empluk_output_inference", 9, 4.8)

format_num <- function(x, digits = 3L) {
  ifelse(is.na(x), "---", formatC(x, format = "f", digits = digits))
}

format_se <- function(x) {
  ifelse(is.na(x), "", paste0("(", format_num(x), ")"))
}

coefficient_lines <- unlist(lapply(seq_len(nrow(coefficient_results)), function(i) {
  d <- coefficient_results[i, ]
  c(
    sprintf(
      "%s & %s & %s & %s \\\\",
      d$label,
      format_num(d$linear_estimate),
      format_num(d$pl_estimate),
      format_num(d$fd_estimate)
    ),
    sprintf(
      " & %s & %s & %s \\\\",
      format_se(d$linear_cluster_se),
      format_se(d$pl_cluster_se),
      format_se(d$fd_cluster_se)
    )
  )
}), use.names = FALSE)

comparison_lines <- vapply(seq_len(nrow(model_comparison)), function(i) {
  d <- model_comparison[i, ]
  sprintf(
    "%s & %.2f & %.2f & %.2f & %.2f & %.4f & %.1f & %.1f \\\\",
    d$model,
    d$wage_edf,
    d$capital_edf,
    d$output_edf,
    d$model_df,
    d$RMSE,
    d$AIC,
    d$BIC
  )
}, character(1))

table_tex <- c(
  "% Generated by empirical_application_empluk/run_manuscript_application.R.",
  "% Do not edit numerical entries by hand.",
  "",
  "\\begin{table}[H]",
  "\\centering",
  "\\caption{Linear coefficients in the EmplUK application}",
  "\\label{tab:empluk_coefficients}",
  "\\begin{tabular}{lccc}",
  "\\toprule",
  "Covariate & Linear FE & Partially linear FE & Partially linear FD \\\\",
  "\\midrule",
  coefficient_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{0.92\\textwidth}",
  sprintf(
    paste0(
      "\\tablenote\\ The dependent variable is log employment in the levels specifications and its first difference in the FD specification. ",
      "All specifications include year effects; the levels models also include firm effects. ",
      "Standard errors clustered by firm are reported in parentheses. In the partially linear models, log output enters through a smooth function and therefore has no scalar coefficient; ",
      "its estimated EDF is %.2f in the factor fixed-effects specification and %.2f in the first-difference specification."
    ),
    smooth_edf(pl_fe, "log_output"), fd_output_edf
  ),
  "\\end{minipage}",
  "\\end{table}",
  "",
  "\\begin{table}[H]",
  "\\centering",
  "\\footnotesize",
  "\\caption{Comparison of specifications in the EmplUK application}",
  "\\label{tab:empluk_model_comparison}",
  "\\begin{tabular}{lrrrrrrr}",
  "\\toprule",
  "Model & Wage EDF & Capital EDF & Output EDF & Model df & RMSE & AIC & BIC \\\\",
  "\\midrule",
  "\\multicolumn{8}{l}{\\textit{Panel A: levels}} \\\\",
  "\\midrule",
  comparison_lines[1:3],
  "\\midrule",
  "\\multicolumn{8}{l}{\\textit{Panel B: first differences}} \\\\",
  "\\midrule",
  comparison_lines[4],
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{0.94\\textwidth}",
  "\\tablenote\\ Linear components are assigned one degree of freedom. Model df is the effective parameter count returned by \\code{logLik(fit)} and used in AIC and BIC; when available, it uses \\code{mgcv}'s corrected EDF and therefore need not equal a simple sum of the displayed component EDFs. It includes year effects, the estimated scale parameter, and firm effects in the levels models. RMSE is calculated for log employment in Panel A and its first difference in Panel B. AIC and BIC are reported descriptively; fREML selects the smoothing parameters. RMSE, AIC, and BIC should not be compared across panels because the dependent variables, estimation samples, and objectives differ.",
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(table_tex, file.path(paper_dir, "empluk_application_results.tex"))

summary_lines <- c(
  paste("Original observations:", nrow(df)),
  paste("Levels observations:", nrow(df)),
  paste("Firms:", nlevels(df$firm_fe)),
  paste("Levels years:", paste(range(df$year), collapse = "--")),
  paste("Partially linear output EDF:", format_num(smooth_edf(pl_fe, "log_output"))),
  paste("First-difference observations:", nrow(df_fd)),
  paste("First-difference output EDF:", format_num(fd_output_edf)),
  paste("Uniform critical value:", format_num(output_uniform_critical)),
  paste("First-difference uniform critical value:", format_num(output_uniform_critical_fd)),
  "",
  "Coefficient results:",
  capture.output(print(coefficient_results, row.names = FALSE)),
  "",
  "Model comparison:",
  capture.output(print(model_comparison, row.names = FALSE)),
  "",
  "Basis sensitivity:",
  capture.output(print(basis_sensitivity, row.names = FALSE)),
  "",
  "Partially linear model summary:",
  capture.output(print(summary(pl_fe))),
  "",
  "First-difference partially linear model summary:",
  capture.output(print(summary(pl_fd))),
  "",
  "Additive GAM summary:",
  capture.output(print(summary(gam_fe))),
  "",
  "Session information:",
  capture.output(print(sessionInfo()))
)
summary_lines <- sub("[[:space:]]+$", "", summary_lines)
writeLines(summary_lines, file.path(output_dir, "manuscript_summary.txt"))

message("Completed manuscript EmplUK application.")
message("Outputs: ", output_dir)
message("Figures: ", figure_dir)
message("LaTeX tables: ", file.path(paper_dir, "empluk_application_results.tex"))
