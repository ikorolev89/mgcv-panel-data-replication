# Exploratory linear and semiparametric fixed-effects models for plm::Snmesp.
#
# Run from the repository root with:
#   Rscript empirical_application_snmesp/run_snmesp_application.R

suppressPackageStartupMessages({
  library(plm)
  library(mgcv)
  library(ggplot2)
  library(sandwich)
  library(lmtest)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
  app_dir <- dirname(script_path)
} else {
  app_dir <- normalizePath(getwd(), mustWork = TRUE)
}

data_dir <- file.path(app_dir, "data")
output_dir <- file.path(app_dir, "outputs")
figure_dir <- file.path(app_dir, "figures")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(182017)

data("Snmesp", package = "plm", envir = environment())
df <- as.data.frame(Snmesp)

required_columns <- c("firm", "year", "n", "w", "y", "i", "k", "f")
if (!identical(names(df), required_columns)) {
  stop("Unexpected Snmesp columns: ", paste(names(df), collapse = ", "))
}
if (anyNA(df)) stop("Snmesp unexpectedly contains missing values.")

df <- df[order(df$firm, df$year), ]
df$firm_fe <- factor(df$firm)
df$year_fe <- factor(df$year)

panel_counts <- table(df$firm_fe)
if (length(unique(panel_counts)) != 1L) stop("The panel is not balanced.")

write.csv(
  df[required_columns],
  file.path(data_dir, "Snmesp.csv"),
  row.names = FALSE
)

covariates <- c("w", "y", "i", "k", "f")
variable_labels <- c(
  w = "Log wages",
  y = "Log output",
  i = "Log intermediate inputs",
  k = "Log capital",
  f = "Real cash flow"
)

# Linear two-way FE benchmark. The reported covariance matrix clusters by firm.
linear_formula <- n ~ w + y + i + k + f + firm_fe + year_fe
linear_fe <- lm(linear_formula, data = df)
linear_vcov_cluster <- sandwich::vcovCL(
  linear_fe,
  cluster = df$firm_fe,
  type = "HC1"
)
linear_test <- lmtest::coeftest(linear_fe, vcov. = linear_vcov_cluster)

linear_coefficients <- data.frame(
  variable = covariates,
  label = unname(variable_labels[covariates]),
  estimate = linear_test[covariates, "Estimate"],
  cluster_se = linear_test[covariates, "Std. Error"],
  statistic = linear_test[covariates, "t value"],
  p_value = linear_test[covariates, "Pr(>|t|)"],
  row.names = NULL
)
linear_coefficients$conf_low <-
  linear_coefficients$estimate - 1.96 * linear_coefficients$cluster_se
linear_coefficients$conf_high <-
  linear_coefficients$estimate + 1.96 * linear_coefficients$cluster_se

write.csv(
  linear_coefficients,
  file.path(output_dir, "linear_fe_coefficients.csv"),
  row.names = FALSE
)

# The modest basis size is deliberate: T is only eight, and this is a first-pass
# shape diagnostic rather than a search for fine local features.
smooth_formula <- n ~
  s(w, bs = "cr", k = 8) +
  s(y, bs = "cr", k = 8) +
  s(i, bs = "cr", k = 8) +
  s(k, bs = "cr", k = 8) +
  s(f, bs = "cr", k = 8) +
  firm_fe + year_fe

semi_fe <- mgcv::bam(
  smooth_formula,
  data = df,
  method = "fREML",
  discrete = TRUE,
  select = TRUE
)

semi_summary <- summary(semi_fe)

rmse <- function(actual, fitted) sqrt(mean((actual - fitted)^2))

linear_loglik <- logLik(linear_fe)
semi_loglik <- logLik(semi_fe)
model_comparison <- data.frame(
  model = c("Linear two-way FE", "Additive semiparametric two-way FE"),
  observations = nrow(df),
  effective_parameters = c(
    attr(linear_loglik, "df"),
    attr(semi_loglik, "df")
  ),
  AIC = c(AIC(linear_fe), AIC(semi_fe)),
  BIC = c(BIC(linear_fe), BIC(semi_fe)),
  in_sample_RMSE = c(
    rmse(df$n, fitted(linear_fe)),
    rmse(df$n, fitted(semi_fe))
  ),
  adjusted_R_squared = c(
    summary(linear_fe)$adj.r.squared,
    semi_summary$r.sq
  ),
  deviance_explained = c(
    1 - deviance(linear_fe) / sum((df$n - mean(df$n))^2),
    semi_summary$dev.expl
  )
)
write.csv(
  model_comparison,
  file.path(output_dir, "model_comparison.csv"),
  row.names = FALSE
)

smooth_table <- as.data.frame(semi_summary$s.table)
smooth_table$term <- rownames(smooth_table)
rownames(smooth_table) <- NULL

term_predictions_observed <- predict(semi_fe, type = "terms")
linear_slopes <- setNames(linear_coefficients$estimate, linear_coefficients$variable)

make_effect_data <- function(variable) {
  term_name <- paste0("s(", variable, ")")
  x <- df[[variable]]
  limits <- as.numeric(quantile(x, probs = c(0.01, 0.99), names = FALSE))
  x_grid <- seq(limits[1], limits[2], length.out = 240L)

  newdata <- data.frame(
    w = median(df$w),
    y = median(df$y),
    i = median(df$i),
    k = median(df$k),
    f = median(df$f),
    firm_fe = factor(levels(df$firm_fe)[1], levels = levels(df$firm_fe)),
    year_fe = factor(levels(df$year_fe)[1], levels = levels(df$year_fe))
  )
  newdata <- newdata[rep(1L, length(x_grid)), ]
  newdata[[variable]] <- x_grid

  term_fit <- predict(
    semi_fe,
    newdata = newdata,
    type = "terms",
    se.fit = TRUE
  )
  smooth_effect <- as.numeric(term_fit$fit[, term_name])
  smooth_se <- as.numeric(term_fit$se.fit[, term_name])

  # This projection isolates curvature in the estimated smooth. It is fitted
  # over the empirical distribution of x, not uniformly over the plotting grid.
  observed_smooth <- as.numeric(term_predictions_observed[, term_name])
  best_linear_projection <- lm(observed_smooth ~ x)
  projected_on_grid <- as.numeric(
    predict(best_linear_projection, newdata = data.frame(x = x_grid))
  )

  linear_fe_effect <- linear_slopes[[variable]] * (x_grid - mean(x))
  nonlinear_deviation <- smooth_effect - projected_on_grid
  central_limits <- as.numeric(
    quantile(x, probs = c(0.05, 0.95), names = FALSE)
  )
  central_index <- x_grid >= central_limits[1] & x_grid <= central_limits[2]

  effect_data <- data.frame(
    variable = variable,
    label = unname(variable_labels[[variable]]),
    x = x_grid,
    smooth = smooth_effect,
    conf_low = smooth_effect - 1.96 * smooth_se,
    conf_high = smooth_effect + 1.96 * smooth_se,
    linear_fe = linear_fe_effect,
    smooth_linear_projection = projected_on_grid,
    nonlinear_deviation = nonlinear_deviation
  )

  residualized_observed <- residuals(best_linear_projection)
  diagnostics <- data.frame(
    variable = variable,
    label = unname(variable_labels[[variable]]),
    plot_p01 = limits[1],
    plot_p99 = limits[2],
    central_p05 = central_limits[1],
    central_p95 = central_limits[2],
    full_minimum = min(x),
    full_maximum = max(x),
    nonlinear_sd_empirical = sd(residualized_observed),
    nonlinear_range_p01_p99 = diff(range(nonlinear_deviation)),
    maximum_absolute_deviation_p01_p99 = max(abs(nonlinear_deviation)),
    nonlinear_range_p05_p95 = diff(range(nonlinear_deviation[central_index])),
    maximum_absolute_deviation_p05_p95 =
      max(abs(nonlinear_deviation[central_index]))
  )

  list(effect_data = effect_data, diagnostics = diagnostics)
}

effect_results <- lapply(covariates, make_effect_data)
effect_data <- do.call(rbind, lapply(effect_results, `[[`, "effect_data"))
curvature_diagnostics <- do.call(
  rbind,
  lapply(effect_results, `[[`, "diagnostics")
)

write.csv(
  effect_data,
  file.path(output_dir, "effect_curves.csv"),
  row.names = FALSE
)

smooth_diagnostics <- merge(
  data.frame(
    variable = sub("^s\\((.*)\\)$", "\\1", smooth_table$term),
    edf = smooth_table$edf,
    reference_df = smooth_table$Ref.df,
    F_statistic = smooth_table$F,
    smooth_zero_test_p_value = smooth_table$`p-value`
  ),
  curvature_diagnostics,
  by = "variable",
  sort = FALSE
)
smooth_diagnostics <- smooth_diagnostics[
  match(covariates, smooth_diagnostics$variable),
]
write.csv(
  smooth_diagnostics,
  file.path(output_dir, "smooth_diagnostics.csv"),
  row.names = FALSE
)

concurvity_result <- mgcv::concurvity(semi_fe, full = FALSE)
concurvity_estimate <- as.data.frame(concurvity_result$estimate)
concurvity_estimate$term_explained <- rownames(concurvity_estimate)
rownames(concurvity_estimate) <- NULL
concurvity_estimate <- concurvity_estimate[
  , c("term_explained", setdiff(names(concurvity_estimate), "term_explained"))
]
write.csv(
  concurvity_estimate,
  file.path(output_dir, "concurvity.csv"),
  row.names = FALSE
)

sample_rug <- function(variable, n = 500L) {
  values <- df[[variable]]
  limits <- as.numeric(quantile(values, probs = c(0.01, 0.99), names = FALSE))
  values <- values[values >= limits[1] & values <= limits[2]]
  if (length(values) > n) values <- sample(values, n)
  data.frame(
    variable = variable,
    label = unname(variable_labels[[variable]]),
    x = values
  )
}
rug_data <- do.call(rbind, lapply(covariates, sample_rug))
label_levels <- unname(variable_labels[covariates])
effect_data$label <- factor(effect_data$label, levels = label_levels)
rug_data$label <- factor(rug_data$label, levels = label_levels)

base_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title.position = "plot",
    legend.position = "bottom"
  )

smooth_plot <- ggplot(effect_data, aes(x = x)) +
  geom_ribbon(
    aes(ymin = conf_low, ymax = conf_high),
    fill = "#77AADD",
    alpha = 0.28
  ) +
  geom_line(aes(y = smooth, color = "Additive FE smooth"), linewidth = 0.9) +
  geom_line(
    aes(y = linear_fe, color = "Linear FE fit"),
    linewidth = 0.8,
    linetype = "22"
  ) +
  geom_rug(
    data = rug_data,
    aes(x = x),
    inherit.aes = FALSE,
    alpha = 0.10,
    sides = "b"
  ) +
  facet_wrap(~label, scales = "free", ncol = 2) +
  scale_color_manual(values = c(
    "Additive FE smooth" = "#1F5A89",
    "Linear FE fit" = "#B24C3D"
  )) +
  labs(
    title = "Snmesp: additive covariate functions versus linear FE fits",
    subtitle = paste(
      "Firm and year effects in both models; shaded intervals are pointwise and model-based.",
      "Axes show the 1st--99th percentiles."
    ),
    x = NULL,
    y = "Partial effect on log employment",
    color = NULL
  ) +
  base_theme

deviation_plot <- ggplot(effect_data, aes(x = x, y = nonlinear_deviation)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.5) +
  geom_line(color = "#6A3D9A", linewidth = 0.9) +
  geom_rug(
    data = rug_data,
    aes(x = x),
    inherit.aes = FALSE,
    alpha = 0.10,
    sides = "b"
  ) +
  facet_wrap(~label, scales = "free", ncol = 2) +
  labs(
    title = "Estimated nonlinear deviations",
    subtitle = paste(
      "Each additive function minus its best linear projection over the empirical covariate distribution.",
      "Axes show the 1st--99th percentiles."
    ),
    x = NULL,
    y = "Deviation from linearity (log employment)"
  ) +
  base_theme

diagnostic_data <- rbind(
  data.frame(
    model = "Linear two-way FE",
    fitted = fitted(linear_fe),
    residual = residuals(linear_fe)
  ),
  data.frame(
    model = "Additive semiparametric two-way FE",
    fitted = fitted(semi_fe),
    residual = residuals(semi_fe)
  )
)

fit_diagnostic_plot <- ggplot(
  diagnostic_data,
  aes(x = fitted, y = residual)
) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.5) +
  geom_point(alpha = 0.16, size = 0.65, color = "#2B6F6D") +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
              color = "#B24C3D", linewidth = 0.8) +
  facet_wrap(~model, ncol = 1) +
  labs(
    title = "Residual patterns under the two specifications",
    x = "Fitted log employment",
    y = "Residual"
  ) +
  base_theme

save_plot_pair <- function(plot, stem, width, height) {
  ggsave(
    file.path(figure_dir, paste0(stem, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = 220,
    bg = "white"
  )
  ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    device = grDevices::pdf,
    bg = "white"
  )
}

save_plot_pair(smooth_plot, "smooth_effects_vs_linear", 10, 10)
save_plot_pair(deviation_plot, "nonlinear_deviations", 10, 10)
save_plot_pair(fit_diagnostic_plot, "model_fit_diagnostics", 9, 8)

summary_path <- file.path(output_dir, "model_summaries.txt")
summary_connection <- file(summary_path, open = "wt")
sink(summary_connection)
cat("Exploratory Snmesp application\n")
cat("==============================\n\n")
cat("Run date:", format(Sys.time()), "\n")
cat("Observations:", nrow(df), "\n")
cat("Firms:", nlevels(df$firm_fe), "\n")
cat("Years:", paste(range(df$year), collapse = "--"), "\n")
cat("Observations per firm:", unique(panel_counts), "\n\n")

cat("Variable summaries\n")
cat("------------------\n")
print(summary(df[required_columns]))

cat("\n\nLinear two-way FE coefficients (firm-clustered covariance)\n")
cat("---------------------------------------------------------\n")
print(linear_coefficients, row.names = FALSE)

cat("\n\nAdditive semiparametric two-way FE summary\n")
cat("-------------------------------------------\n")
print(semi_summary)

cat("\n\nSmooth diagnostics\n")
cat("------------------\n")
print(smooth_diagnostics, row.names = FALSE)

cat("\n\nModel comparison\n")
cat("----------------\n")
print(model_comparison, row.names = FALSE)

cat("\n\nBasis-dimension check\n")
cat("---------------------\n")
print(mgcv::k.check(semi_fe))

cat("\n\nSession information\n")
cat("-------------------\n")
print(sessionInfo())
sink()
close(summary_connection)

message("Completed exploratory Snmesp application.")
message("Outputs: ", output_dir)
message("Figures: ", figure_dir)
