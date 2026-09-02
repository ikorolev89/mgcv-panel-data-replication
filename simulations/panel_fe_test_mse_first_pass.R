# First-pass out-of-sample MSE simulation for fixed-effects panel estimators.
#
# Design:
#   n = 200 units, T periods, p = 1 covariate
#   DGPs: configurable linear, quadratic, and sine regression functions
#   Accuracy target: out-of-sample MSE for g(x) on 1000 fresh test points
#
# Methods:
#   1. Linear FE
#   2. Linear FD
#   3. Fixed-K spline FE in levels
#   4. Fixed-K first-differenced spline
#   5. mgcv::bam with unit fixed effects
#   6. mgcv::bam with unit random effects
#   7. mgcv::bam on first differences with one shared smooth
#   8. optional grf forests: levels, levels + unit, raw differenced, antisym differenced
#
# Usage:
#   Rscript simulations/panel_fe_test_mse_first_pass.R 10 7 500 4 baseline 200 FALSE
#   Rscript simulations/panel_fe_test_mse_first_pass.R 10 7 500 4 sine_complexity 200 FALSE
#   Rscript simulations/panel_fe_test_mse_first_pass.R 10 7 500 4 paper_main 200 FALSE hetero
#   Rscript simulations/panel_fe_test_mse_first_pass.R 10 7 500 4 paper_main 200 FALSE hetero_ar1 0.5
#
# Arguments:
#   1. number of simulation draws
#   2. spline basis dimension K
#   3. number of trees per forest
#   4. number of time periods T
#   5. DGP set: baseline, paper_main, sine_complexity, all, or one DGP name
#   6. number of units n
#   7. include forests: TRUE or FALSE
#   8. error design: iid, hetero, or hetero_ar1; default iid
#   9. AR(1) rho for hetero_ar1; default 0.5

suppressPackageStartupMessages({
  library(mgcv)
  library(ggplot2)
  library(splines)
})

args <- commandArgs(trailingOnly = TRUE)
n_sims <- if (length(args) >= 1) as.integer(args[[1]]) else 10L
spline_df <- if (length(args) >= 2) as.integer(args[[2]]) else 7L
num_trees <- if (length(args) >= 3) as.integer(args[[3]]) else 500L
t_periods <- if (length(args) >= 4) as.integer(args[[4]]) else 4L
dgp_arg <- if (length(args) >= 5) args[[5]] else "baseline"
n_units <- if (length(args) >= 6) as.integer(args[[6]]) else 200L
include_forests <- if (length(args) >= 7) {
  tolower(args[[7]]) %in% c("true", "t", "1", "yes", "y")
} else {
  FALSE
}
error_design <- if (length(args) >= 8) args[[8]] else "iid"
rho <- if (length(args) >= 9) as.numeric(args[[9]]) else 0.5

if (is.na(n_sims) || n_sims < 1L) stop("n_sims must be a positive integer.")
if (is.na(spline_df) || spline_df < 3L) stop("spline_df must be at least 3.")
if (is.na(num_trees) || num_trees < 1L) stop("num_trees must be a positive integer.")
if (is.na(t_periods) || t_periods < 2L) stop("T must be at least 2.")
if (is.na(n_units) || n_units < 2L) stop("n_units must be at least 2.")
if (!error_design %in% c("iid", "hetero", "hetero_ar1")) {
  stop("error design must be iid, hetero, or hetero_ar1.")
}
if (is.na(rho) || abs(rho) >= 1) {
  stop("rho must be numeric with absolute value less than 1.")
}

n_test <- 1000L
sigma_mu <- 1.25
sigma_e <- 0.25

available_dgps <- c("linear", "quadratic", "sin_half_pi", "sin_2pi", "sin_8pi")
dgp_sets <- list(
  baseline = c("linear", "sin_2pi"),
  paper_main = c("linear", "quadratic", "sin_2pi"),
  sine_complexity = c("sin_half_pi", "sin_8pi"),
  all = available_dgps
)
if (dgp_arg %in% names(dgp_sets)) {
  dgp_list <- dgp_sets[[dgp_arg]]
} else if (dgp_arg %in% available_dgps) {
  dgp_list <- dgp_arg
} else if (dgp_arg == "nonlinear") {
  dgp_list <- "sin_2pi"
} else {
  stop(
    "DGP set must be one of ",
    paste(c(names(dgp_sets), available_dgps, "nonlinear"), collapse = ", "),
    "."
  )
}

out_dir <- file.path(getwd(), "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Simulation draws: ", n_sims)
message("n units: ", n_units)
message("T periods: ", t_periods)
message("test points per draw: ", n_test)
message("Spline basis dimension K: ", spline_df)
message("Trees per forest: ", num_trees)
message("DGPs: ", paste(dgp_list, collapse = ", "))
message("Include forests: ", include_forests)
message("Error design: ", error_design)
if (error_design == "hetero_ar1") message("AR(1) rho: ", rho)

if (include_forests && !requireNamespace("grf", quietly = TRUE)) {
  stop("Forest runs use grf::regression_forest. Install it with install.packages('grf').")
}

g_fun <- function(x, dgp) {
  switch(
    dgp,
    linear = 0.5 + 1.4 * x,
    quadratic = 0.5 + 1.4 * x + 0.75 * (x - 0.5)^2,
    sin_half_pi = sin(0.5 * pi * x),
    sin_2pi = sin(2 * pi * x),
    sin_8pi = sin(8 * pi * x)
  )
}

draw_x <- function(n) {
  alpha <- rbeta(n, shape1 = 2.2, shape2 = 2.2)
  pmin(pmax(alpha + rnorm(n, sd = 0.18), 0), 1)
}

normalized_error_variance <- function(index, sigma_e) {
  raw_var <- 1 + 2 * exp(0.75 * index)
  sigma_e^2 * raw_var / mean(raw_var)
}

draw_panel_errors <- function(unit, t_periods, variance_index, sigma_e,
                              error_design, rho) {
  n_obs <- length(unit)
  if (error_design == "iid") {
    return(rnorm(n_obs, sd = sigma_e))
  }

  error_var <- normalized_error_variance(variance_index, sigma_e)
  if (error_design == "hetero") {
    return(rnorm(n_obs, sd = sqrt(error_var)))
  }

  epsilon <- numeric(n_obs)
  for (i in unique(unit)) {
    idx <- which(unit == i)
    ti <- length(idx)
    if (ti != t_periods) stop("Each unit must have T observations.")
    ar_corr <- outer(seq_len(ti), seq_len(ti), function(s, t) rho^abs(s - t))
    sd_i <- sqrt(error_var[idx])
    omega <- ar_corr * outer(sd_i, sd_i)
    chol_omega <- try(chol(omega), silent = TRUE)
    if (inherits(chol_omega, "try-error")) {
      omega <- omega + diag(1e-10 * mean(diag(omega)), ti)
      chol_omega <- chol(omega)
    }
    epsilon[idx] <- as.numeric(t(chol_omega) %*% rnorm(ti))
  }

  epsilon
}

simulate_panel <- function(dgp) {
  unit <- rep(seq_len(n_units), each = t_periods)
  time <- rep(seq_len(t_periods), times = n_units)

  alpha_i <- rbeta(n_units, shape1 = 2.2, shape2 = 2.2)
  mu_i <- 1.1 * (alpha_i - mean(alpha_i)) + rnorm(n_units, sd = sigma_mu)

  x <- numeric(n_units * t_periods)
  for (i in seq_len(n_units)) {
    idx <- unit == i
    x[idx] <- pmin(pmax(alpha_i[i] + rnorm(t_periods, sd = 0.18), 0), 1)
  }

  e <- draw_panel_errors(
    unit = unit,
    t_periods = t_periods,
    variance_index = x,
    sigma_e = sigma_e,
    error_design = error_design,
    rho = rho
  )

  data.frame(
    unit = factor(unit),
    time = time,
    x = x,
    y = g_fun(x, dgp) + mu_i[unit] + e
  )
}

make_first_differences <- function(panel_df) {
  panel_df <- panel_df[order(panel_df$unit, panel_df$time), ]
  split_df <- split(panel_df, panel_df$unit)

  diffs <- lapply(split_df, function(d) {
    data.frame(
      unit = d$unit[-1],
      time = d$time[-1],
      dy = diff(d$y),
      x_now = d$x[-1],
      x_lag = d$x[-nrow(d)]
    )
  })

  do.call(rbind, diffs)
}

make_bs <- function(x) {
  B <- bs(
    x,
    df = spline_df,
    degree = 3,
    intercept = FALSE,
    Boundary.knots = c(0, 1)
  )
  colnames(B) <- paste0("b", seq_len(ncol(B)))
  B
}

predict_bs <- function(B, x_new) {
  B_new <- predict(B, newx = x_new)
  colnames(B_new) <- colnames(B)
  B_new
}

make_unit_dummies <- function(unit, unit_levels) {
  mm <- model.matrix(
    ~ unit - 1,
    data = data.frame(unit = factor(unit, levels = unit_levels))
  )
  colnames(mm) <- paste0("unit_", unit_levels)
  mm
}

fit_grf <- function(X, y, seed) {
  if (!include_forests) stop("fit_grf called when include_forests is FALSE.")
  grf::regression_forest(
    X = as.matrix(X),
    Y = as.numeric(y),
    num.trees = num_trees,
    seed = seed
  )
}

predict_grf <- function(fit, X_new) {
  if (!include_forests) stop("predict_grf called when include_forests is FALSE.")
  as.numeric(predict(fit, newdata = as.matrix(X_new))$predictions)
}

center_to_truth <- function(estimate, truth) {
  estimate - mean(estimate) + mean(truth)
}

mse <- function(estimate, truth) {
  mean((estimate - truth)^2)
}

predict_h_raw <- function(fit, a, b) {
  predict_grf(fit, data.frame(x_a = a, x_b = b))
}

predict_h_sym <- function(fit, a, b) {
  0.5 * (
    predict_h_raw(fit, a, b) -
      predict_h_raw(fit, b, a)
  )
}

symmetrize_differences <- function(diff_df) {
  forward <- data.frame(
    dy = diff_df$dy,
    x_a = diff_df$x_now,
    x_b = diff_df$x_lag
  )

  reverse <- data.frame(
    dy = -diff_df$dy,
    x_a = diff_df$x_lag,
    x_b = diff_df$x_now
  )

  rbind(forward, reverse)
}

fit_one_draw <- function(sim_id, dgp) {
  df <- simulate_panel(dgp)
  df_diff <- make_first_differences(df)

  x_test <- draw_x(n_test)
  truth <- g_fun(x_test, dgp)
  reference_unit <- levels(df$unit)[1]

  B <- make_bs(df$x)
  B_df <- as.data.frame(B)
  spline_cols <- colnames(B_df)
  df_spline <- cbind(df, B_df)
  B_test <- predict_bs(B, x_test)
  test_spline <- cbind(
    data.frame(unit = factor(reference_unit, levels = levels(df$unit))),
    as.data.frame(B_test)
  )

  linear_fe <- lm(y ~ x + unit, data = df)
  linear_hat <- predict(
    linear_fe,
    newdata = data.frame(x = x_test, unit = factor(reference_unit, levels = levels(df$unit)))
  )
  linear_hat <- center_to_truth(linear_hat, truth)

  linear_fd_df <- data.frame(dy = df_diff$dy, dx = df_diff$x_now - df_diff$x_lag)
  linear_fd <- lm(dy ~ dx - 1, data = linear_fd_df)
  linear_fd_hat <- as.numeric(coef(linear_fd)[["dx"]] * x_test)
  linear_fd_hat <- center_to_truth(linear_fd_hat, truth)

  spline_formula <- as.formula(paste("y ~", paste(spline_cols, collapse = " + "), "+ unit"))
  spline_fe <- lm(spline_formula, data = df_spline)
  spline_hat <- predict(spline_fe, newdata = test_spline)
  spline_hat <- center_to_truth(spline_hat, truth)

  B_now <- predict_bs(B, df_diff$x_now)
  B_lag <- predict_bs(B, df_diff$x_lag)
  DB <- B_now - B_lag
  colnames(DB) <- spline_cols
  fd_fit <- lm(dy ~ . - 1, data = cbind(data.frame(dy = df_diff$dy), as.data.frame(DB)))
  fd_coef <- coef(fd_fit)[spline_cols]
  fd_hat <- as.numeric(B_test %*% fd_coef)
  fd_hat <- center_to_truth(fd_hat, truth)

  bam_fe_fit <- bam(
    y ~ s(x, k = 20) + unit,
    data = df,
    method = "fREML",
    select = TRUE,
    discrete = TRUE
  )
  bam_fe_hat <- predict(
    bam_fe_fit,
    newdata = data.frame(x = x_test, unit = factor(reference_unit, levels = levels(df$unit)))
  )
  bam_fe_hat <- center_to_truth(bam_fe_hat, truth)

  bam_re_fit <- bam(
    y ~ s(x, k = 20) + s(unit, bs = "re"),
    data = df,
    method = "fREML",
    select = TRUE,
    discrete = TRUE
  )
  bam_re_hat <- predict(
    bam_re_fit,
    newdata = data.frame(x = x_test, unit = factor(reference_unit, levels = levels(df$unit)))
  )
  bam_re_hat <- center_to_truth(bam_re_hat, truth)

  X_pair <- cbind(df_diff$x_now, df_diff$x_lag)
  L_pair <- cbind(rep(1, nrow(df_diff)), rep(-1, nrow(df_diff)))
  bam_fd_fit <- bam(
    dy ~ s(X_pair, by = L_pair, k = 20),
    data = list(dy = df_diff$dy, X_pair = X_pair, L_pair = L_pair),
    method = "fREML",
    select = TRUE,
    discrete = TRUE
  )
  bam_fd_hat <- predict(
    bam_fd_fit,
    newdata = data.frame(X_pair = x_test, L_pair = 1)
  )
  bam_fd_hat <- center_to_truth(bam_fd_hat, truth)

  base_result <- data.frame(
    sim = sim_id,
    dgp = dgp,
    method = c(
      "linear FE",
      "linear FD",
      "spline FE, fixed K",
      "FD spline, fixed K",
      "mgcv::bam FE",
      "mgcv::bam RE",
      "mgcv::bam FD"
    ),
    mse = c(
      mse(linear_hat, truth),
      mse(linear_fd_hat, truth),
      mse(spline_hat, truth),
      mse(fd_hat, truth),
      mse(bam_fe_hat, truth),
      mse(bam_re_hat, truth),
      mse(bam_fd_hat, truth)
    )
  )

  if (!include_forests) return(base_result)

  forest_levels <- fit_grf(
    X = data.frame(x = df$x),
    y = df$y,
    seed = 100000 + sim_id
  )
  forest_levels_hat <- predict_grf(forest_levels, data.frame(x = x_test))
  forest_levels_hat <- center_to_truth(forest_levels_hat, truth)

  unit_train <- make_unit_dummies(df$unit, levels(df$unit))
  unit_test <- make_unit_dummies(rep(reference_unit, n_test), levels(df$unit))
  forest_unit <- fit_grf(
    X = cbind(x = df$x, unit_train),
    y = df$y,
    seed = 200000 + sim_id
  )
  forest_unit_hat <- predict_grf(forest_unit, cbind(x = x_test, unit_test))
  forest_unit_hat <- center_to_truth(forest_unit_hat, truth)

  diff_raw_df <- data.frame(
    dy = df_diff$dy,
    x_a = df_diff$x_now,
    x_b = df_diff$x_lag
  )
  diff_sym_df <- symmetrize_differences(df_diff)

  forest_diff_raw <- fit_grf(
    X = diff_raw_df[, c("x_a", "x_b")],
    y = diff_raw_df$dy,
    seed = 300000 + sim_id
  )
  forest_diff_sym <- fit_grf(
    X = diff_sym_df[, c("x_a", "x_b")],
    y = diff_sym_df$dy,
    seed = 400000 + sim_id
  )

  x_anchor <- 0.5
  diff_raw_hat <- predict_h_raw(forest_diff_raw, x_test, rep(x_anchor, n_test))
  diff_raw_hat <- center_to_truth(diff_raw_hat, truth)

  diff_sym_hat <- predict_h_sym(forest_diff_sym, x_test, rep(x_anchor, n_test))
  diff_sym_hat <- center_to_truth(diff_sym_hat, truth)

  rbind(
    base_result,
    data.frame(
      sim = sim_id,
      dgp = dgp,
      method = c(
        "grf: levels, no FE",
        "grf: levels + unit",
        "grf: raw differenced",
        "grf: antisym differenced"
      ),
      mse = c(
        mse(forest_levels_hat, truth),
        mse(forest_unit_hat, truth),
        mse(diff_raw_hat, truth),
        mse(diff_sym_hat, truth)
      )
    )
  )
}

set.seed(20260503)
losses_list <- vector("list", n_sims * length(dgp_list))
idx <- 0L

for (dgp in dgp_list) {
  message("\nDGP: ", dgp)

  for (s in seq_len(n_sims)) {
    if (s == 1L || s %% 10L == 0L || s == n_sims) {
      message("Running draw ", s, " / ", n_sims)
    }

    idx <- idx + 1L
    losses_list[[idx]] <- fit_one_draw(s, dgp)
  }
}

losses <- do.call(rbind, losses_list)
losses$rmse <- sqrt(losses$mse)
losses$error_design <- error_design
losses$rho <- if (error_design == "hetero_ar1") rho else NA_real_

method_order <- c(
  "linear FE",
  "linear FD",
  "spline FE, fixed K",
  "FD spline, fixed K",
  "mgcv::bam FE",
  "mgcv::bam RE",
  "mgcv::bam FD",
  "grf: levels, no FE",
  "grf: levels + unit",
  "grf: raw differenced",
  "grf: antisym differenced"
)
method_order <- method_order[method_order %in% unique(as.character(losses$method))]

losses$dgp <- factor(losses$dgp, levels = dgp_list)
losses$method <- factor(losses$method, levels = method_order)

summary_by_method <- do.call(
  rbind,
  lapply(split(losses, list(losses$dgp, losses$method), drop = TRUE), function(d) {
    data.frame(
      dgp = as.character(d$dgp[1]),
      method = as.character(d$method[1]),
      mean_mse = mean(d$mse),
      sd_mse = sd(d$mse),
      se_mse = sd(d$mse) / sqrt(nrow(d)),
      median_mse = median(d$mse),
      mean_rmse = mean(d$rmse)
    )
  })
)
row.names(summary_by_method) <- NULL
summary_by_method$error_design <- error_design
summary_by_method$rho <- if (error_design == "hetero_ar1") rho else NA_real_
summary_by_method$dgp <- factor(summary_by_method$dgp, levels = dgp_list)
summary_by_method$method <- factor(summary_by_method$method, levels = method_order)
summary_by_method <- summary_by_method[order(summary_by_method$dgp, summary_by_method$method), ]

format_suffix_num <- function(x) {
  out <- gsub("\\.", "p", format(x, scientific = FALSE, trim = TRUE))
  gsub("-", "m", out, fixed = TRUE)
}

error_suffix <- if (error_design == "iid") {
  ""
} else if (error_design == "hetero") {
  "_hetero"
} else {
  paste0("_hetero_ar1_rho", format_suffix_num(rho))
}

dgp_suffix <- paste(dgp_list, collapse = "_")
forest_suffix <- if (include_forests) "_forests" else "_no_forests"
suffix <- paste0(
  "M", n_sims,
  "_",
  "n", n_units,
  "_T", t_periods,
  "_K", spline_df,
  "_test", n_test,
  "_", dgp_suffix,
  forest_suffix,
  error_suffix
)
loss_path <- file.path(out_dir, paste0("panel_fe_test_mse_losses_", suffix, ".csv"))
summary_path <- file.path(out_dir, paste0("panel_fe_test_mse_summary_", suffix, ".csv"))
plot_path <- file.path(out_dir, paste0("panel_fe_test_mse_", suffix, ".png"))

write.csv(losses, loss_path, row.names = FALSE)
write.csv(summary_by_method, summary_path, row.names = FALSE)

mse_plot <- ggplot(summary_by_method, aes(x = method, y = mean_mse, fill = method)) +
  geom_col(width = 0.72, show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = pmax(mean_mse - 1.96 * se_mse, 0), ymax = mean_mse + 1.96 * se_mse),
    width = 0.18
  ) +
  facet_wrap(~ dgp, scales = "free_x") +
  coord_flip() +
  labs(
    title = "Out-of-sample MSE for g(x) in FE panel simulations",
    subtitle = paste0(
      n_sims, " draws; n = ", n_units, ", T = ", t_periods,
      ", test points = ", n_test,
      ", errors = ", error_design,
      if (error_design == "hetero_ar1") paste0(" (rho = ", rho, ")") else "",
      ", spline K = ", spline_df
    ),
    x = NULL,
    y = "mean test MSE"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(plot_path, mse_plot, width = 10, height = 6.5, dpi = 180, bg = "white")

message("\nSummary:")
print(summary_by_method, row.names = FALSE)

message("\nWrote:")
message("  ", loss_path)
message("  ", summary_path)
message("  ", plot_path)
