# Diagnostic coverage simulation for confidence intervals/bands for g in panel GAMs.
#
# Models:
#   np: y_it = g(x_it) + mu_i + e_it
#   pl: y_it = beta * x1_it + g(x2_it) + mu_i + e_it
#
# Methods:
#   1. bam FE: factor unit effects
#   2. bam RE: penalized unit effects
#   3. bam FD: first differences with one constrained smooth
#
# Covariance estimators for centered g on an evaluation grid:
#   1. unit-cluster sandwich ignoring penalties
#   2. penalty-adjusted unit-cluster sandwich
#   3. penalty-adjusted cluster sandwich plus V_B - V_F
#   4. mgcv observation-level sandwich with freq = TRUE
#   5. default mgcv observation-level sandwich (includes V_B - V_F)
#
# Usage:
#   Rscript simulations/inference_diagnostics/panel_fe_g_inference_diagnostics.R \
#     1000 200 4 sin_2pi both 50 499 hetero
#   Rscript simulations/inference_diagnostics/panel_fe_g_inference_diagnostics.R \
#     1000 200 4 sin_2pi both 50 499 hetero_ar1 0.5
#
# Arguments:
#   1. Monte Carlo draws M
#   2. number of units n
#   3. number of periods T
#   4. DGP: linear, sin_half_pi, sin_2pi, sin_8pi
#   5. model: np, pl, or both
#   6. number of grid points J
#   7. Gaussian draws for uniform critical values
#   8. error design: iid, hetero, or hetero_ar1; default iid
#   9. AR(1) rho for hetero_ar1; default 0.5
#   10. small-sample df correction: substantive or full; default substantive
#   11. Monte Carlo seed; default 20260505
#   12. nonempty chunk label

suppressPackageStartupMessages({
  library(mgcv)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
M <- if (length(args) >= 1) as.integer(args[[1]]) else 50L
n_units <- if (length(args) >= 2) as.integer(args[[2]]) else 200L
t_periods <- if (length(args) >= 3) as.integer(args[[3]]) else 4L
dgp <- if (length(args) >= 4) args[[4]] else "sin_2pi"
model_arg <- if (length(args) >= 5) args[[5]] else "both"
n_grid <- if (length(args) >= 6) as.integer(args[[6]]) else 50L
n_sup_draws <- if (length(args) >= 7) as.integer(args[[7]]) else 499L
error_design <- if (length(args) >= 8) args[[8]] else "iid"
rho <- if (length(args) >= 9) as.numeric(args[[9]]) else 0.5
df_correction <- if (length(args) >= 10) args[[10]] else "substantive"
mc_seed <- if (length(args) >= 11) as.integer(args[[11]]) else 20260505L
chunk_label <- if (length(args) >= 12) args[[12]] else ""

available_dgps <- c("linear", "sin_half_pi", "sin_2pi", "sin_8pi")
available_models <- c("np", "pl", "both")
available_error_designs <- c("iid", "hetero", "hetero_ar1")

if (is.na(M) || M < 1L) stop("M must be a positive integer.")
if (is.na(n_units) || n_units < 2L) stop("n_units must be at least 2.")
if (is.na(t_periods) || t_periods < 2L) stop("T must be at least 2.")
if (!dgp %in% available_dgps) {
  stop("DGP must be one of ", paste(available_dgps, collapse = ", "), ".")
}
if (!model_arg %in% available_models) {
  stop("model must be np, pl, or both.")
}
if (!error_design %in% available_error_designs) {
  stop("error design must be iid, hetero, or hetero_ar1.")
}
if (is.na(n_grid) || n_grid < 5L) stop("n_grid must be at least 5.")
if (is.na(n_sup_draws) || n_sup_draws < 99L) stop("n_sup_draws must be at least 99.")
if (is.na(rho) || abs(rho) >= 1) {
  stop("rho must be numeric with absolute value less than 1.")
}
if (!df_correction %in% c("substantive", "full")) {
  stop("df correction must be substantive or full.")
}
if (is.na(mc_seed)) stop("Monte Carlo seed must be an integer.")
if (!nzchar(chunk_label)) stop("A nonempty chunk label is required.")

model_list <- if (model_arg == "both") c("np", "pl") else model_arg
methods <- c("bam_fe", "bam_re", "bam_fd")
vcov_types <- c(
  "classic_cluster",
  "penalty_cluster",
  "penalty_cluster_delta",
  "mgcv_sandwich_freq",
  "mgcv_sandwich_default"
)

beta_true <- 1
sigma_mu <- 1.25
sigma_e <- 0.5
k_smooth <- 20
alpha <- 0.05
z_crit <- qnorm(1 - alpha / 2)

x_grid <- seq(0.05, 0.95, length.out = n_grid)
selected_points <- c(0.25, 0.50, 0.75)
selected_grid_idx <- vapply(
  selected_points,
  function(x0) which.min(abs(x_grid - x0)),
  integer(1)
)

## Common Gaussian draws make the uniform critical-value approximation paired
## across covariance estimators without changing the Monte Carlo DGP stream.
set.seed(20260506)
sup_z <- matrix(rnorm(n_sup_draws * n_grid), nrow = n_sup_draws, ncol = n_grid)

out_dir <- file.path(getwd(), "simulations", "inference_diagnostics", "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Monte Carlo draws M: ", M)
message("n units: ", n_units)
message("T periods: ", t_periods)
message("DGP: ", dgp)
message("Models: ", paste(model_list, collapse = ", "))
message("Grid points: ", n_grid)
message("Gaussian sup-t draws: ", n_sup_draws)
message("Error design: ", error_design)
if (error_design == "hetero_ar1") message("AR(1) rho: ", rho)
message("Small-sample df correction: ", df_correction)
message("Monte Carlo seed: ", mc_seed)
message("Chunk label: ", chunk_label)

g_fun <- function(x) {
  switch(
    dgp,
    linear = 0.5 + 1.4 * x,
    sin_half_pi = sin(0.5 * pi * x),
    sin_2pi = sin(2 * pi * x),
    sin_8pi = sin(8 * pi * x)
  )
}

truth_grid <- g_fun(x_grid)
truth_centered <- truth_grid - mean(truth_grid)

normalized_error_variance <- function(index, sigma_e) {
  raw_var <- 1 + 2 * exp(0.75 * index)
  sigma_e^2 * raw_var / mean(raw_var)
}

draw_panel_errors <- function(unit_id, t_periods, variance_index, sigma_e,
                              error_design, rho) {
  n_obs <- length(unit_id)
  if (error_design == "iid") {
    return(rnorm(n_obs, sd = sigma_e))
  }

  error_var <- normalized_error_variance(variance_index, sigma_e)
  if (error_design == "hetero") {
    return(rnorm(n_obs, sd = sqrt(error_var)))
  }

  epsilon <- numeric(n_obs)
  for (i in unique(unit_id)) {
    idx <- which(unit_id == i)
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

simulate_panel <- function() {
  unit_id <- rep(seq_len(n_units), each = t_periods)
  time <- rep(seq_len(t_periods), times = n_units)

  alpha_i <- rbeta(n_units, shape1 = 2.2, shape2 = 2.2)
  mu_i <- 1.1 * (alpha_i - mean(alpha_i)) + rnorm(n_units, sd = sigma_mu)

  x2 <- numeric(n_units * t_periods)
  x1 <- numeric(n_units * t_periods)
  for (i in seq_len(n_units)) {
    idx <- unit_id == i
    x2_i <- pmin(pmax(alpha_i[i] + rnorm(t_periods, sd = 0.18), 0), 1)
    x1_i <- 0.5 * x2_i + 0.5 * alpha_i[i] + rnorm(t_periods, sd = 0.45)
    x2[idx] <- x2_i
    x1[idx] <- x1_i
  }

  e_np <- draw_panel_errors(
    unit_id = unit_id,
    t_periods = t_periods,
    variance_index = x2,
    sigma_e = sigma_e,
    error_design = error_design,
    rho = rho
  )
  e_pl <- draw_panel_errors(
    unit_id = unit_id,
    t_periods = t_periods,
    variance_index = x1 + x2,
    sigma_e = sigma_e,
    error_design = error_design,
    rho = rho
  )

  data.frame(
    unit = factor(unit_id),
    unit_id = unit_id,
    time = time,
    x1 = x1,
    x2 = x2,
    y_np = g_fun(x2) + mu_i[unit_id] + e_np,
    y_pl = beta_true * x1 + g_fun(x2) + mu_i[unit_id] + e_pl
  )
}

make_first_differences <- function(df, model_type) {
  y_name <- if (model_type == "np") "y_np" else "y_pl"
  df <- df[order(df$unit_id, df$time), ]
  split_df <- split(df, df$unit_id)

  diffs <- lapply(split_df, function(d) {
    data.frame(
      unit = d$unit[-1],
      unit_id = d$unit_id[-1],
      time = d$time[-1],
      dy = diff(d[[y_name]]),
      dx1 = diff(d$x1),
      x2_now = d$x2[-1],
      x2_lag = d$x2[-nrow(d)]
    )
  })

  do.call(rbind, diffs)
}

smooth_columns <- function(fit, label) {
  hits <- vapply(fit$smooth, function(sm) identical(sm$label, label), logical(1))
  if (!any(hits)) {
    stop("Could not find smooth label '", label, "'. Labels are: ",
         paste(vapply(fit$smooth, `[[`, character(1), "label"), collapse = ", "))
  }
  sm <- fit$smooth[[which(hits)[1]]]
  sm$first.para:sm$last.para
}

is_unit_effect_coef <- function(coef_names) {
  startsWith(coef_names, "unit") | startsWith(coef_names, "s(unit).")
}

df_correction_edf <- function(fit) {
  if (df_correction == "full") return(sum(fit$edf))

  coef_names <- names(fit$edf)
  keep <- !is_unit_effect_coef(coef_names)
  sum(fit$edf[keep])
}

df_correction_rank <- function(fit, R) {
  if (df_correction == "full") return(qr(R)$rank)

  coef_names <- colnames(R)
  keep <- !is_unit_effect_coef(coef_names)
  qr(R[, keep, drop = FALSE])$rank
}

cluster_small_sample_correction <- function(G, nobs, df_used) {
  (G / (G - 1)) * ((nobs - 1) / max(nobs - df_used, 1))
}

cluster_vcov_penalty <- function(fit, cluster, R = NULL) {
  if (is.null(R)) R <- model.matrix(fit)
  u <- residuals(fit, type = "response")
  scale <- summary(fit)$scale
  bread_inv <- vcov(fit, freq = FALSE) / scale

  score <- R * as.numeric(u)
  score_g <- rowsum(score, group = as.factor(cluster), reorder = FALSE)
  meat <- crossprod(score_g)

  G <- nrow(score_g)
  nobs <- nrow(R)
  edf <- df_correction_edf(fit)
  corr <- cluster_small_sample_correction(G, nobs, edf)

  corr * bread_inv %*% meat %*% bread_inv
}

within_transform <- function(A, cluster) {
  A <- as.matrix(A)
  group <- as.factor(cluster)
  group_sums <- rowsum(A, group = group, reorder = FALSE)
  group_counts <- as.numeric(table(group))
  group_means <- group_sums / group_counts
  A - group_means[as.integer(group), , drop = FALSE]
}

residualize_ignoring_penalty <- function(target, R, target_cols, cluster) {
  coef_names <- colnames(R)
  unit_cols <- is_unit_effect_coef(coef_names)
  intercept_cols <- coef_names == "(Intercept)"
  nuisance_cols <- !(seq_len(ncol(R)) %in% target_cols)

  if (any(unit_cols)) {
    keep <- nuisance_cols & !unit_cols & !intercept_cols
    target_use <- within_transform(target, cluster)
    nuisance <- within_transform(R[, keep, drop = FALSE], cluster)
  } else {
    target_use <- as.matrix(target)
    nuisance <- R[, nuisance_cols, drop = FALSE]
  }

  if (!ncol(nuisance)) return(target_use)
  qr.resid(qr(nuisance), target_use)
}

safe_inverse <- function(A) {
  A <- (A + t(A)) / 2
  ans <- try(solve(A), silent = TRUE)
  if (!inherits(ans, "try-error")) return(ans)

  eig <- eigen(A, symmetric = TRUE)
  tol <- max(dim(A)) * max(abs(eig$values)) * .Machine$double.eps
  inv_values <- ifelse(abs(eig$values) > tol, 1 / eig$values, 0)
  eig$vectors %*% (inv_values * t(eig$vectors))
}

cluster_vcov_classic_smooth <- function(fit, cluster, smooth_cols, R = NULL) {
  if (is.null(R)) R <- model.matrix(fit)
  u <- residuals(fit, type = "response")

  S_tilde <- residualize_ignoring_penalty(
    R[, smooth_cols, drop = FALSE], R, smooth_cols, cluster
  )

  bread <- crossprod(S_tilde)
  score <- S_tilde * as.numeric(u)
  score_g <- rowsum(score, group = as.factor(cluster), reorder = FALSE)
  meat <- crossprod(score_g)

  G <- nrow(score_g)
  nobs <- nrow(R)
  rank_r <- df_correction_rank(fit, R)
  corr <- cluster_small_sample_correction(G, nobs, rank_r)

  bread_inv <- safe_inverse(bread)
  ans <- corr * bread_inv %*% meat %*% bread_inv

  if (M <= 5L && any(is_unit_effect_coef(colnames(R)))) {
    S_full <- qr.resid(
      qr(R[, -smooth_cols, drop = FALSE]),
      R[, smooth_cols, drop = FALSE]
    )
    score_g_full <- rowsum(
      S_full * as.numeric(u),
      group = as.factor(cluster),
      reorder = FALSE
    )
    bread_inv_full <- safe_inverse(crossprod(S_full))
    full_ans <- corr * bread_inv_full %*%
      crossprod(score_g_full) %*% bread_inv_full
    if (max(abs(ans - full_ans)) > 1e-8) {
      stop("Within-transformed classic smooth covariance identity failed.")
    }
  }

  ans
}

symmetrize <- function(A) {
  (A + t(A)) / 2
}

cluster_projected_cov <- function(fit, cluster, L, R = NULL) {
  if (is.null(R)) R <- model.matrix(fit)
  u <- residuals(fit, type = "response")
  scale <- summary(fit)$scale
  bread_inv <- vcov(fit, freq = FALSE) / scale

  score <- R * as.numeric(u)
  score_g <- rowsum(score, group = as.factor(cluster), reorder = FALSE)
  projected_score <- score_g %*% t(L %*% bread_inv)

  G <- nrow(score_g)
  nobs <- nrow(R)
  corr <- cluster_small_sample_correction(G, nobs, df_correction_edf(fit))
  symmetrize(corr * crossprod(projected_score))
}

mgcv_projected_cov <- function(fit, L, freq, R = NULL) {
  if (is.null(R)) R <- model.matrix(fit)
  mu <- fit$fitted.values
  w <- fit$family$mu.eta(fit$linear.predictors) *
    (fit$y - mu) / (fit$sig2 * fit$family$variance(mu))
  inflation <- nrow(R) / (nrow(R) - sum(fit$edf))
  projected_score <- (w * R) %*% t(L %*% fit$Vp)
  B2 <- if (freq) {
    matrix(0, nrow(L), nrow(L))
  } else {
    L %*% (fit$Vp - fit$Ve) %*% t(L)
  }
  symmetrize(inflation * crossprod(projected_score) + B2)
}

delta_hat <- function(fit) {
  symmetrize(vcov(fit, freq = FALSE) - vcov(fit, freq = TRUE))
}

grid_cov_to_supcrit <- function(C, draws = n_sup_draws) {
  C <- (C + t(C)) / 2
  se <- sqrt(pmax(diag(C), 0))
  valid <- is.finite(se) & se > 1e-10
  if (!any(valid)) return(NA_real_)

  C_std <- C[valid, valid, drop = FALSE] / outer(se[valid], se[valid])
  C_std <- (C_std + t(C_std)) / 2
  eig <- eigen(C_std, symmetric = TRUE)
  vals <- pmax(eig$values, 0)
  A <- eig$vectors %*% diag(sqrt(vals), nrow = length(vals))
  z <- sup_z[seq_len(draws), seq_len(length(vals)), drop = FALSE]
  sim <- z %*% t(A)

  as.numeric(quantile(apply(abs(sim), 1, max), probs = 1 - alpha, names = FALSE))
}

make_grid_lpmatrix <- function(fit_info, model_type, method) {
  fit <- fit_info$fit
  if (method %in% c("bam_fe", "bam_re")) {
    nd <- data.frame(
      x1 = 0,
      x2 = x_grid,
      unit = factor(fit_info$reference_unit, levels = fit_info$unit_levels)
    )
    Rg <- predict(fit, newdata = nd, type = "lpmatrix")
    smooth_cols <- smooth_columns(fit, "s(x2)")
  } else {
    nd <- data.frame(
      dx1 = 0,
      X_pair = x_grid,
      L_pair = 1
    )
    Rg <- predict(fit, newdata = nd, type = "lpmatrix")
    smooth_label <- vapply(fit$smooth, `[[`, character(1), "label")[1]
    smooth_cols <- smooth_columns(fit, smooth_label)
  }

  L_full <- Rg
  L_full[, -smooth_cols] <- 0
  L_full <- sweep(L_full, 2, colMeans(L_full), "-")

  L_smooth <- Rg[, smooth_cols, drop = FALSE]
  L_smooth <- sweep(L_smooth, 2, colMeans(L_smooth), "-")

  list(L_full = L_full, L_smooth = L_smooth, smooth_cols = smooth_cols)
}

fit_bam_method <- function(df, model_type, method) {
  y_name <- if (model_type == "np") "y_np" else "y_pl"
  reference_unit <- levels(df$unit)[1]

  if (method == "bam_fe") {
    rhs <- if (model_type == "np") {
      paste0("s(x2, k = ", k_smooth, ") + unit")
    } else {
      paste0("x1 + s(x2, k = ", k_smooth, ") + unit")
    }
    fit <- bam(
      as.formula(paste(y_name, "~", rhs)),
      data = df,
      method = "fREML",
      select = TRUE,
      discrete = TRUE
    )
    return(list(
      fit = fit,
      cluster = df$unit_id,
      reference_unit = reference_unit,
      unit_levels = levels(df$unit)
    ))
  }

  if (method == "bam_re") {
    rhs <- if (model_type == "np") {
      paste0("s(x2, k = ", k_smooth, ") + s(unit, bs = 're')")
    } else {
      paste0("x1 + s(x2, k = ", k_smooth, ") + s(unit, bs = 're')")
    }
    fit <- bam(
      as.formula(paste(y_name, "~", rhs)),
      data = df,
      method = "fREML",
      select = TRUE,
      discrete = TRUE
    )
    return(list(
      fit = fit,
      cluster = df$unit_id,
      reference_unit = reference_unit,
      unit_levels = levels(df$unit)
    ))
  }

  if (method == "bam_fd") {
    df_fd <- make_first_differences(df, model_type)
    X_pair <- cbind(df_fd$x2_now, df_fd$x2_lag)
    L_pair <- cbind(rep(1, nrow(df_fd)), rep(-1, nrow(df_fd)))
    dat <- list(
      dy = df_fd$dy,
      dx1 = df_fd$dx1,
      X_pair = X_pair,
      L_pair = L_pair
    )
    rhs <- if (model_type == "np") {
      paste0("s(X_pair, by = L_pair, k = ", k_smooth, ")")
    } else {
      paste0("dx1 + s(X_pair, by = L_pair, k = ", k_smooth, ")")
    }
    fit <- bam(
      as.formula(paste("dy ~", rhs)),
      data = dat,
      method = "fREML",
      select = TRUE,
      discrete = TRUE
    )
    return(list(
      fit = fit,
      cluster = df_fd$unit_id,
      reference_unit = reference_unit,
      unit_levels = levels(df$unit)
    ))
  }

  stop("Unknown method: ", method)
}

evaluate_fit_all <- function(fit_info, model_type, method) {
  grid <- make_grid_lpmatrix(fit_info, model_type, method)
  fit <- fit_info$fit
  R <- model.matrix(fit)

  estimate <- as.numeric(grid$L_full %*% coef(fit))

  Delta <- delta_hat(fit)
  V_smooth_classic <- cluster_vcov_classic_smooth(
    fit, fit_info$cluster, grid$smooth_cols, R
  )
  C_penalty_cluster <- cluster_projected_cov(
    fit, fit_info$cluster, grid$L_full, R
  )
  C_delta <- symmetrize(grid$L_full %*% Delta %*% t(grid$L_full))
  C_mgcv_freq <- mgcv_projected_cov(
    fit, grid$L_full, freq = TRUE, R = R
  )

  covariance_grid <- list(
    classic_cluster = grid$L_smooth %*% V_smooth_classic %*% t(grid$L_smooth),
    penalty_cluster = C_penalty_cluster,
    penalty_cluster_delta = C_penalty_cluster + C_delta,
    mgcv_sandwich_freq = C_mgcv_freq,
    mgcv_sandwich_default = C_mgcv_freq + C_delta
  )

  if (M <= 5L) {
    V_penalty_full <- cluster_vcov_penalty(fit, fit_info$cluster, R)
    checks <- list(
      covariance_grid$penalty_cluster -
        grid$L_full %*% V_penalty_full %*% t(grid$L_full),
      covariance_grid$mgcv_sandwich_freq - grid$L_full %*%
        vcov(fit, sandwich = TRUE, freq = TRUE) %*% t(grid$L_full),
      covariance_grid$mgcv_sandwich_default - grid$L_full %*%
        vcov(fit, sandwich = TRUE, freq = FALSE) %*% t(grid$L_full)
    )
    if (max(abs(unlist(checks))) > 1e-8) {
      stop("Projected grid covariance identity failed.")
    }
  }

  do.call(rbind, lapply(vcov_types, function(vcov_type) {
    C <- symmetrize(covariance_grid[[vcov_type]])
    se <- sqrt(pmax(diag(C), 0))
    point_lower <- estimate - z_crit * se
    point_upper <- estimate + z_crit * se
    point_covered <- truth_centered >= point_lower & truth_centered <= point_upper

    sup_crit <- grid_cov_to_supcrit(C)
    uniform_lower <- estimate - sup_crit * se
    uniform_upper <- estimate + sup_crit * se
    uniform_covered <- all(
      truth_centered >= uniform_lower & truth_centered <= uniform_upper
    )

    data.frame(
      model = model_type,
      method = method,
      vcov_type = vcov_type,
      x = x_grid,
      truth = truth_centered,
      estimate = estimate,
      se = se,
      point_covered = as.integer(point_covered),
      point_width = point_upper - point_lower,
      sup_crit = sup_crit,
      uniform_covered = as.integer(uniform_covered),
      uniform_width = uniform_upper - uniform_lower
    )
  }))
}

one_draw <- function(sim_id) {
  df <- simulate_panel()
  out <- list()
  idx <- 0L

  for (model_type in model_list) {
    for (method in methods) {
      fit_info <- fit_bam_method(df, model_type, method)
      idx <- idx + 1L
      ans <- evaluate_fit_all(fit_info, model_type, method)
      ans$sim <- sim_id
      out[[idx]] <- ans
    }
  }

  do.call(rbind, out)
}

set.seed(mc_seed)
results <- vector("list", M)

for (m in seq_len(M)) {
  if (m == 1L || m %% 10L == 0L || m == M) {
    message("Monte Carlo draw ", m, " / ", M)
  }

  ans <- try(one_draw(m), silent = TRUE)
  if (inherits(ans, "try-error")) {
    warning("Draw ", m, " failed: ", conditionMessage(attr(ans, "condition")))
    ans <- NULL
  }
  results[[m]] <- ans
}

results <- do.call(rbind, results)
results$error_design <- error_design
results$rho <- if (error_design == "hetero_ar1") rho else NA_real_
results$df_correction <- df_correction
results$M <- M
results$n <- n_units
results$T <- t_periods
results$K <- k_smooth
results$dgp <- dgp
results$mc_seed <- mc_seed
results$chunk <- chunk_label
results <- results[, c("sim", "model", "method", "vcov_type", setdiff(names(results), c("sim", "model", "method", "vcov_type")))]

by_x_mean <- aggregate(
  cbind(point_covered, estimate, se) ~ model + method + vcov_type + x,
  data = results,
  FUN = mean
)
names(by_x_mean)[names(by_x_mean) == "point_covered"] <- "point_coverage"
names(by_x_mean)[names(by_x_mean) == "estimate"] <- "mean_estimate"
names(by_x_mean)[names(by_x_mean) == "se"] <- "mean_reported_se"

by_x_sd <- aggregate(
  estimate ~ model + method + vcov_type + x,
  data = results,
  FUN = sd
)
names(by_x_sd)[names(by_x_sd) == "estimate"] <- "empirical_sd"

coverage_by_x <- merge(
  by_x_mean,
  by_x_sd,
  by = c("model", "method", "vcov_type", "x")
)
coverage_by_x$truth <- truth_centered[match(coverage_by_x$x, x_grid)]
coverage_by_x$bias <- coverage_by_x$mean_estimate - coverage_by_x$truth
coverage_by_x$se_to_sd <- coverage_by_x$mean_reported_se / coverage_by_x$empirical_sd
coverage_by_x$abs_bias <- abs(coverage_by_x$bias)
coverage_by_x$error_design <- error_design
coverage_by_x$rho <- if (error_design == "hetero_ar1") rho else NA_real_
coverage_by_x$df_correction <- df_correction
coverage_by_x$M <- M
coverage_by_x$n <- n_units
coverage_by_x$T <- t_periods
coverage_by_x$K <- k_smooth
coverage_by_x$dgp <- dgp

selected_coverage <- do.call(
  rbind,
  lapply(seq_along(selected_points), function(j) {
    d <- coverage_by_x[abs(coverage_by_x$x - x_grid[selected_grid_idx[j]]) < 1e-12, ]
    d$selected_x <- selected_points[j]
    d$grid_x <- x_grid[selected_grid_idx[j]]
    d[, c("model", "method", "vcov_type", "selected_x", "grid_x", "point_coverage")]
  })
)

selected_wide <- reshape(
  selected_coverage[, c("model", "method", "vcov_type", "selected_x", "point_coverage")],
  idvar = c("model", "method", "vcov_type"),
  timevar = "selected_x",
  direction = "wide"
)

summary_base <- aggregate(
  cbind(point_covered, point_width, uniform_covered, uniform_width, sup_crit) ~ model + method + vcov_type,
  data = results,
  FUN = mean
)
names(summary_base)[names(summary_base) == "point_covered"] <- "avg_pointwise_coverage"
names(summary_base)[names(summary_base) == "point_width"] <- "avg_pointwise_width"
names(summary_base)[names(summary_base) == "uniform_covered"] <- "uniform_coverage"
names(summary_base)[names(summary_base) == "uniform_width"] <- "avg_uniform_width"
names(summary_base)[names(summary_base) == "sup_crit"] <- "avg_sup_crit"

min_coverage <- aggregate(
  point_coverage ~ model + method + vcov_type,
  data = coverage_by_x,
  FUN = min
)
names(min_coverage)[names(min_coverage) == "point_coverage"] <- "min_pointwise_coverage"

diagnostic_means <- aggregate(
  cbind(abs_bias, empirical_sd, mean_reported_se, se_to_sd) ~
    model + method + vcov_type,
  data = coverage_by_x,
  FUN = mean
)
max_abs_bias <- aggregate(
  abs_bias ~ model + method + vcov_type,
  data = coverage_by_x,
  FUN = max
)
names(max_abs_bias)[names(max_abs_bias) == "abs_bias"] <- "max_abs_bias"

n_success <- aggregate(
  sim ~ model + method + vcov_type,
  data = unique(results[, c("sim", "model", "method", "vcov_type")]),
  FUN = length
)
names(n_success)[names(n_success) == "sim"] <- "n_success"

summary_results <- merge(summary_base, min_coverage, by = c("model", "method", "vcov_type"))
summary_results <- merge(summary_results, selected_wide, by = c("model", "method", "vcov_type"), all.x = TRUE)
summary_results <- merge(summary_results, diagnostic_means, by = c("model", "method", "vcov_type"), all.x = TRUE)
summary_results <- merge(summary_results, max_abs_bias, by = c("model", "method", "vcov_type"), all.x = TRUE)
summary_results <- merge(summary_results, n_success, by = c("model", "method", "vcov_type"), all.x = TRUE)
names(summary_results) <- sub("point_coverage\\.", "coverage_at_", names(summary_results))
summary_results <- summary_results[
  order(summary_results$model, summary_results$method, match(summary_results$vcov_type, vcov_types)),
]
row.names(summary_results) <- NULL
summary_results$error_design <- error_design
summary_results$rho <- if (error_design == "hetero_ar1") rho else NA_real_
summary_results$df_correction <- df_correction
summary_results$M <- M
summary_results$n <- n_units
summary_results$T <- t_periods
summary_results$K <- k_smooth
summary_results$dgp <- dgp
summary_results$mc_seed <- mc_seed
summary_results$chunk <- chunk_label

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

file_suffix <- paste0(
  "M", M,
  "_n", n_units,
  "_T", t_periods,
  "_", dgp,
  "_", paste(model_list, collapse = "_"),
  "_J", n_grid,
  error_suffix,
  "_df", df_correction,
  "_chunk", chunk_label
)

results_path <- file.path(out_dir, paste0("panel_fe_g_diagnostics_grid_", file_suffix, ".csv.gz"))
coverage_x_path <- file.path(out_dir, paste0("panel_fe_g_diagnostics_by_x_", file_suffix, ".csv"))
summary_path <- file.path(out_dir, paste0("panel_fe_g_diagnostics_summary_", file_suffix, ".csv"))
plot_path <- file.path(out_dir, paste0("panel_fe_g_diagnostics_coverage_curve_", file_suffix, ".png"))
se_ratio_plot_path <- file.path(out_dir, paste0("panel_fe_g_diagnostics_se_ratio_", file_suffix, ".png"))

write.csv(results, gzfile(results_path), row.names = FALSE)
write.csv(coverage_by_x, coverage_x_path, row.names = FALSE)
write.csv(summary_results, summary_path, row.names = FALSE)

coverage_plot <- ggplot(
  coverage_by_x,
  aes(x = x, y = point_coverage, color = vcov_type)
) +
  geom_hline(yintercept = 0.95, color = "grey45", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  facet_grid(model ~ method) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  labs(
    title = "Pointwise coverage for centered g(x)",
    subtitle = paste0(
      "DGP: ", dgp,
      "; errors: ", error_design,
      if (error_design == "hetero_ar1") paste0(" (rho = ", rho, ")") else "",
      "; M = ", M,
      "; n = ", n_units,
      "; T = ", t_periods
    ),
    x = "x",
    y = "pointwise coverage",
    color = "covariance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title.position = "plot",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(plot_path, coverage_plot, width = 9, height = 6.5, dpi = 180, bg = "white")

se_ratio_plot <- ggplot(
  coverage_by_x,
  aes(x = x, y = se_to_sd, color = vcov_type)
) +
  geom_hline(yintercept = 1, color = "grey45", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  facet_grid(model ~ method) +
  labs(
    title = "Reported standard error relative to Monte Carlo standard deviation",
    subtitle = paste0(
      "DGP: ", dgp,
      "; errors: ", error_design,
      if (error_design == "hetero_ar1") paste0(" (rho = ", rho, ")") else "",
      "; M = ", M,
      "; n = ", n_units,
      "; T = ", t_periods
    ),
    x = "x",
    y = "mean reported SE / empirical SD",
    color = "covariance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title.position = "plot",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(se_ratio_plot_path, se_ratio_plot, width = 9, height = 6.5, dpi = 180, bg = "white")

message("\nSummary:")
print(summary_results, row.names = FALSE)

message("\nWrote:")
message("  ", results_path)
message("  ", coverage_x_path)
message("  ", summary_path)
message("  ", plot_path)
message("  ", se_ratio_plot_path)
