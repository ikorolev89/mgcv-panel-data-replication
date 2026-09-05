# Diagnostic inference simulation for beta in a semiparametric FE panel model:
#
#   y_it = beta * x1_it + g(x2_it) + mu_i + e_it.
#
# Methods:
#   1. bam FE: y ~ x1 + s(x2) + unit
#   2. bam RE: y ~ x1 + s(x2) + s(unit, bs = "re")
#   3. bam FD: dy ~ dx1 + [g(x2_now) - g(x2_lag)]
#
# Covariance estimators:
#   1. unit-cluster sandwich ignoring penalties
#   2. penalty-adjusted unit-cluster sandwich
#   3. penalty-adjusted cluster sandwich plus V_B - V_F
#   4. mgcv observation-level sandwich with freq = TRUE
#   5. default mgcv observation-level sandwich (includes V_B - V_F)
#
# Usage:
#   Rscript simulations/inference_diagnostics/panel_fe_beta_inference_diagnostics.R \
#     1000 0 200 4 baseline 1 1 hetero
#   Rscript simulations/inference_diagnostics/panel_fe_beta_inference_diagnostics.R \
#     1000 0 200 4 baseline 1.075 1 hetero_ar1 0.5
#
# Arguments:
#   1. Monte Carlo draws M
#   2. bootstrap draws B; use 0 to disable bootstrap
#   3. number of units n
#   4. number of periods T
#   5. DGP: baseline, linear_g, or sine_g
#   6. true beta used in the DGP; default 1
#   7. null beta used in tests/bootstrap; default 1
#   8. error design: iid, hetero, or hetero_ar1; default iid
#   9. AR(1) rho for hetero_ar1; default 0.5
#   10. small-sample df correction: substantive or full; default substantive

suppressPackageStartupMessages({
  library(mgcv)
})

args <- commandArgs(trailingOnly = TRUE)
M <- if (length(args) >= 1) as.integer(args[[1]]) else 50L
B <- if (length(args) >= 2) as.integer(args[[2]]) else 0L
n_units <- if (length(args) >= 3) as.integer(args[[3]]) else 200L
t_periods <- if (length(args) >= 4) as.integer(args[[4]]) else 2L
dgp <- if (length(args) >= 5) args[[5]] else "baseline"
beta_true <- if (length(args) >= 6) as.numeric(args[[6]]) else 1
beta_null <- if (length(args) >= 7) as.numeric(args[[7]]) else 1
error_design <- if (length(args) >= 8) args[[8]] else "iid"
rho <- if (length(args) >= 9) as.numeric(args[[9]]) else 0.5
df_correction <- if (length(args) >= 10) args[[10]] else "substantive"

if (is.na(M) || M < 1L) stop("M must be a positive integer.")
if (is.na(B) || B != 0L) stop("This diagnostic script requires B = 0.")
if (is.na(n_units) || n_units < 2L) stop("n_units must be at least 2.")
if (is.na(t_periods) || t_periods < 2L) stop("T must be at least 2.")
if (is.na(beta_true)) stop("beta_true must be numeric.")
if (is.na(beta_null)) stop("beta_null must be numeric.")
if (!dgp %in% c("baseline", "linear_g", "sine_g")) {
  stop("DGP must be baseline, linear_g, or sine_g.")
}
if (!error_design %in% c("iid", "hetero", "hetero_ar1")) {
  stop("error design must be iid, hetero, or hetero_ar1.")
}
if (is.na(rho) || abs(rho) >= 1) {
  stop("rho must be numeric with absolute value less than 1.")
}
if (!df_correction %in% c("substantive", "full")) {
  stop("df correction must be substantive or full.")
}

sigma_mu <- 1.25
sigma_e <- 0.5
k_smooth <- 20
methods <- c("bam_fe", "bam_re", "bam_fd")
vcov_types <- c(
  "classic_cluster",
  "penalty_cluster",
  "penalty_cluster_delta",
  "mgcv_sandwich_freq",
  "mgcv_sandwich_default"
)

out_dir <- file.path(getwd(), "simulations", "inference_diagnostics", "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Monte Carlo draws M: ", M)
message("Bootstrap draws B: ", B)
message("n units: ", n_units)
message("T periods: ", t_periods)
message("DGP: ", dgp)
message("True beta: ", beta_true)
message("Null beta: ", beta_null)
message("Error design: ", error_design)
if (error_design == "hetero_ar1") message("AR(1) rho: ", rho)
message("Small-sample df correction: ", df_correction)

g_fun <- function(x) {
  if (dgp == "linear_g") return(0.5 + 1.2 * x)
  sin(2 * pi * x)
}

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

  e <- draw_panel_errors(
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
    y = beta_true * x1 + g_fun(x2) + mu_i[unit_id] + e
  )
}

make_first_differences <- function(df) {
  df <- df[order(df$unit_id, df$time), ]
  split_df <- split(df, df$unit_id)

  diffs <- lapply(split_df, function(d) {
    data.frame(
      unit = d$unit[-1],
      unit_id = d$unit_id[-1],
      time = d$time[-1],
      dy = diff(d$y),
      dx1 = diff(d$x1),
      x2_now = d$x2[-1],
      x2_lag = d$x2[-nrow(d)]
    )
  })

  do.call(rbind, diffs)
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

cluster_vcov_classic_beta <- function(fit, beta_name, cluster, R = NULL) {
  if (is.null(R)) R <- model.matrix(fit)
  u <- residuals(fit, type = "response")
  beta_idx <- match(beta_name, names(coef(fit)))

  x_tilde <- as.numeric(residualize_ignoring_penalty(
    R[, beta_idx, drop = FALSE], R, beta_idx, cluster
  ))

  bread <- sum(x_tilde^2)
  score_g <- rowsum(x_tilde * as.numeric(u), group = as.factor(cluster), reorder = FALSE)
  meat <- sum(score_g^2)

  G <- nrow(score_g)
  nobs <- nrow(R)
  rank_r <- df_correction_rank(fit, R)
  corr <- cluster_small_sample_correction(G, nobs, rank_r)

  ans <- corr * meat / bread^2

  if (M <= 5L && any(is_unit_effect_coef(colnames(R)))) {
    x_tilde_full <- as.numeric(qr.resid(
      qr(R[, -beta_idx, drop = FALSE]), R[, beta_idx]
    ))
    score_g_full <- rowsum(
      x_tilde_full * as.numeric(u),
      group = as.factor(cluster),
      reorder = FALSE
    )
    full_ans <- corr * sum(score_g_full^2) / sum(x_tilde_full^2)^2
    if (abs(ans - full_ans) > 1e-9) {
      stop("Within-transformed classic beta covariance identity failed.")
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
  symmetrize(
    vcov(fit, sandwich = FALSE, freq = FALSE) -
      vcov(fit, sandwich = FALSE, freq = TRUE)
  )
}

beta_diagnostic_results <- function(fit, beta_name, cluster) {
  beta_hat <- unname(coef(fit)[beta_name])
  beta_idx <- match(beta_name, names(coef(fit)))
  R <- model.matrix(fit)
  L_beta <- matrix(0, nrow = 1L, ncol = length(coef(fit)))
  L_beta[1L, beta_idx] <- 1

  V_penalty_cluster <- cluster_projected_cov(fit, cluster, L_beta, R)
  Delta <- delta_hat(fit)
  Delta_beta <- as.numeric(L_beta %*% Delta %*% t(L_beta))
  V_cluster_delta <- V_penalty_cluster + Delta_beta
  V_mgcv_freq <- mgcv_projected_cov(fit, L_beta, freq = TRUE, R = R)
  V_mgcv_default <- V_mgcv_freq + Delta_beta

  if (M <= 5L) {
    V_penalty_full <- cluster_vcov_penalty(fit, cluster, R)
    checks <- c(
      V_penalty_cluster[1L, 1L] - V_penalty_full[beta_idx, beta_idx],
      V_mgcv_freq[1L, 1L] - vcov(fit, sandwich = TRUE, freq = TRUE)[beta_idx, beta_idx],
      V_mgcv_default[1L, 1L] - vcov(fit, sandwich = TRUE, freq = FALSE)[beta_idx, beta_idx]
    )
    if (max(abs(checks)) > 1e-9) stop("Projected beta covariance identity failed.")
  }

  variances <- c(
    classic_cluster = cluster_vcov_classic_beta(fit, beta_name, cluster, R),
    penalty_cluster = V_penalty_cluster[1L, 1L],
    penalty_cluster_delta = V_cluster_delta[1L, 1L],
    mgcv_sandwich_freq = V_mgcv_freq[1L, 1L],
    mgcv_sandwich_default = V_mgcv_default[1L, 1L]
  )

  se <- sqrt(pmax(variances, 0))
  t_stat <- (beta_hat - beta_null) / se
  covered <- abs(beta_hat - beta_true) <= qnorm(0.975) * se

  data.frame(
    vcov_type = names(variances),
    beta_hat = beta_hat,
    variance = as.numeric(variances),
    se = as.numeric(se),
    t_stat = as.numeric(t_stat),
    reject = as.integer(abs(t_stat) > qnorm(0.975)),
    beta_covered = as.integer(covered),
    delta_beta_variance = Delta_beta
  )
}

fit_unrestricted <- function(df, method, y_name = "y", dy_name = "dy") {
  if (method == "bam_fe") {
    f <- as.formula(paste0(y_name, " ~ x1 + s(x2, k = ", k_smooth, ") + unit"))
    fit <- bam(f, data = df, method = "fREML", select = TRUE, discrete = TRUE)
    return(list(fit = fit, beta_name = "x1", cluster = df$unit_id))
  }

  if (method == "bam_re") {
    f <- as.formula(paste0(y_name, " ~ x1 + s(x2, k = ", k_smooth, ") + s(unit, bs = 're')"))
    fit <- bam(f, data = df, method = "fREML", select = TRUE, discrete = TRUE)
    return(list(fit = fit, beta_name = "x1", cluster = df$unit_id))
  }

  if (method == "bam_fd") {
    X_pair <- cbind(df$x2_now, df$x2_lag)
    L_pair <- cbind(rep(1, nrow(df)), rep(-1, nrow(df)))
    dat <- list(
      dx1 = df$dx1,
      X_pair = X_pair,
      L_pair = L_pair
    )
    dat[[dy_name]] <- df[[dy_name]]
    fit <- bam(
      as.formula(paste0(dy_name, " ~ dx1 + s(X_pair, by = L_pair, k = ", k_smooth, ")")),
      data = dat,
      method = "fREML",
      select = TRUE,
      discrete = TRUE
    )
    return(list(fit = fit, beta_name = "dx1", cluster = df$unit_id))
  }

  stop("Unknown method: ", method)
}

fit_restricted <- function(df, method) {
  if (method == "bam_fe") {
    df$y_null <- df$y - beta_null * df$x1
    fit <- bam(
      y_null ~ s(x2, k = k_smooth) + unit,
      data = df,
      method = "fREML",
      select = TRUE,
      discrete = TRUE
    )
    return(list(fit = fit, fitted_null = beta_null * df$x1 + fitted(fit), resid_null = residuals(fit)))
  }

  if (method == "bam_re") {
    df$y_null <- df$y - beta_null * df$x1
    fit <- bam(
      y_null ~ s(x2, k = k_smooth) + s(unit, bs = "re"),
      data = df,
      method = "fREML",
      select = TRUE,
      discrete = TRUE
    )
    return(list(fit = fit, fitted_null = beta_null * df$x1 + fitted(fit), resid_null = residuals(fit)))
  }

  if (method == "bam_fd") {
    df$dy_null <- df$dy - beta_null * df$dx1
    X_pair <- cbind(df$x2_now, df$x2_lag)
    L_pair <- cbind(rep(1, nrow(df)), rep(-1, nrow(df)))
    fit <- bam(
      dy_null ~ s(X_pair, by = L_pair, k = k_smooth),
      data = list(
        dy_null = df$dy_null,
        X_pair = X_pair,
        L_pair = L_pair
      ),
      method = "fREML",
      select = TRUE,
      discrete = TRUE
    )
    return(list(fit = fit, fitted_null = beta_null * df$dx1 + fitted(fit), resid_null = residuals(fit)))
  }

  stop("Unknown method: ", method)
}

one_method_diagnostics <- function(df, method) {
  dfit <- if (method == "bam_fd") make_first_differences(df) else df
  unrestricted <- fit_unrestricted(dfit, method)

  ans <- beta_diagnostic_results(
    unrestricted$fit,
    unrestricted$beta_name,
    unrestricted$cluster
  )
  ans$method <- method
  ans[, c("method", "vcov_type", setdiff(names(ans), c("method", "vcov_type")))]
}

set.seed(20260504)
results <- vector("list", M * length(methods))
idx <- 0L

for (m in seq_len(M)) {
  if (m == 1L || m %% 10L == 0L || m == M) {
    message("Monte Carlo draw ", m, " / ", M)
  }

  df <- simulate_panel()
  for (method in methods) {
    idx <- idx + 1L
    ans <- try(one_method_diagnostics(df, method), silent = TRUE)
    if (inherits(ans, "try-error")) {
      ans <- data.frame(
        method = method,
        vcov_type = vcov_types,
        beta_hat = NA_real_,
        variance = NA_real_,
        se = NA_real_,
        t_stat = NA_real_,
        reject = NA_integer_,
        beta_covered = NA_integer_,
        delta_beta_variance = NA_real_
      )
    }
    ans$sim <- m
    results[[idx]] <- ans
  }
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
results$beta_true <- beta_true
results$beta_null <- beta_null
results <- results[, c("sim", "method", "vcov_type", setdiff(names(results), c("sim", "method", "vcov_type")))]

summary_results <- do.call(
  rbind,
  lapply(split(results, list(results$method, results$vcov_type), drop = TRUE), function(d) {
    data.frame(
      method = d$method[1],
      vcov_type = d$vcov_type[1],
      n_success = sum(is.finite(d$beta_hat)),
      rejection_frequency = mean(d$reject, na.rm = TRUE),
      beta_coverage = mean(d$beta_covered, na.rm = TRUE),
      mean_beta_hat = mean(d$beta_hat, na.rm = TRUE),
      sd_beta_hat = sd(d$beta_hat, na.rm = TRUE),
      mean_se = mean(d$se, na.rm = TRUE),
      se_to_sd = mean(d$se, na.rm = TRUE) / sd(d$beta_hat, na.rm = TRUE),
      mean_delta_beta_variance = mean(d$delta_beta_variance, na.rm = TRUE)
    )
  })
)
row.names(summary_results) <- NULL
summary_results <- summary_results[order(summary_results$method, match(summary_results$vcov_type, vcov_types)), ]
summary_results$error_design <- error_design
summary_results$rho <- if (error_design == "hetero_ar1") rho else NA_real_
summary_results$df_correction <- df_correction
summary_results$M <- M
summary_results$n <- n_units
summary_results$T <- t_periods
summary_results$K <- k_smooth
summary_results$dgp <- dgp
summary_results$beta_true <- beta_true
summary_results$beta_null <- beta_null

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

suffix <- paste0(
  "M", M,
  "_B", B,
  "_n", n_units,
  "_T", t_periods,
  "_", dgp,
  "_btrue", format_suffix_num(beta_true),
  "_bnull", format_suffix_num(beta_null),
  error_suffix,
  "_df", df_correction
)
results_path <- file.path(out_dir, paste0("panel_fe_beta_diagnostics_results_", suffix, ".csv"))
summary_path <- file.path(out_dir, paste0("panel_fe_beta_diagnostics_summary_", suffix, ".csv"))

write.csv(results, results_path, row.names = FALSE)
write.csv(summary_results, summary_path, row.names = FALSE)

message("\nSummary:")
print(summary_results, row.names = FALSE)

message("\nWrote:")
message("  ", results_path)
message("  ", summary_path)
