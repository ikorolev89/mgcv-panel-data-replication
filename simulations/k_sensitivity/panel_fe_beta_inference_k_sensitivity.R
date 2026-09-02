# Inference-only simulation for beta in a semiparametric FE panel model:
#
#   y_it = beta * x1_it + g(x2_it) + mu_i + e_it.
#
# Methods:
#   1. bam FE: y ~ x1 + s(x2) + unit
#   2. bam RE: y ~ x1 + s(x2) + s(unit, bs = "re")
#   3. bam FD: dy ~ dx1 + [g(x2_now) - g(x2_lag)]
#
# Inference:
#   - Baltagi-Li/An-Hsiao-Li-style unit-cluster sandwich, ignoring penalties
#   - penalty-adjusted unit-cluster sandwich t-test
#   - optional unit wild-cluster bootstrap p-value based on the bootstrap t-statistic
#
# Usage:
#   Rscript simulations/panel_fe_beta_inference_bootstrap.R 100 0 200 2 baseline
#   Rscript simulations/panel_fe_beta_inference_bootstrap.R 50 0 500 4 baseline 1.1 1
#   Rscript simulations/panel_fe_beta_inference_bootstrap.R 50 0 500 4 baseline 1 1 hetero
#   Rscript simulations/k_sensitivity/panel_fe_beta_inference_k_sensitivity.R 50 0 500 4 sine_g 1 1 hetero_ar1 0.5 substantive 20
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
#   11. basis dimension K for the penalized smooth; default 20

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
k_smooth <- if (length(args) >= 11) as.integer(args[[11]]) else 20L

if (is.na(M) || M < 1L) stop("M must be a positive integer.")
if (is.na(B) || B < 0L) stop("B must be a nonnegative integer.")
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
if (is.na(k_smooth) || k_smooth < 3L) stop("smooth K must be at least 3.")

sigma_mu <- 1.25
sigma_e <- 0.5
methods <- c("bam_fe", "bam_re", "bam_fd")

out_dir <- file.path(getwd(), "outputs")
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
message("Penalized smooth basis dimension K: ", k_smooth)

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

cluster_vcov_penalty <- function(fit, cluster) {
  R <- predict(fit, type = "lpmatrix")
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

cluster_vcov_classic_beta <- function(fit, beta_name, cluster) {
  R <- predict(fit, type = "lpmatrix")
  u <- residuals(fit, type = "response")
  beta_idx <- match(beta_name, names(coef(fit)))

  x <- R[, beta_idx]
  Z <- R[, -beta_idx, drop = FALSE]
  x_tilde <- as.numeric(qr.resid(qr(Z), x))

  bread <- sum(x_tilde^2)
  score_g <- rowsum(x_tilde * as.numeric(u), group = as.factor(cluster), reorder = FALSE)
  meat <- sum(score_g^2)

  G <- nrow(score_g)
  nobs <- nrow(R)
  rank_r <- df_correction_rank(fit, R)
  corr <- cluster_small_sample_correction(G, nobs, rank_r)

  corr * meat / bread^2
}

beta_tstat <- function(fit, beta_name, cluster, type = c("penalty", "classic")) {
  type <- match.arg(type)
  beta_hat <- unname(coef(fit)[beta_name])
  if (type == "penalty") {
    V <- cluster_vcov_penalty(fit, cluster)
    idx <- match(beta_name, names(coef(fit)))
    var_beta <- V[idx, idx]
  } else {
    var_beta <- cluster_vcov_classic_beta(fit, beta_name, cluster)
  }
  se <- sqrt(max(var_beta, 0))
  t <- (beta_hat - beta_null) / se

  list(beta = beta_hat, se = se, t = t)
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

wild_cluster_test <- function(df, method) {
  dfit <- if (method == "bam_fd") make_first_differences(df) else df

  unrestricted <- fit_unrestricted(dfit, method)
  obs_penalty <- beta_tstat(
    unrestricted$fit,
    unrestricted$beta_name,
    unrestricted$cluster,
    type = "penalty"
  )
  obs_classic <- beta_tstat(
    unrestricted$fit,
    unrestricted$beta_name,
    unrestricted$cluster,
    type = "classic"
  )
  target_smooth <- unrestricted$fit$smooth[[1L]]
  smooth_cols <- target_smooth$first.para:target_smooth$last.para
  smooth_edf <- sum(unrestricted$fit$edf[smooth_cols])
  edf_ceiling <- length(smooth_cols)
  edf_ratio <- smooth_edf / edf_ceiling

  if (B == 0L) {
    return(data.frame(
      method = method,
      beta_hat = obs_penalty$beta,
      se_classic_cluster = obs_classic$se,
      t_classic_cluster = obs_classic$t,
      reject_classic_cluster = as.integer(abs(obs_classic$t) > qnorm(0.975)),
      se_penalty_cluster = obs_penalty$se,
      t_penalty_cluster = obs_penalty$t,
      reject_penalty_cluster = as.integer(abs(obs_penalty$t) > qnorm(0.975)),
      p_boot = NA_real_,
      reject_wild_boot = NA_integer_,
      boot_success = NA_integer_,
      K_smooth = k_smooth,
      smooth_edf = smooth_edf,
      edf_ceiling = edf_ceiling,
      edf_ratio = edf_ratio,
      binding_90 = as.integer(edf_ratio >= 0.90),
      binding_95 = as.integer(edf_ratio >= 0.95)
    ))
  }

  restricted <- fit_restricted(dfit, method)
  clusters <- sort(unique(dfit$unit_id))
  n_clusters <- length(clusters)
  t_boot <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    w <- sample(c(-1, 1), n_clusters, replace = TRUE)
    names(w) <- clusters
    y_star <- restricted$fitted_null + restricted$resid_null * w[as.character(dfit$unit_id)]

    if (method == "bam_fd") {
      dfit$dy_boot <- y_star
      boot_fit <- try(fit_unrestricted(dfit, method, dy_name = "dy_boot"), silent = TRUE)
    } else {
      dfit$y_boot <- y_star
      boot_fit <- try(fit_unrestricted(dfit, method, y_name = "y_boot"), silent = TRUE)
    }

    if (!inherits(boot_fit, "try-error")) {
      boot_t <- try(
        beta_tstat(boot_fit$fit, boot_fit$beta_name, boot_fit$cluster, type = "penalty")$t,
        silent = TRUE
      )
      if (!inherits(boot_t, "try-error") && is.finite(boot_t)) {
        t_boot[b] <- boot_t
      }
    }
  }

  p_boot <- mean(abs(t_boot) >= abs(obs_penalty$t), na.rm = TRUE)

  data.frame(
    method = method,
    beta_hat = obs_penalty$beta,
    se_classic_cluster = obs_classic$se,
    t_classic_cluster = obs_classic$t,
    reject_classic_cluster = as.integer(abs(obs_classic$t) > qnorm(0.975)),
    se_penalty_cluster = obs_penalty$se,
    t_penalty_cluster = obs_penalty$t,
    reject_penalty_cluster = as.integer(abs(obs_penalty$t) > qnorm(0.975)),
    p_boot = p_boot,
    reject_wild_boot = as.integer(p_boot < 0.05),
    boot_success = sum(is.finite(t_boot)),
    K_smooth = k_smooth,
    smooth_edf = smooth_edf,
    edf_ceiling = edf_ceiling,
    edf_ratio = edf_ratio,
    binding_90 = as.integer(edf_ratio >= 0.90),
    binding_95 = as.integer(edf_ratio >= 0.95)
  )
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
    ans <- try(wild_cluster_test(df, method), silent = TRUE)
    if (inherits(ans, "try-error")) {
      ans <- data.frame(
        method = method,
        beta_hat = NA_real_,
        se_classic_cluster = NA_real_,
        t_classic_cluster = NA_real_,
        reject_classic_cluster = NA_integer_,
        se_penalty_cluster = NA_real_,
        t_penalty_cluster = NA_real_,
        reject_penalty_cluster = NA_integer_,
        p_boot = NA_real_,
        reject_wild_boot = NA_integer_,
        boot_success = NA_integer_,
        K_smooth = k_smooth,
        smooth_edf = NA_real_,
        edf_ceiling = NA_real_,
        edf_ratio = NA_real_,
        binding_90 = NA_integer_,
        binding_95 = NA_integer_
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
results <- results[, c("sim", "method", setdiff(names(results), c("sim", "method")))]

summary_results <- do.call(
  rbind,
  lapply(split(results, results$method), function(d) {
    min_boot_success <- if (all(is.na(d$boot_success))) NA_real_ else min(d$boot_success, na.rm = TRUE)
    mean_boot_success <- if (all(is.na(d$boot_success))) NA_real_ else mean(d$boot_success, na.rm = TRUE)
    data.frame(
      method = d$method[1],
      n_success = sum(is.finite(d$beta_hat)),
      rejection_classic_cluster = mean(d$reject_classic_cluster, na.rm = TRUE),
      rejection_penalty_cluster = mean(d$reject_penalty_cluster, na.rm = TRUE),
      rejection_wild_boot = mean(d$reject_wild_boot, na.rm = TRUE),
      mean_beta_hat = mean(d$beta_hat, na.rm = TRUE),
      sd_beta_hat = sd(d$beta_hat, na.rm = TRUE),
      mean_se_classic_cluster = mean(d$se_classic_cluster, na.rm = TRUE),
      mean_se_penalty_cluster = mean(d$se_penalty_cluster, na.rm = TRUE),
      mean_p_boot = mean(d$p_boot, na.rm = TRUE),
      min_boot_success = min_boot_success,
      mean_boot_success = mean_boot_success,
      K_smooth = k_smooth,
      mean_edf = mean(d$smooth_edf, na.rm = TRUE),
      mean_edf_ceiling = mean(d$edf_ceiling, na.rm = TRUE),
      mean_edf_ratio = mean(d$edf_ratio, na.rm = TRUE),
      fraction_binding_90 = mean(d$binding_90, na.rm = TRUE),
      fraction_binding_95 = mean(d$binding_95, na.rm = TRUE)
    )
  })
)
row.names(summary_results) <- NULL
summary_results$error_design <- error_design
summary_results$rho <- if (error_design == "hetero_ar1") rho else NA_real_
summary_results$df_correction <- df_correction

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
  "_K", k_smooth,
  error_suffix,
  "_df", df_correction
)
results_path <- file.path(out_dir, paste0("panel_fe_beta_inference_results_", suffix, ".csv"))
summary_path <- file.path(out_dir, paste0("panel_fe_beta_inference_summary_", suffix, ".csv"))

write.csv(results, results_path, row.names = FALSE)
write.csv(summary_results, summary_path, row.names = FALSE)

message("\nSummary:")
print(summary_results, row.names = FALSE)

message("\nWrote:")
message("  ", results_path)
message("  ", summary_path)
