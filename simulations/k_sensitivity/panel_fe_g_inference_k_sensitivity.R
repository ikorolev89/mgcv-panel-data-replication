# Coverage simulation for confidence intervals/bands for g in panel GAMs.
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
# Inference for centered g on an evaluation grid:
#   - classic cluster sandwich, ignoring penalties
#   - penalty-adjusted cluster sandwich
#   - pointwise 95% intervals
#   - uniform 95% bands from a Gaussian sup-t approximation
#
# Usage:
#   Rscript simulations/panel_fe_g_inference_coverage.R 50 200 4 sin_2pi both 50 499
#   Rscript simulations/panel_fe_g_inference_coverage.R 50 200 4 sin_2pi both 50 499 hetero
#   Rscript simulations/k_sensitivity/panel_fe_g_inference_k_sensitivity.R 50 200 4 sin_2pi both 50 499 hetero_ar1 0.5 substantive 20
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
#   11. basis dimension K for the penalized smooth; default 20

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
k_smooth <- if (length(args) >= 11) as.integer(args[[11]]) else 20L

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
if (is.na(k_smooth) || k_smooth < 3L) stop("smooth K must be at least 3.")

model_list <- if (model_arg == "both") c("np", "pl") else model_arg
methods <- c("bam_fe", "bam_re", "bam_fd")
se_types <- c("classic", "penalty")

beta_true <- 1
sigma_mu <- 1.25
sigma_e <- 0.5
alpha <- 0.05
z_crit <- qnorm(1 - alpha / 2)

x_grid <- seq(0.05, 0.95, length.out = n_grid)
selected_points <- c(0.25, 0.50, 0.75)
selected_grid_idx <- vapply(
  selected_points,
  function(x0) which.min(abs(x_grid - x0)),
  integer(1)
)

out_dir <- file.path(getwd(), "outputs")
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
message("Penalized smooth basis dimension K: ", k_smooth)

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

safe_inverse <- function(A) {
  A <- (A + t(A)) / 2
  ans <- try(solve(A), silent = TRUE)
  if (!inherits(ans, "try-error")) return(ans)

  eig <- eigen(A, symmetric = TRUE)
  tol <- max(dim(A)) * max(abs(eig$values)) * .Machine$double.eps
  inv_values <- ifelse(abs(eig$values) > tol, 1 / eig$values, 0)
  eig$vectors %*% (inv_values * t(eig$vectors))
}

cluster_vcov_classic_smooth <- function(fit, cluster, smooth_cols) {
  R <- predict(fit, type = "lpmatrix")
  u <- residuals(fit, type = "response")

  S <- R[, smooth_cols, drop = FALSE]
  Z <- R[, -smooth_cols, drop = FALSE]
  S_tilde <- qr.resid(qr(Z), S)

  bread <- crossprod(S_tilde)
  score <- S_tilde * as.numeric(u)
  score_g <- rowsum(score, group = as.factor(cluster), reorder = FALSE)
  meat <- crossprod(score_g)

  G <- nrow(score_g)
  nobs <- nrow(R)
  rank_r <- df_correction_rank(fit, R)
  corr <- cluster_small_sample_correction(G, nobs, rank_r)

  bread_inv <- safe_inverse(bread)
  corr * bread_inv %*% meat %*% bread_inv
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
  z <- matrix(rnorm(draws * length(vals)), nrow = draws)
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

evaluate_fit <- function(fit_info, model_type, method, se_type) {
  grid <- make_grid_lpmatrix(fit_info, model_type, method)
  fit <- fit_info$fit
  target_smooth <- fit$smooth[[1L]]
  target_cols <- target_smooth$first.para:target_smooth$last.para
  smooth_edf <- sum(fit$edf[target_cols])
  edf_ceiling <- length(target_cols)
  edf_ratio <- smooth_edf / edf_ceiling

  estimate <- as.numeric(grid$L_full %*% coef(fit))

  if (se_type == "penalty") {
    V <- cluster_vcov_penalty(fit, fit_info$cluster)
    C <- grid$L_full %*% V %*% t(grid$L_full)
  } else {
    V_smooth <- cluster_vcov_classic_smooth(fit, fit_info$cluster, grid$smooth_cols)
    C <- grid$L_smooth %*% V_smooth %*% t(grid$L_smooth)
  }

  C <- (C + t(C)) / 2
  se <- sqrt(pmax(diag(C), 0))
  point_lower <- estimate - z_crit * se
  point_upper <- estimate + z_crit * se
  point_covered <- truth_centered >= point_lower & truth_centered <= point_upper

  sup_crit <- grid_cov_to_supcrit(C)
  uniform_lower <- estimate - sup_crit * se
  uniform_upper <- estimate + sup_crit * se
  uniform_covered <- all(truth_centered >= uniform_lower & truth_centered <= uniform_upper)

  data.frame(
    model = model_type,
    method = method,
    se_type = se_type,
    x = x_grid,
    truth = truth_centered,
    estimate = estimate,
    se = se,
    point_lower = point_lower,
    point_upper = point_upper,
    point_covered = as.integer(point_covered),
    point_width = point_upper - point_lower,
    sup_crit = sup_crit,
    uniform_lower = uniform_lower,
    uniform_upper = uniform_upper,
    uniform_covered = as.integer(uniform_covered),
    uniform_width = uniform_upper - uniform_lower,
    K_smooth = k_smooth,
    smooth_edf = smooth_edf,
    edf_ceiling = edf_ceiling,
    edf_ratio = edf_ratio,
    binding_90 = as.integer(edf_ratio >= 0.90),
    binding_95 = as.integer(edf_ratio >= 0.95)
  )
}

one_draw <- function(sim_id) {
  df <- simulate_panel()
  out <- list()
  idx <- 0L

  for (model_type in model_list) {
    for (method in methods) {
      fit_info <- fit_bam_method(df, model_type, method)
      for (se_type in se_types) {
        idx <- idx + 1L
        ans <- evaluate_fit(fit_info, model_type, method, se_type)
        ans$sim <- sim_id
        out[[idx]] <- ans
      }
    }
  }

  do.call(rbind, out)
}

set.seed(20260505)
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
results <- results[, c("sim", "model", "method", "se_type", setdiff(names(results), c("sim", "model", "method", "se_type")))]

coverage_by_x <- aggregate(
  point_covered ~ model + method + se_type + x,
  data = results,
  FUN = mean
)
names(coverage_by_x)[names(coverage_by_x) == "point_covered"] <- "point_coverage"
coverage_by_x$df_correction <- df_correction

selected_coverage <- do.call(
  rbind,
  lapply(seq_along(selected_points), function(j) {
    d <- coverage_by_x[abs(coverage_by_x$x - x_grid[selected_grid_idx[j]]) < 1e-12, ]
    d$selected_x <- selected_points[j]
    d$grid_x <- x_grid[selected_grid_idx[j]]
    d[, c("model", "method", "se_type", "selected_x", "grid_x", "point_coverage")]
  })
)

selected_wide <- reshape(
  selected_coverage[, c("model", "method", "se_type", "selected_x", "point_coverage")],
  idvar = c("model", "method", "se_type"),
  timevar = "selected_x",
  direction = "wide"
)

summary_base <- aggregate(
  cbind(point_covered, point_width, uniform_covered, uniform_width, sup_crit,
        smooth_edf, edf_ceiling, edf_ratio, binding_90, binding_95) ~ model + method + se_type,
  data = results,
  FUN = mean
)
names(summary_base)[names(summary_base) == "point_covered"] <- "avg_pointwise_coverage"
names(summary_base)[names(summary_base) == "point_width"] <- "avg_pointwise_width"
names(summary_base)[names(summary_base) == "uniform_covered"] <- "uniform_coverage"
names(summary_base)[names(summary_base) == "uniform_width"] <- "avg_uniform_width"
names(summary_base)[names(summary_base) == "sup_crit"] <- "avg_sup_crit"
names(summary_base)[names(summary_base) == "smooth_edf"] <- "mean_edf"
names(summary_base)[names(summary_base) == "edf_ceiling"] <- "mean_edf_ceiling"
names(summary_base)[names(summary_base) == "edf_ratio"] <- "mean_edf_ratio"
names(summary_base)[names(summary_base) == "binding_90"] <- "fraction_binding_90"
names(summary_base)[names(summary_base) == "binding_95"] <- "fraction_binding_95"

min_coverage <- aggregate(
  point_coverage ~ model + method + se_type,
  data = coverage_by_x,
  FUN = min
)
names(min_coverage)[names(min_coverage) == "point_coverage"] <- "min_pointwise_coverage"

summary_results <- merge(summary_base, min_coverage, by = c("model", "method", "se_type"))
summary_results <- merge(summary_results, selected_wide, by = c("model", "method", "se_type"), all.x = TRUE)
names(summary_results) <- sub("point_coverage\\.", "coverage_at_", names(summary_results))
summary_results <- summary_results[order(summary_results$model, summary_results$method, summary_results$se_type), ]
row.names(summary_results) <- NULL
summary_results$error_design <- error_design
summary_results$rho <- if (error_design == "hetero_ar1") rho else NA_real_
summary_results$df_correction <- df_correction
summary_results$K_smooth <- k_smooth

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
  "_K", k_smooth,
  error_suffix,
  "_df", df_correction
)

results_path <- file.path(out_dir, paste0("panel_fe_g_inference_grid_", file_suffix, ".csv"))
coverage_x_path <- file.path(out_dir, paste0("panel_fe_g_inference_coverage_by_x_", file_suffix, ".csv"))
summary_path <- file.path(out_dir, paste0("panel_fe_g_inference_summary_", file_suffix, ".csv"))
plot_path <- file.path(out_dir, paste0("panel_fe_g_inference_coverage_curve_", file_suffix, ".png"))

write.csv(results, results_path, row.names = FALSE)
write.csv(coverage_by_x, coverage_x_path, row.names = FALSE)
write.csv(summary_results, summary_path, row.names = FALSE)

coverage_plot <- ggplot(
  coverage_by_x,
  aes(x = x, y = point_coverage, color = method, linetype = se_type)
) +
  geom_hline(yintercept = 0.95, color = "grey45", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ model, ncol = 1) +
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
    color = "method",
    linetype = "SE"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title.position = "plot",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(plot_path, coverage_plot, width = 9, height = 6.5, dpi = 180, bg = "white")

message("\nSummary:")
print(summary_results, row.names = FALSE)

message("\nWrote:")
message("  ", results_path)
message("  ", coverage_x_path)
message("  ", summary_path)
message("  ", plot_path)
