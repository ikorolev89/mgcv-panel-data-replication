# Shared numerical implementation and replication-keyed RNG for all g-coverage runs.
# Data, Gaussian integration, and fitting have disjoint L'Ecuyer streams.
# K and covariance-method lists deliberately do not enter the stream key.
g_rng_scheme <- "paired-replication-v1"
g_master_seed <- 20260903L
g_design_seed <- NA_integer_
g_streams <- list()
g_current_fit_state <- NULL
sup_z <- NULL

initialize_g_streams <- function(replications) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  })
  RNGkind("L'Ecuyer-CMRG")
  g_design_seed <<- as.integer((as.double(g_master_seed) +
    1000003 * n_units + 10007 * t_periods +
    1009 * match(error_design, c("iid", "hetero", "hetero_ar1"))) %%
    2147483646 + 1)
  set.seed(g_design_seed)
  state <- .Random.seed
  g_streams <<- lapply(seq_len(replications), function(i) {
    data_state <- state
    gaussian_state <- parallel::nextRNGStream(data_state)
    fit_state <- parallel::nextRNGStream(gaussian_state)
    state <<- parallel::nextRNGStream(fit_state)
    list(data = data_state, gaussian = gaussian_state, fit = fit_state)
  })
}

begin_g_replication <- function(sim_id) {
  streams <- g_streams[[sim_id]]
  assign(".Random.seed", streams$gaussian, envir = .GlobalEnv)
  sup_z <<- matrix(rnorm(n_sup_draws * n_grid), nrow = n_sup_draws)
  g_current_fit_state <<- streams$fit
  assign(".Random.seed", streams$data, envir = .GlobalEnv)
}

set_g_fit_stream <- function(model_type, method) {
  index <- (match(model_type, c("np", "pl")) - 1L) * 3L +
    match(method, c("bam_fe", "bam_re", "bam_fd"))
  state <- g_current_fit_state
  for (i in seq_len(index)) state <- parallel::nextRNGSubStream(state)
  assign(".Random.seed", state, envir = .GlobalEnv)
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
  se <- sqrt(pmax(diag(C), 0))
  valid <- is.finite(se) & se > 1e-10
  if (!any(valid)) return(0)
  correlation <- C[valid, valid, drop = FALSE] / outer(se[valid], se[valid])
  eig <- eigen(symmetrize(correlation), symmetric = TRUE)
  # The symmetric PSD root avoids arbitrary eigenvector-sign coupling of CRNs.
  root <- eig$vectors %*% (sqrt(pmax(eig$values, 0)) * t(eig$vectors))
  z <- sup_z[seq_len(draws), valid, drop = FALSE]
  maxima <- apply(abs(z %*% root), 1L, max)
  as.numeric(quantile(maxima, probs = 1 - alpha, names = FALSE))
}
