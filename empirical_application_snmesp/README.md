# Exploratory `Snmesp` application

This folder is deliberately separate from the manuscript. It contains a
first-pass comparison of a linear two-way fixed-effects model and an additive
semiparametric two-way fixed-effects model for the public `plm::Snmesp` panel.

The common specification uses log employment (`n`) as the outcome and log
wages (`w`), log output (`y`), log intermediate inputs (`i`), log capital (`k`),
and real cash flow (`f`) as contemporaneous covariates. Firm and year effects
are included in both models:

```text
Linear:         n_it = alpha_i + lambda_t + X_it' beta + u_it
Semiparametric: n_it = alpha_i + lambda_t + sum_j g_j(X_jit) + u_it
```

The exercise is exploratory. In particular, it does not resolve the dynamic
and endogeneity issues emphasized in the original applications of these data.
The curves should therefore be read as conditional relationships, not causal
effects or structural labor-demand functions.

## Run

From the repository root:

```r
Rscript empirical_application_snmesp/run_snmesp_application.R
```

The script requires `plm`, `mgcv`, `ggplot2`, `sandwich`, and `lmtest`. It
reloads the data from `plm` and refreshes the local CSV copy on every run.

## Generated files

- `data/Snmesp.csv`: local copy of the public package data used in the run.
- `outputs/linear_fe_coefficients.csv`: linear slopes with firm-clustered
  standard errors.
- `outputs/model_comparison.csv`: fit statistics for the two specifications.
- `outputs/smooth_diagnostics.csv`: effective degrees of freedom and descriptive
  curvature measures for each estimated smooth.
- `outputs/effect_curves.csv`: numerical values underlying the two curve plots.
- `outputs/concurvity.csv`: an additive-model analogue of collinearity.
- `outputs/model_summaries.txt`: complete model and data summaries.
- `figures/smooth_effects_vs_linear.png` and `.pdf`: estimated additive
  functions with pointwise model-based intervals, conditional on the selected
  smoothing parameters, and the linear-FE fits.
- `figures/nonlinear_deviations.png` and `.pdf`: each estimated smooth after
  removing its empirical best linear projection.
- `figures/model_fit_diagnostics.png` and `.pdf`: fitted-value and residual
  comparisons.

The plots restrict each horizontal axis to its 1st--99th percentile. The models
still use every observation. This matters especially for cash flow, which has
several extreme values.
