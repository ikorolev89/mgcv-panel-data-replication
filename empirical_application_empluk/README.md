# `EmplUK` manuscript application and exploratory comparison

`run_manuscript_application.R` generates the revised paper's empirical tables
and figures from the public `plm::EmplUK` panel. It includes linear factor fixed
effects, partially linear factor fixed effects and first differences, and an
additive factor-FE model. `run_empluk_application.R` is an additional exploratory
comparison of linear and additive semiparametric two-way fixed-effects models.

The specifications are

```text
Linear:
  log(emp_it) = alpha_i + lambda_t
                + beta_w log(wage_it)
                + beta_k log(capital_it)
                + beta_y log(output_it) + u_it

Semiparametric:
  log(emp_it) = alpha_i + lambda_t
                + g_w(log(wage_it))
                + g_k(log(capital_it))
                + g_y(log(output_it)) + u_it.
```

The manuscript application also estimates a partially linear specification
that is linear in log wages and log capital and smooth in log output.

The original unbalanced panel is retained. Sector is omitted because it is
time-invariant and therefore absorbed by the firm effects.

## Run

From the repository root:

```r
Rscript empirical_application_empluk/run_empluk_application.R
```

To regenerate the manuscript tables and figures:

```r
Rscript empirical_application_empluk/run_manuscript_application.R
```

The exercise is descriptive. The original Arellano--Bond application treats
employment dynamically and is concerned with predetermined or endogenous
regressors. The curves here should therefore be interpreted as conditional
relationships rather than causal or structural effects.

## Generated files

- `data/EmplUK.csv`: local copy of the package data used in the run.
- `outputs/linear_fe_coefficients.csv`: linear coefficients with standard
  errors clustered by firm.
- `outputs/model_comparison.csv`: descriptive fit statistics.
- `outputs/smooth_diagnostics.csv`: effective degrees of freedom and curvature
  measures.
- `outputs/basis_sensitivity.csv`: fit and smooth complexity for basis
  dimensions 5, 6, 8, 10, and 12.
- `outputs/effect_curves.csv`: numerical values underlying the effect plots.
- `outputs/concurvity.csv`: additive-model concurvity diagnostics.
- `outputs/model_summaries.txt`: full data and model summaries.
- `figures/smooth_effects_vs_linear.png` and `.pdf`: additive functions and
  their linear-FE counterparts.
- `figures/nonlinear_deviations.png` and `.pdf`: estimated smooths after
  subtracting their empirical best linear projections.
- `figures/model_fit_diagnostics.png` and `.pdf`: residual comparisons.
- `manuscript_outputs/`: numerical results and summaries for the paper.
- `manuscript_figures/`: publication-ready comparison and inference figures.
- `../paper/empluk_application_results.tex`: generated LaTeX tables.

The exploratory plots show the 1st--99th percentile of each logged covariate
with pointwise model-based intervals. The manuscript script constructs
penalty-adjusted firm-clustered pointwise and grid-uniform confidence bands,
conditional on the selected smoothing parameters. Models use all observations.
