# Results of the basis-dimension sensitivity exercise

## Scope

The paper's simulation designs were rerun at penalized-smooth basis dimensions
`K = 10, 20, 40`, using `M = 1000` and the three reported sample
configurations: `(n,T) = (200,4), (200,8), (500,4)`. The error design is the
reported heteroskedastic AR(1) process with `rho = 0.5`. The unpenalized spline
benchmarks in the RMSE design remain at `K = 7`.

All 36 jobs completed. All beta summaries have `n_success = 1000`; all nine
coverage files contain the expected 1,000 simulations and 600,000 grid-level
rows. The non-bam RMSE results are exactly equal across the three runs, which
confirms that the K comparisons use paired Monte Carlo draws.
The `K = 20` summaries reproduce the manuscript's reported entries to their
displayed precision.

## EDF and binding definitions

For the substantive smooth, `smooth_edf` is the sum of the coefficient-level
EDF entries in its parameter block. `edf_ceiling` is the number of coefficients
in that block: normally `K - 1` for a centered levels smooth and `K` for the
first-difference linear-functional smooth. The diagnostics use

- `binding_90 = 1{smooth_edf / edf_ceiling >= 0.90}`;
- `binding_95 = 1{smooth_edf / edf_ceiling >= 0.95}`.

Both fractions are reported because 90% is a useful early-warning threshold,
whereas 95% is a stricter definition of a nearly exhausted basis.

## Main findings

### Baseline K=20 versus K=40

The results are very stable when the baseline basis is doubled.

- Across all 27 bam RMSE cells, the largest relative change is 1.18%. For the
  sine DGP, the largest change is 0.73%; for the quadratic DGP it is 0.22%.
- Penalty-adjusted beta rejection frequencies differ by at most 0.2 percentage
  points. Mean beta estimates differ by at most `5.4e-6`.
- Average pointwise coverage differs by at most 0.068 percentage points,
  minimum pointwise coverage by at most 0.4 points, and uniform coverage by at
  most 0.7 points.
- The 90% and 95% binding fractions are zero in every RMSE, beta, and coverage
  cell at both `K = 20` and `K = 40`.

This supports retaining `K = 20` in the reported simulations and documenting
`K = 40` as a sensitivity check.

### K=10 versus K=20

The smaller basis is not uniformly innocuous, although it does not degrade the
reported performance measures.

- Linear and quadratic RMSE are effectively unchanged: the largest relative
  differences are 0.33% and 0.24%, respectively.
- In every sine-DGP RMSE cell, `K = 10` has lower RMSE than `K = 20`. The
  largest relative difference is 12.31%. Thus the smaller basis is acting as
  additional regularization; the result is not invariance to K.
- Beta rejection frequencies differ by at most 0.6 percentage points and mean
  beta estimates by at most `4.4e-5`.
- Average pointwise coverage differs by at most 0.458 percentage points,
  minimum pointwise coverage by at most 1.2 points, and uniform coverage by at
  most 1.3 points.

The EDF diagnostics flag why the RMSE comparison changes. In the sine RMSE
design, the levels FE/RE smooths at `K = 10` cross the 90% threshold in nearly
every draw. They cross the stricter 95% threshold in about 0.6% of draws for
`(200,4)`, and in roughly 98--100% of draws for `(200,8)` and `(500,4)`. The FD
smooth does not bind under either rule. In the beta and coverage designs, the
levels smooths at `K = 10` cross the 90% rule in about 26--29% of draws for
`(200,8)` and about 99% for `(500,4)`, but never cross the 95% rule.

Consequently, a binding diagnostic should be read as evidence that enlarging
the basis can change the fitted complexity, not as evidence that the smaller
basis must have worse finite-sample RMSE.

### EDF can increase with sample size at fixed K

The simulations directly illustrate this point. For the sine RMSE design, the
average levels-FE EDF at fixed `K = 10` rises from 8.44 for `(200,4)` to 8.67
for `(200,8)` and 8.74 for `(500,4)`, against a ceiling of 9. At fixed `K = 20`,
it rises from 11.49 to 12.74 and 13.26, against a ceiling of 19. At fixed
`K = 40`, it rises from 12.05 to 13.64 and 14.34, against a ceiling of 39.

Thus a fixed basis dimension does not imply fixed estimated EDF. It does,
however, impose a fixed upper bound; the binding diagnostics show whether that
bound is empirically active.

## Practical reading

The results support a recommendation to start with a reasonably generous
basis, check EDF relative to its ceiling, and refit with a larger basis. In
these designs, `K = 20` passes that check: doubling it to 40 changes neither
fit accuracy nor inference materially, and its ceiling never binds. The
`K = 10` results also show why a symmetric "halve and double K and expect the
same answer" statement would be too strong: a smaller basis may add useful
regularization and alter finite-sample RMSE even when inference changes little.

## Output map

- `combined/mse_summary_all_K.csv`: all RMSE summaries and EDF diagnostics.
- `combined/beta_summary_all_K.csv`: all beta size/power summaries.
- `combined/coverage_summary_all_K.csv`: all confidence-band summaries.
- `combined/edf_binding_summary.csv`: compact EDF/binding results across blocks.
- `combined/*_stability_ranges.csv`: within-cell ranges across all three K's.
- `combined/mse_by_K.png`, `beta_rejection_by_K.png`, and `coverage_by_K.png`:
  graphical comparisons.
- `outputs/`: raw and per-design results (about 1.6 GB).
- `logs/`, `status/`, and `run_manifest.csv`: reproducibility and completion
  records.

No manuscript file was changed as part of this exercise.
