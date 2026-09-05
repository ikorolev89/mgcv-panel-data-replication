# Inference diagnostic results

The exercise uses `K = 20`, `M = 1000`, the three sample configurations from
Tables 2--5, and two error designs: heteroskedastic serially independent errors
and heteroskedastic AR(1) errors with rho equal to 0.5. Both the nonparametric
and partially linear designs are included for inference on centered `g(x)`.

All requested cells completed with 1,000 successful replications. The final
pooled outputs contain 180 beta-summary rows, 180 function-summary rows, and
9,000 point-specific diagnostic rows. The function designs use the paired,
replication-keyed reruns documented in `../../replication/PAIRED_G_RUNS.md`.

## Beta inference

The table reports empirical size averaged across the six sample/error cells,
followed by the range across those cells.

| Method | Covariance | Mean size | Range |
|---|---|---:|---:|
| FE | Classic cluster | 0.054 | 0.050--0.059 |
| FE | Penalty cluster | 0.056 | 0.050--0.062 |
| FE | Penalty cluster + Delta | 0.056 | 0.050--0.062 |
| FE | mgcv sandwich, frequentist | 0.056 | 0.048--0.068 |
| FE | mgcv sandwich, default | 0.056 | 0.048--0.068 |
| FD | Classic cluster | 0.057 | 0.053--0.061 |
| FD | Penalty cluster | 0.059 | 0.055--0.063 |
| FD | Penalty cluster + Delta | 0.059 | 0.054--0.063 |
| FD | mgcv sandwich, frequentist | 0.090 | 0.080--0.102 |
| FD | mgcv sandwich, default | 0.090 | 0.080--0.102 |

The unit-cluster estimators have size close to 0.05. The observation-level
`mgcv` sandwiches over-reject for FD because they do not account for clustering.
Adding Delta has almost no effect on beta standard errors, size, or power.

## Inference on g(x)

The table averages over both function designs, both error designs, and all three
sample configurations (12 cells per method).

| Method | Covariance | Pointwise coverage | Uniform coverage | Mean SE / MC SD |
|---|---|---:|---:|---:|
| FE | Classic cluster | 0.997 | 1.000 | 1.574 |
| FE | Penalty cluster | 0.944 | 0.929 | 0.981 |
| FE | Penalty cluster + Delta | 0.970 | 0.983 | 1.116 |
| FE | mgcv sandwich, frequentist | 0.946 | 0.937 | 0.984 |
| FE | mgcv sandwich, default | 0.971 | 0.986 | 1.119 |
| FD | Classic cluster | 0.995 | 0.999 | 1.534 |
| FD | Penalty cluster | 0.943 | 0.924 | 0.978 |
| FD | Penalty cluster + Delta | 0.965 | 0.975 | 1.084 |
| FD | mgcv sandwich, frequentist | 0.912 | 0.845 | 0.873 |
| FD | mgcv sandwich, default | 0.947 | 0.945 | 0.990 |

The main patterns are:

- The current penalty-adjusted cluster covariance yields reasonably calibrated
  pointwise FE and FD inference while retaining unit clustering. Its average
  coverage is about 0.943--0.944, and its reported standard errors are close to
  the Monte Carlo standard deviation.
- Adding Delta raises average pointwise coverage to approximately 0.965--0.970
  and uniform coverage to approximately 0.975--0.983. It is therefore a useful
  conservative sensitivity check, but it does not improve pointwise calibration
  relative to the 0.95 target in these designs.
- The classic covariance that ignores penalization is much too conservative:
  its standard errors are roughly 53--57 percent larger than the empirical
  standard deviation and its coverage is essentially one.
- The observation-level frequentist `mgcv` sandwich undercovers for FD. The
  default `mgcv` sandwich adds the same Delta term and performs better, but it
  still does not provide cluster-robust inference.
- RE function estimates display appreciable bias when `n = 500` because the DGP
  correlates unit heterogeneity with the regressors. Their mean absolute bias is
  about 0.015, compared with about 0.002--0.003 for FE and FD. The resulting RE coverage
  should not be interpreted as a pure covariance-estimator comparison.

With 1,000 replications, the Monte Carlo standard error of an individual 0.95
coverage rate is approximately 0.007. This is a rough scale for individual
cells, not the standard error of the averages across grid points and designs.
The persistent 0.965--0.970 pointwise coverage after adding Delta is conservative.

## Practical implication

The simulations support retaining the penalty-adjusted unit-cluster covariance
as the baseline for FE and FD pointwise inference. Reporting the Delta-augmented
version as a diagnostic or robustness check is informative, especially for
simultaneous bands, but the results do not support replacing the baseline with
it automatically. The unpenalized covariance is too conservative, and the
observation-level `mgcv` sandwiches are not substitutes for clustering.
