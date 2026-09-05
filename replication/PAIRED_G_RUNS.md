# Paired smooth-coverage simulations (September 3, 2026)

This run supersedes all earlier smooth-coverage results: main Tables 4--5,
Appendix A3--A4, and the smooth-coverage/EDF entries in A7--A8. Estimation-RMSE
and beta-test simulations are unchanged. The twelve unique cases are the six
K=20 inference diagnostics (three n,T pairs, two error processes), plus K=10
and K=40 for the three AR(1) configurations. Each case contains 1,000
replications, both model classes, and all three estimators. The three main runs
and three K=20 sensitivity runs are derived from the diagnostic replications,
not fitted independently. This yields eighteen sets of saved coverage results.

## Random-number convention

All three entry-point scripts source `simulations/g_coverage_shared.R`.
The scheme is `paired-replication-v1`, master seed 20260903, with
L'Ecuyer-CMRG streams. The recorded design seed is

```
(20260903 + 1000003*n + 10007*T + 1009*error_index) %% 2147483646 + 1
```

where error_index is 1 for iid, 2 for hetero, and 3 for hetero_ar1.
Successive replications reserve three disjoint streams: data generation,
Gaussian integration, and fitting. Fitting has model/estimator-specific
substreams. Basis dimension and covariance-method lists do not enter the key.
Every replication draws a fresh 499-by-50 standard-normal matrix, shared by
the covariance methods and estimators within that replication. A symmetric
positive-semidefinite square root of each estimated correlation matrix maps
those common draws into the corresponding Gaussian process. This avoids
arbitrary eigenvector-sign differences in the common-random-number coupling.

The DGP is independent of the number of Gaussian draws. Replication IDs are
independent of chunk boundaries, worker count, covariance-method order, and
completion order. Common random numbers improve pairing; they do not remove
finite-Gaussian-draw or Monte Carlo approximation error. Earlier use of a single
Gaussian matrix across all diagnostic replications has been discontinued.

## Run and refresh

Run from the repository root on Unix/macOS (six workers is the default):

```sh
Rscript replication/run_paired_g_simulations.R --test
Rscript replication/run_paired_g_simulations.R --run 6
Rscript replication/run_paired_g_simulations.R --summarize
Rscript replication/run_paired_g_simulations.R --promote
Rscript replication/compact_g_results.R
Rscript replication/verify_paired_g.R --record-output-hashes
Rscript simulations/inference_diagnostics/combine_inference_diagnostics.R 1000
Rscript simulations/inference_diagnostics/generate_appendix_tables.R
Rscript simulations/k_sensitivity/summarize_k_sensitivity.R
Rscript simulations/k_sensitivity/generate_appendix_tables.R
Rscript replication/generate_main_tables.R --refresh-coverage-snapshot
Rscript replication/verify_replication.R
```

The preliminary test checks exact overlap up to numerical precision, fresh
integration draws, reproducibility, method-order invariance, and independence
of the DGP and point estimates from the integration-draw count. Production
work is split into 96 resumable chunks of 125 global replication IDs. Cached
chunks are accepted only if their arguments, software versions, and code
fingerprints match. A failed replication stops its chunk rather than silently
reducing the denominator. With the same numerical environment, running an
individual entry-point script produces the same corresponding results.

The driver stages results in `tmp/paired_g_20260903`; promotion preserves prior
matching files in `tmp/g_coverage_before_paired_20260903`. To force a fresh run,
move the staging directory aside. Raw grid traces and chunk RDS files are
local/reproducible and are not included in the public package.

## Provenance and verification

`paired_g_manifest.csv` records the 96 completed chunks, design/replication
keys, integration size, software, start/end times, and chunk checksums.
`paired_g_code_md5.csv` records the scripts actually used. The checksum file
`paired_g_output_md5.csv` binds the 54 published compact, summary, and
point-specific files to this run. It is recorded only after the replication
overlap checks pass. Code fingerprints also prevent reuse of stale chunks.

The fast verifier checks every replication count and summary, all provenance
files, and exact equality of the Classic/Penalty replication-level quantities
across main, diagnostic, and K=20 results. The compact records retain EDF and
binding diagnostics as well as coverage, widths, and critical values.

Earlier September 2 chunk schedules and outputs are retained locally for
historical comparison only. They are not inputs to the current tables or part
of the refreshed public export. Refreshing the local package does not publish
it; updating the public repository is a separate action.
