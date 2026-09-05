# Inference diagnostics

This folder supplies the revised manuscript's appendix comparing five covariance
estimators. It uses separate scripts from the main simulation exercise.

The five estimators are:

1. `classic_cluster`: unit-cluster sandwich that ignores the smoothing penalty.
2. `penalty_cluster`: penalty-adjusted unit-cluster sandwich currently used in the paper.
3. `penalty_cluster_delta`: the preceding estimator plus
   `vcov(fit, sandwich = FALSE, freq = FALSE) -`
   `vcov(fit, sandwich = FALSE, freq = TRUE)`.
4. `mgcv_sandwich_freq`: `vcov(fit, sandwich = TRUE, freq = TRUE)`.
5. `mgcv_sandwich_default`: `vcov(fit, sandwich = TRUE, freq = FALSE)`.

The full exercise uses `K = 20`, `M = 1000`, the three sample configurations in
Tables 2--5, and two error processes: heteroskedastic serially independent errors
and heteroskedastic AR(1) errors with rho equal to 0.5.

Run all designs from the project root with:

```sh
Rscript simulations/inference_diagnostics/run_all_inference_diagnostics.R 1000 6
```

The runner is resumable. A job is skipped when its nonempty summary file already
exists. Results are written to `outputs/`, with one log per job in `logs/`.

After all jobs finish, consolidate the results with:

```sh
Rscript simulations/inference_diagnostics/combine_inference_diagnostics.R 1000
```

The consolidated CSV files and comparison plots are placed in
`outputs/combined/`.

Generate the four appendix tables from the consolidated output with:

```sh
Rscript simulations/inference_diagnostics/generate_appendix_tables.R
```

This writes `paper/inference_method_appendix.tex`, which supplies the numerical
entries in the revised manuscript's inference appendix tables.

The completed-run interpretation is recorded in `RESULTS_SUMMARY.md`. The
machine-readable completion check is `outputs/combined/final_manifest_M1000.csv`.

All function-coverage designs use the paired, replication-keyed scheme described
in `../../replication/PAIRED_G_RUNS.md`: fresh Gaussian draws per replication,
common draws across covariance methods, and separate data and integration RNG
streams. Main and K=20 sensitivity outputs share the matching replications.
The direct script uses the same stream keys as the parallel rerun driver.

Full grid traces are omitted from Git; compact replication records are in
`replication/compact/`. The paired manifest, code fingerprints, and output
checksums are in `replication/paired_g_*.csv`. The earlier seeded chunk schedules
and original interrupted-run manifests are historical and are not used by the
current tables. The general runner above resumes beta results but reruns
function coverage; use the paired driver for validated chunk-level resume.
