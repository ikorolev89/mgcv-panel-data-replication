# Inference diagnostics

This folder supplies the revised manuscript's appendix comparing five covariance
estimators. It uses separate scripts from the main simulation exercise.

The five estimators are:

1. `classic_cluster`: unit-cluster sandwich that ignores the smoothing penalty.
2. `penalty_cluster`: penalty-adjusted unit-cluster sandwich currently used in the paper.
3. `penalty_cluster_delta`: the preceding estimator plus
   `vcov(fit, freq = FALSE) - vcov(fit, freq = TRUE)`.
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

This writes `paper/inference_method_appendix.tex`, which is included by the
revised manuscript.

The completed-run interpretation is recorded in `RESULTS_SUMMARY.md`. The
machine-readable completion check is `outputs/combined/final_manifest_M1000.csv`.

The original runner manifest records an interrupted initial attempt; the final
completion manifest above describes the completed results. The public package
includes only the final manifest. See the public root README for the historical
chunk-seed limitation in the `n=500,T=4` function diagnostics. Full grid traces
are omitted from Git; compact replication records are in `replication/compact/`.
