# Basis-dimension sensitivity simulations

This folder supplies the revised manuscript's basis-dimension appendix. Its
simulation programs rerun the designs reported in the paper at
penalized-smooth basis dimensions `K = 10, 20, 40`:

- out-of-sample accuracy for `g` (Table 1 designs);
- size and power for inference on `beta`;
- pointwise and uniform coverage for `g` in the NP and PL models.

Each design uses `M = 1000`, the paper's three `(n,T)` configurations, and the
heteroskedastic AR(1) error design with `rho = 0.5`. Accuracy and beta inference
retain their original seeds. Function coverage uses the paired scheme in
`../../replication/PAIRED_G_RUNS.md`, shared with the main and inference-method
results. The unpenalized fixed-spline benchmarks remain at `K = 7`; only the
penalized `bam` smooths vary.

## EDF diagnostics

For the substantive smooth, `smooth_edf` is the sum of the coefficient-level
EDF entries in that smooth's parameter block. `edf_ceiling` is the rank of the
corresponding fitted design block: typically `K - 1` for a centered levels
smooth and also `K - 1` for the first-difference linear-functional smooth,
because differencing annihilates the constant basis direction. `edf_ratio` is
their ratio. Two transparent operational definitions of a binding basis are
saved:

- `binding_90`: `edf_ratio >= 0.90`;
- `binding_95`: `edf_ratio >= 0.95`.

The summaries report the average EDF, average ceiling and ratio, and the
fractions binding under both thresholds.

## Running

From this directory:

```sh
Rscript run_all_k_sensitivity.R 9
```

The argument is the number of concurrent workers. Logs, completion markers,
and results go to `logs/`, `status/`, and `outputs/`. Completed non-coverage
jobs are skipped; coverage jobs are rerun because stale completion markers do
not establish the RNG convention. For validated chunk-level coverage resume,
use the paired driver from the repository root.

After the runs, execute from the repository root:

```sh
Rscript simulations/k_sensitivity/summarize_k_sensitivity.R
Rscript simulations/k_sensitivity/generate_appendix_tables.R
```

These commands also regenerate the combined results and appendix table from
the included summaries without rerunning the simulations. Full grid traces
are omitted from Git; `replication/compact/` retains replication-level coverage
records, and `outputs/` retains point-specific coverage summaries.
