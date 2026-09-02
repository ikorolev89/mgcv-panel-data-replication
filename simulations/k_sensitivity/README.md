# Basis-dimension sensitivity simulations

This folder supplies the revised manuscript's basis-dimension appendix. Its
simulation programs rerun the designs reported in the paper at
penalized-smooth basis dimensions `K = 10, 20, 40`:

- out-of-sample accuracy for `g` (Table 1 designs);
- size and power for inference on `beta`;
- pointwise and uniform coverage for `g` in the NP and PL models.

Each design uses `M = 1000`, the paper's three `(n,T)` configurations, the
heteroskedastic AR(1) error design with `rho = 0.5`, and the original random
seeds. The unpenalized fixed-spline benchmarks remain at `K = 7`; only the
penalized `bam` smooths vary.

## EDF diagnostics

For the substantive smooth, `smooth_edf` is the sum of the coefficient-level
EDF entries in that smooth's parameter block. `edf_ceiling` is the number of
coefficients in the block: typically `K - 1` for a centered levels smooth and
`K` for the first-difference linear-functional smooth. `edf_ratio` is their
ratio. Two transparent operational definitions of a binding basis are saved:

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
and results go to `logs/`, `status/`, and `outputs/`. Completed jobs are skipped
when the runner is invoked again.

After the runs, execute from the repository root:

```sh
Rscript simulations/k_sensitivity/summarize_k_sensitivity.R
Rscript simulations/k_sensitivity/generate_appendix_tables.R
```

These commands also regenerate the combined results and appendix table from
the included summaries without rerunning the simulations. Full grid traces
are omitted from Git; `replication/compact/` retains replication-level coverage
records, and `outputs/` retains point-specific coverage summaries.
