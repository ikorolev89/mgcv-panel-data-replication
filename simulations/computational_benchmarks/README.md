# Computational benchmarks

This directory contains the reproducible timing and memory benchmark reported
in the manuscript. It uses the three `bam` estimators and the three panel sizes
from the simulation section, with the sine DGP, heteroskedastic AR(1) errors,
`K = 20`, and ten independently simulated panels per table entry.

Run from the project root:

```r
Rscript simulations/computational_benchmarks/run_benchmarks.R
```

The timing driver requires macOS. To regenerate only the table from the saved
results on any platform, use:

```sh
Rscript simulations/computational_benchmarks/run_benchmarks.R --tables-only
```

The driver launches every method/sample-size combination in a fresh R process
under `/usr/bin/time -l`. It writes numerical results to `outputs/`, diagnostic
logs to `logs/`, and the generated LaTeX table to
`paper/computation_benchmarks.tex`.

The reported fit time is elapsed (wall-clock) time for `mgcv::bam`. The
covariance time is elapsed time for the penalty-adjusted unit-cluster
covariance calculation. Peak RSS includes the R
session and loaded packages and is therefore a process-level memory measure,
not merely the size of the fitted model object.
