# Reproducing the two chunked function-coverage designs

The `n=500,T=4` appendix designs with heteroskedastic independent errors and
heteroskedastic AR(1) errors each use 1,000 Monte Carlo replications. Each design
is divided into four chunks of 250 replications. Every draw estimates both NP
and PL models, all three estimators, and all five covariance methods on the
50-point evaluation grid with 499 Gaussian critical-value draws and `K=20`.

## Seeds and execution

`seeded_g_schedule.csv` is the authoritative schedule. Independent-error chunks
`h1`--`h4` use Monte Carlo seeds 20260902--20260905; AR(1) chunks `a1`--`a4` use
20261902--20261905. Each chunk separately initializes the shared Gaussian draw
matrix with seed 20260506, then initializes its Monte Carlo stream with its
recorded seed. The RNG settings are explicitly Mersenne-Twister, Inversion, and
Rejection.

From the repository root:

```sh
# Both designs; number of concurrent R processes is the second argument.
Rscript simulations/inference_diagnostics/run_seeded_g_diagnostics.R both 8

# A single design, with two workers.
Rscript simulations/inference_diagnostics/run_seeded_g_diagnostics.R hetero 2
Rscript simulations/inference_diagnostics/run_seeded_g_diagnostics.R hetero_ar1 2
```

An optional third argument selects a separate output root, for example
`tmp/seeded-check`. The driver always reruns the requested chunks, checks that
all 30 model/method/covariance cells contain 250 successful replications, then
combines them in chunk order. Combined simulation IDs 1--250, 251--500, 501--750,
and 751--1000 correspond to chunks 1--4. Worker count and finish order do not
affect that mapping. The standard `run_all_inference_diagnostics.R` uses this
same driver for these two full-size designs; it runs inner chunks serially
because its outer design jobs already run concurrently.

After rerunning these two designs in the repository root, refresh their compact
records, consolidate the appendix results, and regenerate the table fragment:

```sh
Rscript replication/compact_g_results.R simulations/inference_diagnostics/outputs/panel_fe_g_diagnostics_grid_M1000_n500_T4_*.csv.gz
Rscript simulations/inference_diagnostics/combine_inference_diagnostics.R 1000
Rscript simulations/inference_diagnostics/generate_appendix_tables.R
```

The first command uses shell expansion to select the two full grid traces.
If an alternate output root was used, copy its completed outputs into the
repository's inference output directory before running these commands.

## Provenance and checks

`outputs/seeded_g_manifest_hetero.csv` and
`outputs/seeded_g_manifest_hetero_ar1.csv` record the exact arguments, seeds, RNG
type, software versions, start and finish times, elapsed time, completion status,
MD5 fingerprints of the chunk script, combiner, and schedule, and MD5 checksums
of the generated chunk grid and summary files. The small chunk summaries are
included in the public package; the full compressed grid traces are reproducible
and remain outside Git.

The combiner rejects missing or duplicate grid records. The main replication
verifier checks the manifest against the schedule and code fingerprints,
validates the chunk-summary checksums, and recomputes every chunk's coverage,
width, and critical-value means from the corresponding block of compact
replication records. The usual checks then validate the combined summaries.

```sh
Rscript replication/verify_replication.R
```

These seeded runs supersede the original two chunked results whose individual
seeds had not been retained. The earlier results remain in private Git history
and a local comparison archive; they are not inputs to the current paper or
replication checks.
