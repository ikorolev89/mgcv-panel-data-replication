# Export validation, 2 September 2026

Validation used R 4.6.0 with the versions recorded in `environment/`.

## Initial export checks

- Parsed every exported R file.
- Ran `verify_replication.R` from a fresh public export: all 66 saved simulation
  runs passed, with 1,000 records per design-method cell and all 225 numerical
  entries in main Tables 1--5 matching the revision.
- Regenerated the main simulation tables, both appendix table files, the saved
  benchmark table, and the EmplUK table file. All five generated LaTeX files
  matched the supplied files byte for byte.
- Refit the manuscript EmplUK application directly from `plm::EmplUK`. All seven
  numerical output CSVs matched the supplied results within `1e-10`; both figures
  were regenerated successfully.
- Ran two-draw execution checks for all nine simulation scripts, covering main
  accuracy/beta/function inference, all three sensitivity scripts, and the beta,
  function, and chunked inference diagnostic scripts. All completed successfully
  in a temporary output directory.
- Audited the explicit export list and verified its SHA-256 file manifest.

## Seeded appendix refresh

The two `n=500,T=4` function-inference designs were rerun on 2 September 2026
using the eight recorded seeds in
`simulations/inference_diagnostics/seeded_g_schedule.csv`. Each error design has
four chunks of 250 draws, with 1,000 successful replications in every
model/estimator/covariance cell. The eight concurrent chunks completed in
23.6 minutes of wall time, excluding combination and document updates.

The recorded commands, seed schedule, mathematical-source fingerprints, and
checksums of all eight raw grids and chunk summaries were checked against the
completed outputs. These results replace the two historical runs whose chunk
seeds had not been retained.

- Refreshed the two compact coverage files and all affected point-specific and
  pooled inference summaries. The remaining 120 pooled function-summary cells
  are unchanged; only the 60 cells belonging to these two designs changed.
- Ran the full verifier in both the private workspace and a fresh public export:
  all 66 saved runs passed, including chunk-to-combined checks, and all 225
  numerical entries in main Tables 1--5 still match the revision.
- Regenerated the four pooled inference CSV files and the updated appendix table
  fragment from the fresh public export; all matched byte for byte.
- Synchronized Appendix Tables A3--A4 and the numerical bias/coverage passages
  in the working manuscript and response. All four document builds succeeded
  without undefined references or overfull boxes. The paper remains 52 pages
  and the response 12 pages; the changed pages were rendered and checked.

Across the 60 refreshed cells, the largest absolute change is 0.00886 in average
pointwise coverage and 0.027 in uniform coverage. The qualitative conclusions
are unchanged.

The full Monte Carlo suite was not rerun. Full grid traces are represented by
compact replication records plus point-specific summaries; the precise
verification scope is documented in the root README. Hardware timings are
retained from the original benchmark.
