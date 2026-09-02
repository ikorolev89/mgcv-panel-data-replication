# Export validation, 2 September 2026

Validation used R 4.6.0 with the versions recorded in `environment/`.

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

These checks do not rerun the full Monte Carlo suite. Full grid traces are
represented by compact replication records plus point-specific summaries; the
precise verification scope and historical chunk-seed limitation are documented
in the root README. Hardware timings are retained from the original benchmark.
