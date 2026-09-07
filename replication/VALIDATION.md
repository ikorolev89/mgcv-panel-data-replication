# Export validation, 3 September 2026

Validation used R 4.6.0 and `mgcv` 1.9-4 in the environment recorded under
`environment/`.

## Paired smooth-coverage refresh

All smooth-coverage results were rerun using the `paired-replication-v1`
convention documented in `PAIRED_G_RUNS.md`. Each Monte Carlo replication
uses separate data-generation, Gaussian-integration, and fitting streams; it
draws a fresh 499-by-50 Gaussian matrix and shares that matrix across covariance
methods. The replication key does not depend on worker count, execution order,
basis dimension, or covariance-method order.

- The preliminary test matched main, K=20, and diagnostic results, checked
  repeated and reversed execution, changed the Gaussian-draw count without
  changing the DGP or point estimates, and reproduced results with one and two
  workers.
- The twelve unique designs completed 1,000 replications each: six K=20
  inference diagnostics (three n,T configurations and two error processes) and
  six K=10/K=40 AR(1) sensitivity designs. This represents 12,000
  case-replications and 72,000 fitted models.
- The 96 recorded chunks of 125 global replication IDs ran from
  2026-09-03 16:24:08 UTC through 17:27:47 UTC. Completed chunks survived the
  change from six to twelve workers. The accepted chunks contain no missing or
  failed replications and total 9.44 process-hours.
- The three main and three K=20 sensitivity outputs were derived from the
  coincident diagnostic replications. Thus the refresh produces eighteen
  reported coverage sets from twelve unique fitting designs.
- The previous matching result files were retained locally under
  `tmp/g_coverage_before_paired_20260903`.

## Numerical verification

The full traces were compacted into replication-level records and checked
against the point-specific and aggregate summaries.

- Main, K=20, and diagnostic Classic/Penalty records agree exactly for every
  shared replication, model, method, and reported quantity.
- The RNG manifest covers all 96 chunks; code fingerprints match the scripts
  used; checksums bind the 54 published compact, summary, and point-specific
  files to the run.
- The full verifier passes all 66 saved simulation runs with 1,000 replications
  per design-method cell: 12 estimation-accuracy, 36 beta-inference, and 18
  function-coverage sets.
- All 225 numerical cells in main Tables 1--5 match the refreshed revision
  snapshot. Tables 4--5, A3--A4, and A7--A8 and their EDF diagnostics were
  regenerated from the checked summaries.
- The combined diagnostic and K-sensitivity outputs and their plots were
  regenerated. The manuscript/response synchronization check passes for all
  inlined tables, bias ranges, and numerical coverage statements.

RMSE, beta inference, computation benchmarks, and the EmplUK application were
not affected by the Gaussian-draw change. Raw grid traces and resumable chunk
files remain local; the public package uses compact replication records plus
point-specific summaries.

## Natural-sample EmplUK refresh

The empirical application was rerun using the natural sample for each
transformation. The levels estimators use all 1,031 observations. Every firm's
observed years form a consecutive sequence, so first differencing yields 891
observations by removing only the first observed year for each of the 140 firms.

- The application script checks that all three levels models use 1,031
  observations and that the first-difference model uses 891 observations.
- The partially linear factor fixed-effects output EDF is 2.87, compared with
  0.96 for first differences. The factor fixed-effects estimates are -0.341 for
  log wages and 0.545 for log capital; the corresponding first-difference
  estimates are -0.430 and 0.397.
- Tables 6--7, Figures 1--2, their generated source files, and the manuscript
  summary were regenerated. Two complete application runs produced identical
  CSV and TeX outputs.

## Public-export and document checks

- A clean public export contained 238 allowlisted files (97.9 MB before
  regeneration). In that export, the full numerical verifier, the paired-RNG
  test, all eight simulation smoke cases, and the natural-sample empirical
  application regeneration passed.
- Regenerating the diagnostic and basis-sensitivity summaries and appendix
  fragments inside the export left the three checked table fragments byte-for-
  byte unchanged. No obsolete seeded-run or chunk file was exported.
- The clean and marked manuscripts compile to 60 pages and the response letter
  to 14 pages. The build logs contain no overfull boxes, unresolved references,
  or multiply defined labels. The affected main and appendix tables, the marked
  replacements, the revised empirical-application pages, and the revised
  response pages were rendered and visually checked for clipping, overlap, and
  page-break problems.


## Submission correction checks, 4 September 2026

The EmplUK application now enforces the known marginal-normal lower bound on
its simulated uniform-band cutoff. The FD cutoff is 1.9599639845 (reported as
1.960), so its nearly linear output function has coincident pointwise and
uniform bands. All 180 grid rows for both FE and FD satisfy interval nesting.
A rank-one covariance check returns the known normal cutoff, while a check with
independent components retains a larger simultaneous cutoff. The coefficient
estimates and factor-FE inference are unchanged.

The manuscript and response distinguish the strong K=10 binding in the sine
estimation-accuracy designs from the weaker binding in the inference designs.
The latter never cross the 95% EDF threshold. The public K-sensitivity findings
also report the corrected FD binding frequencies.

The source-workspace verifier and manuscript/response synchronization check
passed. A fresh 238-file export passed the full 66-run numerical verifier, all
eight simulation smoke cases, and the paired-RNG execution test. Regenerating
the application reproduced all eight CSV files and its LaTeX tables exactly;
regenerating all five simulation/application table fragments reproduced them
byte for byte. The corrected manuscript PDFs each contain 60 pages and the
response contains 14 pages, with no build warnings or unresolved references.

All saved simulation uniform cutoffs exceed 2.696, above the 1.960 marginal
bound. The application correction therefore requires no Monte Carlo rerun and
changes no saved simulation result or simulation-code fingerprint. This release
includes the paired smooth-coverage refresh documented above and removes the
superseded September 2 seeded-chunk files.

The affected code examples, binding discussion, empirical figures, and response
passages were rendered and visually checked in the clean and marked documents.
The marked and clean PDFs have identical extracted content, and the response's
page references remain accurate.

## Final manuscript publication checks, 7 September 2026

The public export now includes the finalized unmarked manuscript PDF, linked
from the root README and covered by the SHA-256 file manifest. The export
contains 239 allowlisted files plus the manifest (240 files, approximately
98.8 MB). Manuscript source, marked revisions, and referee correspondence remain
outside the public export.

The final LaTeX passes completed with no unresolved references or overfull
boxes. The clean and marked manuscripts each contain 60 pages and have
identical extracted text; the response contains 14 pages. The rebuilt PDFs
preserve the author's finalized text. The referee acknowledgment and affected
page breaks were rendered and visually checked, and the response's manuscript
page references remain correct. The numerical manuscript/response
synchronization check passed.

A fresh public export passed verification of all 66 saved simulation runs,
all 225 numerical cells in main Tables 1--5, paired RNG provenance, and the
54 saved-result checksums. Simulation code and results are unchanged from
the previous public release.

The final abstract contains 96 words when hyphenated compounds are counted as
one. The author's locally compiled clean and marked manuscript PDFs were
checked against the source and their title pages visually inspected. Only
page 1 changed; all later pages and reference page numbers remain unchanged.
