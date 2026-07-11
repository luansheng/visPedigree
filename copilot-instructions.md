# visPedigree Copilot Instructions

## Working commands

- First-time environment setup: `make bootstrap`
- Environment consistency check: `make status`
- Full test suite: `make test`
- Single test file or pattern: `make test TEST_FILTER=visped-layout`
  Replace `visped-layout` with another file stem such as
  `tidyped-alignment`, `pedmat-compact`, or `pedanalysis`.
- CRAN-like local check: `make check`
- Pkgdown site build: `make docs`

## High-level architecture

- [`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
  in `R/tidyped.R` is the package hub. Most exported workflows assume a
  `tidyped` object rather than raw pedigree data; the intended workflow
  is “tidy once, reuse many times”.
- `R/utils-s3.R` defines the `tidyped` contract: constructor/recovery
  helpers
  ([`new_tidyped()`](https://luansheng.github.io/visPedigree/reference/new_tidyped.md),
  [`as_tidyped()`](https://luansheng.github.io/visPedigree/reference/as_tidyped.md),
  [`ensure_tidyped()`](https://luansheng.github.io/visPedigree/reference/ensure_tidyped.md),
  [`ensure_complete_tidyped()`](https://luansheng.github.io/visPedigree/reference/ensure_complete_tidyped.md)),
  metadata access
  ([`pedmeta()`](https://luansheng.github.io/visPedigree/reference/pedmeta.md)),
  safe subsetting via `[.tidyped`, and print/summary/plot methods.
- [`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
  has two internal paths. Raw inputs go through validation, loop
  detection, candidate tracing, topological sorting, and generation
  assignment. Existing `tidyped` objects with `cand` use a faster
  integer-indexed C++ path.
- `R/pedanalysis.R` is the main analysis layer (`pedstats`, `pedne`,
  `pedecg`, `pedcontrib`, `pedancestry`, `pedpartial`, `pediv`, etc.).
  `R/inbreed.R` is intentionally thin and delegates the heavy work to
  the matrix/inbreeding backend.
- `R/pedmatrix.R` is the relationship-matrix engine
  ([`pedmat()`](https://luansheng.github.io/visPedigree/reference/pedmat.md),
  [`query_relationship()`](https://luansheng.github.io/visPedigree/reference/query_relationship.md),
  [`expand_pedmat()`](https://luansheng.github.io/visPedigree/reference/expand_pedmat.md),
  [`summary_pedmat()`](https://luansheng.github.io/visPedigree/reference/summary_pedmat.md))
  with compact full-sib support and attribute-rich S3 objects.
- Visualization is split across `R/visped.R` (public API) plus
  `R/visped_graph.R`, `R/visped_layout.R`, `R/visped_style.R`, and
  `R/visped_render.R`; `R/vismat.R` handles matrix heatmaps/histograms
  and `R/vispstat.R` handles `pedstats` plotting.
- `R/splitped.R` returns `splitped` lists of standalone `tidyped`
  groups; isolated `Gen == 0` records are tracked separately instead of
  being treated as regular groups.
- `src/pedigree_calculators.cpp` contains the performance-critical
  kernels (tracing, topological order, generation assignment,
  inbreeding, matrix construction, ancestry/contribution routines).
  `R/RcppExports.R` and `src/RcppExports.cpp` are generated files.
- When changing structural behavior around `tidyped`, consult
  `vignettes/tidyped-structure.Rmd` and
  `.github/prompts/vsPedigree.prompt.md` first.

## Key conventions

- New exported functions should be lowercase without underscores. The
  repo allows underscores for R-idiomatic `is_*`, `as_*`, and `has_*`
  helpers, but legacy names like
  [`query_relationship()`](https://luansheng.github.io/visPedigree/reference/query_relationship.md)
  are exceptions, not templates.
- Core pedigree columns use PascalCase and have stable semantics: `Ind`,
  `Sire`, `Dam`, `Sex`, `Gen`, `IndNum`, `SireNum`, `DamNum`. Many
  algorithms assume `IndNum` is exactly the current row order and that
  `SireNum`/`DamNum` are row pointers.
- `tidyped` metadata lives in `attr(x, "ped_meta")` with fixed fields
  `selfing`, `bisexual_parents`, and `genmethod`. Preserve it when
  restoring or transforming pedigree objects.
- `tidyped` is an S3 class on top of `data.table`. Prefer `data.table`
  operations, preserve class with
  [`new_tidyped()`](https://luansheng.github.io/visPedigree/reference/new_tidyped.md)
  /
  [`as_tidyped()`](https://luansheng.github.io/visPedigree/reference/as_tidyped.md)
  / `ensure_*()` helpers rather than manual `class<-`, and return
  `data.table` results with trailing `[]` so they print visibly after
  `:=`.
- Use the right guard for the job: completeness-sensitive recursive
  analysis and matrix code should go through
  [`ensure_complete_tidyped()`](https://luansheng.github.io/visPedigree/reference/ensure_complete_tidyped.md),
  while visualization and structure-light workflows usually use
  [`validate_tidyped()`](https://luansheng.github.io/visPedigree/reference/validate_tidyped.md)
  or
  [`ensure_tidyped()`](https://luansheng.github.io/visPedigree/reference/ensure_tidyped.md).
- Row subsetting can intentionally degrade a `tidyped` to plain
  `data.table` if parent rows are removed. To extract a valid analysis
  subset, use `tidyped(tp, cand = ids, trace = "up")` instead of naive
  row filtering.
- `pedmat` objects store important metadata in attributes
  (`attr(x, "ped")`, `"method"`, `"call_info"`, `"compact_map"`, etc.),
  not as list elements.
- Large-heatmap behavior in `R/vismat.R` is controlled by `VISMAT_*`
  constants; tests and docs expect those thresholds to stay
  synchronized.
- Prefer extending existing object-centric workflows over adding new
  exported entry points, and push heavy pedigree kernels toward the
  existing Rcpp/OpenMP layer rather than re-implementing them in slow R
  loops.
- Repo-local language preference: internal reasoning may be English,
  final user-facing confirmations should be in Chinese, and code
  comments should be in English.
