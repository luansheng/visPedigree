# Changes in version 1.9.0 released on 11 Jul 2026
## New features
1. **`visped()` SVG output**: `visped()` now writes SVG files when `file` ends in `.svg`; other file names continue to use PDF output by default.
2. **`visped()` custom labels**: Added a `labelvar` argument to display labels from either a user-selected pedigree column or a row-aligned character vector. Compact full-sib family nodes continue to show family size.
3. **`visped()` symbol schemes**: Added a `shapeby` argument for choosing how node shapes are encoded. The default `shapeby = "sex"` uses circles for females, squares for males, diamonds for unknown sex, and hexagons for monoecious individuals. Set `shapeby = "role"` to keep the legacy role-based symbol scheme in which real individuals are circles and compact full-sib family summaries are rectangles.
4. **Matrix-free relationship products**: Added `pedprod()` to compute $Ax$, $AX$, $A^{-1}x$, and $A^{-1}X$ directly from a complete pedigree. The implementation uses the $A = TDT'$ factorization and avoids materializing the dense additive relationship matrix.
5. **`visped()` custom generation labels**: `genlab` now also accepts an unnamed character vector, assigning one user-provided label to each displayed generation from top to bottom. The existing logical values retain their behavior: `TRUE` draws the default `G1`, `G2`, ... labels and `FALSE` omits them.

## Improvements
1. **Matrix-free `pedrel()` summaries**: `pedrel()` now computes grouped mean relationships and coancestries from batched $Ag$ products. It traces the union of selected ancestors once and no longer constructs or scans a dense relationship matrix for each group.
2. **Large-group support in `pedrel()`**: The former dense and compact matrix size guards are no longer needed. The development-only `force`, `max_dense`, and `max_compact` no-op arguments were removed before release, while `compact` remains as a backward-compatible no-op argument. `Status` and `Message` continue to identify groups skipped because fewer than two individuals were selected or groups that failed during tracing.
3. **Matrix-free coancestry analyses**: `pedne(method = "coancestry")` now obtains sampled pair relationships from batched $AX$ products instead of constructing a dense triangular relationship matrix. `pediv()` reuses the same products for $N_e$ and founder genome equivalents, while `pedhalflife()` requests only the $f_g$ statistic it needs.
4. **Matrix-free grouped heatmaps**: `vismat(tidyped_object, by = ...)` now computes grouped additive relationships as $W' A W$ without first constructing the full individual-level relationship matrix.
5. **Clearer compact full-sib family summaries in `visped()`**: Compact full-sib family labels now use the explicit `FS×N` form instead of a bare number, avoiding confusion with an individual ID. Compact family summaries use a green-grey fill with a darker frame, while unknown-sex individuals use a neutral-grey fill with a grey frame.

## Bug fixes
1. **Selfing in `Ainv`**: Corrected the sire-dam diagonal cross-term in Henderson's inverse construction when the sire and dam are the same individual. `pedmat(method = "Ainv")` now remains the numerical inverse of `A` for pedigrees created with `selfing = TRUE`.

## Documentation
1. **Plant pedigree example in `tidy-pedigree` vignette**: Added section 3.9 demonstrating `selfing = TRUE` for monoecious species, including sex inference, inbreeding coefficients under self-fertilization, and the `summary()` output for plant pedigrees.
2. **Plant pedigree visualization in `draw-pedigree` vignette**: Added section 1.1.3 showing a multi-generation self-pollinating pedigree with `visped()`. Monoecious individuals are drawn as hexagons with teal edges, selfing edges are shown in teal, and inbreeding coefficients are displayed with `showf = TRUE`.

# Changes in version 1.8.1 released on 28 Mar 2026
## New features
1. **`visped()` generation labels**: Added a `genlab` argument to `visped()` for drawing generation labels (`G1`, `G2`, ...) on the left margin of pedigree plots. The default is `FALSE`, so existing plots are unchanged unless `genlab = TRUE` is requested.

## Documentation
1. **`draw-pedigree` vignette**: Added an example showing how to display generation labels with `visped(..., genlab = TRUE)`.

# Changes in version 1.8.0 released on 25 Mar 2026
## Performance
1. **Inbreeding calculation**: Replaced the Meuwissen and Luo (1992) linear path-trace algorithm in `cpp_calculate_inbreeding()` with the Sargolzaei and Iwaisaki (2005) LAP (Longest Ancestral Path) bucket method. At N = 1,000,000, the C++ kernel completes in about 0.15 s (previously about 15 s), and the full `inbreed()` call returns in about 0.40 s. The implementation uses O(1) ancestor retrieval via bucket pop, O(1) duplicate suppression via the `L[k] == 0` check, and O(m_i) path-coefficient reset. Results are numerically identical to the previous implementation (maximum difference < 2e-15).
2. **`tidyped()` candidate tracing**: When the input is already a `tidyped` object and `cand` is specified, the fast path now uses three C++ BFS functions (`cpp_trace_ancestors()`, `cpp_trace_descendants()`, and `cpp_topo_order()`) instead of rebuilding an `igraph` object. At N = 1,000,000 with 200 candidates, elapsed time drops from about 0.91 s to about 0.056 s. The fast path now also applies sibling and mate generation alignment for `genmethod = "bottom"`, matching the output of the full path.
3. **`pedgenint()` generation-interval lookup**: Replaced the character-key `data.table` join used to look up each individual's generation number with a pre-computed integer index array. Benchmark on a representative dataset: 4.24 s to 0.49 s.
4. **`pedmat(sparse = TRUE)` matrix coercion**: The conversion from a dense numeric matrix to a `Matrix::dgeMatrix` now bypasses the S4 dispatch overhead of `as(mat, "dgeMatrix")` via a direct `methods::new()` call backed by a package-level class cache. Benchmark on a representative subpedigree: 1.24 s to 0.60 s.
5. **`pedrel(compact = TRUE)` safety guard**: Added an early error when the number of reference individuals exceeds 200,000 in compact mode, preventing unintended O(N^2) matrix allocation.

## Bug fixes
1. **`tidyped()` fast path with `addnum = FALSE` and `cand`**: When the input `tidyped` object was created with `addnum = FALSE`, passing `cand` could raise `"None of the specified candidates were found in the pedigree."` because the fast-path BFS looked up `ped_dt$IndNum`, which was `NULL`. The fix temporarily adds integer index columns for the BFS and removes them from the output when `addnum = FALSE`.
2. **`pediv()` and `pedne()` coancestry outputs when `ECG` was not pre-computed**: Internal `merge(..., by = "Ind", all.x = TRUE)` calls sorted output alphabetically by `Ind`, breaking the `IndNum == row-index` invariant used by the fast-path BFS in `tidyped()`. This could produce incorrect `fg`, `MeanCoan`, and `NeCoancestry` values, such as `fg` near 240 instead of about 19. Fixed by restoring `IndNum` order immediately after the merge.
3. **`summary_pedmat()` density for dense `Matrix` subclasses**: `summary_pedmat()` reported `Density = 100%` for all `A`, `D`, and `AA` matrices returned as `Matrix::dgeMatrix` objects. It now uses `Matrix::nnzero() / (nrow * ncol)` for all `Matrix` subclasses, giving the correct fill ratio for both `dgeMatrix` and `dgCMatrix`.

## Documentation
1. **Inbreeding references**: Updated `inbreed()` and vignette references to cite Sargolzaei and Iwaisaki (2005) instead of Meuwissen and Luo (1992).

## Internal changes
1. **`align_bottom_generations()` helper**: Consolidated the sibling and mate generation-alignment block previously duplicated between the main and fast paths of `tidyped()` into a single internal helper.
2. **C++11 compatibility**: Replaced two C++17 structured-binding usages in the BFS functions with explicit C++11 equivalents for compatibility with GCC 8.
3. **`methods` dependency**: Listed the `methods` package under `Imports` in `DESCRIPTION`, as required when `methods::getClass()` and `methods::new()` are called at runtime.

# Changes in version 1.7.0 released on 23 Mar 2026
## New features
1. **`pediv()` retained genetic diversity (`GeneDiv`)**: `pediv()$summary` gains a `GeneDiv = 1 - MeanCoan` column, the pedigree-based retained genetic diversity of the reference population. Values lie in $[0, 1]$, with larger values indicating more diversity retained relative to an unrelated base population. `print.pediv()` displays it alongside `fg` and `MeanCoan`.
2. **`vismat()` large-pedigree representative view**: When the original pedigree has more than 5,000 individuals, `vismat()` no longer attempts a full N × N matrix expansion. It uses the compact K × K representative-individual matrix directly and adds sibling-count labels of the form `ID (×n)` to each axis tick. When `compact = TRUE` and `by` is supplied, group means are computed algebraically from the K × K matrix without expanding to N × N.

## Internal changes
1. **`vismat()` threshold constants**: Hardcoded values controlling large-pedigree behavior (`5000`, `2000`, and `50`) were refactored into named constants (`VISMAT_EXPAND_MAX`, `VISMAT_LABEL_MAX`, and `VISMAT_WARN_THRESHOLD`) at the top of `R/vismat.R`.

# Changes in version 1.6.2 released on 23 Mar 2026
## New features
1. **`pedrel()` coancestry scale**: Added a `scale` parameter to `pedrel()` supporting `"relationship"` (default, returns mean $a_{ij}$) and `"coancestry"` (returns corrected mean coancestry $\bar{c}$). The coancestry scale uses the diagonal-corrected formula of Caballero and Toro (2000), accounting for self-coancestry within the reference group.

## API changes
1. **`vispstat()` internal-only backend**: `vispstat()` is now the internal backend for `plot.pedstats()`. Users should call the standard S3 method with `plot(stats_obj)`.

## Bug fixes
1. **Spurious subsetting warnings**: Internal group-by slicing in `pedrel()`, `pedne()`, `pediv()`, `pedhalflife()`, and `pedgenint()` no longer triggers false-positive `[.tidyped]` warnings when the subset is only used for ID extraction.
2. **Internal class-restoration messages**: `pedrel(compact = TRUE)` no longer emits class-restoration messages caused by early-return branches in `compact_ped_for_matrix()`. Those branches now preserve the `tidyped` class by returning `data.table::copy(ped)`.
3. **`vispstat()` cleanup**: Removed dead-code variables and added unit tests for the `genint` branch of `vispstat()`.

## Documentation
1. **`vispstat()` generation-interval documentation**: Updated the text to describe mean values accurately and remove the misleading "mean +/- SD" claim.
2. **`pedigree-analysis.Rmd` section 9**: Split the "Average Relationship Trends with `pedrel()`" section into two subsections covering both `scale` options. Added the Caballero and Toro (2000) diagonal-corrected coancestry formula, a worked `scale = "coancestry"` example, and guidance on choosing a scale.
3. **`relationship-matrix.Rmd` updates**: Added a compact-to-`vismat()` direct path, expanded the visualization examples, and replaced incorrect performance thresholds with a reference table.

# Changes in version 1.6.1 released on 21 Mar 2026
## Improvements
1. **Diversity notation**: Output columns in `pedhalflife()$timeseries` are now lowercase (`fe`, `fa`, `fg`, `lnfe`, `lnfa`, `lnfg`, `lnfafe`, and `lnfgfa`) to match population genetics notation.
2. **`plot.pedhalflife()` log view**: In `type = "log"` mode, the plot includes an OLS regression line for total diversity decay ($\ln f_g \sim \text{Time}$) and a vertical reference line for the diversity half-life $T_{1/2}$.
3. **Time-unit labeling**: `plot.pedhalflife()` and `print.pedhalflife()` now use the name of the `timevar` column, such as `Gen` or `Year`, for axis and summary labels.

# Changes in version 1.6.0 released on 20 Mar 2026
## New features
1. **Information-theoretic diversity half-life (`pedhalflife()`)**: Added `pedhalflife()` to track $f_e$, $f_a$, and $f_g$ across time points and fit a log-linear decay model for the rate of genetic diversity loss. The total loss rate $\lambda_{total}$ is decomposed into foundation bottleneck ($\lambda_e$), breeding bottleneck ($\lambda_b$), and genetic drift ($\lambda_d$). The diversity half-life $T_{1/2} = \ln 2 / \lambda_{total}$ is reported in units of the `timevar` column. S3 `print()` and `plot()` methods are provided, with both log-scale (`type = "log"`) and raw-scale (`type = "raw"`) views.

# Changes in version 1.5.0 released on 20 Mar 2026
## New features
1. **Shannon-entropy effective founders and ancestors (`feH`, `faH`)**: `pedcontrib()` and `pediv()` now compute two additional diversity statistics based on the Hill number of order $q = 1$ (Shannon entropy): `feH`, the effective number of founders under equal entropy weighting, and `faH`, the effective number of ancestors under equal entropy weighting. They satisfy $N_{\mathrm{Founder}} \ge f_e^H \ge f_e$ and $N_{\mathrm{Ancestor}} \ge f_a^H \ge f_a$, respectively. Both are computed as $\exp(-\sum p_i \ln p_i)$ and complement the classical quadratic ($q = 2$) effective numbers $f_e$ and $f_a$ (Lacy 1989; Boichard et al. 1997).

# Changes in version 1.4.1 released on 15 Mar 2026
## Bug fixes
1. **Fail-fast incomplete-pedigree analysis**: `inbreed()` and other completeness-sensitive analysis functions now error on row-truncated subsets with missing parent records. This prevents incorrect results, such as zero inbreeding, caused by calculating on partial ancestry data.

# Changes in version 1.4.0 released on 15 Mar 2026
## New features
1. **`tidyped` class redesign**: Refined the internal `tidyped` class architecture around a clearer metadata contract and safer S3/data.table interaction model for repeated downstream analysis and extension.
2. **Safer `tidyped` workflows**: Added `is_tidyped()`, `pedmeta()`, `has_inbreeding()`, and `has_candidates()` to make class checks and metadata inspection explicit.
3. **Fast candidate tracing from existing `tidyped` objects**: `tidyped()` now uses a fast path when the input is already a valid `tidyped` object and `cand` is supplied, avoiding repeated global validation and preprocessing.
4. **Workflow coverage and developer documentation**: Added a workflow vignette, a `tidyped` structure and extension vignette, and regression tests covering safe subsetting, `:=` by-reference behavior, and split workflow semantics.

## Bug fixes
1. **By-reference mutation for `tidyped`**: Replaced class and metadata attachment paths with `data.table::setattr()` so subsequent `:=` operations keep true by-reference behavior instead of writing into shallow copies.
2. **Safe row subsetting**: Added `[.tidyped` interception so incomplete subsets degrade to plain `data.table` objects with a warning, while complete subsets preserve `tidyped` structure and rebuild pedigree indices.
3. **Class recovery**: Core analysis entry points now cooperate with `ensure_tidyped()` and `validate_tidyped()` to recover valid `tidyped` objects after common class-dropping operations.

## Documentation
1. **Pkgdown article navigation**: Reorganized vignette order, restored `draw-pedigree` to the recommended reading sequence, and exposed `tidyped` developer notes through a dedicated pkgdown developer-documentation entry.

# Changes in version 1.3.5 released on 14 Mar 2026
## New features
1. **S3 class protection**: Added `as_tidyped()` and an internal `ensure_tidyped()` mechanism to handle cases where R operations such as `merge()`, `rbind()`, and `dplyr` verbs strip the custom S3 class from `data.table` objects. Major analysis functions can restore the class when the underlying data structure is still valid and inform the user.

## Bug fixes
1. **Analysis entry points**: Updated core analysis functions, including `pedstats()`, `pedne()`, `pediv()`, and `pedrel()`, to use the class-recovery logic.

# Changes in version 1.3.4 released on 14 Mar 2026
## Bug fixes
1. **`data.table` return visibility**: Functions returning `data.table` or `tidyped` objects now explicitly return with `[]` so results auto-print in the R console and knitted documents after internal `data.table` operations such as `:=` and `set*`. Affected functions included `pedancestry()`, `pedpartial()`, `pedne()`, `pedrel()`, and `tidyped()`.
2. **Side-effect prevention**: Updated `calc_ne_demographic()` to operate on a copy of the input pedigree instead of modifying the user's data by reference.

## Documentation
1. **Coding standard**: Added a `data.table` return-visibility rule to `Positron.md`.

# Changes in version 1.3.3 released on 14 Mar 2026
## Bug fixes
1. **Vignette API synchronization**: Replaced outdated `pedinbreed_class()` calls in the pedigree analysis vignette with the current `pedfclass()` interface and aligned examples with the current `reference`, `foundervar`, and `cycle` argument names.

## Documentation
1. **`pedigree-analysis.Rmd` rewrite**: Reorganized the main pedigree analysis vignette into thematic sections covering pedigree overview, pedigree completeness (`pedecg()`), generation intervals (`pedgenint()`), subpopulation structure (`pedsubpop()`), diversity indicators (`pediv()`), effective population size (`pedne()`), average relationship trends (`pedrel()`), inbreeding classification (`pedfclass()`), and ancestry / partial inbreeding diagnostics.
2. **Theory expansion**: Added formulas, interpretation notes, and breeding-use explanations for Equivalent Complete Generations (ECG), generation intervals, effective numbers of founders / ancestors / founder genomes (`f_e`, `f_a`, `f_g`), three effective population size definitions (`N_e` by demographic, inbreeding, and coancestry methods), and average additive relationship (`MeanRel`).
3. **Reference update**: Expanded the vignette bibliography to include Wright (1922, 1931), Lacy (1989), Boichard et al. (1997), Caballero and Toro (2000), Cervantes et al. (2011), and Gutierrez et al. (2008, 2009).

## Internal changes
1. **Analysis regression coverage**: Added unit tests verifying that `pedancestry()` proportions sum to 1 in a multi-line admixture pedigree and that `pedrel()` returns identical results between `compact = TRUE` and `compact = FALSE` on the same pedigree.

# Changes in version 1.3.2 released on 13 Mar 2026
## New features
1. **Core analysis examples**: Added `@examples` to `pedne()`, `pedecg()`, and `pedsubpop()`.
2. **Pedigree connectivity analysis (`pedsubpop()`)**: Enhanced `pedsubpop()` to distinguish pedigree splitting via `splitped()` from grouping and summary analysis. It now reports counts of total individuals, sires, dams, and founders within subgroups or connected components.

## API changes
1. **`pedfclass()` rename**: Renamed `pedinbreedclass()` to `pedfclass()` to align with the package naming guide and provide a shorter API.
2. **`pedfclass()` output refinement**: Renamed the returned class column from `F_Class` to `FClass` and added user-defined inbreeding class breakpoints through `breaks` and `labels`.
3. **`pedgenint()` parameter rename**: Renamed `cycle_length` to `cycle`.
4. **`pedgenint()` and `pedstats()` `unit` parameter**: Removed `"gen"` from `unit` options. The `unit` parameter now accepts `"year"`, `"month"`, `"day"`, or `"hour"`.
5. **`pedgenint()` `timevar` definition**: Clarified `timevar` as a birth-date column. Numeric year inputs are converted to `Date` values using `"YYYY-07-01"` with an informational message. Character date strings are parsed with `as.POSIXct(..., tz = "UTC")`.

## Improvements
1. **`pedsubpop()` documentation**: Refined internal documentation to clarify its use cases alongside `splitped()`.
2. **`.parse_to_numeric_time()` rewrite**: Rewrote the internal time parser to handle `Date`, `POSIXct`, character date strings, and numeric years. `POSIXct` conversions now use `tz = "UTC"` to avoid DST-related artifacts.

## Bug fixes
1. **`vispstat()` pathway filter**: Fixed an issue where the generation-interval bar chart could drop pathways because of an overly broad `factor()` filter. The filter now uses explicit `%in% c("SS", "SD", "DS", "DD")` subsetting.

# Changes in version 1.3.1 released on 12 Mar 2026
## New features
1. **Selfing support in `tidyped()`**: Added the `selfing` argument to support plant and aquaculture pedigrees where an individual can appear as both sire and dam, resolving sex-conflict errors for biologically valid selfing pedigrees (#10).
2. **`pedrel()` ancestral tracing**: `pedrel()` now uses full ancestral tracing via `tidyped(ped, cand = ...)` when calculating subgroup relationships, avoiding underestimation caused by ancestor truncation in deep-inbred populations.

## API changes
1. **`pedancestry()` parameter rename**: Renamed `labelvar` to `foundervar` and `labels` to `target_labels` to make the ancestry-tracing interface more explicit. Old argument names are no longer supported because the function was still under active development.
2. **`pedecg()` parameter cleanup**: Removed the short-lived `reference` argument, which only filtered rows after a full ECG pass and did not define a reference population or prune the pedigree before calculation.

## Improvements
1. **Academic nomenclature alignment**: Updated documentation for `pedrel()` and `pedne()` to distinguish additive genetic relationship ($a_{ij}$) from coancestry ($f_{ij}$). `pedrel()` now states that it returns $a_{ij} = 2f_{ij}$, and `pedne()` documents that its `"coancestry"` method is based on $f_{ij}$.
2. **Monoecious individuals**: Individuals acting as both parents are identified as `"monoecious"` in the `Sex` column.
3. **`visped()` monoecious styling**: `visped()` uses teal (`#26a69a`) for `"monoecious"` individuals.
4. **Role-specific pedigree edges**: Pedigree edges are colored by the parent's role in each mating (sire, dam, or selfing) rather than invariant node sex.
5. **`tidyped` summary methods**: `summary()` and `print()` methods for `tidyped` objects now report the count and percentage of monoecious individuals.
6. **`pedancestry()` initialization**: Optimized initialization on large pedigrees by using vectorized matrix indexing, reducing overhead for pedigrees with more than 25,000 nodes.

## Bug fixes
1. **`pedrel()` deep-inbreeding regression coverage**: Added a unit test verifying relationship calculation in deep-inbreeding scenarios, including Gen 4 relationships reaching 1.0.

# Changes in version 1.3.0 released on 10 Mar 2026
## New features
1. **Founder genome equivalents (`fg`)**: Integrated calculation of founder genome equivalents into `pediv()`. The implementation evaluates mean coancestry with diagonal intra-cohort correction through adaptive scaling, keeping computational costs linear relative to the reference cohort size.
2. **Reproducible parameter inference**: Added a `seed` argument to `pedne()` and `pediv()` for reproducible sampling in effective population size and `fg` calculations that use Monte Carlo approximations.

# Changes in version 1.2.3 released on 08 Mar 2026
## Bug fixes
1. **Trace edge highlighting in `visped()`**: Fixed incorrect edge highlighting when using `trace = "all"`. When a node was highlighted as both an ancestor and a parent of descendants, cross-path edges could be highlighted incorrectly. The fix separates upward and downward trace paths and uses `trace_edges` to control highlighted edges.
2. **Focal-node upward edge in `trace = "down"`**: Fixed an issue where the focal node's upward connection to its parents' family node was highlighted when tracing downward only. Individual-to-family edges are now highlighted only when the individual appears as a child in the traced path.

# Changes in version 1.2.2 released on 08 Mar 2026
## New features
1. **Unified diversity analysis (`pediv()`)**: Added `pediv()` as a single entry point that aggregates founder contributions (`f_e`), ancestor contributions (`f_a`), and three `N_e` estimates (coancestry, inbreeding, and demographic) into a `pediv` S3 object. A `print.pediv()` method provides a formatted summary table.
2. **`complex_ped` dataset**: Added `complex_ped`, a multi-generation pedigree dataset for testing deeper ancestry tracing and cross-generation diversity analyses.

## API changes
1. **Reference-population parameter rename**: Standardized the reference population parameter name across relevant analysis functions: `pedne(..., reference = NULL)`, `pedcontrib(..., reference = NULL)`, and `pedrel(..., reference = NULL)`. These previously used `cand`. Old `cand` arguments are no longer supported in those functions.

## Documentation
1. **`pedigree-analysis.Rmd` rewrite**: Restructured the pedigree analysis vignette with expanded theory explanations for `f_e`, `f_a`, and `N_e`, updated examples using `pediv()` and `reference`, and additional interpretation for breeding decisions.
2. **Workspace organization**: Moved development-only files (`MACOS_OPENMP_FIX.md`, `manuscript.md`, and analysis scripts) into `sandbox/`, with corresponding `.gitignore` and `.Rbuildignore` rules.

# Changes in version 1.2.1 released on 07 Mar 2026
## New features
1. **Ancestral analysis (`pedcontrib()`)**: Added algorithms for assessing genetic diversity through gene-origin probabilities. `pedcontrib()` computes effective number of founders (`f_e`) through recursive gene derivation and effective number of ancestors (`f_a`) through Boichard's iterative algorithm.
2. **Missing-parent conservation in `pedcontrib()`**: Single missing parents are augmented with temporary phantom parents before processing so probability mass is conserved for half-founder records.
3. **Ancestry proportions (`pedancestry()`)**: Added `pedancestry()` to trace line origins and monitor surviving gene proportions from specified historic founder lines or strains to descendants.
4. **Partial inbreeding (`pedpartial()`)**: Added `pedpartial()` to decompose the overall inbreeding coefficient into fractions attributed to targeted ancestors, following Meuwissen and Luo (1992).
5. **`half_founder_ped` dataset**: Added `half_founder_ped`, an ENDOG-derived dataset containing records with a single known parent for testing phantom-parent corrections.

## Performance
1. **Peeling core engine**: Rebuilt the C++ array engine backing the `f_a` and `f_e` calculations. The bounded O(K × N) array-state implementation avoids excessive latency on deep pedigrees with more than 180,000 nodes.

## Documentation
1. **Reference index**: Expanded `_pkgdown.yml` mappings to expose pedigree statistical functions such as `pedancestry()`, `pedcontrib()`, `pedpartial()`, and `pedecg()` in the reference documentation.
2. **Analysis vignette**: Updated `vignettes/pedigree-analysis.Rmd` with Boichard-style genetic bottleneck interpretations (`f_e` versus `f_a`) and examples for targeted lineage flow.

# Changes in version 1.2.0 released on 04 Mar 2026
## New features
1. **Effective population size methods in `pedne()`**: Expanded `pedne()` to support three methods for estimating `N_e`: `"coancestry"` based on the rate of coancestry ($\Delta c$), `"inbreeding"` based on the individual rate of inbreeding ($\Delta F$), and `"demographic"` based on the numbers of breeding males ($N_m$) and females ($N_f$). The coancestry method is now the default and may provide earlier signals of diversity loss in selected populations than inbreeding-based estimates.
2. **Parallel and sampled coancestry estimation**: Added OpenMP multi-threading for the `method = "coancestry"` path via `ncores`, and added `nsamples` for sampled estimation on large pedigrees.

## Performance
1. **Coancestry backend**: Implemented the C++ backend `cpp_calculate_sampled_coancestry_delta()` using `RcppArmadillo`, replacing previous R-level coancestry calculations for larger datasets.

## Documentation
1. **`pedne()` parameter scopes**: Documented that `ncores` and `nsamples` apply to `method = "coancestry"`.
2. **`pedne()` method descriptions**: Expanded documentation for the three `N_e` calculation methods to support method selection.

# Changes in version 1.1.1 released on 02 Mar 2026
## New features
1. **`pedgenint()` sex-independent pathways**: Added `SO` (sire-to-offspring) and `DO` (dam-to-offspring) generation intervals alongside the standard four pathways. These pathways support settings such as aquaculture or early-stage screening where offspring sex may be unknown.

## API changes
1. **`pedne()` interface standardization**: Renamed `timevar` to `by` and `cohort` to `cand`; removed unused or misleading parameters (`unit`, `cycle_length`, and `maxgen`). Old `timevar` and `cohort` arguments are retained with deprecation warnings.
2. **`vismat()` parameter alignment**: Renamed `grouping` to `by` for grouping consistency. The old `grouping` argument is retained with a deprecation warning.

## Bug fixes
1. **`pedrel()` correctness**: Fixed an error where mean average relatedness was calculated by summing the full relationship matrix, including traced ancestors, and dividing by only the target subgroup size. The function now subsets the relationship matrix and handles `NUsed < 2`. Output columns `N` and `MeanRel` were replaced with `NTotal`, `NUsed`, and `MeanRel`.
2. **`pedgenint()` aggregation**: `pedgenint()` now outputs the appropriate unweighted mixture standard deviation for generation intervals alongside the unweighted four-pathway average.
3. **`pedgenint()` sample size**: Fixed underestimation of the `Average` pathway `N` by evaluating all parent-offspring pairs through `calc_all_pathway()`.
4. **`pedcontrib()` effective numbers**: Standardized effective founder (`Ne_f`) and effective ancestor (`Ne_a`) calculations so they use the full untruncated cohort before reporting the top-ranked rows. The result list now includes total and reported count variables.
5. **`pedcontrib()` deep-pedigree latency**: Replaced a string-named-vector backward pass with a pure integer-indexed backward pass for large and deep pedigrees, including cases with more than 200,000 records.
6. **`pedpartial()` and `pedancestry()` input compatibility**: Missing numeric identifier columns in incoming pedigrees, such as objects created with `addnum = FALSE`, no longer break `pedpartial()` or `pedancestry()`. The pedigree propagation loop in `pedancestry()` was also simplified to use direct vector lookup.
7. **`pedne()` performance bottleneck**: Removed obsolete O(N^2) individual traversal in `calc_ancestral_f()`, using the direct formula by Gutierrez et al.

# Changes in version 1.1.0 released on 01 Mar 2026
## New features
1. **Pedigree analysis module**: Added pedigree analysis and population-genetics functions: `pedstats()` for summary and demographic statistics; `pedrel()` for average relatedness within groups; `pedgenint()` for sire-sire, sire-dam, dam-sire, dam-dam, and overall generation intervals; `pedcontrib()` for founder and ancestor genetic contributions; `pedancestry()` for ancestral-lineage proportions; and `pedpartial()` for partial inbreeding attributed to specified ancestors.
2. **Pedigree analysis visualization**: Added `vispstat()` for plotting generation-interval bar charts and depth-related distributions such as Equivalent Complete Generations.

# Changes in version 1.0.1 released on 31 Jan 2026
## Bug fixes
1. **Compact matrix correctness**: Fixed incorrect `A`, `D`, and `AA` relationship values in `compact = TRUE` mode for parent-offspring and avuncular pairs caused by merging parent individuals with their non-parent siblings.
2. **Pedigree compression strategy**: Updated compaction logic to preserve the genetic identity of any individual that appears as a sire or dam, ensuring parents have unique entries in the relationship matrix.
3. **Sibling row and column expansion**: Fixed `expand_pedmat()` so sibling off-diagonal elements are calculated from parent kinship rather than by duplicating representative diagonal values.
4. **Generation alignment**: Fixed `tidyped(..., genmethod = "bottom")` so sibling consistency is prioritized over mate alignment, ensuring full siblings are aligned to the same generation.
5. **`visped()` edge highlighting**: Fixed edge highlighting so relationship edges are only emphasized when `trace` is used.
6. **Shared-parent and shared-child paths**: Corrected edge highlighting for cases where a parent has multiple families or a family has multiple children.
7. **`visped()` layout with `showf = TRUE`**: Fixed a layout optimization failure by using immutable individual IDs.

# Changes in version 1.0.0 released on 24 Jan 2026
## Breaking changes
1. **Matrix API names**: Renamed `pedmatrix()` to `pedmat()`, `expand_pedmatrix()` to `expand_pedmat()`, and `summary_pedmatrix()` to `summary_pedmat()`. Legacy function names were removed.
2. **`pedmat()` default method**: Changed the default `method` of `pedmat()` from `"f"` to `"A"`.
3. **Thread argument name**: Standardized `n_threads` to `threads` across functions.

## New features
1. **Family assignment and summary**: `tidyped()` now assigns a `Family` column identifying full-sib groups. `summary.tidyped()` reports family counts, sizes, largest families, and offspring summaries.
2. **Pedigree splitting (`splitped()`)**: Added `splitped()` to detect disconnected pedigree components, exclude isolated individuals, and return re-indexed `tidyped` objects for separate analysis or visualization.
3. **Relationship matrix support**: `pedmat()` supports additive (`A`, `Ainv`), dominance (`D`, `Dinv`), and additive-by-additive epistatic (`AA`, `AAinv`) relationship matrices.
4. **Relationship matrix visualization (`vismat()`)**: Added `vismat()` for heatmaps and histograms of `pedmat` objects, `tidyped` objects, and standard matrices. Heatmaps can be annotated with family groups when a pedigree is provided.

## Documentation
1. **CRAN preparation**: Updated vignette generation, S3 method dispatch, and documentation for CRAN submission.

## Internal changes
1. **Portable compilation**: Standardized `src/Makevars` for cross-platform compatibility.
2. **Dependencies**: Moved `RcppArmadillo` to `LinkingTo`.

# Changes in version 0.7.3 released on 13 Jan 2026
## Breaking changes
1. **`pedmatrix()` return value and method selection**: `pedmatrix()` now requires a single `method` value, such as `method = "A"`. It returns the requested matrix or vector directly instead of a named list. Requesting multiple methods in one call now raises an error; use repeated calls for multiple outputs.

## New features
1. **Rcpp relationship calculations**: Added `pedmatrix()` with Rcpp implementations for the additive relationship matrix (`A`), sparse inverse additive matrix (`Ainv`) using Henderson's rules, dominance matrix (`D`), and inbreeding coefficients (`f`) using the Meuwissen and Luo (1992) path-tracing algorithm.

## Improvements
1. **Default inbreeding backend**: `inbreed()` now uses the native Rcpp implementation by default, moving `nadiv` to `Suggests`.
2. **Documentation and website**: Updated package documentation and vignettes. The package website is available at [https://luansheng.github.io/visPedigree/](https://luansheng.github.io/visPedigree/).

# Changes in version 0.7.2 released on 12 Jan 2026
## New features
1. **Flexible generation assignment**: Added the `genmethod` parameter to `tidyped()`. Users can choose `"top"` (top-aligned, default) or `"bottom"` (bottom-aligned) generation inference. The `"top"` method aligns founders at generation 1; the `"bottom"` method aligns terminal nodes at the bottom for visualization of unrelated introduced parents.

## Improvements
1. **Default generation assignment**: Changed the default generation assignment method to `"top"`.
2. **Pkgdown website**: Generated and published the package website at [https://luansheng.github.io/visPedigree/](https://luansheng.github.io/visPedigree/).
3. **Automated documentation deployment**: Added a GitHub Actions workflow for documentation updates and GitHub Pages deployment.

# Changes in version 0.7.1 released on 11 Jan 2026
## Improvements
1. **User feedback**: Standardized filtering notifications. The message `"Note: Removed N isolated individuals..."` now appears consistently for all pedigree sizes when Gen 0 individuals are present.
2. **Trace semantics**: Corrected `trace = "all"` in `tidyped()` and `visped()` so it returns the union of ancestors and descendants (`"up"` plus `"down"`) instead of the entire connected component.

## Performance
1. **Large-pedigree `visped()` rendering**: Optimized `visped()` for large pedigrees through attribute handling and vectorized rendering, reducing repeated `igraph` attribute lookups for 100,000+ individuals.
2. **Vectorized tracing**: Refactored `trace_ped_candidates()` in `tidyped()` to use vectorized `igraph::neighborhood()` calls. In one benchmark, 37,000 candidates in a 178,000-individual pedigree were traced in about 1.2 s.
3. **Early isolated-individual filtering**: Implemented early filtering of isolated Gen 0 individuals in `prepare_ped_graph()` to streamline graph conversion and layout.

# Changes in version 0.7.0 released on 10 Jan 2026
## Breaking changes
1. **Graph-based `tidyped()` core**: Reimplemented the pedigree tidying engine around a directed acyclic graph representation, with cycle detection and generation inference based on topological sorting.
2. **Modular `visped()` architecture**: Split the previous monolithic `visped.R` implementation into `visped_layout.R`, `visped_graph.R`, `visped_style.R`, and `visped_render.R`.

## New features
1. **`visped()` layout parameters**: Added `pagewidth` to control PDF page width and `symbolsize` to scale node sizes relative to label dimensions.
2. **Two-pass rendering**: Added a two-pass strategy in `plot_ped_igraph()` so edges connect at node centers in vector PDF outputs.
3. **Highlight tracing**: Added ancestry and descendant highlighting through the `trace` parameter in `visped()`.

## Bug fixes
1. **`outline = TRUE` rendering**: Fixed rendering failure by correcting graph attribute indexing.
2. **Coordinate precision**: Improved coordinate calculation precision to reduce overlap in high-density generations.

## Internal changes
1. **Testing framework**: Modernized the unit testing suite to `testthat` 3rd edition and removed legacy `context()` warnings.

# Changes in version 0.6.2 released on 01 Jan 2026
## New features
1. **`summary.tidyped()`**: Added a `summary()` method for `tidyped` objects to report counts of individuals, founders, sex distribution, and related pedigree statistics.

## Bug fixes
1. **`tidyped(..., inbreed = TRUE)`**: Fixed failure caused by incorrect class assignment order.
2. **`visped(..., showf = TRUE)`**: Added a warning when the `f` column is missing instead of erroring.
3. **Vignette navigation**: Fixed broken internal navigation links in package vignettes.

# Changes in version 0.6.1 released on 30 Dec 2025
## New features
1. **Highlight styling**: Implemented opaque highlighting effects for clearer visualization.
2. **`visped()` trace direction**: Added the `trace` option to `visped()` to control ancestry tracing direction.

# Changes in version 0.6.0 released on 28 Dec 2025
## New features
1. **Strict `tidyped` S3 class structure**: Added the `new_tidyped()` constructor and `validate_tidyped()` validator for `tidyped` objects.

# Changes in version 0.5.0 released on 26 Dec 2025
## New features
1. **`visped()` highlighting**: Added the `highlight` parameter to `visped()` for highlighting specific individuals with a character vector of IDs or a list of custom colors.
2. **`visped()` inbreeding display**: Added the `showf` parameter to `visped()` to display inbreeding coefficients on the pedigree graph.
3. **`tidyped()` inbreeding calculation**: Added the `inbreed` parameter to `tidyped()` to calculate inbreeding coefficients using `nadiv`.
4. **Standalone `inbreed()`**: Refactored `inbreed()` as a standalone tool that operates on `tidyped` objects.

## Performance
1. **`repeloverlap()`**: Optimized `repeloverlap()` with `data.table`.

## Bug fixes
1. **`visped()` compact-highlight-showf crash**: Fixed a crash when combining `compact = TRUE`, `highlight`, and `showf = TRUE` by refactoring `ped2igraph()` to delay label modification until after layout calculation.
2. **Documentation grammar**: Fixed grammar and phrasing across function documentation for CRAN compliance.
3. **`R CMD check` notes**: Fixed `data.table` non-standard evaluation notes by adding `R/globals.R`.

# Changes in version 0.4.1 released on 25 Dec 2025
## Bug fixes
1. **Edge overlap detection**: Fixed overlapping edge detection for small pedigree graphs.
2. **Compact-mode coloring**: Improved coloring consistency for compact mode.

# Changes in version 0.2.6 released on 31 Mar 2020
## Bug fixes
1. **Candidate generation tracing**: Fixed a bug where candidates were traced to generation `n + 1` when `tracegen = n`. Reported by Mianyu Liu.

# Changes in version 0.2.5 released on 25 Feb 2020
## Bug fixes
1. **`trace = "all"` in `tidyped()`**: Fixed cases where `tidyped()` did not work with `trace = "all"` ([issue comment](https://github.com/luansheng/visPedigree/issues/2#issue-568599008)).

# Changes in version 0.2.4.1 released on 24 Feb 2020
## Bug fixes
1. **Repeated `tidyped()` calls**: Fixed an unexpected column named `NA` when a `tidyped` object was tidied again with `tidyped()`.

# Changes in version 0.2.4 released on 12 Jun 2019
## Bug fixes
1. **Input mutation**: Fixed a bug where a `data.table` supplied as `ped` could be changed by `tidyped()` or `visped()`.

# Changes in version 0.2.3 released on 05 Mar 2019
## Bug fixes
1. **Generation inference**: Fixed incorrect generation-number inference for individuals.

# Changes in version 0.2.2 released on 28 Jan 2019
## Bug fixes
1. **Candidate retention**: Fixed a bug where a tidied pedigree did not include candidates absent from the original `Ind` column when `cand` was not `NULL`.

# Changes in version 0.2.1 released on 17 Nov 2018
## Bug fixes
1. **Node overlap**: Repelled overlapping nodes caused by very small differences (digits > 7) among node x positions.
