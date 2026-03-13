# Changes in version 1.3.2 released on 13 Mar 2026
## New features
1. **Added Comprehensive Examples**: Added `@examples` to core analysis functions `pedne()`, `pedecg()`, and `pedsubpop()` to improve documentation completeness and provide immediate value to users.
2. **Pedigree Connectivity Analysis (`pedsubpop`)**: Enhanced `pedsubpop()` to better distinguish between pedigree splitting (via `splitped`) and grouping/summary analysis. It now provides clear counts of total individuals, sires, dams, and founders within subgroups or connected components.

## API Changes
1. **`pedfclass()` rename**: Renamed the inbreeding-class summary helper from `pedinbreedclass()` to `pedfclass()` to align with the package naming guide and provide a shorter, clearer user-facing API.
2. **`pedfclass()` output refinement**: Renamed the returned class column from `F_Class` to `FClass`, and added support for user-defined inbreeding class breakpoints through the `breaks` and `labels` arguments.
3. **`pedgenint()` parameter rename**: Renamed `cycle_length` to `cycle` for consistency with the package naming guide.
4. **`pedgenint()` / `pedstats()` `unit` parameter**: Removed `"gen"` from `unit` options. The `unit` parameter now only accepts `"year"`, `"month"`, `"day"`, or `"hour"`. The previous `"gen"` option produced incorrect results when combined with date inputs.
5. **`pedgenint()` `timevar` definition**: Clarified `timevar` as a **birth date** column. Numeric year inputs (e.g., `2020`) are now automatically converted to `Date` (`"YYYY-07-01"`) with an informational message. Character date strings are parsed via `as.POSIXct()` with `tz = "UTC"` to avoid DST artifacts.

## Minor improvements and bug fixes
1. **Documentation Audit**: Refined internal documentation for `pedsubpop()` to clarify its use cases alongside `splitped()`.
2. **`vispstat()` pathway filter fix**: Fixed an issue where the generation-interval bar chart could silently drop pathways due to an overly broad `factor()` filter. Now uses explicit `%in% c("SS", "SD", "DS", "DD")` subsetting.
3. **`.parse_to_numeric_time()` rewrite**: Completely rewrote the internal time parser to handle `Date`, `POSIXct`, character date strings, and numeric years robustly. All `POSIXct` conversions now use `tz = "UTC"` to prevent DST-related artifacts.

# Changes in version 1.3.1 released on 12 Mar 2026
## New features
1. Added the `selfing` argument to `tidyped()` to support plant and aquaculture pedigrees where an individual can appear as both Sire and Dam, resolving biologically impossible sex conflict errors (#10).
2. **`pedrel()` logic upgrade**: Modified `pedrel()` to use full ancestral tracing via `tidyped(ped, cand = ...)` when calculating sub-group relationships. This ensures that relationships in deep-inbred populations (e.g., full-sib mating over multiple generations) are calculated correctly rather than being underestimated due to ancestor truncation.

## API Changes
1. **`pedancestry()` parameter rename**: Renamed `labelvar` to `foundervar` and `labels` to `target_labels` to align with the package naming guide and make the ancestry-tracing interface more explicit. Old argument names are no longer supported because this function is still under active development.
2. **`pedecg()` parameter cleanup**: Removed the short-lived `reference` argument. It only filtered rows after a full ECG pass and did not define a true reference population or prune the pedigree before calculation. Users should subset the returned table directly if needed.

## Minor improvements and bug fixes
1. **Academic nomenclature alignment**: Updated documentation for `pedrel()` and `pedne()` to explicitly distinguish between **Additive Genetic Relationship ($a_{ij}$)** and **Coancestry ($f_{ij}$)**. 
    - `pedrel()` now clearly states it returns $a_{ij} = 2f_{ij}$.
    - `pedne()` documentation now specifies that its `"coancestry"` method is based on $f_{ij}$.
2. Individuals acting as both parents are now identified as `"monoecious"` in the `Sex` column.
3. `visped()` now uses a distinct teal color (`#26a69a`) to render `"monoecious"` individuals, ensuring clear visual separation from males, females, and highlighted nodes.
4. Pedigree edges are now colored based on the parent's role in a specific mating (Sire blue, Dam gold, Selfing teal) rather than invariant node sex, allowing monoecious individuals to display role-appropriate connection colors.
5. `summary()` and `print()` methods for `tidyped` objects now accurately report the count and percentage of monoecious individuals.
6. **Efficiency Optimization**: Optimized `pedancestry()` initialization on large pedigrees by using vectorized matrix indexing, significantly reducing overhead for pedigrees with >25k nodes.
7. Added a new unit test for `pedrel()` to verify correct relationship calculation in deep-inbreeding scenarios (Gen 4 relationships reaching 1.0).

# Changes in version 1.3.0 released on 10 Mar 2026

## New Features
1. **Founder Genome Equivalents ($f_g$)**: Integrated the robust calculation of Founder Genome Equivalents into `pediv()`. It directly evaluates the mean coancestry while properly correcting for diagonal intra-cohort elements through adaptive scaling, keeping computational costs linear relative to the reference cohort size.
2. **Reproducible Parameter Inference**: Added a `seed` argument to both `pedne()` and `pediv()` functions functions enabling precise reproducible sampling for effective population size estimations (Ne) and $f_g$ computations using Monte Carlo approximations.

# Changes in version 1.2.3 released on 08 Mar 2026

## Bug Fixes
1. **Trace edge highlighting in `visped()`**: Fixed incorrect edge highlighting when using `trace = "all"`. Previously, when a node was highlighted as both an ancestor (via upward tracing) and a parent of descendants (via downward tracing), the cross-path edges were incorrectly highlighted. For example, `visped(tp, highlight = "X", trace = "all")` would incorrectly highlight the edge from N to Z1/Z2, even though that parent-child relationship is not on X's trace path. The fix separates the up and down trace paths and uses `trace_edges` to precisely control which edges are highlighted.
2. **Focal node upward edge in `trace = "down"`**: Fixed an issue where the focal node's upward connection to its parents' family node was incorrectly shown as highlighted when tracing downward only. For example, `visped(tp, highlight = "X", trace = "down")` would show X's edge to the UxV family node in solid black, even though X's ancestors are not part of the downward trace. Now, `individual → family` edges are only highlighted when the individual appears as a child in the traced path.

# Changes in version 1.2.2 released on 08 Mar 2026

## New Features
1. **Unified Diversity Analysis (`pediv`)**: Added `pediv()` as a single entry-point wrapper that aggregates founder contributions ($f_e$), ancestor contributions ($f_a$), and all three Ne estimates (coancestry, inbreeding, demographic) into one consolidated `pediv` S3 object. A dedicated `print.pediv()` method provides a formatted summary table.
2. **New Dataset (`complex_ped`)**: Added `complex_ped`, a multi-generation pedigree dataset suitable for testing deeper ancestry tracing and cross-generation diversity analyses.

## API Changes
1. **Parameter Rename (`cand` → `reference`)**: Standardized the reference population parameter name across all relevant analysis functions:
    - `pedne(... , reference = NULL)` (previously `cand`)
    - `pedcontrib(... , reference = NULL)` (previously `cand`)
    - `pedrel(... , reference = NULL)` (previously `cand`)
    
    The new name better reflects population genetics terminology, where the target group is the **reference population** rather than a set of "candidates". *(Note: Old `cand` argument is no longer supported; please update existing scripts.)*

## Documentation
1. **Vignette Rewrite (`pedigree-analysis.Rmd`)**: Completely restructured the pedigree analysis vignette with expanded theory explanations for $f_e$, $f_a$, and Ne, updated code examples using `pediv()` and the new `reference` parameter, and improved narrative linking the statistical outputs to practical breeding decisions.
2. **Workspace Reorganization**: Moved development-only files (`MACOS_OPENMP_FIX.md`, `manuscript.md`, analysis scripts) into `sandbox/` to keep the package root clean. Added corresponding `.gitignore` and `.Rbuildignore` rules.

# Changes in version 1.2.1 released on 07 Mar 2026

## New Features
1. **Ancestral Analysis (`pedcontrib`)**: Added robust algorithms for assessing genetic diversity through gene origin probabilities. Computes the **effective number of founders ($f_e$)** via recursive gene derivation and the **effective number of ancestors ($f_a$)** via Boichard's iterative algorithm.
2. **Missing Parent Conservation ("Phantom Parents")**: Implemented correct probability mass conservation. In `pedcontrib`, single missing parents (half-founders) are seamlessly augmented with temporary "phantom parents" before processing, overcoming the critical issue of gene probability leakage found in earlier tools.
3. **Ancestry Proportions (`pedancestry`)**: Added `pedancestry()` function to trace line origins and monitor the surviving proportion of genes from specified historic founder lines or strains down to modern descendants.
4. **Partial Inbreeding (`pedpartial`)**: Engineered the Meuwissen & Luo (1992) based partial inbreeding decomposition `pedpartial()`. Enables breaking down the overall inbreeding coefficient into discrete fractions attributed to specifically targeted ancestors.
5. **New Dataset (`half_founder_ped`)**: Added empirical ENDOG dataset containing instances of strictly missing single parents (sire known, dam unknown, etc.) specifically engineered to test and validate phantom-parent corrections.

## Performance
1. **Peeling Core Engine**: Rebuilt the C++ core array engine backing the $f_a$ and $f_e$ calculations. Execution latency for incredibly massive and deep graphs (>180,000 nodes) was resolved, avoiding hanging scenarios by limiting computational bounds optimally inside $O(K \times N)$ array states.

## Documentation
1. **Analysis Indexing**: Expanded `_pkgdown.yml` configuration mapping to fully expose all newly engineered high-level pedigree statistical functions (`pedancestry`, `pedcontrib`, `pedpartial`, `pedecg`, etc.) within the main Reference documentation.
2. **Analysis Vignettes**: Updated `vignettes/pedigree-analysis.Rmd` carefully illustrating Boichard's genetic bottleneck interpretations ($f_e$ vs $f_a$) alongside new code examples for tracing targeted lineage flows.

# Changes in version 1.2.0 released on 04 Mar 2026

## New Features
1. **Enhanced Effective Population Size (Ne) Calculation**:
    - The `pedne()` function has been significantly expanded and now supports three robust methods for estimating Ne in breeding populations:
        - **`coancestry`** (New Default): Based on the rate of coancestry ($\Delta c$). This method captures the loss of genetic potential and is considered the "gold standard" for populations under selection pressure. It typically yields a smaller, more conservative Ne estimate than inbreeding-based methods, providing a better early warning signal for genetic diversity loss.
        - **`inbreeding`**: Based on the individual rate of inbreeding ($\Delta F$). This method reflects the realized inbreeding but may overestimate Ne in managed populations where mating between relatives is actively avoided.
        - **`demographic`**: A census-based method using the number of breeding males ($N_m$) and females ($N_f$).
2. **Parallel Processing Support**:
    - Introduced **OpenMP multi-threading** for the computationally intensive `coancestry` method. Users can now specify `ncores` to speed up large-scale matrix calculations.
    - Added a `nsamples` parameter to allow efficient estimation on massive pedigrees by sampling subsets of each cohort.

## Performance
1. **C++ Optimization**:
    - Implemented a high-performance C++ backend (`cpp_calculate_sampled_coancestry_delta`) using `RcppArmadillo`. This replaces the previous R-based logic for coancestry calculations, enabling the analysis of much larger datasets.

## Documentation
1. **Clarified Parameter Scopes**: Updated documentation for `pedne()` to explicitly state that `ncores` and `nsamples` parameters are specific to the `method = "coancestry"` calculation path.
2. **Method Descriptions**: Expanded details on the three Ne calculation methods to help users choose the most appropriate metric for their breeding program.

# Changes in version 1.1.1 released on 02 Mar 2026

## New Features
1. **pedgenint Sex-Independent Pathways**: Added evaluation of `SO` (Sire-to-Offspring) and `DO` (Dam-to-Offspring) generation intervals alongside the standard 4 pathways. This is especially useful for aquatic species (like shrimp) or early-stage screening where offspring sex might remain unknown.

## API Changes and Refactoring
1. **pedne Interface Standardization**: 
    - Renamed arguments `timevar` to `by`, and `cohort` to `cand` to harmonize parameter naming conventions across the package.
    - Removed unused and misleading parameters (`unit`, `cycle_length`, `maxgen`). The effective population size Ne calculation innately depends on Equivalent Complete Generations (ECG), making it independent of scalar temporal units.
2. **vismat Parameter Alignment**: Renamed `grouping` argument to `by` to maintain grouping consistency.  
*(Note: Old arguments `timevar`, `cohort` in `pedne()` and `grouping` in `vismat()` are retained for backward compatibility but will display a deprecation warning.)*

## Bug fixes
1. **pedrel Correctness**: Fixed a critical calculation bug in `pedrel()` where the mean average relatedness calculation erroneously divided the sum of the full relationship matrix (including all traced ancestors) by only the size of the target subgroup. It now cleanly subsets the relationship matrix, and correctly handles boundary limits (`NUsed < 2`). The output columns `N` and `MeanRel` behavior has been replaced with `NTotal`, `NUsed`, and `MeanRel`.
2. **pedgenint Aggregation**: Fixed `pedgenint()` to output appropriate unweighted mixture standard deviation for generating generation intervals alongside its unweighted 4-pathway average interval estimate.
3. **pedgenint Sample Size (N)**: Fixed an issue where the `Average` pathway N was severely underestimated. It now accurately evaluates all parent-offspring pairs via `calc_all_pathway()`.
4. **pedcontrib Accuracy**: Standardized effective founders (`Ne_f`) and effective ancestors (`Ne_a`) calculation in `pedcontrib()` to ensure they are calculated based upon the full un-truncated cohort before outputting strictly the `top` n-ranked figures. Results list has been augmented with variables tracking the `_total` and `_reported` count values.
5. **pedcontrib Deep Pedigree Latency**: Replaced a string-named vector backward pass with a pure integer-indexed backward pass, resolving instances where evaluating contributions on deep, large pedigrees (e.g., > 200,000 records) would hang indefinitely due to scaling constraints.
6. **pedpartial / pedancestry Input Compatibility**: Ensured missing numeric identifiers in incoming pedigrees (e.g. `addnum = FALSE`) do not break `pedpartial()` or `pedancestry()`. Increased performance of the pedigree propagation loop in `pedancestry` by dropping an internal array linear probe algorithm with an immediate linear vector lookup.
7. **pedne Performance bottleneck**: Removed an obsolete `O(N^2)` individual traversal evaluation (`calc_ancestral_f()`), streamlining calculation purely around the efficient direct formula by Gutiérrez et al.

# Changes in version 1.1.0 released on 01 Mar 2026

## New Features
1. **Pedigree Analysis Module**: Introduced a comprehensive suite of pedigree analysis and population genetics tools.
    - `pedstats()`: Calculate holistic and demographic statistics.
    - `pedrel()`: Formulate average relatedness within specific population groupings.
    - `pedgenint()`: Compute distinct breeding pathways (SS, SD, DS, DD) and overall population generation intervals.
    - `pedcontrib()`: Determine genetic contributions from founders (`Ne_f`) and prominent ancestors (`Ne_a`) utilizing iterative gene flow derivations.
    - `pedancestry()`: Establish proportionality of ancestral lineages on subsequent descendants.
    - `pedpartial()`: Decompose inbreeding mechanisms to detect fractional/partial origins from core ancestors.
2. **Pedigree Analysis Visualization**: Added `vispstat()` to intuitively render bar charts of generation intervals and histogram distributions detailing depth tracking factors (like Equivalent Complete Generations).

# Changes in version 1.0.1 released on 31 Jan 2026
## Bug fixes
1. **Compact Matrix Correctness**: Fixed a critical data integrity bug in `compact = TRUE` mode where relationship values (A, D, AA) were incorrect for parent-offspring and avuncular pairs due to improper merging of parent individuals with their non-parent siblings.
2. **Pedigree Compression Strategy**: Updated compaction logic to preserve original genetic identity of any individual that appears as a sire or dam, ensuring parents always have unique entries in the relationship matrix.
3. **Sibling Row/Column Expansion**: Fixed `expand_pedmat()` to correctly handle sibling off-diagonal elements by dynamically calculating relationship values based on parent kinship, rather than simply duplicating representative diagonal values.
4. **Generation Alignment Logic**: Fixed `tidyped(..., genmethod = "bottom")` to prioritize **Sibling Consistency** (P1) over **Mate Alignment** (P2). This ensures that full siblings are always aligned to the same generation.
5. **`visped()` edge highlighting**: Fixed edge highlighting logic so relationship edges are only emphasized when `trace` is used.
6. **Shared-parent/shared-child paths**: Corrected edge highlighting for cases where a parent has multiple families or a family has multiple children.
7. **`visped()` layout**: Fixed layout optimization failure when `showf = TRUE`. The layout algorithm now correctly uses immutable individual IDs.

# Changes in version 1.0.0 released on 24 Jan 2026
## API Standardization (BREAKING)
To provide a clean and intuitive API for v1.0.0, core function names and behaviors have been standardized:
- **`pedmatrix`** is renamed to **`pedmat`**.
- **`pedmat` default `method` is now `"A"`** (Additive Relationship Matrix). Previously it was `"f"` (Inbreeding Coefficients).
- **`expand_pedmatrix`** is renamed to **`expand_pedmat`**.
- **`summary_pedmatrix`** is renamed to **`summary_pedmat`**.
- The parameter **`n_threads`** is standardized to **`threads`** across all functions.
- Legacy function names (`pedmatrix`, etc.) have been removed. Please use `pedmat()` directly.

## New Features
1. **Family Assignment and Summary**: 
    - `tidyped()` now automatically assigns and includes a `Family` column, identifying full-sib groups.
    - `summary.tidyped()` has been updated to provide family statistics (count, sizes, top largest families) and richer offspring analysis.
2. **Pedigree Splitting (`splitped`)**: Added `splitped()` function to detect and split disconnected pedigree components. It efficiently identifies independent sub-populations (connected components) using graph theory, excludes isolated individuals, and returns a list of re-indexed `tidyped` objects ready for separate analysis or visualization.
3. **Comprehensive Matrix Support**: `pedmat()` (formerly `pedmatrix`) now fully supports 6 types of genetic relationship matrices: Additive (A, Ainv), Dominance (D, Dinv), and Additive-by-Additive Epistatic (AA, AAinv).
4. **Relationship Matrix Visualization (`vismat`)**: Added `vismat()` function for visualizing relationship matrices (A, D, AA, etc.) with heatmaps and histograms. It supports `pedmat` objects, `tidyped` objects (auto-calculates A matrix), and standard matrices. Heatmaps can be annotated with family groups when a pedigree is provided.

## CRAN Submission & Internal Improvements
This release marks the first stable version 1.0.0, polished for CRAN.

1.  **Portable Compilation**: Standardized `src/Makevars` for cross-platform compatibility (removed GNU/platform-specific extensions).
2.  **Dependencies**: Moved `RcppArmadillo` to `LinkingTo` to optimize package structure.
3.  **Documentation & S3**: Fixed `vignette` generation, resolved `diag` S3 method dispatch, and cleaned up documentation for CRAN compliance.

# Changes in version 0.7.3 released on 13 Jan 2026
## New behavior (BREAKING)
1. **Simplified `pedmatrix()` return and single-method enforcement**: `pedmatrix()` now requires a single `method` argument (e.g., `method = "A").` When a single method is requested, the function returns the corresponding matrix or vector directly (not a named list). Requesting multiple methods in one call will now raise an error. Use repeated calls for multiple outputs.

## New features
1. **High-Performance Genetic Relationship Calculations**: Introduced `pedmatrix()` function implemented in Rcpp for efficient computation of:
    - Additive relationship matrix (A) using the tabular recursive algorithm.
    - Sparse inverse additive matrix (A-Inverse) using Henderson's rules.
    - Dominance matrix (D) using the tabular approach.
    - Inbreeding coefficients (f) using the Meuwissen & Luo (1992) path-tracing algorithm.

## Improvements
1. **Default Inbreeding Calculation Method**: The `inbreed()` function now uses the native Rcpp implementation by default, moving the `nadiv` package to `Suggests`.
2. **Documentation and Website**: Updated package documentation and vignettes to reflect new features and improvements. The official package website is available at [https://luansheng.github.io/visPedigree/](https://luansheng.github.io/visPedigree/).

# Changes in version 0.7.2 released on 12 Jan 2026
## New features
1. **Flexible Generation Assignment**: Added `genmethod` parameter to `tidyped()`. Users can now choose between `"top"` (top-aligned, default) and `"bottom"` (bottom-aligned) methods for generation inference. 
    - The `"top"` method aligns founders at Generation 1, which is more appropriate for biological pedigrees and prevents "founder drift" in pedigrees with varying depths.
    - The `"bottom"` method aligns terminal nodes at the bottom, useful for visualizing introductions of unrelated exogenous parents.

## Improvements
1. **Default Logic Change**: Switched the default generation assignment method to `"top"` (top-down) for more intuitive biological visualization.
2. **Pkgdown Documentation**: Generated and published the official package website at [https://luansheng.github.io/visPedigree/](https://luansheng.github.io/visPedigree/).
3. **Automated CI/CD**: Added GitHub Actions workflow for automatic documentation updates and deployment via GitHub Pages.

# Changes in version 0.7.1 released on 11 Jan 2026
## Performance optimizations
1. **Large Pedigree Performance**: Optimized `visped` performance for displaying large pedigrees through efficient attribute handling and vectorized rendering. Computation time for 100k+ individuals reduced significantly by avoiding redundant `igraph` attribute lookups.
2. **Vectorized Tracing**: Refactored `trace_ped_candidates` in `tidyped` to use vectorized `igraph::neighborhood` calls, achieving ~150x speedup for large candidate lists (e.g., 37k candidates in a 178k individual pedigree traced in ~1.2s).
3. **Early Filtering**: Implemented unified early filtering of isolated individuals (Gen 0) in `prepare_ped_graph` to streamline downstream graph conversion and layout algorithms.

## Improvements
1. **User Feedback**: Standardized filtering notifications. The message "Note: Removed N isolated individuals..." now appears consistently for all pedigree sizes when Gen 0 individuals are present.
2. **Refined Tracing**: Corrected `trace = "all"` logic in both `tidyped` and `visped`. It now correctly retrieves the union of ancestors and descendants ("up" + "down") instead of the entire connected component (undirected search).

# Changes in version 0.7.0 released on 10 Jan 2026
## Breaking changes & Major Refactoring
1. **Graph-based `tidyped` Core**: Reimplemented the pedigree tidying engine using formal graph theory principles (Directed Acyclic Graphs). Improved loop detection and generation inference accuracy using topological sorting.
2. **Modular Architecture**: Split the monolithic `visped.R` into functional modules: `visped_layout.R`, `visped_graph.R`, `visped_style.R`, and `visped_render.R` for better maintainability.

## New features
1. **New Parameters in `visped()`**: 
    - `pagewidth`: Allows users to specify the PDF page width (default 200 inches) to accommodate different pedigree scales.
    - `symbolsize`: A scaling factor (default 1) to adjust node sizes relative to label dimensions, providing finer control over whitespace.
2. **Two-Pass Rendering Engine**: Introduced a two-pass strategy in `plot_ped_igraph()` to ensure edges connect exactly at node centers, eliminating visual gaps in vector PDF outputs.
3. **Enhanced Highlighting**: Added support for real-time ancestry and descendant highlighting via the `trace` parameter in `visped()`.

## Bug fixes
1. Fixed rendering failure in `outline = TRUE` mode by correcting attribute indexing in the graph object.
2. modernized the unit testing suite to `testthat` 3rd edition, removing all legacy `context()` warnings.
3. Improved coordinate calculation precision to prevent overlapping in high-density generations.

# Changes in version 0.6.2 released on 01 Jan 2026
## New features
1. Added `summary()` method for `tidyped` objects to provide quick pedigree statistics (number of individuals, founders, sex distribution, etc.).

## Bug fixes
1. Fixed an issue where `tidyped(..., inbreed=TRUE)` failed due to incorrect class assignment order.
2. Fixed `visped(..., showf=TRUE)` to gracefully handle missing `f` columns by warning the user instead of erroring.
3. Fixed broken internal navigation links in package vignettes.

# Changes in version 0.6.1 released on 30 Dec 2025
## New features
1. Implemented opaque highlighting effects for better visualization clarity.
2. Added `trace` option to `visped()` to control ancestry tracing direction.

# Changes in version 0.6.0 released on 28 Dec 2025
## New features
1. Implemented strict S3 class structure for `tidyped` objects with `new_tidyped()` constructor and `validate_tidyped()` validator to ensure data integrity.

# Changes in version 0.5.0 released on 26 Dec 2025
## New features
1. Added `highlight` parameter to `visped()` function. Users can now highlight specific individuals using a character vector of IDs or a list for custom colors.
2. Added `showf` parameter to `visped()` function to display inbreeding coefficients on the pedigree graph.
3. Added `inbreed` parameter to `tidyped()` function to calculate inbreeding coefficients using the `nadiv` package.
4. Refactored `inbreed()` function as a standalone tool that operates on `tidyped` objects.
5. Optimized `repeloverlap()` function using `data.table` for significantly better performance.

## Bug fixes
1. Fixed a critical crash in `visped()` when combining `compact = TRUE`, `highlight`, and `showf = TRUE` by refactoring `ped2igraph()` to delay label modification until after layout calculation.
2. Fixed documentation grammar and phrasing across all functions for CRAN compliance.
3. Fixed `R CMD check` notes related to `data.table` non-standard evaluation by adding `R/globals.R`.

# Changes in version 0.4.1 released on 25 Dec 2025
## Bug fixes
1. Fixed overlapping edge detection for small pedigree graphs.
2. Improved coloring consistency for compact mode.

# Changes in version 0.2.6 released on 31 Mar 2020
## Bug fixes
1. Fixed a bug that the number of generations for candidates would be traced to n+1 when tracegen=n. This bug is found by Mianyu Liu.

# Changes in version 0.2.5 released on 25 Feb 2020
## Bug fixes
1. The tidyped() does not work with trace='all' in [certain cases](https://github.com/luansheng/visPedigree/issues/2#issue-568599008)

# Changes in version 0.2.4.1 released on 24 Feb 2020
## Bug fixes
1. An unexpected column with the name as NA occured when a tidyped object is tidyed again using the tidyped()

# Changes in version 0.2.4 released on 12 June 2019
## Bug fixes
1. The data.table used as the input parameter 'ped' may be changed in tidyped() and visped().

# Changes in version 0.2.3 released on 05 Mar 2019
## Bug fixes
1. The generation number of individuals is not inferred rightly.

# Changes in version 0.2.2 released on 28 Jan 2019
## Bug fixes
1. The tidied pedigree will not include the candidates which are not in the Ind column of the origin pedigree when the cand parameter is not NULL.

# Changes in version 0.2.1 released on 17 Nov 2018
## Bug fixes
1. Repel the overlapping nodes due to very small differences (digits > 7) among x positions of nodes

