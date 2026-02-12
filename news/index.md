# Changelog

## Changes in version 1.0.1 released on 31 Jan 2026

### Bug fixes

1.  **Compact Matrix Correctness**: Fixed a critical data integrity bug
    in `compact = TRUE` mode where relationship values (A, D, AA) were
    incorrect for parent-offspring and avuncular pairs due to improper
    merging of parent individuals with their non-parent siblings.
2.  **Pedigree Compression Strategy**: Updated compaction logic to
    preserve original genetic identity of any individual that appears as
    a sire or dam, ensuring parents always have unique entries in the
    relationship matrix.
3.  **Sibling Row/Column Expansion**: Fixed
    [`expand_pedmat()`](https://luansheng.github.io/visPedigree/reference/expand_pedmat.md)
    to correctly handle sibling off-diagonal elements by dynamically
    calculating relationship values based on parent kinship, rather than
    simply duplicating representative diagonal values.
4.  **Generation Alignment Logic**: Fixed
    `tidyped(..., genmethod = "bottom")` to prioritize **Sibling
    Consistency** (P1) over **Mate Alignment** (P2). This ensures that
    full siblings are always aligned to the same generation.
5.  **[`visped()`](https://luansheng.github.io/visPedigree/reference/visped.md)
    edge highlighting**: Fixed edge highlighting logic so relationship
    edges are only emphasized when `trace` is used.
6.  **Shared-parent/shared-child paths**: Corrected edge highlighting
    for cases where a parent has multiple families or a family has
    multiple children.
7.  **[`visped()`](https://luansheng.github.io/visPedigree/reference/visped.md)
    layout**: Fixed layout optimization failure when `showf = TRUE`. The
    layout algorithm now correctly uses immutable individual IDs.

## Changes in version 1.0.0 released on 24 Jan 2026

### API Standardization (BREAKING)

To provide a clean and intuitive API for v1.0.0, core function names and
behaviors have been standardized: - **`pedmatrix`** is renamed to
**`pedmat`**. - **`pedmat` default `method` is now `"A"`** (Additive
Relationship Matrix). Previously it was `"f"` (Inbreeding
Coefficients). - **`expand_pedmatrix`** is renamed to
**`expand_pedmat`**. - **`summary_pedmatrix`** is renamed to
**`summary_pedmat`**. - The parameter **`n_threads`** is standardized to
**`threads`** across all functions. - Legacy function names
(`pedmatrix`, etc.) are preserved as deprecated wrappers to ensure
backward compatibility.

### New Features

1.  **Family Assignment and Summary**:
    - [`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
      now automatically assigns and includes a `Family` column,
      identifying full-sib groups.
    - [`summary.tidyped()`](https://luansheng.github.io/visPedigree/reference/summary.tidyped.md)
      has been updated to provide family statistics (count, sizes, top
      largest families) and richer offspring analysis.
2.  **Pedigree Splitting (`splitped`)**: Added
    [`splitped()`](https://luansheng.github.io/visPedigree/reference/splitped.md)
    function to detect and split disconnected pedigree components. It
    efficiently identifies independent sub-populations (connected
    components) using graph theory, excludes isolated individuals, and
    returns a list of re-indexed `tidyped` objects ready for separate
    analysis or visualization.
3.  **Comprehensive Matrix Support**:
    [`pedmat()`](https://luansheng.github.io/visPedigree/reference/pedmat.md)
    (formerly `pedmatrix`) now fully supports 6 types of genetic
    relationship matrices: Additive (A, Ainv), Dominance (D, Dinv), and
    Additive-by-Additive Epistatic (AA, AAinv).
4.  **Relationship Matrix Visualization (`vismat`)**: Added
    [`vismat()`](https://luansheng.github.io/visPedigree/reference/vismat.md)
    function for visualizing relationship matrices (A, D, AA, etc.) with
    heatmaps and histograms. It supports `pedmat` objects, `tidyped`
    objects (auto-calculates A matrix), and standard matrices. Heatmaps
    can be annotated with family groups when a pedigree is provided.

### CRAN Submission & Internal Improvements

This release marks the first stable version 1.0.0, polished for CRAN.

1.  **Portable Compilation**: Standardized `src/Makevars` for
    cross-platform compatibility (removed GNU/platform-specific
    extensions).
2.  **Dependencies**: Moved `RcppArmadillo` to `LinkingTo` to optimize
    package structure.
3.  **Documentation & S3**: Fixed `vignette` generation, resolved `diag`
    S3 method dispatch, and cleaned up documentation for CRAN
    compliance.

## Changes in version 0.7.3 released on 13 Jan 2026

### New behavior (BREAKING)

1.  **Simplified `pedmatrix()` return and single-method enforcement**:
    `pedmatrix()` now requires a single `method` argument (e.g.,
    `method = "A").` When a single method is requested, the function
    returns the corresponding matrix or vector directly (not a named
    list). Requesting multiple methods in one call will now raise an
    error. Use repeated calls for multiple outputs.

### New features

1.  **High-Performance Genetic Relationship Calculations**: Introduced
    `pedmatrix()` function implemented in Rcpp for efficient computation
    of:
    - Additive relationship matrix (A) using the tabular recursive
      algorithm.
    - Sparse inverse additive matrix (A-Inverse) using Henderson’s
      rules.
    - Dominance matrix (D) using the tabular approach.
    - Inbreeding coefficients (f) using the Meuwissen & Luo (1992)
      path-tracing algorithm.

### Improvements

1.  **Default Inbreeding Calculation Method**: The
    [`inbreed()`](https://luansheng.github.io/visPedigree/reference/inbreed.md)
    function now uses the native Rcpp implementation by default, moving
    the `nadiv` package to `Suggests`.
2.  **Documentation and Website**: Updated package documentation and
    vignettes to reflect new features and improvements. The official
    package website is available at
    <https://luansheng.github.io/visPedigree/>.

## Changes in version 0.7.2 released on 12 Jan 2026

### New features

1.  **Flexible Generation Assignment**: Added `genmethod` parameter to
    [`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md).
    Users can now choose between `"top"` (top-aligned, default) and
    `"bottom"` (bottom-aligned) methods for generation inference.
    - The `"top"` method aligns founders at Generation 1, which is more
      appropriate for biological pedigrees and prevents “founder drift”
      in pedigrees with varying depths.
    - The `"bottom"` method aligns terminal nodes at the bottom, useful
      for visualizing introductions of unrelated exogenous parents.

### Improvements

1.  **Default Logic Change**: Switched the default generation assignment
    method to `"top"` (top-down) for more intuitive biological
    visualization.
2.  **Pkgdown Documentation**: Generated and published the official
    package website at <https://luansheng.github.io/visPedigree/>.
3.  **Automated CI/CD**: Added GitHub Actions workflow for automatic
    documentation updates and deployment via GitHub Pages.

## Changes in version 0.7.1 released on 11 Jan 2026

CRAN release: 2026-01-21

### Performance optimizations

1.  **Large Pedigree Performance**: Optimized `visped` performance for
    displaying large pedigrees through efficient attribute handling and
    vectorized rendering. Computation time for 100k+ individuals reduced
    significantly by avoiding redundant `igraph` attribute lookups.
2.  **Vectorized Tracing**: Refactored `trace_ped_candidates` in
    `tidyped` to use vectorized
    [`igraph::neighborhood`](https://r.igraph.org/reference/ego.html)
    calls, achieving ~150x speedup for large candidate lists (e.g., 37k
    candidates in a 178k individual pedigree traced in ~1.2s).
3.  **Early Filtering**: Implemented unified early filtering of isolated
    individuals (Gen 0) in `prepare_ped_graph` to streamline downstream
    graph conversion and layout algorithms.

### Improvements

1.  **User Feedback**: Standardized filtering notifications. The message
    “Note: Removed N isolated individuals…” now appears consistently for
    all pedigree sizes when Gen 0 individuals are present.
2.  **Refined Tracing**: Corrected `trace = "all"` logic in both
    `tidyped` and `visped`. It now correctly retrieves the union of
    ancestors and descendants (“up” + “down”) instead of the entire
    connected component (undirected search).

## Changes in version 0.7.0 released on 10 Jan 2026

### Breaking changes & Major Refactoring

1.  **Graph-based `tidyped` Core**: Reimplemented the pedigree tidying
    engine using formal graph theory principles (Directed Acyclic
    Graphs). Improved loop detection and generation inference accuracy
    using topological sorting.
2.  **Modular Architecture**: Split the monolithic `visped.R` into
    functional modules: `visped_layout.R`, `visped_graph.R`,
    `visped_style.R`, and `visped_render.R` for better maintainability.

### New features

1.  **New Parameters in
    [`visped()`](https://luansheng.github.io/visPedigree/reference/visped.md)**:
    - `pagewidth`: Allows users to specify the PDF page width (default
      200 inches) to accommodate different pedigree scales.
    - `symbolsize`: A scaling factor (default 1) to adjust node sizes
      relative to label dimensions, providing finer control over
      whitespace.
2.  **Two-Pass Rendering Engine**: Introduced a two-pass strategy in
    [`plot_ped_igraph()`](https://luansheng.github.io/visPedigree/reference/plot_ped_igraph.md)
    to ensure edges connect exactly at node centers, eliminating visual
    gaps in vector PDF outputs.
3.  **Enhanced Highlighting**: Added support for real-time ancestry and
    descendant highlighting via the `trace` parameter in
    [`visped()`](https://luansheng.github.io/visPedigree/reference/visped.md).

### Bug fixes

1.  Fixed rendering failure in `outline = TRUE` mode by correcting
    attribute indexing in the graph object.
2.  modernized the unit testing suite to `testthat` 3rd edition,
    removing all legacy `context()` warnings.
3.  Improved coordinate calculation precision to prevent overlapping in
    high-density generations.

## Changes in version 0.6.2 released on 01 Jan 2026

### New features

1.  Added [`summary()`](https://rdrr.io/r/base/summary.html) method for
    `tidyped` objects to provide quick pedigree statistics (number of
    individuals, founders, sex distribution, etc.).

### Bug fixes

1.  Fixed an issue where `tidyped(..., inbreed=TRUE)` failed due to
    incorrect class assignment order.
2.  Fixed `visped(..., showf=TRUE)` to gracefully handle missing `f`
    columns by warning the user instead of erroring.
3.  Fixed broken internal navigation links in package vignettes.

## Changes in version 0.6.1 released on 30 Dec 2025

### New features

1.  Implemented opaque highlighting effects for better visualization
    clarity.
2.  Added `trace` option to
    [`visped()`](https://luansheng.github.io/visPedigree/reference/visped.md)
    to control ancestry tracing direction.

## Changes in version 0.6.0 released on 28 Dec 2025

### New features

1.  Implemented strict S3 class structure for `tidyped` objects with
    [`new_tidyped()`](https://luansheng.github.io/visPedigree/reference/new_tidyped.md)
    constructor and
    [`validate_tidyped()`](https://luansheng.github.io/visPedigree/reference/validate_tidyped.md)
    validator to ensure data integrity.

## Changes in version 0.5.0 released on 26 Dec 2025

### New features

1.  Added `highlight` parameter to
    [`visped()`](https://luansheng.github.io/visPedigree/reference/visped.md)
    function. Users can now highlight specific individuals using a
    character vector of IDs or a list for custom colors.
2.  Added `showf` parameter to
    [`visped()`](https://luansheng.github.io/visPedigree/reference/visped.md)
    function to display inbreeding coefficients on the pedigree graph.
3.  Added `inbreed` parameter to
    [`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
    function to calculate inbreeding coefficients using the `nadiv`
    package.
4.  Refactored
    [`inbreed()`](https://luansheng.github.io/visPedigree/reference/inbreed.md)
    function as a standalone tool that operates on `tidyped` objects.
5.  Optimized
    [`repeloverlap()`](https://luansheng.github.io/visPedigree/reference/repeloverlap.md)
    function using `data.table` for significantly better performance.

### Bug fixes

1.  Fixed a critical crash in
    [`visped()`](https://luansheng.github.io/visPedigree/reference/visped.md)
    when combining `compact = TRUE`, `highlight`, and `showf = TRUE` by
    refactoring
    [`ped2igraph()`](https://luansheng.github.io/visPedigree/reference/ped2igraph.md)
    to delay label modification until after layout calculation.
2.  Fixed documentation grammar and phrasing across all functions for
    CRAN compliance.
3.  Fixed `R CMD check` notes related to `data.table` non-standard
    evaluation by adding `R/globals.R`.

## Changes in version 0.4.1 released on 25 Dec 2025

## Changes in version 0.2.6 released on 31 Mar 2020

### New features

### Bug fixes

1.  Fixed a bug that the number of generations for candidates would be
    traced to n+1 when tracegen=n. This bug is found by Mianyu Liu.

## Changes in version 0.2.5 released on 25 Feb 2020

### New features

### Bug fixes

1.  The tidyped() does not work with trace=‘all’ in [certain
    cases](https://github.com/luansheng/visPedigree/issues/2#issue-568599008)

## Changes in version 0.2.4.1 released on 24 Feb 2020

### New features

### Bug fixes

1.  An unexpected column with the name as NA occured when a tidyped
    object is tidyed again using the tidyped()

## Changes in version 0.2.4 released on 12 June 2019

### New features

### Bug fixes

1.  The data.table used as the input parameter ‘ped’ may be changed in
    tidyped() and visped().

## Changes in version 0.2.3 released on 05 Mar 2019

### New features

### Bug fixes

1.  The generation number of individuals is not inferred rightly.

## Changes in version 0.2.2 released on 28 Jan 2019

### New features

### Bug fixes

1.  The tidied pedigree will not include the candidates which are not in
    the Ind column of the origin pedigree when the cand parameter is not
    NULL.

## Changes in version 0.2.1 released on 17 Nov 2018

### New features

### Bug fixes

1.  Repel the overlapping nodes due to very small differences (digits
    \> 7) among x positions of nodes
