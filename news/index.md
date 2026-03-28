# Changelog

## Changes in version 1.8.1 released on 28 Mar 2026

### New features

1.  **[`visped()`](https://luansheng.github.io/visPedigree/reference/visped.md)
    generation labels**: Added a new `genlab` argument to
    [`visped()`](https://luansheng.github.io/visPedigree/reference/visped.md)
    for drawing generation labels (`G1`, `G2`, …) on the left margin of
    pedigree plots. This helps identify each row in deep pedigrees. The
    default is `FALSE`, so existing plots are unchanged unless
    `genlab = TRUE` is requested.

### Documentation

1.  **`draw-pedigree` vignette updated**: Added a new example showing
    how to display generation labels with `visped(..., genlab = TRUE)`.

## Changes in version 1.8.0 released on 25 Mar 2026

### Bug Fixes

1.  **[`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
    fast path with `addnum = FALSE` + `cand`**: When the input `tidyped`
    object was created with `addnum = FALSE` (no `IndNum`/`SireNum`/
    `DamNum` columns), passing a `cand` argument would always raise
    `"None of the specified candidates were found in the pedigree."`
    because the fast-path BFS looked up `ped_dt$IndNum` which was
    `NULL`. The fix temporarily adds integer index columns for the BFS
    and removes them from the output when `addnum = FALSE`.

2.  **[`pediv()`](https://luansheng.github.io/visPedigree/reference/pediv.md)
    /
    [`pedne()`](https://luansheng.github.io/visPedigree/reference/pedne.md)
    incorrect `fg`, `MeanCoan`, and `NeCoancestry` values when `ECG` was
    not pre-computed**: Inside
    [`pediv()`](https://luansheng.github.io/visPedigree/reference/pediv.md)
    and
    [`pedne()`](https://luansheng.github.io/visPedigree/reference/pedne.md),
    a `merge(..., by = "Ind", all.x = TRUE)` was used to attach `ECG`
    values.
    [`data.table::merge()`](https://rdrr.io/pkg/data.table/man/merge.html)
    sorts output by the key column (`Ind`, alphabetically), breaking the
    `IndNum == row-index` invariant that the fast-path BFS in
    [`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
    depends on. The result was that the traced pedigree used for
    coancestry calculation was incorrect, producing wildly wrong `fg`
    and related values (e.g., `fg ≈ 240` instead of `fg ≈ 19`). Fixed by
    calling `setorder(ped_dt, IndNum)` immediately after the merge.

3.  **[`summary_pedmat()`](https://luansheng.github.io/visPedigree/reference/summary_pedmat.md)
    reported `Density = 100%` for all `A`, `D`, and `AA` matrices**: The
    speed branch changed the return type of `pedmat(sparse = TRUE)` for
    these methods from `dgCMatrix` (sparse) to `dgeMatrix` (dense) to
    avoid the O(N²) zero-scan performed by
    `Matrix::Matrix(..., sparse = TRUE)`. However,
    [`summary_pedmat()`](https://luansheng.github.io/visPedigree/reference/summary_pedmat.md)
    contained a hard-coded branch that returned `sparsity = 1.0` for any
    non-`sparseMatrix` `Matrix` subclass (including `dgeMatrix`), always
    printing `Density = 100%` instead of the true fill ratio. Fixed by
    using `Matrix::nnzero() / (nrow * ncol)` for all `Matrix`
    subclasses; `nnzero()` works correctly for both `dgeMatrix` and
    `dgCMatrix`.

### Performance

1.  **~100× faster inbreeding calculation**: Replaced the Meuwissen &
    Luo (1992) linear path-trace algorithm in `cpp_calculate_inbreeding`
    with the Sargolzaei & Iwaisaki (2005) LAP (Longest Ancestral Path)
    bucket method. At N = 1,000,000, the C++ kernel now completes in
    ~0.15 s (previously ~15 s), and the full
    [`inbreed()`](https://luansheng.github.io/visPedigree/reference/inbreed.md)
    call returns in ~0.40 s. The key improvements are O(1) ancestor
    retrieval via bucket pop (vs. O(gap) linear scan), O(1) duplicate
    suppression via `L[k] == 0` check, and O(m_i) path-coefficient reset
    (vs. O(k) full-array scan). All results are numerically identical to
    the previous implementation (max difference \< 2 × 10⁻¹⁵).

2.  **~16× faster
    [`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
    candidate tracing**: When the input is already a `tidyped` object
    and `cand` is specified, the fast path now uses three new pure C++
    BFS functions (`cpp_trace_ancestors`, `cpp_trace_descendants`,
    `cpp_topo_order`) instead of rebuilding an igraph object. At N =
    1,000,000 with 200 candidates, elapsed time drops from ~0.91 s to
    ~0.056 s. The fast path now also correctly applies sibling and mate
    generation alignment for `genmethod = "bottom"`, matching the output
    of the full path.

3.  **~9× faster
    [`pedgenint()`](https://luansheng.github.io/visPedigree/reference/pedgenint.md)
    generation-interval lookup**: Replaced the character-key
    `data.table` join used to look up each individual’s generation
    number with a pre-computed integer index array. Benchmark on a
    representative dataset: 4.24 s → 0.49 s.

4.  **~2× faster `pedmat(sparse = TRUE)` matrix coercion**: The
    conversion from a dense numeric matrix to a `Matrix::dgeMatrix` now
    bypasses the S4 dispatch overhead of `as(mat, "dgeMatrix")` via a
    direct [`methods::new()`](https://rdrr.io/r/methods/new.html) call
    backed by a package-level class cache. The cached class object is
    looked up once per session instead of on every call. Benchmark: 1.24
    s → 0.60 s for a representative subpedigree.

5.  **`pedrel(compact = TRUE)` 200 k safety guard**: Added an early
    error when the number of reference individuals exceeds 200,000 in
    compact mode, preventing inadvertent multi-hour runs from silent
    O(N²) matrix allocation.

### Internal changes

1.  **`align_bottom_generations()` extracted as shared helper**: The
    ~60-line sibling/mate generation-alignment block previously
    duplicated between the main path and the fast path of
    [`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
    has been consolidated into a single internal helper. Both paths now
    call it consistently.

2.  **C++11 compatibility**: Two C++17 structured-binding usages
    (`auto [cur, d] = q.front()`) in the BFS functions have been
    replaced with explicit C++11 equivalents, ensuring compatibility
    with GCC 8 (CRAN’s minimum required compiler).

3.  **`methods` declared in `Imports`**: The `methods` package is now
    explicitly listed under `Imports` in `DESCRIPTION`, as required by
    CRAN policy when
    [`methods::getClass()`](https://rdrr.io/r/methods/getClass.html) and
    [`methods::new()`](https://rdrr.io/r/methods/new.html) are called at
    runtime.

### Documentation

1.  Updated
    [`inbreed()`](https://luansheng.github.io/visPedigree/reference/inbreed.md)
    and vignette references to cite Sargolzaei & Iwaisaki (2005) instead
    of Meuwissen & Luo (1992).

## Changes in version 1.7.0 released on 23 Mar 2026

### New features

1.  **[`pediv()`](https://luansheng.github.io/visPedigree/reference/pediv.md)
    retained genetic diversity (`GeneDiv`)**: `pediv()$summary` gains a
    new column `GeneDiv = 1 - MeanCoan`, the pedigree-based retained
    genetic diversity of the reference population. Values lie in
    $\lbrack 0,1\rbrack$; higher values indicate more diversity retained
    relative to an unrelated base population.
    [`print.pediv()`](https://luansheng.github.io/visPedigree/reference/print.pediv.md)
    displays it alongside `fg` and `MeanCoan`. This dimensionless
    complement to `fg` is easier to communicate to non-specialist
    stakeholders.

2.  **[`vismat()`](https://luansheng.github.io/visPedigree/reference/vismat.md)
    large-pedigree representative view**: When the original pedigree has
    more than 5,000 individuals,
    [`vismat()`](https://luansheng.github.io/visPedigree/reference/vismat.md)
    no longer attempts a full N × N matrix expansion. Instead, it uses
    the compact matrix (K × K representative individuals) directly and
    adds sibling-count labels of the form `ID (×n)` to each axis tick.
    When `compact = TRUE` and `by` is supplied, the function now
    computes group means algebraically from the K × K matrix without
    expanding to N × N, eliminating memory overflow for very large
    pedigrees (e.g., N ≈ 178,000 in `big_family_size_ped`).

### Internal changes

1.  **[`vismat()`](https://luansheng.github.io/visPedigree/reference/vismat.md)
    named threshold constants**: Hardcoded magic numbers controlling
    large-pedigree behavior (`5000`, `2000`, `50`) have been refactored
    into named constants (`VISMAT_EXPAND_MAX`, `VISMAT_LABEL_MAX`,
    `VISMAT_WARN_THRESHOLD`) at the top of `R/vismat.R` for easier
    maintenance.

## Changes in version 1.6.2 released on 23 Mar 2026

### New features

1.  **[`pedrel()`](https://luansheng.github.io/visPedigree/reference/pedrel.md)
    coancestry scale**: Added a `scale` parameter to
    [`pedrel()`](https://luansheng.github.io/visPedigree/reference/pedrel.md)
    supporting `"relationship"` (default, returns mean $a_{ij}$) and
    `"coancestry"` (returns corrected mean coancestry $\bar{c}$). The
    coancestry scale uses the diagonal-corrected formula of Caballero &
    Toro (2000), properly accounting for self-coancestry within the
    reference group.

### API Changes

1.  **[`vispstat()`](https://luansheng.github.io/visPedigree/reference/vispstat.md)
    internal downgrade**: The
    [`vispstat()`](https://luansheng.github.io/visPedigree/reference/vispstat.md)
    function has been downgraded to an internal-only function.
    Specifically, it is now the internal backend for
    [`plot.pedstats()`](https://luansheng.github.io/visPedigree/reference/vispstat.md).
    Users should use the standard `plot(stats_obj)` S3 method instead.
2.  **Standardized
    [`vispstat()`](https://luansheng.github.io/visPedigree/reference/vispstat.md)
    documentation**: Updated the documentation for generation intervals
    in
    [`vispstat()`](https://luansheng.github.io/visPedigree/reference/vispstat.md)
    to accurately reflect the visualization of mean values (removing the
    misleading “± SD” claim).

### Documentation

1.  **`pedigree-analysis.Rmd` §9 expanded**: Section 9 (“Average
    Relationship Trends with
    [`pedrel()`](https://luansheng.github.io/visPedigree/reference/pedrel.md)”)
    split into two sub-sections to cover both `scale` options. Added
    §9.2 with the Caballero & Toro (2000) diagonal-corrected coancestry
    formula, a worked example using `scale = "coancestry"` (returning
    `MeanCoan`), and guidance on when to prefer each scale.
2.  **`relationship-matrix.Rmd` updated**: Added §3.2
    compact-to-[`vismat()`](https://luansheng.github.io/visPedigree/reference/vismat.md)
    direct path and expanded §4.1 with five sub-sections
    (`reorder = FALSE`, `ids`, `by = "Gen"`, `by = "Family"`, compact
    auto-expand). Replaced incorrect performance thresholds in §5 with
    an accurate reference table.

### Bug fixes

1.  **Eliminated spurious subsetting warnings**: Internal group-by
    slicing in
    [`pedrel()`](https://luansheng.github.io/visPedigree/reference/pedrel.md),
    [`pedne()`](https://luansheng.github.io/visPedigree/reference/pedne.md),
    [`pediv()`](https://luansheng.github.io/visPedigree/reference/pediv.md),
    [`pedhalflife()`](https://luansheng.github.io/visPedigree/reference/pedhalflife.md),
    and
    [`pedgenint()`](https://luansheng.github.io/visPedigree/reference/pedgenint.md)
    previously triggered false-positive `[.tidyped]` warnings
    (“Subsetting removed parent records…”). Fixed by using
    [`as.data.table()`](https://rdrr.io/pkg/data.table/man/as.data.table.html)
    to bypass the completeness guard when the subset is only used for ID
    extraction, not for pedigree computation.
2.  **Eliminated internal class-restoration messages**:
    `pedrel(compact = TRUE)` previously emitted “Note: ‘ped’ lost its
    tidyped class … Restoring automatically.” messages. Root cause was
    [`compact_ped_for_matrix()`](https://luansheng.github.io/visPedigree/reference/compact_ped_for_matrix.md)
    returning a plain `data.table` in its early-return branches (no
    full-sib families or no compactable members). Fixed by returning
    `data.table::copy(ped)` to preserve the `tidyped` class.
3.  **Internal cleanup**: Removed dead code variables and added
    comprehensive unit tests for the `genint` branch of
    [`vispstat()`](https://luansheng.github.io/visPedigree/reference/vispstat.md).

## Changes in version 1.6.1 released on 21 Mar 2026

### Improvements

1.  **Standardized diversity notation**: Output columns in
    `pedhalflife()$timeseries` are now lowercase (`fe`, `fa`, `fg`,
    `lnfe`, `lnfa`, `lnfg`, `lnfafe`, `lnfgfa`) to maintain consistency
    with standard population genetics nomenclature.
2.  **Enhanced
    [`plot.pedhalflife()`](https://luansheng.github.io/visPedigree/reference/pedhalflife.md)
    visualization**: In `type = "log"` mode, the plot now includes an
    OLS regression line for total diversity decay
    ($\ln f_{g} \sim \text{Time}$) and a vertical reference line for the
    diversity half-life $T_{1/2}$.
3.  **Improved unit labeling**:
    [`plot.pedhalflife()`](https://luansheng.github.io/visPedigree/reference/pedhalflife.md)
    and
    [`print.pedhalflife()`](https://luansheng.github.io/visPedigree/reference/pedhalflife.md)
    now automatically use the name of the `timevar` column (e.g., “Gen”,
    “Year”) for axis and summary labels.

## Changes in version 1.6.0 released on 20 Mar 2026

### New features

1.  **Information-theoretic diversity half-life
    ([`pedhalflife()`](https://luansheng.github.io/visPedigree/reference/pedhalflife.md))**:
    New function that tracks $f_{e}$, $f_{a}$, and $f_{g}$ across time
    points and fits a log-linear decay model to quantify the rate of
    genetic diversity loss. The total loss rate $\lambda_{total}$ is
    decomposed into three additive components: foundation bottleneck
    ($\lambda_{e}$, unequal founder contributions), breeding bottleneck
    ($\lambda_{b}$, overuse of key ancestors), and genetic drift
    ($\lambda_{d}$, random sampling loss). The diversity half-life
    $T_{1/2} = \ln 2/\lambda_{total}$ is reported in units of the
    `timevar` column (e.g., generations, years). S3 `print` and `plot`
    methods are provided. The
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
    supports both log-scale (`type = "log"`) and raw-scale
    (`type = "raw"`) views of the decay trajectory.

## Changes in version 1.5.0 released on 20 Mar 2026

### New features

1.  **Shannon-entropy effective founders and ancestors (`feH`, `faH`)**:
    [`pedcontrib()`](https://luansheng.github.io/visPedigree/reference/pedcontrib.md)
    and
    [`pediv()`](https://luansheng.github.io/visPedigree/reference/pediv.md)
    now compute two additional diversity statistics based on the Hill
    number of order $q = 1$ (Shannon entropy):
    - `feH` — effective number of founders under equal entropy
      weighting. Satisfies the inequality
      $N_{Founder} \geq f_{e}^{H} \geq f_{e}$.
    - `faH` — effective number of ancestors under equal entropy
      weighting. Satisfies the inequality
      $N_{Ancestor} \geq f_{a}^{H} \geq f_{a}$. Both are computed from
      the vector of genetic contributions using the formula
      $\exp\left( - \sum p_{i}\ln p_{i} \right)$ and complement the
      classical quadratic ($q = 2$) effective numbers $f_{e}$ and
      $f_{a}$ (Lacy 1989; Boichard et al. 1997). The new columns appear
      in `pediv()$summary` and `pedcontrib()$summary` alongside the
      existing `fe`, `fa`, and `fg` columns.

## Changes in version 1.4.1 released on 15 Mar 2026

### Bug fixes

1.  **Fail-fast incomplete pedigree analysis**:
    [`inbreed()`](https://luansheng.github.io/visPedigree/reference/inbreed.md)
    and other completeness-sensitive analysis functions now error on
    row-truncated subsets with missing parent records. This prevents
    incorrect results (e.g., zero inbreeding) caused by calculating on
    partial ancestry data.

## Changes in version 1.4.0 released on 15 Mar 2026

### New features

1.  **`tidyped` class redesign and optimization**: Refined the internal
    `tidyped` class architecture around a clearer metadata contract and
    safer S3/data.table interaction model, making the object more robust
    for repeated downstream analysis and extension.
2.  **Safer `tidyped` object workflows**: Added
    [`is_tidyped()`](https://luansheng.github.io/visPedigree/reference/is_tidyped.md),
    [`pedmeta()`](https://luansheng.github.io/visPedigree/reference/pedmeta.md),
    [`has_inbreeding()`](https://luansheng.github.io/visPedigree/reference/has_inbreeding.md),
    and
    [`has_candidates()`](https://luansheng.github.io/visPedigree/reference/has_candidates.md)
    to make class checks and metadata inspection explicit and
    user-facing.
3.  **Fast candidate tracing from existing `tidyped` objects**:
    [`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
    now uses a fast path when the input is already a valid `tidyped`
    object and `cand` is supplied, avoiding repeated global validation
    and preprocessing.
4.  **Workflow coverage and developer documentation**: Added a new
    workflow vignette, a `tidyped` structure and extension vignette, and
    focused regression tests covering safe subsetting, `:=` by-reference
    behavior, and split workflow semantics.

### Bug fixes

1.  **Stable by-reference mutation for `tidyped`**: Replaced class and
    metadata attachment paths with
    [`data.table::setattr()`](https://rdrr.io/pkg/data.table/man/setattr.html)
    so subsequent `:=` operations keep true by-reference behavior
    instead of silently writing into shallow copies.
2.  **Safe row subsetting**: Added `[.tidyped` interception so
    incomplete subsets degrade to plain `data.table` objects with a
    warning, while complete subsets preserve `tidyped` structure and
    rebuild pedigree indices correctly.
3.  **More robust class recovery**: Core analysis entry points now
    cooperate with
    [`ensure_tidyped()`](https://luansheng.github.io/visPedigree/reference/ensure_tidyped.md)
    /
    [`validate_tidyped()`](https://luansheng.github.io/visPedigree/reference/validate_tidyped.md)
    to recover valid `tidyped` objects after common class-dropping
    operations.

### Documentation

1.  **Pkgdown article navigation**: Reorganized vignette order, restored
    `draw-pedigree` to the recommended reading sequence, and exposed
    `tidyped` developer notes through a dedicated pkgdown
    developer-documentation entry.

## Changes in version 1.3.5 released on 14 Mar 2026

### New features

1.  **S3 Class Protection**: Added
    [`as_tidyped()`](https://luansheng.github.io/visPedigree/reference/as_tidyped.md)
    and an internal
    [`ensure_tidyped()`](https://luansheng.github.io/visPedigree/reference/ensure_tidyped.md)
    mechanism to robustly handle the “silent class loss” bug. Standard R
    operations like [`merge()`](https://rdrr.io/r/base/merge.html),
    [`rbind()`](https://rdrr.io/r/base/cbind.html), and `dplyr` verbs
    often strip custom S3 classes from `data.table` objects. Major
    analysis functions now automatically detect if the `tidyped` class
    is missing and restore it if the underlying data structure is
    intact, providing an informational message to the user instead of
    erroring.

### Bug fixes

1.  **Robust Analysis Entry**: Updated all 11 core analysis functions
    (including
    [`pedstats()`](https://luansheng.github.io/visPedigree/reference/pedstats.md),
    [`pedne()`](https://luansheng.github.io/visPedigree/reference/pedne.md),
    [`pediv()`](https://luansheng.github.io/visPedigree/reference/pediv.md),
    [`pedrel()`](https://luansheng.github.io/visPedigree/reference/pedrel.md),
    etc.) to use the new auto-recovery logic. This ensures that analysis
    remains user-friendly and reliable even after manual data
    manipulation by the user.

## Changes in version 1.3.4 released on 14 Mar 2026

### Bug fixes

1.  **`data.table` invisibility**: Fixed a subtle but pervasive issue
    where functions returning `data.table` or `tidyped` objects (which
    are based on `data.table`) were returning them invisibly. This
    occurred because internal `data.table` operations like `:=` and
    `set*` set an internal “invisible” flag. Affected functions included
    [`pedancestry()`](https://luansheng.github.io/visPedigree/reference/pedancestry.md),
    [`pedpartial()`](https://luansheng.github.io/visPedigree/reference/pedpartial.md),
    [`pedne()`](https://luansheng.github.io/visPedigree/reference/pedne.md),
    [`pedrel()`](https://luansheng.github.io/visPedigree/reference/pedrel.md),
    [`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md),
    and many others. All relevant functions now explicitly return the
    object using the `[]` syntax to ensure they auto-print correctly in
    the R console and knitted documents.
2.  **Side-effect prevention**: Updated `calc_ne_demographic()` to
    operate on a copy of the input pedigree instead of modifying the
    user’s data by reference.

### Documentation

1.  **Positron Guide**: Added a new coding standard to `Positron.md`
    regarding `data.table` return visibility to prevent regressive
    “invisible output” bugs in future development.

## Changes in version 1.3.3 released on 14 Mar 2026

### Documentation

1.  **`pedigree-analysis.Rmd` rewrite**: Reorganized the main pedigree
    analysis vignette into clearer thematic sections covering pedigree
    overview, pedigree completeness
    ([`pedecg()`](https://luansheng.github.io/visPedigree/reference/pedecg.md)),
    generation intervals
    ([`pedgenint()`](https://luansheng.github.io/visPedigree/reference/pedgenint.md)),
    subpopulation structure
    ([`pedsubpop()`](https://luansheng.github.io/visPedigree/reference/pedsubpop.md)),
    diversity indicators
    ([`pediv()`](https://luansheng.github.io/visPedigree/reference/pediv.md)),
    effective population size
    ([`pedne()`](https://luansheng.github.io/visPedigree/reference/pedne.md)),
    average relationship trends
    ([`pedrel()`](https://luansheng.github.io/visPedigree/reference/pedrel.md)),
    inbreeding classification
    ([`pedfclass()`](https://luansheng.github.io/visPedigree/reference/pedfclass.md)),
    and ancestry / partial inbreeding diagnostics.
2.  **Theory expansion**: Added core formulas, interpretation notes, and
    breeding-use explanations for Equivalent Complete Generations (ECG),
    generation intervals, effective numbers of founders / ancestors /
    founder genomes (`f_e`, `f_a`, `f_g`), three effective population
    size definitions (`N_e` by demographic, inbreeding, and coancestry
    methods), and average additive relationship (`MeanRel`).
3.  **Reference update**: Expanded the vignette bibliography to include
    the key classical references underlying the package’s diversity and
    effective population size metrics, including Wright (1922, 1931),
    Lacy (1989), Boichard et al. (1997), Caballero & Toro (2000),
    Cervantes et al. (2011), and Gutiérrez et al. (2008, 2009).

### Bug fixes

1.  **Vignette API synchronization**: Replaced outdated
    `pedinbreed_class()` calls in the pedigree analysis vignette with
    the current
    [`pedfclass()`](https://luansheng.github.io/visPedigree/reference/pedfclass.md)
    interface and aligned all examples with the current `reference`,
    `foundervar`, and `cycle` argument names.

### Testing

1.  **Analysis regression coverage**: Added focused unit tests to verify
    that
    [`pedancestry()`](https://luansheng.github.io/visPedigree/reference/pedancestry.md)
    proportions sum to 1 in a multi-line admixture pedigree and that
    [`pedrel()`](https://luansheng.github.io/visPedigree/reference/pedrel.md)
    returns identical results between `compact = TRUE` and
    `compact = FALSE` modes on the same pedigree.

## Changes in version 1.3.2 released on 13 Mar 2026

### New features

1.  **Added Comprehensive Examples**: Added `@examples` to core analysis
    functions
    [`pedne()`](https://luansheng.github.io/visPedigree/reference/pedne.md),
    [`pedecg()`](https://luansheng.github.io/visPedigree/reference/pedecg.md),
    and
    [`pedsubpop()`](https://luansheng.github.io/visPedigree/reference/pedsubpop.md)
    to improve documentation completeness and provide immediate value to
    users.
2.  **Pedigree Connectivity Analysis (`pedsubpop`)**: Enhanced
    [`pedsubpop()`](https://luansheng.github.io/visPedigree/reference/pedsubpop.md)
    to better distinguish between pedigree splitting (via `splitped`)
    and grouping/summary analysis. It now provides clear counts of total
    individuals, sires, dams, and founders within subgroups or connected
    components.

### API Changes

1.  **[`pedfclass()`](https://luansheng.github.io/visPedigree/reference/pedfclass.md)
    rename**: Renamed the inbreeding-class summary helper from
    `pedinbreedclass()` to
    [`pedfclass()`](https://luansheng.github.io/visPedigree/reference/pedfclass.md)
    to align with the package naming guide and provide a shorter,
    clearer user-facing API.
2.  **[`pedfclass()`](https://luansheng.github.io/visPedigree/reference/pedfclass.md)
    output refinement**: Renamed the returned class column from
    `F_Class` to `FClass`, and added support for user-defined inbreeding
    class breakpoints through the `breaks` and `labels` arguments.
3.  **[`pedgenint()`](https://luansheng.github.io/visPedigree/reference/pedgenint.md)
    parameter rename**: Renamed `cycle_length` to `cycle` for
    consistency with the package naming guide.
4.  **[`pedgenint()`](https://luansheng.github.io/visPedigree/reference/pedgenint.md)
    /
    [`pedstats()`](https://luansheng.github.io/visPedigree/reference/pedstats.md)
    `unit` parameter**: Removed `"gen"` from `unit` options. The `unit`
    parameter now only accepts `"year"`, `"month"`, `"day"`, or
    `"hour"`. The previous `"gen"` option produced incorrect results
    when combined with date inputs.
5.  **[`pedgenint()`](https://luansheng.github.io/visPedigree/reference/pedgenint.md)
    `timevar` definition**: Clarified `timevar` as a **birth date**
    column. Numeric year inputs (e.g., `2020`) are now automatically
    converted to `Date` (`"YYYY-07-01"`) with an informational message.
    Character date strings are parsed via
    [`as.POSIXct()`](https://rdrr.io/r/base/as.POSIXlt.html) with
    `tz = "UTC"` to avoid DST artifacts.

### Minor improvements and bug fixes

1.  **Documentation Audit**: Refined internal documentation for
    [`pedsubpop()`](https://luansheng.github.io/visPedigree/reference/pedsubpop.md)
    to clarify its use cases alongside
    [`splitped()`](https://luansheng.github.io/visPedigree/reference/splitped.md).
2.  **[`vispstat()`](https://luansheng.github.io/visPedigree/reference/vispstat.md)
    pathway filter fix**: Fixed an issue where the generation-interval
    bar chart could silently drop pathways due to an overly broad
    [`factor()`](https://rdrr.io/r/base/factor.html) filter. Now uses
    explicit `%in% c("SS", "SD", "DS", "DD")` subsetting.
3.  **`.parse_to_numeric_time()` rewrite**: Completely rewrote the
    internal time parser to handle `Date`, `POSIXct`, character date
    strings, and numeric years robustly. All `POSIXct` conversions now
    use `tz = "UTC"` to prevent DST-related artifacts.

## Changes in version 1.3.1 released on 12 Mar 2026

### New features

1.  Added the `selfing` argument to
    [`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
    to support plant and aquaculture pedigrees where an individual can
    appear as both Sire and Dam, resolving biologically impossible sex
    conflict errors
    ([\#10](https://github.com/luansheng/visPedigree/issues/10)).
2.  **[`pedrel()`](https://luansheng.github.io/visPedigree/reference/pedrel.md)
    logic upgrade**: Modified
    [`pedrel()`](https://luansheng.github.io/visPedigree/reference/pedrel.md)
    to use full ancestral tracing via `tidyped(ped, cand = ...)` when
    calculating sub-group relationships. This ensures that relationships
    in deep-inbred populations (e.g., full-sib mating over multiple
    generations) are calculated correctly rather than being
    underestimated due to ancestor truncation.

### API Changes

1.  **[`pedancestry()`](https://luansheng.github.io/visPedigree/reference/pedancestry.md)
    parameter rename**: Renamed `labelvar` to `foundervar` and `labels`
    to `target_labels` to align with the package naming guide and make
    the ancestry-tracing interface more explicit. Old argument names are
    no longer supported because this function is still under active
    development.
2.  **[`pedecg()`](https://luansheng.github.io/visPedigree/reference/pedecg.md)
    parameter cleanup**: Removed the short-lived `reference` argument.
    It only filtered rows after a full ECG pass and did not define a
    true reference population or prune the pedigree before calculation.
    Users should subset the returned table directly if needed.

### Minor improvements and bug fixes

1.  **Academic nomenclature alignment**: Updated documentation for
    [`pedrel()`](https://luansheng.github.io/visPedigree/reference/pedrel.md)
    and
    [`pedne()`](https://luansheng.github.io/visPedigree/reference/pedne.md)
    to explicitly distinguish between **Additive Genetic Relationship
    ($a_{ij}$)** and **Coancestry ($f_{ij}$)**.
    - [`pedrel()`](https://luansheng.github.io/visPedigree/reference/pedrel.md)
      now clearly states it returns $a_{ij} = 2f_{ij}$.
    - [`pedne()`](https://luansheng.github.io/visPedigree/reference/pedne.md)
      documentation now specifies that its `"coancestry"` method is
      based on $f_{ij}$.
2.  Individuals acting as both parents are now identified as
    `"monoecious"` in the `Sex` column.
3.  [`visped()`](https://luansheng.github.io/visPedigree/reference/visped.md)
    now uses a distinct teal color (`#26a69a`) to render `"monoecious"`
    individuals, ensuring clear visual separation from males, females,
    and highlighted nodes.
4.  Pedigree edges are now colored based on the parent’s role in a
    specific mating (Sire blue, Dam gold, Selfing teal) rather than
    invariant node sex, allowing monoecious individuals to display
    role-appropriate connection colors.
5.  [`summary()`](https://rdrr.io/r/base/summary.html) and
    [`print()`](https://rdrr.io/r/base/print.html) methods for `tidyped`
    objects now accurately report the count and percentage of monoecious
    individuals.
6.  **Efficiency Optimization**: Optimized
    [`pedancestry()`](https://luansheng.github.io/visPedigree/reference/pedancestry.md)
    initialization on large pedigrees by using vectorized matrix
    indexing, significantly reducing overhead for pedigrees with \>25k
    nodes.
7.  Added a new unit test for
    [`pedrel()`](https://luansheng.github.io/visPedigree/reference/pedrel.md)
    to verify correct relationship calculation in deep-inbreeding
    scenarios (Gen 4 relationships reaching 1.0).

## Changes in version 1.3.0 released on 10 Mar 2026

### New Features

1.  **Founder Genome Equivalents ($f_{g}$)**: Integrated the robust
    calculation of Founder Genome Equivalents into
    [`pediv()`](https://luansheng.github.io/visPedigree/reference/pediv.md).
    It directly evaluates the mean coancestry while properly correcting
    for diagonal intra-cohort elements through adaptive scaling, keeping
    computational costs linear relative to the reference cohort size.
2.  **Reproducible Parameter Inference**: Added a `seed` argument to
    both
    [`pedne()`](https://luansheng.github.io/visPedigree/reference/pedne.md)
    and
    [`pediv()`](https://luansheng.github.io/visPedigree/reference/pediv.md)
    functions functions enabling precise reproducible sampling for
    effective population size estimations (Ne) and $f_{g}$ computations
    using Monte Carlo approximations.

## Changes in version 1.2.3 released on 08 Mar 2026

### Bug Fixes

1.  **Trace edge highlighting in
    [`visped()`](https://luansheng.github.io/visPedigree/reference/visped.md)**:
    Fixed incorrect edge highlighting when using `trace = "all"`.
    Previously, when a node was highlighted as both an ancestor (via
    upward tracing) and a parent of descendants (via downward tracing),
    the cross-path edges were incorrectly highlighted. For example,
    `visped(tp, highlight = "X", trace = "all")` would incorrectly
    highlight the edge from N to Z1/Z2, even though that parent-child
    relationship is not on X’s trace path. The fix separates the up and
    down trace paths and uses `trace_edges` to precisely control which
    edges are highlighted.
2.  **Focal node upward edge in `trace = "down"`**: Fixed an issue where
    the focal node’s upward connection to its parents’ family node was
    incorrectly shown as highlighted when tracing downward only. For
    example, `visped(tp, highlight = "X", trace = "down")` would show
    X’s edge to the UxV family node in solid black, even though X’s
    ancestors are not part of the downward trace. Now,
    `individual → family` edges are only highlighted when the individual
    appears as a child in the traced path.

## Changes in version 1.2.2 released on 08 Mar 2026

### New Features

1.  **Unified Diversity Analysis (`pediv`)**: Added
    [`pediv()`](https://luansheng.github.io/visPedigree/reference/pediv.md)
    as a single entry-point wrapper that aggregates founder
    contributions ($f_{e}$), ancestor contributions ($f_{a}$), and all
    three Ne estimates (coancestry, inbreeding, demographic) into one
    consolidated `pediv` S3 object. A dedicated
    [`print.pediv()`](https://luansheng.github.io/visPedigree/reference/print.pediv.md)
    method provides a formatted summary table.
2.  **New Dataset (`complex_ped`)**: Added `complex_ped`, a
    multi-generation pedigree dataset suitable for testing deeper
    ancestry tracing and cross-generation diversity analyses.

### API Changes

1.  **Parameter Rename (`cand` → `reference`)**: Standardized the
    reference population parameter name across all relevant analysis
    functions:
    - `pedne(... , reference = NULL)` (previously `cand`)
    - `pedcontrib(... , reference = NULL)` (previously `cand`)
    - `pedrel(... , reference = NULL)` (previously `cand`)

    The new name better reflects population genetics terminology, where
    the target group is the **reference population** rather than a set
    of “candidates”. *(Note: Old `cand` argument is no longer supported;
    please update existing scripts.)*

### Documentation

1.  **Vignette Rewrite (`pedigree-analysis.Rmd`)**: Completely
    restructured the pedigree analysis vignette with expanded theory
    explanations for $f_{e}$, $f_{a}$, and Ne, updated code examples
    using
    [`pediv()`](https://luansheng.github.io/visPedigree/reference/pediv.md)
    and the new `reference` parameter, and improved narrative linking
    the statistical outputs to practical breeding decisions.
2.  **Workspace Reorganization**: Moved development-only files
    (`MACOS_OPENMP_FIX.md`, `manuscript.md`, analysis scripts) into
    `sandbox/` to keep the package root clean. Added corresponding
    `.gitignore` and `.Rbuildignore` rules.

## Changes in version 1.2.1 released on 07 Mar 2026

### New Features

1.  **Ancestral Analysis (`pedcontrib`)**: Added robust algorithms for
    assessing genetic diversity through gene origin probabilities.
    Computes the **effective number of founders ($f_{e}$)** via
    recursive gene derivation and the **effective number of ancestors
    ($f_{a}$)** via Boichard’s iterative algorithm.
2.  **Missing Parent Conservation (“Phantom Parents”)**: Implemented
    correct probability mass conservation. In `pedcontrib`, single
    missing parents (half-founders) are seamlessly augmented with
    temporary “phantom parents” before processing, overcoming the
    critical issue of gene probability leakage found in earlier tools.
3.  **Ancestry Proportions (`pedancestry`)**: Added
    [`pedancestry()`](https://luansheng.github.io/visPedigree/reference/pedancestry.md)
    function to trace line origins and monitor the surviving proportion
    of genes from specified historic founder lines or strains down to
    modern descendants.
4.  **Partial Inbreeding (`pedpartial`)**: Engineered the Meuwissen &
    Luo (1992) based partial inbreeding decomposition
    [`pedpartial()`](https://luansheng.github.io/visPedigree/reference/pedpartial.md).
    Enables breaking down the overall inbreeding coefficient into
    discrete fractions attributed to specifically targeted ancestors.
5.  **New Dataset (`half_founder_ped`)**: Added empirical ENDOG dataset
    containing instances of strictly missing single parents (sire known,
    dam unknown, etc.) specifically engineered to test and validate
    phantom-parent corrections.

### Performance

1.  **Peeling Core Engine**: Rebuilt the C++ core array engine backing
    the $f_{a}$ and $f_{e}$ calculations. Execution latency for
    incredibly massive and deep graphs (\>180,000 nodes) was resolved,
    avoiding hanging scenarios by limiting computational bounds
    optimally inside $O(K \times N)$ array states.

### Documentation

1.  **Analysis Indexing**: Expanded `_pkgdown.yml` configuration mapping
    to fully expose all newly engineered high-level pedigree statistical
    functions (`pedancestry`, `pedcontrib`, `pedpartial`, `pedecg`,
    etc.) within the main Reference documentation.
2.  **Analysis Vignettes**: Updated `vignettes/pedigree-analysis.Rmd`
    carefully illustrating Boichard’s genetic bottleneck interpretations
    ($f_{e}$ vs $f_{a}$) alongside new code examples for tracing
    targeted lineage flows.

## Changes in version 1.2.0 released on 04 Mar 2026

### New Features

1.  **Enhanced Effective Population Size (Ne) Calculation**:
    - The
      [`pedne()`](https://luansheng.github.io/visPedigree/reference/pedne.md)
      function has been significantly expanded and now supports three
      robust methods for estimating Ne in breeding populations:
      - **`coancestry`** (New Default): Based on the rate of coancestry
        ($\Delta c$). This method captures the loss of genetic potential
        and is considered the “gold standard” for populations under
        selection pressure. It typically yields a smaller, more
        conservative Ne estimate than inbreeding-based methods,
        providing a better early warning signal for genetic diversity
        loss.
      - **`inbreeding`**: Based on the individual rate of inbreeding
        ($\Delta F$). This method reflects the realized inbreeding but
        may overestimate Ne in managed populations where mating between
        relatives is actively avoided.
      - **`demographic`**: A census-based method using the number of
        breeding males ($N_{m}$) and females ($N_{f}$).
2.  **Parallel Processing Support**:
    - Introduced **OpenMP multi-threading** for the computationally
      intensive `coancestry` method. Users can now specify `ncores` to
      speed up large-scale matrix calculations.
    - Added a `nsamples` parameter to allow efficient estimation on
      massive pedigrees by sampling subsets of each cohort.

### Performance

1.  **C++ Optimization**:
    - Implemented a high-performance C++ backend
      (`cpp_calculate_sampled_coancestry_delta`) using `RcppArmadillo`.
      This replaces the previous R-based logic for coancestry
      calculations, enabling the analysis of much larger datasets.

### Documentation

1.  **Clarified Parameter Scopes**: Updated documentation for
    [`pedne()`](https://luansheng.github.io/visPedigree/reference/pedne.md)
    to explicitly state that `ncores` and `nsamples` parameters are
    specific to the `method = "coancestry"` calculation path.
2.  **Method Descriptions**: Expanded details on the three Ne
    calculation methods to help users choose the most appropriate metric
    for their breeding program.

## Changes in version 1.1.1 released on 02 Mar 2026

### New Features

1.  **pedgenint Sex-Independent Pathways**: Added evaluation of `SO`
    (Sire-to-Offspring) and `DO` (Dam-to-Offspring) generation intervals
    alongside the standard 4 pathways. This is especially useful for
    aquatic species (like shrimp) or early-stage screening where
    offspring sex might remain unknown.

### API Changes and Refactoring

1.  **pedne Interface Standardization**:
    - Renamed arguments `timevar` to `by`, and `cohort` to `cand` to
      harmonize parameter naming conventions across the package.
    - Removed unused and misleading parameters (`unit`, `cycle_length`,
      `maxgen`). The effective population size Ne calculation innately
      depends on Equivalent Complete Generations (ECG), making it
      independent of scalar temporal units.
2.  **vismat Parameter Alignment**: Renamed `grouping` argument to `by`
    to maintain grouping consistency.  
    *(Note: Old arguments `timevar`, `cohort` in
    [`pedne()`](https://luansheng.github.io/visPedigree/reference/pedne.md)
    and `grouping` in
    [`vismat()`](https://luansheng.github.io/visPedigree/reference/vismat.md)
    are retained for backward compatibility but will display a
    deprecation warning.)*

### Bug fixes

1.  **pedrel Correctness**: Fixed a critical calculation bug in
    [`pedrel()`](https://luansheng.github.io/visPedigree/reference/pedrel.md)
    where the mean average relatedness calculation erroneously divided
    the sum of the full relationship matrix (including all traced
    ancestors) by only the size of the target subgroup. It now cleanly
    subsets the relationship matrix, and correctly handles boundary
    limits (`NUsed < 2`). The output columns `N` and `MeanRel` behavior
    has been replaced with `NTotal`, `NUsed`, and `MeanRel`.
2.  **pedgenint Aggregation**: Fixed
    [`pedgenint()`](https://luansheng.github.io/visPedigree/reference/pedgenint.md)
    to output appropriate unweighted mixture standard deviation for
    generating generation intervals alongside its unweighted 4-pathway
    average interval estimate.
3.  **pedgenint Sample Size (N)**: Fixed an issue where the `Average`
    pathway N was severely underestimated. It now accurately evaluates
    all parent-offspring pairs via `calc_all_pathway()`.
4.  **pedcontrib Accuracy**: Standardized effective founders (`Ne_f`)
    and effective ancestors (`Ne_a`) calculation in
    [`pedcontrib()`](https://luansheng.github.io/visPedigree/reference/pedcontrib.md)
    to ensure they are calculated based upon the full un-truncated
    cohort before outputting strictly the `top` n-ranked figures.
    Results list has been augmented with variables tracking the `_total`
    and `_reported` count values.
5.  **pedcontrib Deep Pedigree Latency**: Replaced a string-named vector
    backward pass with a pure integer-indexed backward pass, resolving
    instances where evaluating contributions on deep, large pedigrees
    (e.g., \> 200,000 records) would hang indefinitely due to scaling
    constraints.
6.  **pedpartial / pedancestry Input Compatibility**: Ensured missing
    numeric identifiers in incoming pedigrees (e.g. `addnum = FALSE`) do
    not break
    [`pedpartial()`](https://luansheng.github.io/visPedigree/reference/pedpartial.md)
    or
    [`pedancestry()`](https://luansheng.github.io/visPedigree/reference/pedancestry.md).
    Increased performance of the pedigree propagation loop in
    `pedancestry` by dropping an internal array linear probe algorithm
    with an immediate linear vector lookup.
7.  **pedne Performance bottleneck**: Removed an obsolete `O(N^2)`
    individual traversal evaluation (`calc_ancestral_f()`), streamlining
    calculation purely around the efficient direct formula by Gutiérrez
    et al.

## Changes in version 1.1.0 released on 01 Mar 2026

### New Features

1.  **Pedigree Analysis Module**: Introduced a comprehensive suite of
    pedigree analysis and population genetics tools.
    - [`pedstats()`](https://luansheng.github.io/visPedigree/reference/pedstats.md):
      Calculate holistic and demographic statistics.
    - [`pedrel()`](https://luansheng.github.io/visPedigree/reference/pedrel.md):
      Formulate average relatedness within specific population
      groupings.
    - [`pedgenint()`](https://luansheng.github.io/visPedigree/reference/pedgenint.md):
      Compute distinct breeding pathways (SS, SD, DS, DD) and overall
      population generation intervals.
    - [`pedcontrib()`](https://luansheng.github.io/visPedigree/reference/pedcontrib.md):
      Determine genetic contributions from founders (`Ne_f`) and
      prominent ancestors (`Ne_a`) utilizing iterative gene flow
      derivations.
    - [`pedancestry()`](https://luansheng.github.io/visPedigree/reference/pedancestry.md):
      Establish proportionality of ancestral lineages on subsequent
      descendants.
    - [`pedpartial()`](https://luansheng.github.io/visPedigree/reference/pedpartial.md):
      Decompose inbreeding mechanisms to detect fractional/partial
      origins from core ancestors.
2.  **Pedigree Analysis Visualization**: Added
    [`vispstat()`](https://luansheng.github.io/visPedigree/reference/vispstat.md)
    to intuitively render bar charts of generation intervals and
    histogram distributions detailing depth tracking factors (like
    Equivalent Complete Generations).

## Changes in version 1.0.1 released on 31 Jan 2026

CRAN release: 2026-02-23

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
(`pedmatrix`, etc.) have been removed. Please use
[`pedmat()`](https://luansheng.github.io/visPedigree/reference/pedmat.md)
directly.

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

### Bug fixes

1.  Fixed overlapping edge detection for small pedigree graphs.
2.  Improved coloring consistency for compact mode.

## Changes in version 0.2.6 released on 31 Mar 2020

### Bug fixes

1.  Fixed a bug that the number of generations for candidates would be
    traced to n+1 when tracegen=n. This bug is found by Mianyu Liu.

## Changes in version 0.2.5 released on 25 Feb 2020

### Bug fixes

1.  The tidyped() does not work with trace=‘all’ in [certain
    cases](https://github.com/luansheng/visPedigree/issues/2#issue-568599008)

## Changes in version 0.2.4.1 released on 24 Feb 2020

### Bug fixes

1.  An unexpected column with the name as NA occured when a tidyped
    object is tidyed again using the tidyped()

## Changes in version 0.2.4 released on 12 June 2019

### Bug fixes

1.  The data.table used as the input parameter ‘ped’ may be changed in
    tidyped() and visped().

## Changes in version 0.2.3 released on 05 Mar 2019

### Bug fixes

1.  The generation number of individuals is not inferred rightly.

## Changes in version 0.2.2 released on 28 Jan 2019

### Bug fixes

1.  The tidied pedigree will not include the candidates which are not in
    the Ind column of the origin pedigree when the cand parameter is not
    NULL.

## Changes in version 0.2.1 released on 17 Nov 2018

### Bug fixes

1.  Repel the overlapping nodes due to very small differences (digits
    \> 7) among x positions of nodes
