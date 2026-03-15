# visPedigree Development Guide

This document outlines the coding standards, naming conventions, and design principles for the `visPedigree` package. All contributors (human and AI) are expected to follow these guidelines to ensure consistency.

## 1. Naming Conventions

### 1.1 Exported Functions (User-Facing)
*   **Style**: All lowercase, **NO underscores**. 
*   **Format**: `[prefix][noun/verb]` or `[noun]`.
*   **Prefixes**:
    *   `ped...`: For core pedigree processing and statistical analysis functions (e.g., `pedstats`, `pedecg`).
    *   `vis...`: For visualization functions (e.g., `visped`, `visstats`).
*   **Examples**:
    *   ✅ `tidyped()`
    *   ✅ `pedstats()`
    *   ✅ `pedmat()`
    *   ❌ `ped_stats()`
    *   ❌ `calculatePedigree()`
*   **Exemptions**: R-idiomatic predicates and coercers (`is_*`, `as_*`, `has_*`) are allowed with underscores, because the R community expects this pattern. Current exemptions: `is_tidyped()`, `as_tidyped()`, `has_inbreeding()`, `has_candidates()`.
*   **Known legacy deviations**: `query_relationship()`, `expand_pedmat()`, `summary_pedmat()` use underscores for readability and backward compatibility. These are considered "grandfathered" and should not be used as templates for new functions. New analytical helpers must follow the no-underscore rule (e.g., `pedstats`, `pedrel`).
*   **Rule of Thumb**: If you are unsure, avoid underscores. Only use them for the `is`/`as`/`has` patterns.
*   **Style**: Snake case (`noun_verb` or `verb_noun`) is allowed to distinguish them from exported API.
*   **Examples**: `validate_and_prepare_ped`, `check_ped_loops`.

### 1.3 Classes and S3 Methods
*   **S3 Class Names**: All lowercase, must match the constructor function name.
    *   Example: `tidyped`, `pedstats`.
*   **S3 Methods**: Follow R standard `method.class`.
    *   Example: `summary.tidyped`, `plot.pedstats`.

### 1.4 List Components (Return Values)
When a function returns a list or an S3 object (which is a list), the names of the list elements should be:
*   **Style**: Snake case (lowercase + underscores) allowed for readability.
*   **Examples**:
    *   `x$summary`
    *   `x$gen_interval`
    *   `x$founders`

### 1.5 Data Columns (data.table)
Columns in `data.table` or `data.frame` objects returned by the package must use **PascalCase** (First letter capitalized). This matches the core `tidyped` object structure.
*   **Reserved Core Columns**: `Ind`, `Sire`, `Dam`, `Sex`, `Gen`, `IndNum`, `SireNum`, `DamNum`, `Family`.
*   **Statistical Columns**:
    *   ✅ `MeanF`, `DeltaF`, `Ne`
    *   ✅ `ECG`, `Completeness`
    *   ✅ `Group`, `Year`
    *   ❌ `mean_f`, `delta_F`, `ne`

### 1.6 Arguments
*   **Style**: Lowercase preferred. Underscores allowed for long arguments if necessary for readability, but short lowercase is better.
*   **Examples**: `ped`, `cand`, `trace`, `by`, `timevar`, `cohort`.

## 2. API Design Principles

*   **Object-Centric**: Statistics and analysis should revolve around the `tidyped` object. 
*   **Tidy Output**: Functions typically return `data.table` objects or lists of `data.table`s.
*   **Performance via `data.table`**: For R-level data manipulation, always prioritize `data.table` syntax and modify-in-place operations (`:=`) over base R to guarantee high performance and minimal memory overhead on large pedigree datasets.
*   **Layered Computation**:
    *   **Level 1 (Basic)**: Pure R implementation.
    *   **Level 2/3 (Heavy)**: Use `Rcpp` / `RcppArmadillo` for computationally intensive tasks (e.g., matrix inversion, iterative gene contributions), leveraging the `IndNum` integer indices.
*   **`data.table` Return Visibility**: Any function that returns a `data.table` object **must** use `return(x[])` (with trailing `[]`) instead of `return(x)`. The empty `[]` clears the invisible flag set by `:=` and `set*()` operations, ensuring the object auto-prints correctly in console and knitr.

## 3. Anti-Bloat Governance Checklist

The package is allowed to grow in lines of code, but it must not grow in uncontrolled complexity. The goal is to keep `visPedigree` extensible without turning it into a fragile, exception-driven codebase.

Before merging any non-trivial change, contributors should check the following:

### 3.1 Public API Gate
*   Do **not** add a new exported function unless an existing exported function cannot reasonably absorb the feature.
*   Prefer extending an existing object-centric workflow over creating a parallel entry point.
*   Every new exported argument should have a clear long-term meaning, not just solve one temporary edge case.
*   New exported helpers must justify their maintenance cost in documentation, tests, and examples.

### 3.2 Delete-or-Consolidate Rule
*   Every new helper should trigger a quick search for older helpers that can now be removed or merged.
*   If the same validation or preprocessing logic appears twice, centralize it.
*   AI-generated code should not merely wrap old logic in a new layer; it should simplify the system where possible.

### 3.3 Prefer Explicit Failure Over Silent Recovery
*   When pedigree structure is invalid for an analysis, prefer a clear error over an implicit repair that may change scientific meaning.
*   Auto-recovery is acceptable only when the underlying semantics are preserved (for example, restoring a dropped class when the pedigree structure is still intact).
*   Never silently compute pedigree-derived quantities on row-truncated ancestry.

## 4. Architectural Boundaries

### 4.1 Core `tidyped` Invariants
The following rules define the package's most important internal contract:

*   Core pedigree columns must remain semantically stable: `Ind`, `Sire`, `Dam`, `Sex`, `Gen`, `IndNum`, `SireNum`, `DamNum`.
*   `IndNum` must correspond to row order when algorithms assume 1..N integer indexing.
*   `SireNum` and `DamNum` must either be 0 / missing-parent encodings or valid row-index references into `IndNum`.
*   A **complete pedigree** means every non-missing `Sire` and `Dam` identifier is present in `Ind`.
*   Row-truncated subsets are not valid `tidyped` pedigrees for completeness-sensitive analyses.
*   Pedigree metadata must be stored consistently in `attr(x, "ped_meta")`.
*   `ped_meta` is a named list with three fixed fields:
    *   `selfing` (logical): whether selfing/monoecious mode was used.
    *   `bisexual_parents` (character vector): IDs appearing as both sire and dam.
    *   `genmethod` (character): generation assignment method (`"top"` or `"bottom"`).
*   All `tidyped` objects must be created through `new_tidyped()`, which attaches the class via `setattr` (no copy) and clears data.table's invisible flag via `x[]`.

### 4.1.1 `[.tidyped` Subsetting Contract
The `[.tidyped` method is the primary guardian of structural integrity:

*   **Complete subset** (all referenced parents still present): class is preserved, `IndNum`/`SireNum`/`DamNum` are rebuilt in-place, metadata is retained.
*   **Incomplete subset** (parent records removed): result is degraded to plain `data.table` with a warning. The user is guided to use `tidyped(tp, cand = ids, trace = "up")`.
*   **`:=` operations**: modify by reference, class and metadata are preserved via `setattr`.
*   **Column-only selections** (core columns missing): returned as plain `data.table` without warning.

### 4.2 Function Classes by Structural Requirement
All exported analysis functions should mentally fall into one of these groups:

1. **Complete-pedigree required** — guarded by `ensure_complete_tidyped()`.  
    Examples: `inbreed()`, `pedecg()`, `pedrel()`, `pedgenint()`, `pedcontrib()`, `pedancestry()`, `pedpartial()`, `pediv()`, `pedmat()`.

2. **Structure-light summaries** — guarded by `ensure_tidyped()`.  
    Examples: `pedsubpop()`, `splitped()`, `pedne(method="demographic")`.

3. **Post-computation summaries on existing results** — guarded by `ensure_tidyped()`, with conditional escalation to `ensure_complete_tidyped()` only when new pedigree recursion is needed.  
    Examples: `pedfclass()` (Class 2 when `f` column already exists; Class 1 when it must compute `f`).

4. **Conditional — parameter-dependent** — guard level depends on which code path the user selects.  
    Examples: `pedstats(ecg=TRUE)` requires Class 1; `pedstats(ecg=FALSE, genint=FALSE)` only needs Class 2. `pedne(method="coancestry")` requires Class 1; `pedne(method="demographic")` only needs Class 2.

5. **Visualization** — guarded by `validate_tidyped()` (auto-recovery without completeness requirement).  
    Examples: `visped()`, `vismat()`, `vispstat()`.

Any new function should explicitly choose one class above before implementation.

### 4.3 Layering Rules
Use this layering model consistently:

*   **User layer**: exported functions, user-facing errors, stable API.
*   **Workflow layer**: orchestration helpers that sequence validation, preprocessing, and aggregation.
*   **Computation layer**: pure pedigree/statistical logic.
*   **Infrastructure layer**: S3 methods, validators, metadata accessors, class restoration, C++ bridges.

Cross-layer leakage should be minimized. In particular, exported analysis functions should not each reinvent class restoration, completeness checks, or index rebuilding.

## 5. File and Function Size Heuristics

These are not rigid limits, but they are review triggers:

*   Exported function larger than ~120 lines: consider splitting workflow and computation.
*   Internal function larger than ~80 lines: consider extracting helpers.
*   File larger than ~800-1000 lines: review whether it contains multiple domains and should be split.
*   More than 2 repeated code blocks with the same shape: refactor to a helper.
*   More than 3 nested control levels: simplify.

Large files are acceptable only when they still have a single coherent responsibility.

### 5.1 Current Status (v1.4.1)
Known files exceeding the threshold — pending refactoring:

*   `R/pedanalysis.R` (~1940 lines): contains 12+ functions spanning generation intervals, ECG, Ne, founder/ancestor contributions, inbreeding classification, partial inbreeding, and diversity analysis. Candidate split: by analytical domain (depth, diversity, inbreeding).
*   `R/pedmatrix.R` (~1460 lines): contains `pedmat()` core, compact mode, `query_relationship()`, `expand_pedmat()`, and all pedmat S3 methods. Candidate split: core vs. compact vs. S3 methods.

## 6. Testing Rules That Protect Architecture

Tests should protect not only numeric outputs, but also package design constraints.

Every important feature should consider coverage in at least one of these categories:

*   **Invariant tests**: class retention, metadata retention, index integrity, complete-pedigree rules.
*   **Regression tests**: previously fixed bugs must stay fixed.
*   **Failure-path tests**: invalid inputs should fail clearly and for the right reason.
*   **Round-trip tests**: object survives common workflows such as subsetting, merging, tracing, and plotting.

For pedigree analysis, tests that lock down structural semantics are often more valuable than adding yet another numeric snapshot.

## 7. AI-Assisted Development Rules

AI contributors are useful for implementation speed, but they have a natural tendency to add wrappers, exceptions, and duplicate helpers. To prevent architectural drift:

*   Do not ask AI only to “make it work”; ask it to preserve invariants and minimize new surface area.
*   Prefer AI changes that simplify existing logic over changes that add compatibility layers.
*   Require AI to explain whether a change introduces a new concept, a new API, or only an internal consolidation.
*   If AI proposes a new helper, review whether an existing helper should be replaced.
*   If AI proposes silent fallback behavior, review whether fail-fast is scientifically safer.

## 8. Refactoring Triggers

Schedule refactoring when one or more of the following happens:

*   The same bug pattern appears in multiple functions.
*   A new feature requires touching many unrelated files.
*   A single file begins to mix validation, orchestration, computation, and printing.
*   New contributors cannot easily tell which helper they should call.
*   You feel reluctant to modify a function because breakage risk is unclear.

Refactoring should be treated as routine maintenance, not as an admission of failure.

## 9. Lightweight Review Checklist for Each PR or Commit

Before merging, ask:

1. Did this change add any new exported API? If yes, was that really necessary?
2. Did this change centralize logic, or duplicate it?
3. Did this change preserve `tidyped` invariants?
4. Did this change make scientific failure modes clearer?
5. Did this change add or update regression tests?
6. Could any old code now be deleted?
7. Is the package simpler, or only more capable?

If capability increased but simplicity decreased, review again before merging.
