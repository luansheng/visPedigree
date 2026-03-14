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

### 1.2 Internal Functions (Private)
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
