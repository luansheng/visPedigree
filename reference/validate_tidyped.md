# Internal validator for tidyped class

Validates a tidyped object. If the object has lost its class but retains
the required columns, it is automatically restored via
[`ensure_tidyped()`](https://luansheng.github.io/visPedigree/reference/ensure_tidyped.md).
Fatal structural problems (e.g., missing core columns) raise an error.

## Usage

``` r
validate_tidyped(x)
```

## Arguments

- x:

  A tidyped object

## Value

The (possibly restored) object if valid, otherwise an error
