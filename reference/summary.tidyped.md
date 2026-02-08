# Summary method for tidyped objects

Summary method for tidyped objects

## Usage

``` r
# S3 method for class 'tidyped'
summary(object, ...)
```

## Arguments

- object:

  A tidyped object.

- ...:

  Additional arguments (ignored).

## Value

A summary.tidyped object (list) containing:

- `n_ind`: Total number of individuals.

- `n_male`, `n_female`, `n_unknown_sex`: Sex composition counts.

- `n_founders`: Number of individuals with no known parents.

- `n_both_parents`: Count of individuals with complete parentage.

- `max_gen`, `gen_dist`: (Optional) Maximum generation and its
  distribution.

- `n_families`, `family_sizes`, `top_families`: (Optional) Family
  statistics.

- `f_stats`, `n_inbred`: (Optional) Inbreeding coefficient statistics.

- `n_cand`, `cand_f_stats`: (Optional) Candidate-specific statistics.
