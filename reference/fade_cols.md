# Fade colors by appending a reduced alpha value

Converts any R color specification to \`#RRGGBB4D\` form. Handles hex
colors (\`#RRGGBB\`, \`#RRGGBBAA\`) and named colors (e.g. \`"red"\`).

## Usage

``` r
fade_cols(x)
```

## Arguments

- x:

  Character vector of colors.

## Value

Character vector of faded hex colors.
