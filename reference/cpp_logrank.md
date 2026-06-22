# Internal fast log-rank test implementation

This is an internal implementation of right-censored logrank test for
use within
[ic_logrank](https://igrave.github.io/icsp2/reference/ic_logrank.md)
where many tests need to be performed efficiently for variance
calculation.

## Usage

``` r
cpp_logrank(time, event, group)
```

## Arguments

- time:

  Numeric vector of survival/censoring times.

- event:

  Logical vector of event indicators (TRUE = event, FALSE = censored).

- group:

  Integer vector of group labels, coded as consecutive integers (i.e.
  `1, 2, ...`). This is not checked.

## Value

An object of class `"cpp_logrank"`, a list with components:

- observed:

  Numeric vector of observed events per group.

- expected:

  Numeric vector of expected events per group.

- variance:

  Numeric variance matrix

- n_by_group:

  Numeric vector of sample sizes per group.

## Details

Compared to Survival::survdiff, this implementation is faster but does
not have the same flexibility or checks.

## Examples

``` r
set.seed(1992)
time  <- rexp(100)
event <- sample(c(TRUE, FALSE), 100, replace = TRUE)
group <- sample(1:2, 100, replace = TRUE)
cpp_logrank(time, event, group)
#> 
#> Log-rank test
#> Sample sizes by group:
#>  1  2 
#> 47 53 
#> 
#>          N Observed Expected     O-E
#> Group 1 47       27  25.6869  1.3131
#> Group 2 53       25  26.3131 -1.3131
#> 
#> Variance-covariance matrix:
#>          1        2
#> 1  12.5806 -12.5806
#> 2 -12.5806  12.5806
```
