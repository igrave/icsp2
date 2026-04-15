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

  Integer vector of event indicators (1 = event, 0 = censored).

- group:

  Integer vector of group labels, coded as consecutive integers starting
  from 1 (i.e. `1, 2, ...`).

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

## Examples

``` r
set.seed(1992)
time  <- rexp(100)
event <- sample(0:1, 100, replace = TRUE)
group <- sample(1:2, 100, replace = TRUE)
cpp_logrank(time, event, group)
#> 
#> Log-rank test
#> Sample sizes by group:
#>  1  2 
#> 47 53 
#> 
#>          N Observed Expected     O-E
#> Group 1 47       20  23.5486 -3.5486
#> Group 2 53       28  24.4514  3.5486
#> 
#> Variance-covariance matrix:
#>          1        2
#> 1  11.8621 -11.8621
#> 2 -11.8621  11.8621
```
