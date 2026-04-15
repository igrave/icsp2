# Control ic_sp2 fitting options

Control ic_sp2 fitting options

## Usage

``` r
ic_sp_control(
  useGA = TRUE,
  maxIter = 10000,
  baseUpdates = 5,
  regStart = NULL,
  derivMethod = c(12, 1),
  updateCovars = TRUE
)
```

## Arguments

- useGA:

  Should constrained gradient ascent step be used?

- maxIter:

  Maximum iterations

- baseUpdates:

  number of baseline updates (ICM + GA) per iteration

- regStart:

  Initial values for regression parameters

- derivMethod:

  Method for derivative calculations.

- updateCovars:

  Should covariates be updated during fitting?

  @description Creates the control options for
  [`ic_sp2()`](https://igrave.github.io/icsp2/reference/ic_sp2.md).
  Defaults not intended to be changed for use in standard analyses.

## Details

The constrained gradient step, controlled by `useGA = TRUE`, is a step
that was added to improve the convergence in a special case. The option
to turn it off is only in place to help demonstrate its utility.

`regStart` also for seeding of initial value of regression parameters.
Intended for use in “warm start" for bootstrap samples and providing
fixed regression parameters.
