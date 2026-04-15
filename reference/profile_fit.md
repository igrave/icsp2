# Refit the model with fixed regression parameters

Refit the model with fixed regression parameters

## Usage

``` r
profile_fit(object, beta = object$coefficients)
```

## Arguments

- object:

  Fitted model object from `ic_sp`

- beta:

  Vector of regression parameters to fix in the refit. Default is the
  original fitted regression parameters.

## Value

Raw fitted model object with fixed regression parameters.
