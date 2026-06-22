# Fit a semi-parametric model for interval-censored data

This function fits a semi-parametric model for interval-censored data
using either the proportional hazards or proportional odds model.
Aliases `ic_sp_ph()` and `ic_sp_po()` are provided for convenience.

## Usage

``` r
ic_sp2(
  formula,
  data,
  weights,
  subset,
  na.action,
  B = c(0, 1),
  control = ic_sp_control(...),
  model = c("ph", "po"),
  profile_ci = NULL,
  ...
)

ic_sp_ph(
  formula,
  data,
  weights,
  subset,
  na.action,
  B = c(0, 1),
  control = ic_sp_control(...),
  model = c("ph", "po"),
  profile_ci = NULL,
  ...
)

ic_sp_po(
  formula,
  data,
  weights,
  subset,
  na.action,
  B = c(0, 1),
  control = ic_sp_control(...),
  model = c("ph", "po"),
  profile_ci = NULL,
  ...
)
```

## Arguments

- formula:

  A model formula with `Surv(l, u, type = 'interval2')` response and
  covariates on the right-hand side. May also contain
  [`strata()`](https://rdrr.io/pkg/survival/man/strata.html) terms.

- data:

  A data frame containing the variables in the formula, including strata
  terms.

- weights:

  Optional vector of weights for each observation, or the name of a
  variable in `data` containing the weights.

- subset:

  Optional expression indicating a subset of the rows of `data` to be
  used in the fit.

- na.action:

  Optional function to handle missing data. Default is `na.omit`.

- B:

  A vector of length 2 giving the lower and upper bounds for the
  observation times. Default is c(0, 1).

- control:

  A list of control settings, with defaults created by
  [`ic_sp_control()`](https://igrave.github.io/icsp2/reference/ic_sp_control.md).

- model:

  Type of model to fit. Choices are `"ph"` for proportional hazards and
  `"po"` for proportional odds. Default is `"ph"`. This is normally
  determined by the function aliases `ic_sp_ph` and `ic_sp_po`.

- profile_ci:

  Confidence level for profile likelihood confidence intervals based on
  the likelihood ratio. Set to `NULL` to skip this profile likelihood
  confidence interval calculations.

- ...:

  Additional arguments passed to control.

## Value

An object of class `"ic_sp2"`, a list containing the fitted model
information, including coefficients, variance-covariance matrix, and
other details.

## Examples

``` r
fit_ph <- ic_sp_ph(Surv(l, u, type = "interval2") ~ grp, data = miceData, profile_ci = 0.95)
summary(fit_ph)
#> Call:
#> ic_sp_ph(formula = Surv(l, u, type = "interval2") ~ grp, data = miceData, 
#>     profile_ci = 0.95)
#> 
#> Coefficients:
#>       Estimate Std.Error z.value p.value
#> grpge   0.6785    0.3537  1.9184   0.055
#> 
#> Log-likelihood: -76.56894 
#> Iterations: 28 

# The CI for the regression parameter(s) based on the likelihood ratio:
fit_ph$profile_ci
#>            lower    upper
#> grpge -0.1590267 1.366948

confint(fit_ph, level = 0.95)
#>              2.5%    97.5%
#> grpge -0.01470717 1.371635
```
