# Simulates interval censored data from regression model with a Weibull baseline

Simulates interval censored data from a regression model with a Weibull
baseline distribution. Used for demonstration

## Usage

``` r
simIC_weib(
  n = 100,
  b1 = 0.5,
  b2 = -0.5,
  model = "ph",
  shape = 2,
  scale = 2,
  inspections = 2,
  inspectLength = 2.5,
  rndDigits = NULL,
  prob_cen = 1
)
```

## Arguments

- n:

  Number of samples simulated

- b1:

  Value of first regression coefficient

- b2:

  Value of second regression coefficient

- model:

  Type of regression model. Options are 'po' (prop. odds) and 'ph' (Cox
  PH)

- shape:

  shape parameter of baseline distribution

- scale:

  scale parameter of baseline distribution

- inspections:

  number of inspections times of censoring process

- inspectLength:

  max length of inspection interval

- rndDigits:

  number of digits to which the inspection time is rounded to, creating
  a discrete inspection time. If `rndDigits = NULL`, the inspection time
  is not rounded, resulting in a continuous inspection time

- prob_cen:

  probability event being censored. If event is uncensored, l == u

## Details

Exact event times are simulated according to regression model: covariate
`x1` is distributed `rnorm(n)` and covariate `x2` is distributed
`1 - 2 * rbinom(n, 1, 0.5)`. Event times are then censored with a case
II interval censoring mechanism with `inspections` different inspection
times. Time between inspections is distributed as
`runif(min = 0, max = inspectLength)`. Note that the user should be
careful in simulation studies not to simulate data where nearly all the
data is right censored (or more over, all the data with x2 = 1 or -1) or
this can result in degenerate solutions!

## Author

Clifford Anderson-Bergman

## Examples

``` r
set.seed(1)
sim_data <- simIC_weib(n = 500, b1 = .3, b2 = -.3, model = 'ph',
                      shape = 2, scale = 2, inspections = 6,
                      inspectLength = 1)
#simulates data from a cox-ph with beta Weibull distribution.

ic_sp_ph(Surv(l, u, type = 'interval2') ~ x1 + x2, data = sim_data)
#> Call:
#> ic_sp_ph(formula = Surv(l, u, type = "interval2") ~ x1 + x2, 
#>     data = sim_data)
#> 
#> Coefficients:
#>         x1         x2 
#>  0.0410343 -0.3506595 
#> 
#> Log-likelihood: -753.3359
#> Number of iterations: 10
ic_sp_po(Surv(l, u, type = 'interval2') ~ x1 + x2, data = sim_data)
#> Call:
#> ic_sp_po(formula = Surv(l, u, type = "interval2") ~ x1 + x2, 
#>     data = sim_data)
#> 
#> Coefficients:
#>          x1          x2 
#> -0.08955382  0.47309564 
#> 
#> Log-likelihood: -759.4884
#> Number of iterations: 8

#'ph' fit looks better than 'po'; the difference between the transformed survival
#function looks more constant
```
