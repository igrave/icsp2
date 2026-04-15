# Sun's Log-Rank Test for Interval-Censored Data

Sun's Log-Rank Test for Interval-Censored Data

## Usage

``` r
ic_logrank(
  formula,
  data,
  subset,
  na.action,
  B = c(0, 1),
  n_samples = 1000,
  type = c("sas", "hly"),
  ...
)
```

## Arguments

- formula:

  A formula with `Surv(l, u, type = 'interval2')` response and a single
  grouping variable on the right-hand side. May also contain
  [`strata()`](https://rdrr.io/pkg/survival/man/strata.html) terms for
  stratified analysis.

- data:

  A data frame containing the variables in the formula.

- subset:

  Optional expression indicating which subset of rows to use.

- na.action:

  Function to handle missing values.

- B:

  A vector of length 2 giving bounds for observation times. Default is
  c(0, 1).

- n_samples:

  The number of "imputation" samples for the variance calculation.
  Default is 1000.

- type:

  The method for calculating the test statistic and variance. Options
  are `"sas"` (default) or `"hly"`.

- ...:

  Additional arguments (currently unused).

## Value

An object of class `ic_logrank` containing:

- logrank:

  The log-rank statistics "observed - expected" for all groups and
  strata

- logrank_overall:

  The log-rank statistics "observed - expected" for all groups

- statistic:

  The overall chi-squared test statistic based on imputation

- df:

  Degrees of freedom

- p.value:

  P-value from chi-squared distribution

- var:

  Variance-covariance matrix

- groups:

  Group levels being compared

- n:

  Sample sizes per group

- method:

  Description of test method

- data.name:

  Name of the data

- call:

  The matched call

## Details

Performs Sun's (1996) non-parametric log-rank test for comparing
survival distributions across groups with interval-censored data. This
test compares group-specific NPMLE survival curves without assuming
proportional hazards.

The test statistic is \\Q = U(0)' I(0)^{-1} U(0)\\, which follows a
chi-squared distribution with k-1 degrees of freedom under the null
hypothesis, where k is the number of groups.

The default test type is `"sas"`, which uses Sun's test statistic
combined with variance estimated based on the Huang, Lee and Yu (2008)
procedure sampling exact observation times from the Turnbull intervals.

Alternatively, the `"hly"` type calculates the test statistic and
variance using the multiple imputation approach of Huang, Lee and Yu
(2008) directly.

Stratified tests are constructed by calculating the \\U\\ and \\V\\
matrices for each stratum separately and then summing the
stratum-specific matrices to give a global test statistic \\\sum\bar{U}'
(\sum\hat{V})^{-1} \sum\bar{U}\\. This is the procedure described in the
SAS documentation for PROC ICLIFETEST.

## References

Sun, J. (1996). A non-parametric test for interval-censored failure time
data with application to AIDS studies. *Statistics in Medicine*, 15(13),
1387-1395.

Huang, J., Lee, C., and Yu, Q. (2008). A Generalized Log-Rank Test for
Interval-Censored Failure Time Data via Multiple Imputation. *Statistics
in Medicine 27:3217–3226*. <http://dx.doi.org/10.1002/sim.3211>

SAS Institute Inc. (2026). *SAS/STAT® 26.03 User's Guide: The ICLIFETEST
Procedure*.
<https://documentation.sas.com/doc/en/statug/latest/statug_iclifetest_details01.htm>
Accessed 14 April 2026.

## Examples

``` r
# Simple two-group comparison
data(miceData)
ic_logrank(Surv(l, u, type = "interval2") ~ grp, data = miceData)
#> 
#> Sun's log-rank test for interval-censored data 
#> ============================================== 
#> 
#> Call:
#> ic_logrank(formula = Surv(l, u, type = "interval2") ~ grp, data = miceData)
#> 
#> Sample sizes by group:
#> Group
#> ce ge 
#> 96 48 
#> 
#> Log-rank Statistics:
#>       ce       ge 
#> -3.40709  3.40709 
#> 
#>   Chi-squared statistic: Q = 0.9351
#>   Degrees of freedom:    df = 1
#>   P-value:               p = 0.3335
#> 
#> Variance-covariance matrix:
#>          ce       ge
#> ce  12.4135 -12.4135
#> ge -12.4135  12.4135
#> 
#> Calculated with  1000  samples
```
