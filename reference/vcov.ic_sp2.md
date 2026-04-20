# Profile Likelihood Covariance for Semi-Parametric Models

Profile Likelihood Covariance for Semi-Parametric Models

## Usage

``` r
# S3 method for class 'ic_sp2'
vcov(object, type = "oim_curvature", fixed = 5, typical = 1, large = 2, ...)
```

## Arguments

- object:

  Fitted model object from `ic_sp`

- type:

  One of `"oim_curvature"` (default), `"oim_fixed"`, or `"opg_fixed"`.
  See details for explanation.

- fixed:

  A fixed factor to multiply by `n^(-1/2)` to determin the perturbation
  size for fixed types.

- typical:

  A typical value for the regression parameters, used to determine the
  scale of `h_n`. Default is 1. This is required for the
  `"oim_curvature"` type.

- large:

  A large value for the regression parameters, used to determine the
  scale of `h_n`. Default is 2. This is required for the
  `"oim_curvature"` type.

- ...:

  Unused.

## Value

Variance-covariance matrix of the regression parameters.

## Details

The covariance matrix is calculated using the profile likelihood
approach. (Murphey and Vand Der Vaart 2000). This method involves
perturbing the regression parameters, updating the baseline hazard
estimates using the `profile_fit` function with the perturbed
parameters, and calculating the change in log-likelihood for each
perturbation, which is then used to compute the covariance matrix. We
borrowing the naming convention from the Stata `stintcox` manual
<https://www.stata.com/manuals/ststintcox.pdf>.

Type `"oim_curvature"` (Boruvka and Cook 2015) uses the curvature of the
log-likelihood function to determine the perturbation size for the
profile likelihood calculations, as described in Boruvka and Cook
(2014). The `typical` and `large` parameters are used to determine the
scale of the perturbations.

Type `"oim_fixed"` (Zeng et al 2016) uses a fixed perturbation size
based on the `fixed` parameter and the sample size `n`.

Type `"opg_fixed"` (Zeng et al 2017) uses a fixed perturbation size
based on the `fixed` parameter and the sample size `n`, but uses the
outer product of gradients instead of the observed information matrix
for the covariance calculation.

For larged values of `fixed` the model fitting may fail to converge.

## References

Murphy, S. A., & Van Der Vaart, A. W. (2000). On Profile Likelihood.
Journal of the American Statistical Association, 95(450), 449–465.
https://doi.org/10.1080/01621459.2000.10474219

Boruvka, A., and Cook, R. J. (2015), A Cox-Aalen Model for
Interval-censored Data. Scand J Statist, 42, 414–426. doi:
10.1111/sjos.12113.

Donglin Zeng, Lu Mao, D. Y. Lin, Maximum likelihood estimation for
semiparametric transformation models with interval-censored data,
Biometrika, Volume 103, Issue 2, June 2016, Pages 253–271,
https://doi.org/10.1093/biomet/asw013

Donglin Zeng, Fei Gao, D. Y. Lin, Maximum likelihood estimation for
semiparametric regression models with multivariate interval-censored
data, Biometrika, Volume 104, Issue 3, September 2017, Pages 505–525,
https://doi.org/10.1093/biomet/asx029
