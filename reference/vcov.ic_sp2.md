# Profile Likelihood Covariance for Semi-Parametric Models

Profile Likelihood Covariance for Semi-Parametric Models

## Usage

``` r
# S3 method for class 'ic_sp2'
vcov(object, typical = 1, large = 2, ...)
```

## Arguments

- object:

  Fitted model object from `ic_sp`

- typical:

  A typical value for the regression parameters, used to determine the
  scale of `h_n`. Default is 1.

- large:

  A large value for the regression parameters, used to determine the
  scale of `h_n`. Default is 2.

- ...:

  Unused.

## Value

Variance-covariance matrix of the regression parameters.

## Details

The covariance matrix is calculated using the profile likelihood
approach described in Boruvka and Cook (2014). This method involves
perturbing the regression parameters and refitting the model to estimate
the curvature of the log-likelihood function, which is then used to
compute the covariance matrix. The `typical` and `large` parameters are
used to determine the scale of the perturbations.

## References

Boruvka, A., and Cook, R. J. (2015), A Cox-Aalen Model for
Interval-censored Data. Scand J Statist, 42, 414–426. doi:
10.1111/sjos.12113.
