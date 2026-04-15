# Lung Tumor Interval Censored Data from Hoel and Walburg 1972

RFM mice were sacrificed and examined for lung tumors. This resulted in
current status interval censored data: if the tumor was present, this
implied left censoring and if no tumor was present this implied right
censoring. Mice were placed in two different groups: conventional
environment or germ free environment.

## Usage

``` r
miceData
```

## Format

An object of class `data.frame` with 144 rows and 3 columns.

## Fields

- `l`:

  left side of observation interval

- `u`:

  right side of observation interval

- `grp`:

  Group for mouse. Either `"ce"` (conventional environment) or `"ge"`
  (germ-free environment)

## References

    Hoel D. and Walburg, H.,(1972), Statistical analysis of survival experiments, \emph{The Annals of Statistics},
    18, 1259-1294

## Examples

``` r
head(miceData)
#>   l   u grp
#> 1 0 381  ce
#> 2 0 477  ce
#> 3 0 485  ce
#> 4 0 515  ce
#> 5 0 539  ce
#> 6 0 563  ce
mice_ph_fit <- ic_sp_ph(Surv(l, u, type = 'interval2') ~ grp, data = miceData)

summary(mice_ph_fit)
#> Call:
#> ic_sp_ph(formula = Surv(l, u, type = "interval2") ~ grp, data = miceData)
#> 
#> Coefficients:
#>       Estimate Std.Error z.value p.value
#> grpge   0.6785    0.3537  1.9184   0.055
#> 
#> Log-likelihood: -76.56894 
#> Iterations: 28 
```
