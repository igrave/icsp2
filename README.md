# icsp2
Semiparametric interval censored survival regression with stratification

```r
# Install 'icsp2' from R-universe:
install.packages('icsp2', repos = c('https://igrave.r-universe.dev', 'https://cloud.r-project.org'))
```


This is a cut down package to calculate Cox models for interval
censored data and perform a log-rank type test. It is especially designed
to mirror the typical clinical trial survival analysis methods that would
be done for right censored data. Therefore it could be used as a 
sensitivity analysis in the presence of interval censoring.

This package borrows heavily from [`icenReg`](https://cran.r-project.org/package=icenReg),
specifically the `ic_sp()` implementation. The main change is the addition of 
support for stratified baseline hazards and analytical derivatives.
We also implement a profile likelihood variance estimation procedure for the
regression coefficients.

The logrank test is our own and aims to replace a 
[PROC ICLIFETEST](https://documentation.sas.com/doc/en/statug/latest/statug_iclifetest_overview.htm)
test using Sun type weights and imputed variance.
