# Lines method for ic_sp2 objects

Lines method for ic_sp2 objects

## Usage

``` r
# S3 method for class 'ic_sp2'
lines(x, y, type = "s", strata, ...)
```

## Arguments

- x:

  an object of class 'ic_sp2'

- y:

  ignored

- type:

  the type of plot to produce, default is "l" for lines

- strata:

  Optional vector of indices to subset the strata to be plotted.

- ...:

  additional arguments passed to the lines function

## Details

Due to the interval censoring, the survival curves are not unique, so
two lines are drawn for each curve: one for the lower bound and one for
the upper bound of the interval. A pair of lines is drawn for each
stratum, if strata are present in the model. Therefore to unambiguously
specify graphics parameters such as colour or line type, the user should
specify them as a vector of length equal to the number of strata (for
strata specific pars) or twice the number of strata (for line specific
pars, e.g. strata 1 lower, strata 1 upper, strata 2 lower, strata 2
upper).
