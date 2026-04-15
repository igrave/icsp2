# Interval censored time from diabetes onset to diabetic nephronpathy

Data set contains interval censored survival time for time from onset of
diabetes to to diabetic nephronpathy. Identical to the `diabetes`
dataset found in the package `glrt`.

## Usage

``` r
IR_diabetes
```

## Format

An object of class `data.frame` with 731 rows and 3 columns.

## Fields

- `left`:

  left side of observation interval

- `right`:

  right side of observation interval

- `gender`:

  gender of subject

## References

Borch-Johnsens, K, Andersen, P and Decker, T (1985). "The effect of
proteinuria on relative mortality in Type I (insulin-dependent) diabetes
mellitus." Diabetologia, 28, 590-596.

## Examples

``` r
 head(IR_diabetes)
#>   left right gender
#> 1   24    27   male
#> 2   22    22 female
#> 3   37    39   male
#> 4   20    20   male
#> 5    1    16   male
#> 6    8    20 female
 diabetes_fit <- ic_sp_po(
     Surv(left, right, type = 'interval2') ~ gender,
     data = IR_diabetes
     )
```
