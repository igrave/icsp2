#' Lines method for ic_sp2 objects
#' @param x an object of class 'ic_sp2'
#' @param y ignored
#' @param type the type of plot to produce, default is "l" for lines
#' @param ... additional arguments passed to the lines function
#' @param strata Optional vector of indicies to subset the strata to be plotted.
#' @export
#'
#' @details
#' Due to the interval censoring, the survival curves are not unique,
#' so two lines are drawn for each curve: one for the lower bound and
#' one for the upper bound of the interval.
#' A pair of lines is drawn for each stratum, if strata are present in the model.
#' Therefore to unambiguously specify graphics parmameters such as colour or
#' line type, the user should specify them as a vector of length equal to the
#' number of strata (for strata specific pars) or twice the number of strata
#' (for line specific pars, ie strata 1 lower, strata 1 upper, strata 2 lower,
#'  strata 2 upper).
#'
lines.ic_sp2 <- function(x, y, type = "s", strata, ...) {
  if (missing(strata)) {
    strata <- seq_along(x$s)
  } else {
    strata <- as.integer(strata)
    stopifnot(strata >= 1 & strata <= length(x$s))
  }
  n_strata <- length(strata)
  pars_dots <- lapply(list(..., type = type), function(par) {
    if (length(par) == 1) {
      rep(par, n_strata * 2)
    } else if (length(par) == n_strata) {
      rep(par, each = 2)
    } else if (length(par) == n_strata * 2) {
      par
    } else {
      rep(par, length = n_strata * 2)
    }
  })

  pars_dots$x <- unlist(
    lapply(x$T_bull_Intervals[strata], function(intervals) {
      list(
        c(intervals[1, 1], intervals[1, ]),
        c(intervals[1, 1], intervals[2, ])
      )
    }),
    recursive = FALSE
  )
  pars_dots$y <- unlist(
    lapply(x$s[strata], function(s) {
      list(s, s)
    }),
    recursive = FALSE
  )
  browser()
  .mapply(lines, pars_dots, list())
  invisible(NULL)
}
