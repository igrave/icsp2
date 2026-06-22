#' Internal fast log-rank test implementation
#'
#' This is an internal implementation of right-censored logrank test
#' for use within [ic_logrank] where many tests need to be performed
#' efficiently for variance calculation.
#' 
#' Compared to Survival::survdiff, this implementation is faster but does
#' not have the same flexibility or checks.
#'
#' @param time Numeric vector of survival/censoring times.
#' @param event Logical vector of event indicators (TRUE = event, FALSE = censored).
#' @param group Integer vector of group labels, coded as consecutive integers
#'   (i.e. `1, 2, ...`). This is not checked.
#'
#' @return An object of class `"cpp_logrank"`, a list with components:
#' \describe{
#'   \item{observed}{Numeric vector of observed events per group.}
#'   \item{expected}{Numeric vector of expected events per group.}
#'   \item{variance}{Numeric variance matrix}
#'   \item{n_by_group}{Numeric vector of sample sizes per group.}
#' }
#' @keywords internal
#' @examples
#' set.seed(1992)
#' time  <- rexp(100)
#' event <- sample(c(TRUE, FALSE), 100, replace = TRUE)
#' group <- sample(1:2, 100, replace = TRUE)
#' cpp_logrank(time, event, group)
#'
#' @export
cpp_logrank <- function(time, event, group) {
  n <- length(time)
  if (length(event) != n || length(group) != n) {
    stop("`time`, `event`, and `group` must have the same length.")
  }
  if (!is.numeric(time)) {
    stop("`time` must be numeric.")
  }

  if (anyNA(event)) {
    stop("`event` must contain only TRUE/FALSE.")
  }
  event <- as.integer(event)

  group <- as.integer(group)

  result <- fast_logrank(time, event, group)
  class(result) <- "cpp_logrank"
  result
}

#' @exportS3Method print cpp_logrank
print.cpp_logrank <- function(x, ..., digits = 4) {
  K <- length(x$observed)
  cat("\nLog-rank test\n")

  cat("Sample sizes by group:\n")
  names(x$n_by_group) <- seq_len(K)
  print(x$n_by_group)
  cat("\n")

  oe <- data.frame(
    N = x$n_by_group,
    Observed = x$observed,
    Expected = round(x$expected, digits),
    `O-E` = round(x$observed - x$expected, digits),
    check.names = FALSE
  )
  rownames(oe) <- paste0("Group ", seq_len(K))
  print(oe)
  cat("\n")

  cat("Variance-covariance matrix:\n")
  dimnames(x$variance) <- list(seq_len(K), seq_len(K))
  print(round(x$variance, digits))

  invisible(x)
}
