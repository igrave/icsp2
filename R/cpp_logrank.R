#' Fast log-rank test via C++
#'
#' This is an internal implementation for use within [ic_logrank] where
#' many tests need to be performed efficiently for variance calculation.
#'
#' @param time Numeric vector of survival/censoring times.
#' @param event Integer vector of event indicators (1 = event, 0 = censored).
#' @param group Integer vector of group labels, coded as consecutive integers
#'   starting from 1 (i.e. `1, 2, ...`).
#'
#' @return An object of class `"cpp_logrank"`, a list with components:
#' \describe{
#'   \item{observed}{Numeric vector of observed events per group.}
#'   \item{expected}{Numeric vector of expected events per group.}
#'   \item{variance}{Numeric variance matrix}
#'   \item{n_by_group}{Numeric vector of sample sizes per group.}
#' }
#'
#' @examples
#' set.seed(1992)
#' time  <- rexp(100)
#' event <- sample(0:1, 100, replace = TRUE)
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

  event <- as.integer(event)
  if (!all(event %in% c(0L, 1L))) {
    stop("`event` must contain only 0 and 1.")
  }

  group <- as.integer(group)
  K <- max(group)
  if (!identical(sort(unique(group)), seq_len(K))) {
    stop("`group` must be consecutive integers starting from 1.")
  }

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
