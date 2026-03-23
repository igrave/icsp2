#' @exportS3Method print summary_ic_sp2
print.summary_ic_sp2 <- function(
  x,
  digits = max(3L, getOption("digits") - 3L),
  ...
) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  printCoefmat(x$coef_matrix, digits = digits)
  cat("\nLog-likelihood:", x$llk, "\n")
  cat("Iterations:", x$iterations, "\n")
  invisible(x)
}

#' @exportS3Method summary ic_sp2
summary.ic_sp2 <- function(object, ...) {
  object$var <- vcov(object)
  object$coef_matrix <- cbind(
    Estimate = object$coefficients,
    Std.Error = sqrt(diag(object$var)),
    z.value = object$coefficients / sqrt(diag(object$var)),
    p.value = 2 * pnorm(-abs(object$coefficients / sqrt(diag(object$var))))
  )

  class(object) <- c("summary_ic_sp2", class(object))
  object
}
