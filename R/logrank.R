#' Sun's Log-Rank Test for Interval-Censored Data
#'
#' @param formula A formula with \code{Surv(l, u, type = 'interval2')} response
#'   and a single grouping variable on the right-hand side. May also contain
#'   \code{strata()} terms for stratified analysis.
#' @param data A data frame containing the variables in the formula.
#' @param subset Optional expression indicating which subset of rows to use.
#' @param na.action Function to handle missing values.
#' @param B A vector of length 2 giving bounds for observation times. Default is c(0, 1).
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{ic_logrank} containing:
#'   \item{statistic}{The test statistic Q}
#'   \item{parameter}{Degrees of freedom}
#'   \item{p.value}{P-value from chi-squared distribution}
#'   \item{score}{Score vector U(0)}
#'   \item{information}{Information matrix I(0)}
#'   \item{var}{Variance-covariance matrix I(0)^{-1}}
#'   \item{groups}{Group levels being compared}
#'   \item{n}{Sample sizes per group}
#'   \item{method}{Description of test method}
#'   \item{data.name}{Name of the data}
#'   \item{call}{The matched call}
#'
#' @details
#' Performs Sun's (1996) log-rank test for comparing survival distributions
#' across groups with interval-censored data. This is a score test from the
#' proportional hazards model under the null hypothesis that all groups have
#' the same survival distribution (β = 0).
#'
#' The test assumes proportional hazards across groups. If \code{strata()}
#' terms are included, separate baseline hazards are estimated for each
#' stratum, but the group comparison is made across all strata assuming
#' a common proportional hazards effect.
#'
#' The test statistic is Q = U(0)' I(0)^{-1} U(0), which follows a
#' chi-squared distribution with k-1 degrees of freedom under the null
#' hypothesis, where k is the number of groups.
#'
#' @references
#' Sun, J. (1996). A non-parametric test for interval-censored failure time
#' data with application to AIDS studies. \emph{Statistics in Medicine},
#' 15(13), 1387-1395.
#'
#' @examples
#' \dontrun{
#' # Simple two-group comparison
#' data(miceData)
#' ic_logrank(Surv(l, r, type = "interval2") ~ grp, data = miceData)
#'
#' # Stratified analysis
#' ic_logrank(Surv(l, r, type = "interval2") ~ trt + strata(site), data = mydata)
#' }
#'
#' @export
ic_logrank <- function(
  formula,
  data,
  subset,
  na.action,
  B = c(0, 1),
  ...
) {
  call <- match.call()

  # Check that formula has the correct structure
  if (length(formula) != 3) {
    stop("Formula must have left and right hand sides")
  }

  # Set up model frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  # Extract terms and check for Surv response
  terms_obj <- terms(
    formula,
    data = data,
    specials = c("strata", "cluster", "Surv")
  )

  if (
    is.null(attr(terms_obj, "specials")$Surv) ||
      attr(terms_obj, "specials")$Surv != 1
  ) {
    stop("Response must be a Surv object")
  }

  if (attr(mf[[1]], "type") != "interval") {
    stop("Response must be Surv(time1, time2, type = 'interval2')")
  }

  # Extract response
  response_mat <- as.matrix(model.response(mf))
  response_mat[response_mat[, 3] == 0, 2] <- Inf
  response_mat <- adjust_intervals(B, response_mat[, 1:2])

  # Check for cluster terms (not supported)
  if (!is.null(attr(terms_obj, "specials")$cluster)) {
    stop("ic_logrank() does not support cluster terms")
  }

  # Extract strata if present
  strata_special <- attr(terms_obj, "specials")$strata
  if (!is.null(strata_special)) {
    strata_vars <- mf[, strata_special, drop = FALSE]
    strata <- interaction(strata_vars, drop = TRUE)
    mf[, strata_special] <- NULL
    terms_obj <- drop.terms(
      terms_obj,
      dropx = strata_special - 1,
      keep.response = TRUE
    )
  } else {
    strata <- factor(rep(1, nrow(mf)))
  }

  # Extract covariates (should be group variable)
  attr(terms_obj, "intercept") <- 1
  x <- model.matrix(terms_obj, data = mf)[, -1, drop = FALSE]

  if (ncol(x) == 0) {
    stop("Formula must include a grouping variable on the right-hand side")
  }

  # Get original group variable (before dummy coding)
  # After strata terms are removed, we should have just the group variable
  # in the model matrix
  remaining_terms <- attr(terms_obj, "term.labels")

  if (length(remaining_terms) != 1) {
    stop("Formula should have exactly one grouping variable (not in strata)")
  }

  group_var_name <- remaining_terms[1]

  # Get the original group variable from the model frame
  # It should be the first non-response column after strata removal
  group_var <- mf[[group_var_name]]
  if (!is.factor(group_var)) {
    group_var <- as.factor(group_var)
  }

  group_levels <- levels(group_var)
  n_groups <- length(group_levels)

  # Get sample sizes per group
  n_per_group <- table(group_var)

  # Check weights
  weights <- model.weights(mf)
  if (is.null(weights)) {
    weights <- rep(1, nrow(mf))
  }

  # For log-rank test, fit at null hypothesis (β = 0) with optimized baseline
  # The C++ code now computes the score U(0) with the optimized baseline
  k <- ncol(x) # Number of dummy variables (n_groups - 1)

  other_info <- list(
    useGA = TRUE,
    maxIter = 10000,
    baselineUpdates = 5,
    useFullHess = TRUE,
    updateCovars = FALSE, # Keep beta fixed at 0
    recenterCovars = FALSE, # Don't recenter for score test
    regStart = rep(0, k),
    derivMethod = c(12, 1)
  )

  # Fit at null to get score vector U(0)
  null_fit <- .fit_ic_sp(
    x = x,
    y = response_mat,
    model_type = "ic_ph",
    weights = weights,
    strata = strata,
    other_info = other_info
  )

  # Extract score vector U(0) and information matrix I(0) from C++
  # Both computed analytically via calcAnalyticRegDervs()
  score <- null_fit$score
  information <- -null_fit$hessian # I(0) = -H(0) where H is the Hessian

  # Variance-covariance matrix is inverse of information
  var_matrix <- solve(information)

  # Compute test statistic: Q = U' I^{-1} U
  test_stat <- as.numeric(t(score) %*% var_matrix %*% score)

  # Degrees of freedom
  df <- k

  # P-value from chi-squared distribution
  p_value <- pchisq(test_stat, df = df, lower.tail = FALSE)

  # Create result object
  result <- list(
    statistic = c(Q = test_stat),
    parameter = c(df = df),
    p.value = p_value,
    score = score,
    information = information,
    var = var_matrix,
    groups = group_levels,
    n = n_per_group,
    strata = if (nlevels(strata) > 1) levels(strata) else NULL,
    method = if (nlevels(strata) > 1) {
      "Sun's stratified log-rank test for interval-censored data"
    } else {
      "Sun's log-rank test for interval-censored data"
    },
    data.name = deparse(substitute(data)),
    call = call
  )

  class(result) <- "ic_logrank"
  result
}


#' @exportS3Method print ic_logrank
print.ic_logrank <- function(x, digits = 4, ...) {
  cat("\n")
  cat(x$method, "\n\n")
  cat("data: ", x$data.name, "\n")

  if (!is.null(x$strata)) {
    cat("strata: ", paste(x$strata, collapse = ", "), "\n")
  }

  cat("\n")
  cat("Groups:\n")
  print(
    data.frame(
      Group = names(x$n),
      N = as.numeric(x$n)
    ),
    row.names = FALSE
  )

  cat("\n")
  cat("Test statistic:\n")
  cat(sprintf(
    "  Q = %.4f, df = %d, p-value = %.4g\n",
    x$statistic,
    x$parameter,
    x$p.value
  ))
  cat("\n")

  invisible(x)
}


#' @exportS3Method summary ic_logrank
summary.ic_logrank <- function(object, ...) {
  cat("\n")
  cat(object$method, "\n")
  cat(strrep("=", nchar(object$method)), "\n\n")

  cat("Call:\n")
  print(object$call)
  cat("\n")

  if (!is.null(object$strata)) {
    cat(
      "Stratification variables:",
      paste(object$strata, collapse = ", "),
      "\n\n"
    )
  }

  cat("Sample sizes by group:\n")
  print(
    data.frame(
      Group = names(object$n),
      N = as.numeric(object$n),
      Proportion = sprintf("%.1f%%", 100 * as.numeric(object$n) / sum(object$n))
    ),
    row.names = FALSE
  )

  cat("\n")
  cat("Test Results:\n")
  cat(sprintf("  Chi-squared statistic: Q = %.4f\n", object$statistic))
  cat(sprintf("  Degrees of freedom:    df = %d\n", object$parameter))
  cat(sprintf("  P-value:               p = %.4g\n", object$p.value))

  if (object$p.value < 0.001) {
    cat("  ***\n")
  } else if (object$p.value < 0.01) {
    cat("  **\n")
  } else if (object$p.value < 0.05) {
    cat("  *\n")
  } else if (object$p.value < 0.1) {
    cat("  .\n")
  }

  cat("\n")
  cat("Score vector U(0):\n")
  print(round(object$score, 4))

  cat("\n")
  cat("Information matrix I(0):\n")
  print(round(object$information, 4))

  cat("\n")
  cat("Variance-covariance matrix:\n")
  print(round(object$var, 4))

  cat("\n")
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")

  invisible(object)
}
