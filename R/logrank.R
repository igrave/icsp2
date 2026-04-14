#' Compute Log-Rank Test Statistic
#'
#' @param npmle_fit Output of .fit_ic_sp() with pooled NPMLE (one entry per stratum)
#' @param response_mat Matrix of response intervals (n x 2)
#' @param group_var Factor of group assignments (same order as response_mat)
#' @param strata Factor of stratum assignments (one level per stratum in npmle_fit)
#' @param n_samples Number of samples for variance imputation
#'
#' @return List with statistics for each stratum and overall
#' @noRd
compute_statistic <- function(
  npmle_fit,
  response_mat,
  group_var,
  strata,
  n_samples = 1000,
  type = c("sas", "hly")
) {
  type <- match.arg(type, c("sas", "hly"))

  mi_list <- by(
    response_mat,
    strata,
    function(y_s) find_maximal_intersections(y_s[, 1], y_s[, 2])
  )
  group_levels <- levels(group_var)
  n_groups <- nlevels(group_var)
  k_groups <- n_groups - 1L
  strata_levels <- levels(strata)

  T_s <- list()
  U_s <- V_s <- list()
  for (s in seq_along(strata_levels)) {
    theta <- npmle_fit$p_hat[[s]]
    m <- length(theta)
    n <- length(mi_list[[s]]$l_inds)
    P_ij <- matrix(0, nrow = n, ncol = m)
    for (i in seq_len(n)) {
      j_i <- seq.int(mi_list[[s]]$l_inds[i] + 1L, mi_list[[s]]$r_inds[i] + 1L)
      P_ij[i, j_i] <- theta[j_i]
    }
    group_var_s <- group_var[strata == strata_levels[s]]

    if (type == "sas") {
      d_ij <- P_ij / rowSums(P_ij)

      n_jl <- d_jl <- matrix(0, nrow = n_groups, ncol = m)

      for (g in seq.int(n_groups)) {
        this_group <- group_levels[g]
        d_jl[g, ] <- colSums(d_ij[group_var_s == this_group, , drop = FALSE])
        n_jl[g, ] <- rev(cumsum(rev(d_jl[g, ])))
      }

      d_j <- colSums(d_jl)
      n_j <- colSums(n_jl)
      T_s[[s]] <- rowSums(d_jl - n_jl * rep(d_j / n_j, each = n_groups))
    }

    H <- n_samples
    U <- matrix(0, nrow = n_groups, ncol = H)
    V <- array(0, dim = c(n_groups, n_groups, H))
    k_subset <- seq.int(k_groups)

    S_ij <- matrix(0, nrow = nrow(P_ij), ncol = ncol(P_ij))
    P_row_sums <- rowSums(P_ij)
    for (i in seq_len(nrow(P_ij))) {
      S_ij[i, ] <- cumsum(P_ij[i, ]) / P_row_sums[i]
    }

    for (h in seq.int(H)) {
      q <- runif(nrow(S_ij))
      j <- rowSums(S_ij <= q) + 1L
      imp <- mi_list[[s]]$mi_r[j] # Does _l or _r make a difference?

      svdf <- survdiff(Surv(imp, rep(1, length(imp))) ~ group_var_s)
      U[, h] <- svdf$obs - svdf$exp
      V[,, h] <- svdf$var
    }

    U_s[[s]] <- rowMeans(U)

    V_s[[s]] <- rowMeans(V, dims = 2) -
      ((U - U_s[[s]]) %*% t(U - U_s[[s]])) / (H - 1)
  }
  if (type == "sas") {
    U_s <- T_s
  }
  U_all <- Reduce("+", U_s)
  V_all <- Reduce("+", V_s)
  statistic <- as.numeric(
    t(U_all[k_subset]) %*% solve(V_all[k_subset, k_subset]) %*% U_all[k_subset]
  )

  U_strata <- do.call(cbind, U_s)
  rownames(U_strata) <- group_levels
  colnames(U_strata) <- strata_levels
  colnames(V_all) <- rownames(V_all) <- names(U_all) <- group_levels
  list(
    U_strata = U_strata,
    statistic = statistic,
    var = V_all,
    U = U_all
  )
}


#' Sun's Log-Rank Test for Interval-Censored Data
#'
#' @param formula A formula with \code{Surv(l, u, type = 'interval2')} response
#'   and a single grouping variable on the right-hand side. May also contain
#'   \code{strata()} terms for stratified analysis.
#' @param data A data frame containing the variables in the formula.
#' @param subset Optional expression indicating which subset of rows to use.
#' @param na.action Function to handle missing values.
#' @param B A vector of length 2 giving bounds for observation times. Default is c(0, 1).
#' @param n_samples The number of "imputation" samples for the variance calculation. Default is 1000.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{ic_logrank} containing:
#'   \item{logrank}{The log-rank statistics "observed - expected" for all groups and strata}
#'   \item{logrank_overall}{The log-rank statistics "observed - expected" for all groups}
#'   \item{statistic}{The overall chi-squared test statistic based on imputation}
#'   \item{df}{Degrees of freedom}
#'   \item{p.value}{P-value from chi-squared distribution}
#'   \item{var}{Variance-covariance matrix I(0)^{-1}}
#'   \item{groups}{Group levels being compared}
#'   \item{n}{Sample sizes per group}
#'   \item{method}{Description of test method}
#'   \item{data.name}{Name of the data}
#'   \item{call}{The matched call}
#'
#' @details
#' Performs Sun's (1996) non-parametric log-rank test for comparing survival
#' distributions across groups with interval-censored data. This test compares
#' group-specific NPMLE survival curves without assuming proportional hazards.
#'
#' The test statistic is \eqn{Q = U(0)' I(0)^{-1} U(0)}, which follows a
#' chi-squared distribution with k-1 degrees of freedom under the null
#' hypothesis, where k is the number of groups.
#'
#' The default test type is "sas", which uses Sun's test statistic combined
#' with variance estimated based on the Huang, Lee and Yu (2008) procedure
#' sampling exact observation times from the Turnball intervals.
#'
#' Alternatively, the "hly" type calculates the test statistic and variance
#' using the multiple imputation approach of Huang, Lee and Yu (2008) directly.
#'
#' Stratified tests are constructed by calculating the \eqn{U} and \eqn{V}
#' matrices for each stratum separately and then summing the stratum-specific
#' matrices to give a global test statistic
#' \eqn{\sum\bar{U}' (\sum\hat{V})^{-1} \sum\bar{U}}.
#' This is the procedure described in the SAS documentation for PROC ICLIFETEST.
#'
#'
#' @references
#' Sun, J. (1996). A non-parametric test for interval-censored failure time
#' data with application to AIDS studies. \emph{Statistics in Medicine},
#' 15(13), 1387-1395.
#'
#' Huang, J., Lee, C., and Yu, Q. (2008). A Generalized Log-Rank Test
#' for Interval-Censored Failure Time Data via Multiple Imputation.
#' \emph{Statistics in Medicine 27:3217–3226}. http://dx.doi.org/10.1002/sim.3211
#'
#' SAS Institute Inc. (2026). \emph{SAS/STAT® 26.03 User's Guide: The ICLIFETEST Procedure}.
#' https://documentation.sas.com/doc/en/statug/latest/statug_iclifetest_details01.htm
#' Accessed 14 April 2026.
#'
#' @examples
#' # Simple two-group comparison
#' data(miceData)
#' ic_logrank(Surv(l, u, type = "interval2") ~ grp, data = miceData)
#'
#' @export
ic_logrank <- function(
  formula,
  data,
  subset,
  na.action,
  B = c(0, 1),
  n_samples = 1000,
  type = c("sas", "hly"),
  ...
) {
  call <- match.call()
  type <- match.arg(type, c("sas", "hly"))

  if (is.numeric(n_samples) && n_samples > 1 && is.finite(n_samples)) {
    n_samples <- as.integer(n_samples)
  } else {
    stop("n_samples should be a positive integer")
  }

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
    original_strata <- interaction(strata_vars, drop = TRUE)
    mf[, strata_special] <- NULL
    terms_obj <- drop.terms(
      terms_obj,
      dropx = strata_special - 1,
      keep.response = TRUE
    )
  } else {
    original_strata <- NULL
  }

  # Extract grouping variable
  if (length(attr(terms_obj, "term.labels")) != 1) {
    stop("Formula should have exactly one grouping variable (not in strata)")
  }

  group_var <- interaction(
    mf[all.vars(delete.response(terms_obj))],
    drop = TRUE
  )

  group_levels <- levels(group_var)
  n_groups <- length(group_levels)

  # Get sample sizes per group
  if (!is.null(original_strata)) {
    n_per_group <- table(Group = group_var, Strata = original_strata)
  } else {
    n_per_group <- table(Group = group_var)
  }
  if (any(n_per_group == 0)) {
    print(n_per_group)
    stop("All groups must have at least one observation")
  }
  # Check weights
  weights <- model.weights(mf)
  if (is.null(weights)) {
    weights <- rep(1, nrow(mf))
  }

  # Strata for pooled NPMLE: use original_strata only (not groups)
  # Under H0 we estimate a single pooled survival curve per stratum
  npmle_strata <- if (!is.null(original_strata)) {
    original_strata
  } else {
    factor(rep(1L, nrow(response_mat)))
  }

  # Fit pooled NPMLE (no covariates) separately within each stratum
  other_info <- list(
    useGA = TRUE,
    maxIter = 10000,
    baselineUpdates = 5,
    useFullHess = FALSE,
    updateCovars = FALSE,
    recenterCovars = FALSE,
    regStart = numeric(0),
    derivMethod = c(12, 1)
  )

  npmle_fit <- .fit_ic_sp(
    x = matrix(nrow = nrow(response_mat), ncol = 0),
    y = response_mat,
    model_type = "ic_ph",
    weights = weights,
    strata = npmle_strata,
    other_info = other_info
  )

  sun_result <- compute_statistic(
    npmle_fit = npmle_fit,
    response_mat = response_mat,
    group_var = group_var,
    strata = npmle_strata,
    n_samples = n_samples,
    type = type
  )

  # Degrees of freedom
  df <- n_groups - 1

  # P-value from chi-squared distribution
  p_value <- sapply(sun_result$statistic, function(stat) {
    pchisq(as.numeric(stat), df = df, lower.tail = FALSE)
  })

  # Create result object
  result <- list(
    logrank = sun_result$U_strata,
    logrank_overall = rowSums(sun_result$U_strata),
    statistic = c(Q = as.numeric(sun_result$statistic)),
    df = c(df = df),
    p.value = p_value,
    var = sun_result$var,
    groups = group_levels,
    n = n_per_group,
    strata = if (!is.null(original_strata)) levels(original_strata) else NULL,
    method = if (!is.null(original_strata)) {
      "Stratified Sun's log-rank test for interval-censored data"
    } else {
      "Sun's log-rank test for interval-censored data"
    },
    data.name = deparse(substitute(data)),
    n_samples = n_samples,
    call = call
  )

  class(result) <- "ic_logrank"
  result
}


#' @exportS3Method print ic_logrank
print.ic_logrank <- function(object, digits = 4, ...) {
  cat("\n")
  cat(object$method, "\n")
  cat(strrep("=", nchar(object$method)), "\n\n")

  cat("Call:\n")
  print(object$call)
  cat("\n")

  cat("Sample sizes by group:\n")
  print(object$n)
  cat("\n")

  cat("Log-rank Statistics:\n")
  print(object$logrank_overall)
  cat("\n")

  cat(sprintf(
    "  Chi-squared statistic: Q = %s\n",
    format(object$statistic, digits = digits)
  ))
  cat(sprintf("  Degrees of freedom:    df = %d\n", object$df))
  cat(sprintf(
    "  P-value:               p = %s\n",
    format.pval(object$p.value, digits = digits)
  ))

  cat("\n")
  cat("Variance-covariance matrix:\n")
  print(round(object$var, 4))
  cat("\n")
  cat("Calculated with ", object$n_samples, " samples\n")
  invisible(object)
}
