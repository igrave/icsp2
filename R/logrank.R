#' Compute Sun's (1996) Log-Rank Test Statistic
#'
#' @param npmle_fit Output of .fit_ic_sp() with pooled NPMLE (one entry per stratum)
#' @param mi_list Output of find_maximal_intersections() per stratum
#' @param group_var Factor of group assignments (same order as response_mat)
#' @param group_levels Character vector of group level names
#' @param strata Factor of stratum assignments (one level per stratum in npmle_fit)
#'
#' @return List with statistic, score, information (V), var (V^-1)
#' @keywords internal
compute_sun_statistic <- function(
  npmle_fit,
  response_mat,
  group_var,
  strata,
  n_samples = 1000
) {
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
    P_ij <- a_ij <- matrix(0, nrow = n, ncol = m)
    for (i in seq_len(n)) {
      j_i <- seq.int(mi_list[[s]]$l_inds[i] + 1L, mi_list[[s]]$r_inds[i] + 1L)
      a_ij[i, j_i] <- 1
      P_ij[i, j_i] <- theta[j_i]
    }
    d_ij <- P_ij / rowSums(P_ij)

    n_jl <- d_jl <- matrix(0, nrow = n_groups, ncol = m)
    group_var_s <- group_var[strata == strata_levels[s]]
    for (g in seq.int(n_groups)) {
      this_group <- group_levels[g]
      d_jl[g, ] <- colSums(d_ij[group_var_s == this_group, , drop = FALSE])
      n_jl[g, ] <- rev(cumsum(rev(d_jl[g, ])))
    }

    d_j <- colSums(d_jl)
    n_j <- colSums(n_jl)
    T_s[[s]] <- rowSums(d_jl - n_jl * rep(d_j / n_j, each = n_groups))

    H <- n_samples
    # k_groups <- k_groups + 1
    U <- matrix(0, nrow = k_groups, ncol = H)
    V <- array(0, dim = c(k_groups, k_groups, H))
    k_subset <- seq.int(k_groups)

    # Calculate this once (TODO check for stratum)
    S_ij <- matrix(0, nrow = nrow(P_ij), ncol = ncol(P_ij))
    P_row_sums <- rowSums(P_ij)
    for (i in seq_len(nrow(P_ij))) {
      S_ij[i, ] <- cumsum(P_ij[i, ]) / P_row_sums[i]
    }

    for (h in seq.int(H)) {
      q <- runif(nrow(S_ij))
      j <- rowSums(S_ij <= q) + 1L
      imp <- mi_list[[s]]$mi_r[j] #(mi_list[[s]]$mi_l[j] + mi_list[[s]]$mi_r[j]) / 2

      svdf <- survdiff(Surv(imp, rep(1, length(imp))) ~ group_var_s)
      U[, h] <- svdf$obs[k_subset] - svdf$exp[k_subset]
      V[,, h] <- svdf$var[k_subset, k_subset]
    }

    U_s[[s]] <- rowMeans(U) # dims=1 (default): average over H columns -> k-vector

    V_s[[s]] <- rowMeans(V, dims = 2) -
      ((U - U_s[[s]]) %*% t(U - U_s[[s]])) / (H - 1)
  }

  U_all <- Reduce("+", U_s)
  V_all <- Reduce("+", V_s)
  statistic <- as.numeric(t(U_all) %*% solve(V_all) %*% U_all)
  list(T_s = T_s, statistic = statistic, var = V_all, U = U_all, U_strata = U_s)
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
#' Performs Sun's (1996) non-parametric log-rank test for comparing survival
#' distributions across groups with interval-censored data. This test compares
#' group-specific NPMLE survival curves without assuming proportional hazards.
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

  # Extract group variable
  remaining_terms <- attr(terms_obj, "term.labels")

  if (length(remaining_terms) != 1) {
    stop("Formula should have exactly one grouping variable (not in strata)")
  }

  group_var_name <- remaining_terms[1]
  group_var <- mf[[group_var_name]]
  if (!is.factor(group_var)) {
    group_var <- as.factor(group_var)
  }

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

  # Compute Sun's score statistic
  sun_result <- compute_sun_statistic(
    npmle_fit = npmle_fit,
    response_mat = response_mat,
    group_var = group_var,
    strata = npmle_strata,
    n_samples = n_samples
  )

  # Degrees of freedom
  df <- n_groups - 1

  # P-value from chi-squared distribution
  p_value <- sapply(sun_result$statistic, function(stat) {
    pchisq(as.numeric(stat), df = df, lower.tail = FALSE)
  })

  # Create result object
  result <- list(
    logrank = unlist(sun_result$T_s),
    statistic = c(Q = as.numeric(sun_result$statistic)),
    parameter = c(df = df),
    p.value = p_value,
    score = sun_result$score,
    information = sun_result$information,
    var = sun_result$var,
    groups = group_levels,
    n = n_per_group,
    strata = if (!is.null(original_strata)) levels(original_strata) else NULL,
    method = if (!is.null(original_strata)) {
      "Stratified Sun's non-parametric log-rank test for interval-censored data"
    } else {
      "Sun's non-parametric log-rank test for interval-censored data"
    },
    data.name = deparse(substitute(data)),
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

  cat("Log-rank scores by group:\n")
  print(object$logrank)

  cat("\n")
  cat(sprintf(
    "  Chi-squared statistic: Q = %s\n",
    format(object$statistic, digits = digits)
  ))
  cat(sprintf("  Degrees of freedom:    df = %d\n", object$parameter))
  cat(sprintf(
    "  P-value:               p = %s\n",
    format.pval(object$p.value, digits = digits)
  ))

  cat("\n")
  cat("Variance-covariance matrix:\n")
  print(round(object$var, 4))

  invisible(object)
}
