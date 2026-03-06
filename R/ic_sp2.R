#' Control ic_sp2 fitting options
#'
#' @param useGA Should constrained gradient ascent step be used?
#' @param maxIter Maximum iterations
#' @param baseUpdates number of baseline updates (ICM + GA) per iteration
#' @param regStart Initial values for regression parameters
#' @param derivMethod Method for derivative calculations.
#' @param updateCovars Should covariates be updated during fitting?
#'
#'  @description
#' Creates the control options for [ic_sp2()].
#' Defaults not intended to be changed for use in standard analyses.
#'
#' @details
#' The constrained gradient step, actived by \code{useGA = TRUE},
#' is a step that was added to improve the convergence in a special case.
#' The option to turn it off is only in place to help demonstrate it's utility.
#'
#'  \code{regStart} also for seeding of initial value of regression parameters.
#'  Intended for use in ``warm start" for bootstrap samples
#'  and providing fixed regression parameters when calculating fit in qq-plots.
#'
#' @export
ic_sp_control <- function(
  useGA = TRUE,
  maxIter = 10000,
  baseUpdates = 5,
  regStart = NULL,
  derivMethod = c(12, 1),
  updateCovars = TRUE
) {
  ans <- list(
    useGA = useGA,
    maxIter = maxIter,
    baseUpdates = baseUpdates,
    regStart = regStart,
    updateCovars = updateCovars,
    derivMethod = derivMethod
  )
  return(ans)
}


#' Fit a semi-parametric model for interval-censored data
#' @param formula A model formula with Surv(l, u, type = 'interval2') response and covariates on the right-hand side. May also contain `strata()` terms.
#' @param data A data frame containing the variables in the formula, including strata terms.
#' @param weights Optional vector of weights for each observation, or the name of a variable in `data` containing the weights.
#' @param B A vector of length 2 giving the lower and upper bounds for the observation times. Default is c(0, 1).
#' @param control A list of control settings, with defaults created by [ic_sp_control()].
#' @return A list containing the fitted model information, including coefficients, variance-covariance matrix, and other details.
#' @export
ic_sp2 <- function(
  formula,
  data,
  weights,
  subset,
  na.action,
  B = c(0, 1),
  control = ic_sp_control(...),
  model = c("ph", "po"),
  ...
) {
  # Information about orginal call to function. Useful for expanding X in predict(fit, newdata)
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(
    c("formula", "data", "subset", "weights", "na.action"),
    names(mf),
    0L
  )
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  terms <- terms(
    formula,
    data = data,
    specials = c("strata", "cluster", "Surv")
  )

  if (attr(terms, "specials")$Surv != 1) {
    stop("Response must be a Surv object.")
  } else if (attr(mf[[1]], "type") != "interval") {
    stop("Response must Surv(time1, time2, type = \"interval2\").")
  } else {
    response_mat <- as.matrix(model.response(mf))
    response_mat[response_mat[, 3] == 0, 2] <- Inf
    response_mat <- adjust_intervals(B, response_mat[, 1:2])
  }

  if (!is.null(attr(terms, "specials")$cluster)) {
    stop("ic_sp_ph() does not support cluster terms.")
  }
  strata_special <- attr(terms, "specials")$strata
  if (!is.null(strata_special)) {
    strata_vars <- mf[, strata_special, drop = FALSE]
    strata <- interaction(strata_vars, drop = TRUE)
    mf[, strata_special] <- NULL
    terms <- drop.terms(terms, dropx = strata_special - 1, keep.response = TRUE)
  } else {
    strata <- factor(rep(1, nrow(mf)))
  }

  attr(terms, "intercept") <- 1
  x <- model.matrix(terms, data = mf)[, -1, drop = FALSE]

  call_info = readingCall(call) # TODO do we need this?

  x_names = colnames(x)
  check_matrix(x)

  model_type <- if (call[[1]] == "ic_sp_ph") {
    "ph"
  } else if (call[[1]] == "ic_sp_po") {
    "po"
  } else {
    match.arg(model, c("ph", "po"), several.ok = FALSE)
  }

  weights <- check_weights(mf)

  if (!is.null(control$regStart)) {
    regStart <- control$regStart
  } else {
    regStart <- rep(0, length(x_names))
  }
  if (length(regStart) != length(x_names)) {
    stop("length of provided regression parameters wrong length")
  }

  other_info <- list(
    useGA = control$useGA,
    maxIter = control$maxIter,
    baselineUpdates = control$baseUpdates,
    useFullHess = TRUE,
    updateCovars = control$updateCovars,
    recenterCovars = (length(x_names) > 0),
    regStart = regStart,
    derivMethod = control$derivMethod
  )

  # Recentering covariates
  covarOffset <- colMeans(x)
  x <- x - rep(1, nrow(x)) %*% t(covarOffset)

  result <- .fit_ic_sp(
    x = x,
    y = response_mat,
    model_type = paste0("ic_", model_type),
    weights = weights,
    strata = strata,
    other_info = other_info
  )

  dataEnv <- list()
  dataEnv[["x"]] <- as.matrix(x, nrow = nrow(response_mat))
  if (ncol(dataEnv$x) == 1) {
    colnames(dataEnv[["x"]]) <- x_names
  }
  dataEnv[["y"]] <- response_mat
  dataEnv[["strata"]] <- strata
  dataEnv[["weights"]] <- weights

  bsMat <- NULL
  covar <- NULL

  names(result$coefficients) <- x_names
  result$covarOffset <- matrix(covarOffset, nrow = 1)
  result$var <- covar
  result$call <- call
  result$formula <- formula
  result$.dataEnv <- new.env()
  if (!missing(data)) {
    result$.dataEnv$data = data
  }
  list2env(dataEnv, envir = result$.dataEnv)
  result$par = "semi-parametric"
  result$model <- model_type
  result$reg_pars <- result$coefficients
  result$terms <- terms
  result$xlevels <- .getXlevels(terms, mf)
  if (result$iterations == control$maxIter) {
    warning("Maximum iterations reached.")
  }
  result$other_info <- other_info
  class(result) <- c(paste0("ic_sp_", model_type), "ic_sp2")
  result
}

#' @rdname ic_sp2
#' @export
ic_sp_ph <- ic_sp2

#' @rdname ic_sp2
#' @export
ic_sp_po <- ic_sp2

#' Fitter function for interval censored semi-parametric models
#'
#' @param y Response matrix from `Surv(..., type = "interval2")`
#' @param x Covariate matrix
#' @param weights Observation weights
#' @param strata Strata factor
#' @param model_type Model type: "ic_ph" or "ic_po"
#' @param other_info List of other fitting options
#' @return Fitted model object
#' @export
#' @details
#' For advanced use only. This function is called internally by `ic_sp_ph`
#'  and `ic_sp_po`.
#'
.fit_ic_sp <- function(
  x,
  y,
  weights,
  strata,
  model_type,
  other_info
) {
  # obsMat,
  # covars,
  # callText = 'ic_ph',
  # weights,
  # strata,
  # other_info

  if (any(y[, 1] > y[, 2])) {
    stop(
      "left side of response interval greater than right side. This is impossible."
    )
  }
  useGA <- other_info$useGA
  maxIter <- other_info$maxIter
  baselineUpdates <- other_info$baselineUpdates
  useFullHess <- other_info$useFullHess
  updateCovars <- other_info$updateCovars
  regStart <- other_info$regStart
  derivMethod <- other_info$derivMethod

  mi_info <- by(y, strata, function(y_s) {
    find_maximal_intersections(y_s[, 1], y_s[, 2])
  })

  covars_list <- lapply(split(seq_len(nrow(y)), strata), function(i) {
    x[i, , drop = FALSE]
  })

  weights <- split(as.numeric(weights), strata)

  model_type_int <- switch(
    model_type,
    ic_ph = 1L,
    ic_po = 2L,
    stop("model_type not recognized in fit_ic_sp")
  )

  linds <- lapply(mi_info, function(x) x$l_inds)
  rinds <- lapply(mi_info, function(x) x$r_inds)

  c_ans <- .Call(
    'ic_sp_ch',
    linds,
    rinds,
    covars_list, # list covariates of each strata
    model_type_int,
    weights, # list of weights
    nlevels(strata), # number of strata
    useGA,
    as.integer(maxIter),
    as.integer(baselineUpdates),
    as.logical(useFullHess),
    as.logical(updateCovars),
    as.double(regStart),
    as.integer(derivMethod)
  )
  names(c_ans) <- c('p_hat', 'coefficients', 'llk', 'iterations', 'score')
  result <- list() #new(model_type)
  result$p_hat <- lapply(c_ans$p_hat, function(p) p / sum(p))
  result$s <- lapply(result$p_hat, function(p) 1 - c(0, cumsum(p)))
  result$coefficients <- c_ans$coefficients
  result$llk <- c_ans$llk
  result$iterations <- c_ans$iterations
  result$score <- c_ans$score
  result[['intervals']] <- lapply(mi_info, function(mi) {
    rbind(mi[['mi_l']], mi[['mi_r']])
  })

  return(result)
}


#' Profile Likelihood Covariance for Semi-Parametric Models
#' @param fit Fitted model object from \code{ic_sp}
#' @param constant Multiplier for the constant `h_n` in the profile likelihood.
#'   Either length 1 or length equal to number of regression parameters.
#' @param ... Additional arguments.
#' @exportS3Method vcov ic_sp2
vcov.ic_sp2 <- function(object, constant = 1, ...) {
  fit <- object
  if (!inherits(fit, "ic_sp_ph") && !inherits(fit, "ic_sp_po")) {
    stop("Fit must be an object of class ic_sp_ph or ic_sp_po.")
  }

  n <- nrow(fit$.dataEnv$data)
  k <- length(fit$coefficients)
  if (length(constant) == 1) {
    # constant <- rep(constant, k)[1]
  } else if (length(constant) != k) {
    stop(
      'Constant must be length 1 or length equal to number of regression parameters'
    )
  }
  llk_beta <- fit$llk
  llk_beta_k <- numeric(k)
  llk_beta_j_k <- matrix(NA, nrow = k, ncol = k)

  h <- constant * sqrt(1 / n)

  call_args <- list(
    other_info = fit$other_info,
    x = fit$.dataEnv$x,
    y = fit$.dataEnv$y,
    model_type = paste0("ic_", fit$model),
    weights = fit$.dataEnv$weights,
    strata = fit$.dataEnv$strata
  )
  call_args$other_info$updateCovars <- FALSE

  for (i in seq_len(k)) {
    # fit the model with beta_i + h
    beta <- fit$coefficients
    beta[i] <- beta[i] + h
    call_args$other_info$regStart <- beta
    new_fit <- do.call(.fit_ic_sp, call_args)
    llk_beta_k[i] <- new_fit$llk
  }

  for (i in seq_len(k)) {
    for (j in seq(from = i, to = k)) {
      # fit the model with beta_i + h_i and beta_j + h_j
      beta <- fit$coefficients
      beta[i] <- beta[i] + h
      beta[j] <- beta[j] + h
      call_args$other_info$regStart <- beta
      new_fit <- do.call(.fit_ic_sp, call_args)
      llk_beta_j_k[i, j] <- new_fit$llk
      if (i != j) {
        llk_beta_j_k[j, i] <- new_fit$llk
      }
    }
  }
  inv_cov <- matrix(NA, nrow = k, ncol = k)
  for (i in seq_len(k)) {
    for (j in seq(from = i, to = k)) {
      inv_cov[i, j] <- (llk_beta -
        llk_beta_k[i] -
        llk_beta_k[j] +
        llk_beta_j_k[i, j]) /
        (h^2)
      if (i != j) {
        inv_cov[j, i] <- inv_cov[i, j]
      }
    }
  }
  result <- list()
  result$inv_cov <- inv_cov
  result <- -solve(inv_cov)

  result
}

profile_fit <- function(object, beta = object$coefficients) {
  call_args <- list(
    other_info = object$other_info,
    x = object$.dataEnv$x,
    y = object$.dataEnv$y,
    model_type = paste0("ic_", object$model),
    weights = object$.dataEnv$weights,
    strata = object$.dataEnv$strata
  )
  call_args$other_info$updateCovars <- FALSE
  call_args$other_info$regStart <- beta
  new_fit <- do.call(.fit_ic_sp, call_args)
  return(new_fit)
}

#' Print method for ic_sp2 objects
#' @param x Fitted model object from \code{ic_sp}
#' @param ... Additional arguments.
#' @exportS3Method print ic_sp2
print.ic_sp2 <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Coefficients:\n")
  # printCoefmat(x$coefficients)
  print(x$coefficients)
  cat("\n")
  cat(paste0("Log-likelihood: ", round(x$llk, 4), "\n"))
  cat(paste0("Number of iterations: ", x$iterations, "\n"))
}
