#' Settings for ic_sp_ph
#'
#' @param useGA Should constrained gradient ascent step be used?
#' @param maxIter Maximum iterations
#' @param baseUpdates number of baseline updates (ICM + GA) per iteration
#' @param regStart Initial values for regression parameters
#' @param derivMethod Method for derivative calculations.
#' @param updateCovars Should covariates be updated during fitting?
#'
#'  @description
#' Creates the control options for the \code{ic_sp} function.
#' Defaults not intended to be changed for use in standard analyses.
#'
#' @details
#' The constrained gradient step, actived by \code{useGA = T},
#' is a step that was added to improve the convergence in a special case.
#' The option to turn it off is only in place to help demonstrate it's utility.
#'
#'  \code{regStart} also for seeding of initial value of regression parameters. Intended for use in ``warm start" for bootstrap samples
#'  and providing fixed regression parameters when calculating fit in qq-plots.
#'
#' @export
ic_sp_settings <- function(
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
#' @param settings A list of control settings created by `ic_sp_settings()`.
#' @return A list containing the fitted model information, including coefficients, variance-covariance matrix, and other details.
#' @export
ic_sp2 <- function(
  formula,
  data,
  weights,
  subset,
  na.action,
  B = c(0, 1),
  settings = ic_sp_settings(),
  model = c("ph", "po")
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
    terms <- terms[attr(terms, "factors")[strata_special, ] == 0]
  } else {
    strata <- factor(rep(1, nrow(mf)))
  }
  attr(terms, "intercept") <- 0
  x <- model.matrix(terms, data = mf)

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

  if (!is.null(settings$regStart)) {
    regStart <- settings$regStart
  } else {
    regStart <- rep(0, length(x_names))
  }
  if (length(regStart) != length(x_names)) {
    stop("length of provided regression parameters wrong length")
  }

  other_info <- list(
    useGA = settings$useGA,
    maxIter = settings$maxIter,
    baselineUpdates = settings$baseUpdates,
    useFullHess = TRUE,
    updateCovars = settings$updateCovars,
    recenterCovars = (length(x_names) > 0),
    regStart = regStart,
    derivMethod = settings$derivMethod
  )

  # Recentering covariates
  covarOffset <- colMeans(x)
  x <- x - rep(1, nrow(x)) %*% t(covarOffset)

  fitInfo <- .fit_ic_sp(
    x = x,
    y = response_mat,
    model_type = paste0("ic_", model_type),
    weights = weights,
    strata = strata,
    other_info = other_info
  )

  dataEnv <- list()
  dataEnv[['x']] <- as.matrix(x, nrow = nrow(response_mat))
  if (ncol(dataEnv$x) == 1) {
    colnames(dataEnv[['x']]) <- x_names
  }
  dataEnv[['y']] <- response_mat
  dataEnv[['strata']] <- strata
  dataEnv[['weights']] <- weights

  bsMat <- NULL
  covar <- NULL

  names(fitInfo$coefficients) <- x_names
  fitInfo$covarOffset <- matrix(covarOffset, nrow = 1)
  fitInfo$var <- covar
  fitInfo$call = call
  fitInfo$formula = formula
  fitInfo$.dataEnv <- new.env()
  if (!missing(data)) {
    fitInfo$.dataEnv$data = data
  }
  list2env(dataEnv, envir = fitInfo$.dataEnv)
  fitInfo$par = 'semi-parametric'
  fitInfo$model <- model_type
  fitInfo$reg_pars <- fitInfo$coefficients
  fitInfo$terms <- call_info$mt
  fitInfo$xlevels <- .getXlevels(call_info$mt, call_info$mf)
  if (fitInfo$iterations == settings$maxIter) {
    warning('Maximum iterations reached in ic_sp.')
  }
  fitInfo$other_info <- other_info

  return(fitInfo)
}


ic_sp_ph <- ic_sp_po <- ic_sp2


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
    findMaximalIntersections(y_s[, 1], y_s[, 2])
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
  myFit <- new(model_type)
  myFit$p_hat <- c_ans$p_hat
  myFit$coefficients <- c_ans$coefficients
  myFit$llk <- c_ans$llk
  myFit$iterations <- c_ans$iterations
  myFit$score <- c_ans$score
  myFit[['T_bull_Intervals']] <- lapply(mi_info, function(mi) {
    rbind(mi[['mi_l']], mi[['mi_r']])
  })
  myFit$p_hat <- lapply(myFit$p_hat, function(p) p / sum(p))

  return(myFit)
}


#' Profile Likelihood Covariance for Semi-Parametric Models
#' @param fit Fitted model object from \code{ic_sp}
#' @param constant Multiplier for the constant `h_n` in the profile likelihood.
#' @param ... Additional arguments.
#' @exportS3Method vcov ic_po
vcov.ic_po <- function(object, constant = 1, ...) {
  fit <- object
  if (!inherits(fit, 'ic_ph') && !inherits(fit, 'ic_po')) {
    stop('fit must be an object of class ic_ph or ic_po')
  }

  n <- nrow(fit$.dataEnv$data)
  k <- length(fit$coefficients)
  llk_beta <- fit$llk
  llk_beta_k <- numeric(k)
  llk_beta_j_k <- matrix(NA, nrow = k, ncol = k)

  h <- constant * sqrt(1 / n)
  for (i in seq_len(k)) {
    beta <- fit$coefficients
    beta[i] <- beta[i] + h
    # fit the model with beta_i + h
    other_info <- fit$other_info
    other_info$updateCovars <- FALSE
    other_info$regStart <- beta
    new_fit <- .fit_ic_sp(
      x = fit$.dataEnv$x,
      y = fit$.dataEnv$y,
      model_type = class(fit)[1],
      weights = fit$.dataEnv$weights,
      strata = fit$.dataEnv$strata,
      other_info = other_info
    )
    llk_beta_k[i] <- new_fit$llk
  }

  for (i in seq_len(k)) {
    for (j in seq(from = i, to = k)) {
      beta <- fit$coefficients
      beta[i] <- beta[i] + h
      beta[j] <- beta[j] + h
      # fit the model with beta_i + h and beta_j + h
      other_info <- fit$other_info
      other_info$updateCovars <- FALSE
      other_info$regStart <- beta
      new_fit <- .fit_ic_sp(
        x = fit$.dataEnv$x,
        y = fit$.dataEnv$y,
        model_type = class(fit)[1],
        weights = fit$.dataEnv$weights,
        strata = fit$.dataEnv$strata,
        other_info = other_info
      )
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
