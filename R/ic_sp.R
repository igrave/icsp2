# fit_ICPH <- function(
#   obsMat,
#   covars,
#   callText = 'ic_ph',
#   weights,
#   strata,
#   other_info
# ) {
#   if (any(obsMat[, 1] > obsMat[, 2])) {
#     stop(
#       "left side of response interval greater than right side. This is impossible."
#     )
#   }
#   useGA <- other_info$useGA
#   maxIter <- other_info$maxIter
#   baselineUpdates <- other_info$baselineUpdates
#   useFullHess <- other_info$useFullHess
#   updateCovars <- other_info$updateCovars
#   regStart <- other_info$regStart
#   derivMethod <- other_info$derivMethod
#   # recenterCovars = FALSE
#   # if(getNumCovars(covars) == 0)	recenterCovars <- FALSE
#   mi_info <- by(obsMat, strata, function(x) {
#     findMaximalIntersections(x[, 1], x[, 2])
#   })
#   #mi_info <- findMaximalIntersections(obsMat[,1], obsMat[,2])
#   # k = length(mi_info[['mi_l']])
#   covars_list <- lapply(split(seq_len(nrow(obsMat)), strata), function(i) {
#     covars[i, , drop = FALSE]
#   })
#   #covars <- by(as.matrix(covars), strata, I)
#   weights <- split(as.numeric(weights), strata)

#   if (callText == 'ic_ph') {
#     fitType = as.integer(1)
#   } else if (callText == 'ic_po') {
#     fitType = as.integer(2)
#   } else {
#     stop('callText not recognized in fit_ICPH')
#   }

#   # if(recenterCovars){
#   #   pca_info <- prcomp(covars, scale. = TRUE)
#   #   covars <- as.matrix(pca_info$x)
#   #   regStart <- solve(pca_info$rotation, (regStart * pca_info$scale) )
#   # }

#   linds <- lapply(mi_info, function(x) x$l_inds)
#   rinds <- lapply(mi_info, function(x) x$r_inds)

#   c_ans <- .Call(
#     'ic_sp_ch',
#     linds,
#     rinds,
#     covars_list, # list covariates of each strata
#     fitType,
#     weights, # list of weights
#     nlevels(strata), # number of strata
#     useGA,
#     as.integer(maxIter),
#     as.integer(baselineUpdates),
#     as.logical(useFullHess),
#     as.logical(updateCovars),
#     as.double(regStart),
#     as.integer(derivMethod)
#   )
#   names(c_ans) <- c('p_hat', 'coefficients', 'llk', 'iterations', 'score')
#   myFit <- new(callText)
#   myFit$p_hat <- c_ans$p_hat
#   myFit$coefficients <- c_ans$coefficients
#   myFit$llk <- c_ans$llk
#   myFit$iterations <- c_ans$iterations
#   myFit$score <- c_ans$score
#   myFit[['T_bull_Intervals']] <- lapply(mi_info, function(mi) {
#     rbind(mi[['mi_l']], mi[['mi_r']])
#   })
#   myFit$p_hat <- lapply(myFit$p_hat, function(p) p / sum(p))
#   # if(recenterCovars == TRUE){
#   #   myFit$pca_coefs <- myFit$coefficients
#   #   myFit$pca_info <- pca_info
#   #   myFit$coefficients <- as.numeric( myFit$pca_info$rotation %*% myFit$coefficients) / myFit$pca_info$scale
#   #   myFit$baseOffset = as.numeric(myFit$coefficients %*% myFit$pca_info$center)
#   # }
#   return(myFit)
# }

#' Profile Likelihood Covariance for Semi-Parametric Models
#' @param fit Fitted model object from \code{ic_sp}
#' @param constant Multiplier for the constant `h_n` in the profile likelihood.
#' @param ... Additional arguments.
#'
profile_likelihood_covariance <- function(fit, constant = 1, ...) {
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
    new_fit <- fit_ICPH(
      fit$.dataEnv$y,
      fit$.dataEnv$x,
      callText = class(fit)[1],
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
      new_fit <- fit_ICPH(
        fit$.dataEnv$y,
        fit$.dataEnv$x,
        callText = class(fit)[1],
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
  result$cov <- -solve(inv_cov)
  result
}
