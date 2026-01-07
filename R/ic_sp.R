#' Semi-Parametric models for Interval Censored Data
#' 
#' @param formula regression formula. Response must be a \code{Surv} object of type \code{'interval2'}or \code{cbind}. See details.
#' @param data dataset
#' @param model What type of model to fit. Current choices are "\code{ph}" (Cox PH) or "\code{po}" (proportional odds)
#' @param weights Vector of case weights. Not standardized; see details
#' @param bs_samples Number of bootstrap samples used for estimation of standard errors 
#' @param useMCores Should multiple cores be used for bootstrap sample? Does not register cluster (see example)
#' @param B Should intervals be open or closed? See details.
#' @param controls Advanced control options 
#' 
#' @description  	Fits a semi-parametric model for interval censored data. 
#' Can fit either a Cox-PH model or a proportional odds model.  
#'
#' The covariance matrix for the regression coefficients is estimated via bootstrapping. 
#' For large datasets, this can become slow so parallel processing can be used to take advantage of multiple cores via the \code{foreach} package. 
#'
#'@details
#'	Response variable should either be of the form \code{cbind(l, u)} or 
#'	\code{Surv(l, u, type = 'interval2')}, where \code{l} and \code{u} are the lower 
#'	and upper ends of the interval known to contain the event of interest. 
#'	Uncensored data can be included by setting \code{l == u}, right censored data 
#'	can be included by setting \code{u == Inf} or \code{u == NA} and left censored 
#'	data can be included by setting \code{l == 0}.
#'
#' The argument \code{B} determines whether the intervals should be open or closed, 
#' i.e. \code{B = c(0,1)} implies that the event occurs within the interval \code{(l,u]}.
#'  The exception is that if \code{l == u}, it is assumed that the event is uncensored, 
#'  regardless of \code{B}. 
#'
#' In regards to weights, they are not standardized. 
#' This means that if `weight[i] = 2`, this is the equivalent to having two 
#' observations with the same values as subject i. 
#'
#' The algorithm used is inspired by the extended ICM algorithm from Wei Pan 1999.
#' However, it uses a conditional Newton Raphson step (for the regression parameters) 
#' and an ICM step (for the baseline survival parameters), rather than one single
#' ICM step (for both sets). In addition, a gradient ascent can also be used
#' to update the baseline parameters. This step is necessary if the
#' data contains many uncensored observations, very similar to how 
#' the EM algorithm greatly accelerates the ICM algorithm for the NPMLE 
#' (gradient ascent is used rather than the EM, as the M step is not 
#' in closed form for semi-parametric models). 
#'
#' Earlier versions of icenReg used an active set algorithm, which was not
#'  as fast for large datasets.
#'
#' @examples
#' set.seed(1)
#'
#' sim_data <- simIC_weib(n = 100, inspections = 5, inspectLength = 1)
#' ph_fit <- ic_sp(Surv(l, u, type = 'interval2') ~ x1 + x2, 
#'                 data = sim_data)	
#' # Default fits a Cox-PH model
#' 
#' summary(ph_fit)		
#' # Regression estimates close to true 0.5 and -0.5 values
#'
#'
#' new_data <- data.frame(x1 = c(0,1), x2 = c(1, 1) )
#' rownames(new_data) <- c('group 1', 'group 2')
#' plot.icenReg_fit(ph_fit, new_data)
#' # plotting the estimated survival curves
#' 
#' po_fit <- ic_sp(Surv(l, u, type = 'interval2') ~ x1 + x2, 
#'                 data = sim_data, model = 'po')
#' # fits a proportional odds model
#' 
#' summary(po_fit)
#'
#' # Not run: how to set up multiple cores
#' # library(doParallel)
#' # myCluster <- makeCluster(2) 
#' # registerDoParallel(myCluster)
#' # fit <- ic_sp(Surv(l, u, type = 'interval2') ~ x1 + x2,
#' #              data = sim_data, useMCores = TRUE
#' #              bs_samples = 500)
#' # stopCluster(myCluster)
#'
#'
#' @author Clifford Anderson-Bergman
#' @references 
#' Pan, W., (1999), Extending the iterative convex minorant algorithm to the Cox model for interval-censored data, \emph{Journal of Computational and Graphical Statistics}, Vol 8(1), pp109-120
#'
#' Wellner, J. A., and Zhan, Y. (1997) A hybrid algorithm for computation of the maximum likelihood estimator from censored data, \emph{Journal of the  American Statistical Association}, Vol 92, pp945-959
#' 
#' Anderson-Bergman, C. (preprint) Revisiting the iterative convex minorant algorithm for interval censored survival regression models
#' @export
ic_sp <- function(
  formula,
  data,
  model = 'ph', 
  weights = NULL,
  strata = NULL,
  bs_samples = 0,
  useMCores = F, 
  B = c(0,1), 
  controls = makeCtrls_icsp() 
){
  recenterCovars = TRUE
  useFullHess = TRUE  
  
  # Information about orginal call to function. Useful for expanding X in predict(fit, newdata)
  call_base = match.call(expand.dots = FALSE)
  call_info = readingCall(call_base)
  
  if(missing(data)) data <- environment(formula)
  checkFor_cluster(formula)
  
  reg_items = make_xy(formula, data)
  yMat <- reg_items$y
  x <- reg_items$x
  xNames = reg_items$xNames
  
  if(length(xNames) == 0 & bs_samples > 0){
    cat('no covariates included, so bootstrapping is not useful. Setting bs_samples = 0')
    bs_samples = 0
  }
  # For semi-parametric model, need to handle case when l == u
  yMat <- adjustIntervals(B, yMat)
  
  checkMatrix(x)

  if(model == 'ph')	callText = 'ic_ph'
  else if(model == 'po')	callText = 'ic_po'
  else stop('invalid choice of model. Current optios are "ph" (cox ph) or "po" (proportional odds)')
  
  weights <- checkWeights(weights, yMat)	
  strata <- checkStrata(strata, yMat)
  if(length(x) == 0) recenterCovars = FALSE
  
  if(!is.null(controls$regStart)) regStart <- controls$regStart
  else                            regStart <- rep(0, length(xNames)) 
  if(length(regStart) != length(xNames)){
    stop("length of provided regression parameters wrong length")
  }
  
  other_info <- list(useGA = controls$useGA, maxIter = controls$maxIter, 
    baselineUpdates = controls$baseUpdates, 
    useFullHess = useFullHess, 
    updateCovars = controls$updateReg,
    recenterCovars = recenterCovars, 
    regStart = regStart,
    derivMethod = controls$derivMethod)  
    
  # Recentering covariates
  covarOffset <- icColMeans(x)
  x <- t(t(x) - covarOffset)
    
  fitInfo <- fit_ICPH(yMat, x, callText, weights, strata, other_info)
  dataEnv <- list()
  dataEnv[['x']] <- as.matrix(x, nrow = nrow(yMat))
  if(ncol(dataEnv$x) == 1) colnames(dataEnv[['x']]) <- xNames
  dataEnv[['y']] <- yMat
  dataEnv[['strata']] <- strata
  dataEnv[['weights']] <- weights
  seeds = as.integer( runif(bs_samples, 0, 2^31) )
  bsMat <- numeric()
  if(useMCores) `%mydo%` <- `%dopar%`
  else          `%mydo%` <- `%do%`
  i <- NULL  #this is only to trick R CMD check, 
  #as it does not recognize the foreach syntax
  if(bs_samples > 0){
    bsMat <- foreach(i = seeds, .combine = 'rbind') %mydo%{
      set.seed(i)
      sampDataEnv <- bs_sampleData(dataEnv, weights)
      ans <- getBS_coef(sampDataEnv, callText = callText, other_info = other_info)
      rm(sampDataEnv)
      return(ans)
    }
  }
      
  if(bs_samples > 0){
    names(fitInfo$coefficients) <- xNames
    colnames(bsMat) <- xNames
    incompleteIndicator <- is.na(bsMat[,1])
    numNA <- sum(incompleteIndicator)
    if(numNA > 0){
      if(numNA / length(incompleteIndicator) >= 0.1)
        cat('warning: ', numNA,
          ' bootstrap samples (out of ', bs_samples, 
          ') were dropped due to singular covariate matrix.',
          'Likely due to very sparse covariate. Be wary of these results.\n', sep = '')
      bsMat <- bsMat[!incompleteIndicator,,drop = F]
    }
    covar <- cov(bsMat)
  } else { 
    bsMat <- NULL
    covar <- NULL
  }
      
  names(fitInfo$coefficients) <- xNames
  fitInfo$covarOffset <- matrix(covarOffset, nrow = 1)
  fitInfo$bsMat <- bsMat
  fitInfo$var <- covar
  fitInfo$call = call_base
  fitInfo$formula = formula
  fitInfo$.dataEnv <- new.env()
  if(!missing(data)){ fitInfo$.dataEnv$data = data }
  list2env(dataEnv, envir = fitInfo$.dataEnv)
  fitInfo$par = 'semi-parametric'
  fitInfo$model = model
  fitInfo$reg_pars <- fitInfo$coefficients
  fitInfo$terms <- call_info$mt
  fitInfo$xlevels <- .getXlevels(call_info$mt, call_info$mf)
  if(fitInfo$iterations == controls$maxIter){
    warning('Maximum iterations reached in ic_sp.')
  }
  fitInfo$other_info <- other_info
  return(fitInfo)
}
    
#' Control Parameters for ic_sp
#' 
#' @param useGA Should constrained gradient ascent step be used?
#' @param maxIter Maximum iterations
#' @param baseUpdates number of baseline updates (ICM + GA) per iteration
#' @param regStart Initial values for regression parameters
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
#' @author Clifford Anderson-Bergman
#' @export
makeCtrls_icsp <- function(useGA = T, maxIter = 10000, baseUpdates = 5,
      regStart = NULL, derivMethod = 1){
  ans <- list(useGA = useGA, maxIter = maxIter, 
  baseUpdates = baseUpdates, 
  regStart = regStart, updateReg = TRUE, derivMethod = derivMethod)
  return(ans)
}
        
        
fit_ICPH <- function(obsMat, covars, callText = 'ic_ph', weights, strata, other_info){
  if(any(obsMat[,1] > obsMat[,2])) 
  stop("left side of response interval greater than right side. This is impossible.")
  useGA <- other_info$useGA
  maxIter <- other_info$maxIter
  baselineUpdates <- other_info$baselineUpdates
  useFullHess <- other_info$useFullHess
  updateCovars <- other_info$updateCovars
  regStart <- other_info$regStart
  derivMethod <- other_info$derivMethod
  # recenterCovars = FALSE
  # if(getNumCovars(covars) == 0)	recenterCovars <- FALSE
  mi_info <- by(obsMat, strata, function(x) findMaximalIntersections(x[,1], x[,2]))
  #mi_info <- findMaximalIntersections(obsMat[,1], obsMat[,2])
  # k = length(mi_info[['mi_l']])
  covars_list <- lapply(split(seq_len(nrow(obsMat)), strata), function(i) covars[i, , drop = FALSE])
  #covars <- by(as.matrix(covars), strata, I)
  weights <- split(as.numeric(weights), strata)
  
  if(callText == 'ic_ph'){fitType = as.integer(1)}
  else if(callText == 'ic_po'){fitType = as.integer(2)}
  else {stop('callText not recognized in fit_ICPH')}
  
  # if(recenterCovars){
  #   pca_info <- prcomp(covars, scale. = TRUE)
  #   covars <- as.matrix(pca_info$x)
  #   regStart <- solve(pca_info$rotation, (regStart * pca_info$scale) )
  # }
  
  linds <- lapply(mi_info, function(x) x$l_inds)
  rinds <- lapply(mi_info, function(x) x$r_inds)

  c_ans <- .Call(
    'ic_sp_ch',
    linds,
    rinds,
    covars_list, # list covariates of each strata
    fitType,
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
  myFit <- new(callText)
  myFit$p_hat <- c_ans$p_hat
  myFit$coefficients <- c_ans$coefficients
  myFit$llk <- c_ans$llk
  myFit$iterations <- c_ans$iterations
  myFit$score <- c_ans$score
  myFit[['T_bull_Intervals']] <- lapply(mi_info, function(mi) rbind(mi[['mi_l']], mi[['mi_r']]))
  myFit$p_hat <- lapply(myFit$p_hat, function(p) p / sum(p)) 
  # if(recenterCovars == TRUE){
  #   myFit$pca_coefs <- myFit$coefficients
  #   myFit$pca_info <- pca_info
  #   myFit$coefficients <- as.numeric( myFit$pca_info$rotation %*% myFit$coefficients) / myFit$pca_info$scale		
  #   myFit$baseOffset = as.numeric(myFit$coefficients %*% myFit$pca_info$center)
  # }
  return(myFit)
}


#' Profile Likelihood Covariance for Semi-Parametric Models
#' @param fit Fitted model object from \code{ic_sp}
#' @param constant Multiplier for the constant `h_n` in the profile likelihood.
#' @param ... Additional arguments.
#' 
profile_likelihood_covariance <- function(fit, constant = 1, ...){
  if(!inherits(fit, 'ic_ph') && !inherits(fit, 'ic_po')){
    stop('fit must be an object of class ic_ph or ic_po')
  }
  
  n <- nrow(fit$.dataEnv$data)
  k <- length(fit$coefficients)
  llk_beta <- fit$llk
  llk_beta_k <- numeric(k)
  llk_beta_j_k <- matrix(NA, nrow = k, ncol = k)

  h <- sqrt(constant / n)
  for (i in seq_len(k)){
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
    for (j in seq(from = i, to = k)){
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
      llk_beta_j_k[i,j] <- new_fit$llk
      if (i != j){
        llk_beta_j_k[j,i] <- new_fit$llk
      }
    }
  }
  inv_cov <- matrix(NA, nrow = k, ncol = k)
  for (i in seq_len(k)){
    for (j in seq(from = i, to = k)){
      inv_cov[i,j] <- (llk_beta - llk_beta_k[i] - llk_beta_k[j] + llk_beta_j_k[i,j]) / (h^2)
      if (i != j){
        inv_cov[j,i] <- inv_cov[i,j]
      }
    }
  }
  result <- list()
  result$inv_cov <- inv_cov
  result$cov <- -solve(inv_cov)
  result
}