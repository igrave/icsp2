#' Get Survival Curve Estimates from icenReg Model
#'
#' @param fit      model fit with \code{ic_par} or \code{ic_sp}
#' @param newdata  \code{data.frame} containing covariates
#' @param p        Percentiles
#' @param q        Quantiles
#' @description
#' Gets probability or quantile estimates from a \code{ic_sp}, \code{ic_par} or \code{ic_bayes} object.
#' Provided estimates conditional on regression parameters found in \code{newdata}.
#' @details
#' For the \code{ic_sp} and \code{ic_par}, the MLE estimate is returned. For \code{ic_bayes},
#' the MAP estimate is returned. To compute the posterior means, use \code{sampleSurv}.
#'
#' If \code{newdata} is left blank, baseline estimates will be returned (i.e. all covariates = 0).
#' If \code{p} is provided, will return the estimated \code{F^{-1}(p | x)}. If \code{q} is provided,
#' will return the estimated \code{F(q | x)}. If neither \code{p} nor \code{q} are provided,
#' the estimated conditional median is returned.
#'
#' In the case of \code{ic_sp}, the MLE of the baseline survival is not necessarily unique,
#' as probability mass is assigned to disjoint Turnbull intervals, but the likelihood function is
#' indifferent to how probability mass is assigned within these intervals. In order to have a well
#' defined estimate returned, we assume probability is assigned uniformly in these intervals.
#' In otherwords, we return *a* maximum likelihood estimate, but don't attempt to characterize *all* maximum
#' likelihood estimates with this function. If that is desired, all the information needed can be
#' extracted with \code{getSCurves}.
#' @examples
#' simdata <- simIC_weib(n = 500, b1 = .3, b2 = -.3,
#' inspections = 6, inspectLength = 1)
#' fit <- ic_par(Surv(l, u, type = 'interval2') ~ x1 + x2,
#'              data = simdata)
#' new_data <- data.frame(x1 = c(1,2), x2 = c(-1,1))
#' rownames(new_data) <- c('grp1', 'grp2')
#'
#' estQ <- getFitEsts(fit, new_data, p = c(.25, .75))
#'
#' estP <- getFitEsts(fit, q = 400)
#' @author Clifford Anderson-Bergman
#' @export
getFitEsts <- function(fit, newdata = NULL, p, q) {
  if (is.null(newdata)) {} else {
    if (!identical(newdata, "midValues")) {
      if (!is(newdata, "data.frame")) {
        stop("newdata should be a data.frame")
      }
    }
  }
  etas <- get_etas(fit, newdata)

  if (missing(p)) {
    p <- NULL
  }
  if (missing(q)) {
    q <- NULL
  }
  if (!is.null(q)) {
    xs <- q
    type = 'q'
  } else {
    type = 'p'
    if (is.null(p)) {
      xs <- 0.5
    } else {
      xs <- p
    }
  }

  if (length(etas) == 1) {
    etas <- rep(etas, length(xs))
  }
  if (length(xs) == 1) {
    xs <- rep(xs, length(etas))
  }
  if (length(etas) != length(xs)) {
    stop('length of p or q must match nrow(newdata) OR be 1')
  }

  regMod <- fit$model

  if (inherits(fit, 'sp_fit')) {
    scurves <- getSCurves(fit, newdata = NULL)
    if (length(scurves) > 1) {
      stop("Models with stratification factors are not supported")
    }
    baselineInfo <- list(
      tb_ints = scurves[[1]]$Tbull_ints,
      s = scurves[[1]]$S_curves$baseline
    )
    baseMod = 'sp'
  }
  if (inherits(fit, 'par_fit') | inherits(fit, 'bayes_fit')) {
    baseMod <- fit$par
    baselineInfo <- fit$baseline
  }

  if (fit$model == 'po' | fit$model == 'ph' | fit$model == 'none') {
    surv_etas <- etas
    scale_etas <- rep(1, length(etas))
  } else if (fit$model == 'aft') {
    scale_etas <- etas
    surv_etas <- rep(1, length(etas))
  } else {
    stop('model not recognized in getFitEsts')
  }

  if (type == 'q') {
    ans <- getSurvProbs(
      xs / scale_etas,
      surv_etas,
      baselineInfo = baselineInfo,
      regMod = regMod,
      baseMod = baseMod
    )
    #    ans <- ans * scale_etas
    return(ans)
  } else if (type == 'p') {
    #    xs <- xs / scale_etas
    ans <- getSurvTimes(
      xs,
      surv_etas,
      baselineInfo = baselineInfo,
      regMod = regMod,
      baseMod = baseMod
    )
    ans <- ans * scale_etas
    return(ans)
  }
}


#' Predictions from icenReg Regression Model
#'
#' @param object   Model fit with \code{ic_par} or \code{ic_sp}
#' @param type     type of prediction. Options include \code{"lp", "response"}
#' @param newdata  \code{data.frame} containing covariates
#' @param ...      other arguments (will be ignored)
#'
#' @description
#' Gets various estimates from an \code{ic_np}, \code{ic_sp} or \code{ic_par} object.
#' @details
#' If \code{newdata} is left blank, will provide estimates for original data set.
#'
#' For the argument \code{type}, there are two options. \code{"lp"} provides the
#' linear predictor for each subject (i.e. in a proportional hazards model,
#' this is the log-hazards ratio, in proportional odds, the log proporitonal odds),
#' \code{"response"} provides the median response value for each subject,
#' *conditional on that subject's covariates, but ignoring their actual response interval*.
#' Use \code{imputeCens} to impute the censored values.
#' @examples
#' simdata <- simIC_weib(n = 500, b1 = .3, b2 = -.3,
#'                       inspections = 6,
#'                       inspectLength = 1)
#'
#' fit <- ic_par(cbind(l, u) ~ x1 + x2,
#'               data = simdata)
#'
#' imputedValues <- predict(fit)
#' @author Clifford Anderson-Bergman
#' @export
predict.icsp2 <- function(
  object,
  type = 'response',
  newdata = NULL,
  ...
) {
  stop("unimplemented")
  #imputeOptions = fullSample, fixedParSample, median
  if (is.null(newdata)) {
    newdata <- object$getRawData()
  }
  if (type == 'lp') {
    return(log(get_etas(object, newdata = newdata)))
  }
  if (type == 'response') {
    return(getFitEsts(fit = object, newdata = newdata, p = 0.5))
  }
  stop('"type" not recognized: options are "lp", "response" and "impute"')
}

#' Impute Interval Censored Data from icenReg Regression Model
#'
#' @param fit         icenReg model fit
#' @param newdata     \code{data.frame} containing covariates and censored intervals. If blank, will use data from model
#' @param imputeType  type of imputation. See details for options
#' @param samples  Number of imputations (ignored if \code{imputeType = "median"})
#'
#' @description
#' Imputes censored responses from data.
#' @details
#'  If \code{newdata} is left blank, will provide estimates for original data set.
#'
#'  There are several options for how to impute. \code{imputeType = 'median'}
#'  imputes the median time, conditional on the response interval, covariates and
#'  regression parameters at the MLE. To get random imputations without accounting
#'  for error in the estimated parameters \code{imputeType ='fixedParSample'} takes a
#'  random sample of the response variable, conditional on the response interval,
#'  covariates and estimated parameters at the MLE. Finally,
#'  \code{imputeType = 'fullSample'} first takes a random sample of the coefficients,
#'  (assuming asymptotic normality) and then takes a random sample
#'  of the response variable, conditional on the response interval,
#'  covariates, and the random sample of the coefficients.
#'
#' @examples
#' simdata <- simIC_weib(n = 500)
#'
#' fit <- ic_par(cbind(l, u) ~ x1 + x2,
#'               data = simdata)
#'
#' imputedValues <- imputeCens(fit)
#' @author Clifford Anderson-Bergman
#' @export
imputeCens <- function(
  fit,
  newdata = NULL,
  imputeType = 'fullSample',
  samples = 5
) {
  if (is.null(newdata)) {
    newdata <- fit$getRawData()
  }
  yMat <- expandY(fit$formula, newdata, fit)
  p1 <- getFitEsts(fit, newdata, q = as.numeric(yMat[, 1]))
  p2 <- getFitEsts(fit, newdata, q = as.numeric(yMat[, 2]))
  ans <- matrix(nrow = length(p1), ncol = samples)
  storage.mode(ans) <- 'double'
  if (imputeType == 'median') {
    isBayes <- is(fit, 'bayes_fit')
    if (isBayes) {
      orgCoefs <- getSamplablePars(fit)
      map_ests <- c(fit$MAP_baseline, fit$MAP_reg_pars)
      setSamplablePars(fit, map_ests)
    }
    p_med <- (p1 + p2) / 2
    ans <- getFitEsts(fit, newdata, p = p_med)
    isLow <- ans < yMat[, 1]
    ans[isLow] <- yMat[isLow, 1]
    isHi <- ans > yMat[, 2]
    ans[isHi] <- yMat[isHi]
    if (isBayes) {
      setSamplablePars(fit, orgCoefs)
    }
    #    rownames(ans) <- rownames(newdata)
    return(ans)
  }
  if (imputeType == 'fixedParSample') {
    isBayes <- is(fit, 'bayes_fit')
    if (isBayes) {
      orgCoefs <- getSamplablePars(fit)
      map_ests <- c(fit$MAP_baseline, fit$MAP_reg_pars)
      setSamplablePars(fit, map_ests)
    }
    for (i in 1:samples) {
      p_samp <- runif(length(p1), p1, p2)
      theseImputes <- getFitEsts(fit, newdata, p = p_samp)
      isLow <- theseImputes < yMat[, 1]
      theseImputes[isLow] <- yMat[isLow, 1]
      isHi <- theseImputes > yMat[, 2]
      theseImputes[isHi] <- yMat[isHi, 2]
      rownames(ans) <- rownames(newdata)
      ans <- fastMatrixInsert(theseImputes, ans, colNum = i)
    }
    if (isBayes) {
      setSamplablePars(fit, orgCoefs)
    }
    rownames(ans) <- rownames(newdata)
    return(ans)
  }
  if (imputeType == 'fullSample') {
    isSP <- is(fit, 'sp_fit')
    isBayes <- is(fit, 'bayes_fit')
    for (i in 1:samples) {
      orgCoefs <- getSamplablePars(fit)
      if (isBayes) {
        sampledCoefs <- sampBayesPar(fit)
      } else if (!isSP) {
        coefVar <- getSamplableVar(fit)
        sampledCoefs <- sampPars(orgCoefs, coefVar)
      } else {
        sampledCoefs <- getBSParSample(fit)
      }
      setSamplablePars(fit, sampledCoefs)
      p1 <- getFitEsts(fit, newdata, q = as.numeric(yMat[, 1]))
      p2 <- getFitEsts(fit, newdata, q = as.numeric(yMat[, 2]))
      p_samp <- runif(length(p1), p1, p2)
      theseImputes <- getFitEsts(fit, newdata, p = p_samp)
      isLow <- theseImputes < yMat[, 1]
      theseImputes[isLow] <- yMat[isLow, 1]
      isHi <- theseImputes > yMat[, 2]
      theseImputes[isHi] <- yMat[isHi, 2]
      fastMatrixInsert(theseImputes, ans, colNum = i)
      setSamplablePars(fit, orgCoefs)
    }
    rownames(ans) <- rownames(newdata)
    return(ans)
  }
  stop('imputeType type not recognized.')
}

#' Draw samples from an icenReg model
#'
#' @param fit         icenReg model fit
#' @param newdata     \code{data.frame} containing covariates. If blank, will use data from model
#' @param sampleType  type of samples See details for options
#' @param samples     Number of samples
#'
#' @description
#' Samples response values from an icenReg fit conditional on covariates, but not
#' censoring intervals. To draw response samples conditional on covariates and
#' restrained to intervals, see \code{imputeCens}.
#'
#'
#' @details
#'  Returns a matrix of samples. Each row of the matrix corresponds with a subject with the
#'  covariates of the corresponding row of \code{newdata}. For each column of the matrix,
#'  the same sampled parameters are used to sample response variables.
#'
#'  If \code{newdata} is left blank, will provide estimates for original data set.
#'
#'  There are several options for how to sample. To get random samples without accounting
#'  for error in the estimated parameters \code{imputeType ='fixedParSample'} takes a
#'  random sample of the response variable, conditional on the response interval,
#'  covariates and estimated parameters at the MLE. Alternatively,
#'  \code{imputeType = 'fullSample'} first takes a random sample of the coefficients,
#'  (assuming asymptotic normality for the ic_par) and then takes a random sample
#'  of the response variable, conditional on the response interval,
#'  covariates, and the random sample of the coefficients.
#'
#' @examples
#' simdata <- simIC_weib(n = 500)
#'
#' fit <- ic_par(cbind(l, u) ~ x1 + x2,
#'               data = simdata)
#'
#' newdata = data.frame(x1 = c(0, 1), x2 = c(1,1))
#'
#' sampleResponses <- ic_sample(fit, newdata = newdata, samples = 100)
#' @author Clifford Anderson-Bergman
#' @export
ic_sample <- function(
  fit,
  newdata = NULL,
  sampleType = 'fullSample',
  samples = 5
) {
  if (is.null(newdata)) {
    newdata <- fit$getRawData()
  }
  yMat <- cbind(rep(-Inf, nrow(newdata)), rep(Inf, nrow(newdata)))
  #  p1 <- getFitEsts(fit, newdata, q = as.numeric(yMat[,1]) )
  #  p2 <- getFitEsts(fit, newdata, q = as.numeric(yMat[,2]) )
  nRow = nrow(newdata)
  p1 = rep(0, nrow(newdata))
  p2 = rep(1, nrow(newdata))
  ans <- matrix(nrow = length(p1), ncol = samples)
  storage.mode(ans) <- 'double'
  isSP <- is(fit, 'sp_fit')
  isBayes <- is(fit, 'bayes_fit')
  if (sampleType == 'fixedParSample') {
    if (isBayes) {
      orgCoefs <- getSamplablePars(fit)
      map_ests <- c(fit$MAP_baseline, fit$MAP_reg_pars)
      setSamplablePars(fit, map_ests)
    }
    for (i in 1:samples) {
      p_samp <- runif(length(p1), p1, p2)
      theseImputes <- getFitEsts(fit, newdata, p = p_samp)
      isLow <- theseImputes < yMat[, 1]
      theseImputes[isLow] <- yMat[isLow, 1]
      isHi <- theseImputes > yMat[, 2]
      theseImputes[isHi] <- yMat[isHi, 2]
      ans[, i] <- fastMatrixInsert(theseImputes, ans, colNum = i)
    }
    if (isBayes) {
      setSamplablePars(fit, orgCoefs)
    }
    rownames(ans) <- rownames(newdata)
    return(ans)
  }
  if (sampleType == 'fullSample') {
    for (i in 1:samples) {
      orgCoefs <- getSamplablePars(fit)
      if (isBayes) {
        sampledCoefs <- sampBayesPar(fit)
      } else if (!isSP) {
        coefVar <- getSamplableVar(fit)
        sampledCoefs <- sampPars(orgCoefs, coefVar)
      } else {
        sampledCoefs <- getBSParSample(fit)
      }
      setSamplablePars(fit, sampledCoefs)
      #      p1 <- getFitEsts(fit, newdata, q = as.numeric(yMat[,1]) )
      #      p2 <- getFitEsts(fit, newdata, q = as.numeric(yMat[,2]) )
      p_samp <- runif(length(p1), p1, p2)
      theseImputes <- getFitEsts(fit, newdata, p = p_samp)
      isLow <- theseImputes < yMat[, 1]
      theseImputes[isLow] <- yMat[isLow, 1]
      isHi <- theseImputes > yMat[, 2]
      theseImputes[isHi] <- yMat[isHi, 2]
      ans[, i] <- theseImputes
      #      fastMatrixInsert(theseImputes, ans, colNum = i)
      setSamplablePars(fit, orgCoefs)
    }
    rownames(ans) <- rownames(newdata)
    return(ans)
  }
  stop('sampleType type not recognized.')
}
