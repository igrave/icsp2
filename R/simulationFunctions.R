#' Simulates Event Time from Survival Regression Model
#' 
#' @description Simulates event times from proportional hazards,
#' proportional odds and accelerated failure time models. 
#' 
#' @param linPred    A numeric vector of linear predictors for
#' simulated data
#' @param model      Regression model type. Options are \code{"ph"}, 
#' \code{"po"} or \code{"aft"} 
#' @param dist       The baseline distrubtion \code{q} function
#' @param paramList  A list of parameters to be passed to the baseline
#' distribution function 
#' @noRd
simEventTime <- function(linPred = 0, model = 'ph', 
                         dist = qweibull, 
                         paramList = list(shape = 1, scale = 1)){
  n = length(linPred)
  rawP <- runif(n)
  nu <- exp(linPred)
  if(model == 'aft'){
    paramList$p <- rawP
    rawTimes <- do.call(dist, paramList)
    ans <- rawTimes * nu
    return(ans)
    
  }
  if (model == "ph") 
    adjFun <- function(x, nu) {
      1 - x^(1/nu)
    }
  else if (model == "po") 
    adjFun <- function(x, nu) {
      1 - x * (1/nu)/(x * 1/nu - x + 1)
    }
  adjP <- adjFun(rawP, nu)
  paramList$p <- adjP
  ans <- do.call(dist, paramList)
  return(ans)
}

#' Simulates interval censored data from regression model with a Weibull baseline
#' 
#' @param n Number of samples simulated
#' @param b1 Value of first regression coefficient
#' @param b2 Value of second regression coefficient
#' @param model Type of regression model. Options are 'po' (prop. odds) and 'ph' (Cox PH)
#' @param shape shape parameter of baseline distribution
#' @param scale scale parameter of baseline distribution
#' @param inspections number of inspections times of censoring process
#' @param inspectLength max length of inspection interval
#' @param rndDigits number of digits to which the inspection time is rounded to, 
#' creating a discrete inspection time. If \code{rndDigits = NULL}, the inspection time is not rounded, 
#' resulting in a continuous inspection time
#' @param prob_cen probability event being censored. If event is uncensored, l == u
#'
#' @description
#' Simulates interval censored data from a regression model with a weibull baseline distribution. Used for demonstration
#' @details 
#' Exact event times are simulated according to regression model: covariate \code{x1} 
#' is distributed \code{rnorm(n)} and covariate \code{x2} is distributed
#' \code{1 - 2 * rbinom(n, 1, 0.5)}. Event times are then censored with a 
#' case II interval censoring mechanism with \code{inspections} different inspection times. 
#' Time between inspections is distributed as \code{runif(min = 0, max = inspectLength)}. 
#' Note that the user should be careful in simulation studies not to simulate data 
#' where nearly all the data is right censored (or more over, all the data with x2 = 1 or -1) 
#' or this can result in degenerate solutions!
#' 
#' @examples 
#' set.seed(1)
#' sim_data <- simIC_weib(n = 500, b1 = .3, b2 = -.3, model = 'ph', 
#'                       shape = 2, scale = 2, inspections = 6, 
#'                       inspectLength = 1)
#' #simulates data from a cox-ph with beta weibull distribution.
#'
#' diag_covar(Surv(l, u, type = 'interval2') ~ x1 + x2, 
#'            data = sim_data, model = 'po')
#' diag_covar(Surv(l, u, type = 'interval2') ~ x1 + x2,
#'            data = sim_data, model = 'ph')
#'
#' #'ph' fit looks better than 'po'; the difference between the transformed survival
#' #function looks more constant
#' @author Clifford Anderson-Bergman
#' @export
simIC_weib <- function (n = 100, b1 = 0.5, b2 = -0.5, model = "ph", shape = 2, 
          scale = 2, inspections = 2, inspectLength = 2.5, rndDigits = NULL, 
          prob_cen = 1) 
{
  x1 <- runif(n, -1, 1)
  x2 <- 1 - 2 * rbinom(n, 1, 0.5)
  linPred <- x1 * b1 + x2 * b2
  trueTimes <- simEventTime(linPred, model = model, dist = qweibull, 
                            paramList = list(shape = shape, scale = scale))
  obsTimes <- runif(n = n, max = inspectLength)
  if (!is.null(rndDigits)) 
    obsTimes <- round(obsTimes, rndDigits)
  l <- rep(0, n)
  u <- rep(0, n)
  caught <- trueTimes < obsTimes
  u[caught] <- obsTimes[caught]
  l[!caught] <- obsTimes[!caught]
  if (inspections > 1) {
    for (i in 2:inspections) {
      oldObsTimes <- obsTimes
      obsTimes <- oldObsTimes + runif(n, max = inspectLength)
      if (!is.null(rndDigits)) 
        obsTimes <- round(obsTimes, rndDigits)
      caught <- trueTimes >= oldObsTimes & trueTimes < 
        obsTimes
      needsCatch <- trueTimes > obsTimes
      u[caught] <- obsTimes[caught]
      l[needsCatch] <- obsTimes[needsCatch]
    }
  }
  else {
    needsCatch <- !caught
  }
  u[needsCatch] <- Inf
  if (sum(l > u) > 0) 
    stop("warning: l > u! Bug in code")
  isCensored <- rbinom(n = n, size = 1, prob = prob_cen) == 
    1
  l[!isCensored] <- trueTimes[!isCensored]
  u[!isCensored] <- trueTimes[!isCensored]
  if (sum(l == Inf) > 0) {
    allTimes <- c(l, u)
    allFiniteTimes <- allTimes[allTimes < Inf]
    maxFiniteTime <- max(allFiniteTimes)
    l[l == Inf] <- maxFiniteTime
  }
  return(data.frame(l = l, u = u, x1 = x1, x2 = x2))
}


#' Simulate Current Status Data
#'
#' @description Simulates current status data from a survival regression model 
#' with a Weibull baseline distribution. 
#'
#' @param n       Number of observations
#' @param b1      Regression coefficient 1
#' @param b2      Regression coefficient 2
#' @param model   Regression model to use. Choices are \code{"ph"}, \code{"po"} or \code{"aft"}
#' @param shape   Baseline shape parameter
#' @param scale   Baseline scale parameter
#' 
#' @details Exact event times are simulated according to the given survival regression model.
#' Two covariates are used; \code{x1 = rnorm(n), x2 = 1 - 2 * rbinom(n, 1, .5)}. After
#' event times are simulated, current status inspection times are simulated following the 
#' exact same conditional distribution as event time (so each event time necessarily has 
#' probability 0.5 of being right censored). 
#' 
#' Returns data in current status format, i.e. inspection time and event indicator. 
#' Use \code{cs2ic} to convert to interval censored format (see example).
#' 
#' @examples 
#' simData <- simCS_weib()
#' fit <- ic_par(cs2ic(time, event) ~ x1 + x2, data = simData)
#' @export
simCS_weib <- function (n = 100, b1 = 0.5, b2 = -0.5, model = "ph", shape = 2, 
                        scale = 2) 
{
  x1 <- runif(n, -1, 1)
  x2 <- 1 - 2 * rbinom(n, 1, 0.5)
  linPred <- x1 * b1 + x2 * b2
  eventTimes <- simEventTime(linPred, model = model, dist = qweibull, 
                            paramList = list(shape = shape, scale = scale))
  inspectTimes <- simEventTime(linPred, model = model, dist = qweibull, 
                            paramList = list(shape = shape, scale = scale))
  event <- eventTimes < inspectTimes
  ans <- data.frame(time = inspectTimes, event = event, x1 = x1, x2 = x2)
  return(ans)
}

#' Simulate Doubly Censored Data
#'
#' @description Simulates doubly censored data from a survival regression model 
#' with a Weibull baseline distribution. 
#'
#' @param n       Number of observations
#' @param b1          Regression coefficient 1
#' @param b2          Regression coefficient 2
#' @param model       Regression model to use. Choices are \code{"ph"}, \code{"po"} or \code{"aft"}
#' @param shape       Baseline shape parameter
#' @param scale       Baseline scale parameter
#' @param lowerLimit  Lower censoring threshold
#' @param upperLimit  Upper censoring threshold
#' @details Exact event times are simulated according to the given survival regression model.
#' Two covariates are used; \code{x1 = rnorm(n), x2 = 1 - 2 * rbinom(n, 1, .5)}. After
#' event times are simulated, all values less than \code{lowerLimit} are left censored
#' and all values less than \code{upperLimit} are right censored.  
#' @examples 
#' simData <- simDC_weib()
#' fit <- ic_par(cbind(l, u) ~ x1 + x2, data = simData)
#' @export
simDC_weib <- function (n = 100, b1 = 0.5, b2 = -0.5, model = "ph", shape = 2, 
                                scale = 2, lowerLimit = 0.75, upperLimit = 2){
  x1 <- runif(n, -1, 1)
  x2 <- 1 - 2 * rbinom(n, 1, 0.5)
  linPred <- x1 * b1 + x2 * b2
  trueTimes <- simEventTime(linPred, model = model, dist = qweibull, 
                          paramList = list(shape = shape, scale = scale))
  l <- numeric()
  u <- numeric()
  isLeftCen <- trueTimes < lowerLimit
  isRightCen <- trueTimes > upperLimit
  isUnCen <- !(isLeftCen | isRightCen)
  l[isLeftCen] <- 0
  u[isLeftCen] <- lowerLimit
  l[isRightCen] <- upperLimit
  u[isRightCen] <- Inf
  l[isUnCen] <- trueTimes[isUnCen]
  u[isUnCen] <- trueTimes[isUnCen]
  ans <- data.frame(l = l, u = u, x1 = x1, x2 = x2)
  return(ans)
}