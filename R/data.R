#' Lung Tumor Interval Censored Data from Hoel and Walburg 1972
#' @name miceData
#' @docType data
#' @field l    left side of observation interval
#' @field u    right side of observation interval
#' @field grp  Group for mouse. Either `"ce"` (conventional environment) or `"ge"` (germ-free environment)
#'
#' @description
#' RFM mice were sacrificed and examined for lung tumors. This resulted in current status interval censored data:
#' if the tumor was present, this implied left censoring and if no tumor was present this implied right censoring.
#' Mice were placed in two different groups: conventional environment or germ free environment.
#' @references
#' 	Hoel D. and Walburg, H.,(1972), Statistical analysis of survival experiments, \emph{The Annals of Statistics},
#' 	18, 1259-1294
#' @examples
#' head(miceData)
#' mice_ph_fit <- ic_sp_ph(Surv(l, u, type = 'interval2') ~ grp, data = miceData)
#'
#' summary(mice_ph_fit)
"miceData"

#' Interval censored time from diabetes onset to diabetic nephronpathy
#' @name IR_diabetes
#' @docType data
#'
#' @field left   left side of observation interval
#' @field right  right side of observation interval
#' @field gender gender of subject
#'
#' @description
#' Data set contains interval censored survival time for time from onset of
#' diabetes to to diabetic nephronpathy. Identical to the \code{diabetes}
#' dataset found in the package \code{glrt}.
#' @examples
#'  head(IR_diabetes)
#'  diabetes_fit <- ic_sp_po(
#'      Surv(left, right, type = 'interval2') ~ gender,
#'      data = IR_diabetes
#'      )
#' @references
#'  Borch-Johnsens, K, Andersen, P and Decker, T (1985).
#'  "The effect of proteinuria on relative mortality in Type I (insulin-dependent) diabetes mellitus."
#'  Diabetologia, 28, 590-596.
"IR_diabetes"
