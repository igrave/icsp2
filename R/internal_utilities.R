###			SEMIPARAMETRIC UTILITIES

adjust_intervals <- function(B = c(0, 1), surv_matrix, eps = 10^-10) {
  is_censored <- surv_matrix[, 2] - surv_matrix[, 1] > (2 * eps)
  if (B[1] == 0) {
    surv_matrix[is_censored, 1] = surv_matrix[is_censored, 1] + eps
  }
  if (B[2] == 0) {
    surv_matrix[is_censored, 2] = surv_matrix[is_censored, 2] - eps
  }
  surv_matrix
}

find_maximal_intersections <- function(lower, upper) {
  all_vals <- sort(unique(c(lower, upper)))
  is_left <- all_vals %in% lower
  is_right <- all_vals %in% upper
  mi_list <- findMI(all_vals, is_left, is_right, lower, upper)
  names(mi_list) <- c('l_inds', 'r_inds', 'mi_l', 'mi_r')
  mi_list
}


check_weights <- function(model.frame) {
  w <- model.weights(model.frame)
  if (is.null(w)) {
    w <- rep(1, nrow(model.frame))
  } else {
    w <- as.numeric(w)
  }
  if (any(is.na(w))) {
    stop('NAs not allowed in weights')
  }
  if (any(w < 0)) {
    stop('negative weights not allowed')
  }
  return(w)
}


check_matrix <- function(x) {
  if (qr(cbind(x, 1))$rank < ncol(x) + 1) {
    stop("Covariate matrix is rank deficient. Check covariates")
  }
}
