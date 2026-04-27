profile_ci <- function(
  object, 
  level = 0.95, 
  tol = 1e-3, 
  search_width = 2.3 / sqrt(-diag(object$hessian)),
  ...
) {
  cutoff <- qchisq(level, df = 1) / 2
  target <- object$llk - cutoff
  k <- length(object$coefficients)


  #for (j in seq_len(k)) {
  ci_list <- mirai_map(seq_len(k), ..., .f = function(j){
    beta_hat <- object$coefficients[j]

    pl_j <- function(val) {
      beta <- object$coefficients
      beta[j] <- val
      icsp2::profile_fit(object, beta)$llk - target
    }

    
    lower <- uniroot(pl_j,
      lower = beta_hat - search_width[j], upper = beta_hat,
      tol = tol
    )$root

    upper <- uniroot(pl_j,
      lower = beta_hat, upper = beta_hat + search_width[j],
      tol = tol
    )$root
    c(lower = lower, upper = upper)
  })
  ci <- do.call(rbind, ci_list[])
  rownames(ci) <- names(object$coefficients)
  ci
  
}
