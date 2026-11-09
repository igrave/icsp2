test_that("PH model works for sim data with exact obs", {
  set.seed(1951)
  n <- 500
  sim_data <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = 1 - 2 * rbinom(n, 1, 0.5),
    ic = rbinom(n, 1, 0.5)
  )
  sim_data$time <- rexp(
    n,
    rate = exp(0.3 * sim_data$x1 - 0.3 * sim_data$x2) / 2
  ) +
    1
  sim_data$l <- ifelse(sim_data$ic == 0, sim_data$time, floor(sim_data$time))
  sim_data$u <- ifelse(sim_data$ic == 0, sim_data$time, ceiling(sim_data$time))
  n

  result_1 <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data,
    control = ic_sp_control(derivMethod = 1)
  )
  result_12 <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data
  )

  result_icr <- icenReg::ic_sp(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data
  )

  expect_equal(result_1$coefficients, result_12$coefficients, tolerance = 1e-6)
  expect_equal(result_1$llk, result_icr$llk, tolerance = 1e-7)
  expect_equal(result_1$intervals[[1]], result_icr$T_bull_Intervals)

  expect_equal(result_1$coefficients, result_icr$coefficients, tolerance = 1e-6)
  expect_equal(result_1$llk, result_icr$llk, tolerance = 1e-7)
  expect_equal(result_1$intervals[[1]], result_icr$T_bull_Intervals)
})


test_that("PH model works for sim data", {
  set.seed(1951)
  sim_data <- simIC_weib(n = 500, inspections = 3, inspectLength = 1)

  result <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data,
    control = ic_sp_control(derivMethod = 1)
  )

  icr_result <- icenReg::ic_sp(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data
  )

  expect_equal(result$coefficients, icr_result$coefficients, tolerance = 1e-6)
  expect_equal(result$llk, icr_result$llk, tolerance = 1e-7)
  expect_equal(result$intervals[[1]], icr_result$T_bull_Intervals)
})


test_that("PH model works for miceData", {
  data(miceData)

  result <- ic_sp_ph(Surv(l, u, type = 'interval2') ~ grp, data = miceData)

  icr_result <- icenReg::ic_sp(
    Surv(l, u, type = 'interval2') ~ grp,
    data = miceData
  )
  expect_equal(result$coefficients, icr_result$coefficients, tolerance = 1e-6)

  expect_equal(result$intervals[[1]], icr_result$T_bull_Intervals)
})


test_that("PH model works for miceData with numerical derivatives", {
  data(miceData)

  result <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ grp,
    data = miceData,
    control = ic_sp_control(derivMethod = 1)
  )

  icr_result <- icenReg::ic_sp(
    Surv(l, u, type = 'interval2') ~ grp,
    data = miceData
  )
  expect_equal(result$coefficients, icr_result$coefficients, tolerance = 1e-7)
  expect_equal(result$intervals[[1]], icr_result$T_bull_Intervals)
})


test_that("PH model works with analytical and numerical derivatives", {
  set.seed(1281)
  n <- 4700
  sim_data <- simIC_weib(n = n, inspections = 5, inspectLength = .4)

  d1 <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data,
    control = ic_sp_control(derivMethod = 1)
  )

  d12 <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data,
    control = ic_sp_control(derivMethod = 12)
  )

  d121 <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data,
    control = ic_sp_control(derivMethod = c(12, 1))
  )

  expect_equal(d1$coefficients, d12$coefficients, tolerance = 1e-4)
  expect_equal(d1$llk, d12$llk)
  expect_equal(d1$coefficients, d121$coefficients, tolerance = 1e-4)
  expect_equal(d1$llk, d121$llk)
})


test_that("PH model works for miceData as stratified", {
  data(miceData)

  result <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ strata(grp),
    data = miceData
  )

  icr_result <- icenReg::ic_np(
    Surv(l, u, type = 'interval2') ~ grp,
    data = miceData
  )
  expect_equal(
    result$llk,
    icr_result$fitList[[1]]$llk + icr_result$fitList[[2]]$llk
  )
  expect_equal(
    result$intervals[["ce"]],
    icr_result$fitList$ce$T_bull_Intervals
  )
  expect_equal(
    result$intervals[["ge"]],
    icr_result$fitList$ge$T_bull_Intervals
  )
})

test_that("PH model works for sim data with variance estimation", {
  set.seed(1951)
  sim_data <- simIC_weib(n = 500, inspections = 3, inspectLength = 1)

  result <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data,
    control = ic_sp_control(derivMethod = c(12, 1))
  )
  result_cov <- vcov(result)

  icr_result <- icenReg::ic_sp(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data,
    model = "ph",
    bs_samples = 1000
  )
  icr_result_cov <- vcov(icr_result)
  expect_true(
    abs(sqrt(diag(result_cov)[1]) - sqrt(diag(icr_result_cov)[1])) < 0.05,
  )
  expect_true(
    abs(sqrt(diag(result_cov)[2]) - sqrt(diag(icr_result_cov)[2])) < 0.05,
  )
})


test_that("PO model works for sim data", {
  sim_data <- simIC_weib(n = 500, inspections = 3, inspectLength = 1)

  result <- ic_sp_po(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data,
    control = ic_sp_control(derivMethod = c(12, 1))
  )

  icr_result <- icenReg::ic_sp(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data,
    model = "po"
  )
  expect_equal(result$coefficients, icr_result$coefficients, tolerance = 1e-7)

  expect_equal(result$intervals[[1]], icr_result$T_bull_Intervals)
})


test_that("PO model works for sim data with variance estimation", {
  set.seed(1977)
  sim_data <- simIC_weib(n = 500, inspections = 3, inspectLength = 1)

  result <- ic_sp_po(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data,
    control = ic_sp_control(derivMethod = c(12, 1))
  )
  result_cov <- vcov(result, "oim_curvature")

  icr_result <- icenReg::ic_sp(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data,
    model = "po",
    bs_samples = 1000
  )
  icr_result_cov <- vcov(icr_result)
  expect_equal(result_cov, icr_result_cov, tolerance = 0.05)
  expect_equal(
    sqrt(diag(result_cov)),
    sqrt(diag(icr_result_cov)),
    tolerance = 0.1
  )
})


test_that("print works for ic_sp2", {
  set.seed(1951)
  sim_data <- simIC_weib(n = 100, inspections = 3, inspectLength = 1)

  result <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data,
    control = ic_sp_control(derivMethod = 1)
  )

  expect_snapshot(
    print(result)
  )
})

test_that("print works for stratified fit", {
  data(miceData)
  md <- miceData
  md$strata <- rep(c("A", "B"), length.out = nrow(md))
  result <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ grp + strata(strata),
    data = md
  )

  expect_snapshot(
    print(result)
  )
  summary_res <- summary(result)
  print(summary_res)
})


test_that("profile_fit works as expected", {
  data(miceData)
  md <- miceData
  md$strata <- rep(c("A", "B"), length.out = nrow(md))
  result <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ grp,
    data = md
  )

  prof_result <- profile_fit(result)

  expect_equal(prof_result$coefficients[[1]], result$coefficients[[1]])
  expect_equal(prof_result$llk, result$llk)

  prof_result_07 <- profile_fit(result, 0.7)
  expect_equal(prof_result_07$coefficients[[1]], 0.7)
  expect_equal(prof_result_07$llk, -76.570288)
})

test_that("profile_fit is stable under coefficient perturbation", {
  set.seed(1951)
  sim_data <- simIC_weib(n = 500, inspections = 3, inspectLength = 1)

  fit <- ic_sp_ph(
    Surv(l, u, type = "interval2") ~ x1 + x2,
    data = sim_data,
    control = ic_sp_control(derivMethod = c(12, 1))
  )

  beta_pert <- fit$coefficients
  beta_pert[1] <- beta_pert[1] + 0.5
  beta_pert[2] <- beta_pert[2] - 0.5

  prof <- expect_no_error(profile_fit(fit, beta_pert))

  expect_true(is.finite(prof$llk))
  expect_true(all(is.finite(prof$coefficients)))
  expect_gt(prof$iterations, 0L)
  expect_lt(prof$iterations, fit$other_info$maxIter)
})

test_that("numerical-derivative fit converges before maxIter in flatter settings", {
  set.seed(1998)
  n <- 80
  sim_data <- simIC_weib(
    n = n,
    inspections = 3,
    inspectLength = 0.6,
    b1 = 0.95,
    b2 = -0.5
  )

  z <- as.data.frame(matrix(rnorm(n * 6), nrow = n, ncol = 6))
  colnames(z) <- paste0("z", 1:6)
  sim_data <- cbind(sim_data, z)

  fit <- expect_no_error(
    ic_sp_ph(
      Surv(l, u, type = "interval2") ~ x1 + x2 + z1 + z2 + z3 + z4 + z5 + z6,
      data = sim_data,
      control = ic_sp_control(derivMethod = 1, maxIter = 250)
    )
  )

  expect_true(is.finite(fit$llk))
  expect_true(all(is.finite(fit$coefficients)))
  expect_gt(fit$iterations, 0)
  expect_lt(fit$iterations, fit$other_info$maxIter)
})


test_that("PH model works with profile_ci", {
  set.seed(1952)
  sim_data <- simIC_weib(n = 500, inspections = 3, inspectLength = 1)

  result <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data,
    profile_ci = 0.95
  )
  ci_res_1 <- result$profile_ci

  ci_expected_1 <- matrix(
    c(
      0.1278399833,
      -0.6858324346,
      0.6228696709,
      -0.3968694990
    ),
    nrow = 2,
    dimnames = list(c("x1", "x2"), c("lower", "upper"))
  )
  expect_equal(ci_res_1, ci_expected_1, tolerance = 1e-6)

  ci_res_2 <- confint(result)
  ci_expected_2 <- matrix(
    c(0.125367883, -0.680544642, 0.623030116, -0.398065529),
    nrow = 2,
    dimnames = list(c("x1", "x2"), c("2.5%", "97.5%"))
  )
  expect_equal(ci_res_2, ci_expected_2, tolerance = 1e-7)
  expect_equal(unname(ci_res_1), unname(ci_res_2), tolerance = 1e-2)
})



test_that("profile_ci lower and upper bounds are consistent", {
  set.seed(1998)
  n <- 25
  sim_data <- simIC_weib(
    n = n, inspections = 3, inspectLength = 0.6, b1 = 0.95, b2 = -0.5
  )
  z <- as.data.frame(matrix(rnorm(n * 6), nrow = n, ncol = 6))
  colnames(z) <- paste0("z", 1:6)
  sim_data <- cbind(sim_data, z)

  fit <- ic_sp_ph(
    Surv(l, u, type = "interval2") ~ x1 + x2 + z1 + z2 + z3 + z4 + z5 + z6,
    data = sim_data,
    profile_ci = 0.95
  )
  ci <- fit$profile_ci
  b <- coef(fit)

  expect_true(all(ci[, "lower"] < b))
  expect_true(all(ci[, "upper"] > b))
  expect_true(all(ci[, "upper"] - ci[, "lower"] > 0))
})


test_that("confint works as expected", {
  set.seed(1950)
  sim_data <- simIC_weib(n = 500, inspections = 3, inspectLength = 1)

  result <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data,
    control = ic_sp_control(derivMethod = c(12, 1))
  )

  ci_res <- confint.ic_sp2(result, type = "oim_curvature", typical = 0.5, large = 1.5)
  
  vcov_res <- vcov(result, type = "oim_curvature", typical = 0.5, large = 1.5)
  ci_manual <- coef(result) + matrix(
    rep(c(-qnorm(0.975), qnorm(0.975)), each = 2) * sqrt(diag(vcov_res)),
    nrow = 2, dimnames = dimnames(ci_res)
  )

  expect_equal(ci_res, ci_manual, tolerance = 1e-6)
  ci_expected <- matrix(
    c(
      0.1701936,
      -0.6513221,
      0.6332324,
      -0.3672596
    ),
    nrow = 2,
    dimnames = list(c("x1", "x2"), c("2.5%", "97.5%"))
  )
  expect_equal(ci_res, ci_expected, tolerance = 1e-6)
})
test_that("invert_information warns and returns NA for bad matrices", {
  expect_warning(
    v <- invert_information(matrix(c(1, NaN, NaN, 1), 2)),
    "non-finite"
  )
  expect_true(all(is.na(v)))
  expect_warning(
    v <- invert_information(matrix(c(1, 2, 2, 1), 2)),
    "not positive definite"
  )
  expect_true(all(is.na(v)))
  m <- matrix(c(2, 0.5, 0.5, 1), 2)
  expect_equal(invert_information(m), solve(m))
})

test_that("vcov central differences agree with a fine central-difference reference", {
  set.seed(5309)
  sim_data <- simIC_weib(n = 500, inspections = 5, inspectLength = 1)
  fit <- ic_sp_ph(Surv(l, u, type = "interval2") ~ x1 + x2, data = sim_data)

  pl <- function(d) profile_fit(fit, fit$coefficients + d)$llk
  h <- 0.02
  E <- diag(h, 2)
  info <- matrix(0, 2, 2)
  for (i in 1:2) {
    info[i, i] <- -(pl(E[i, ]) - 2 * fit$llk + pl(-E[i, ])) / h^2
  }
  info[1, 2] <- info[2, 1] <- -(pl(E[1, ] + E[2, ]) - pl(E[1, ] - E[2, ]) -
    pl(-E[1, ] + E[2, ]) + pl(-E[1, ] - E[2, ])) / (4 * h^2)
  se_ref <- sqrt(diag(solve(info)))

  for (type in c("oim_curvature", "oim_fixed")) {
    for (difference in c("forward", "central")) {
      v <- vcov(fit, type = type, difference = difference)
      expect_true(all(abs(sqrt(diag(v)) / se_ref - 1) < 0.05))
    }
  }
  v <- vcov(fit, difference = "central")
  expect_equal(rownames(v), c("x1", "x2"))
})

test_that("h_max caps the step and a large h_max reproduces the uncapped result", {
  set.seed(4471)
  sim_data <- simIC_weib(n = 50, inspections = 5, inspectLength = 1)
  fit <- ic_sp_ph(Surv(l, u, type = "interval2") ~ x1 + x2, data = sim_data)
  # fixed / sqrt(n) = 0.71 > h_max, so the cap changes the result
  capped <- vcov(fit, type = "oim_fixed", fixed = 5)
  uncapped <- suppressWarnings(
    vcov(fit, type = "oim_fixed", fixed = 5, h_max = Inf)
  )
  smaller <- vcov(fit, type = "oim_fixed", fixed = 0.5 * sqrt(50))
  expect_equal(capped, smaller)
  expect_false(isTRUE(all.equal(capped, uncapped)))
})

test_that("vcov never returns NaN silently", {
  set.seed(8812)
  n <- 25
  sim_data <- simIC_weib(n = n, inspections = 5, inspectLength = 1)
  sim_data[paste0("z", 1:6)] <- matrix(rnorm(n * 6), n)
  fit <- suppressWarnings(ic_sp_ph(
    Surv(l, u, type = "interval2") ~ x1 + x2 + z1 + z2 + z3 + z4 + z5 + z6,
    data = sim_data,
    control = ic_sp_control(maxIter = 250)
  ))
  for (fixed in c(1, 5, 10)) {
    v <- withCallingHandlers(
      vcov(fit, type = "oim_fixed", fixed = fixed, h_max = Inf),
      warning = function(w) invokeRestart("muffleWarning")
    )
    expect_false(any(is.nan(v)))
    expect_true(all(is.na(v)) || all(diag(v) > 0))
  }
})

test_that("integer weights match replicated-row fits", {
  set.seed(7311)
  n <- 80
  dat <- data.frame(x = rnorm(n), w = sample(1:4, n, TRUE))
  t_true <- rexp(n, exp(0.5 * dat$x))
  dat$l <- floor(t_true * 4) / 4
  dat$u <- ifelse(runif(n) < 0.2, Inf, dat$l + 0.25)
  rep_dat <- dat[rep(seq_len(n), dat$w), ]
  fml <- Surv(l, u, type = "interval2") ~ x

  for (fitter in list(ic_sp_ph, ic_sp_po)) {
    for (dm in c(12, 1)) {
      ctl <- ic_sp_control(derivMethod = c(dm, 1))
      fw <- fitter(fml, data = dat, weights = w, control = ctl, profile_ci = NULL)
      fr <- fitter(fml, data = rep_dat, control = ctl, profile_ci = NULL)
      expect_equal(fw$llk, fr$llk, tolerance = 1e-6)
      expect_equal(unname(coef(fw)), unname(coef(fr)), tolerance = 1e-4)
    }
  }
})
