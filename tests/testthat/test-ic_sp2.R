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
    dimnames = list(c("x1", "x2"), c("2.5 %", "97.5 %"))
  )
  expect_equal(ci_res_2, ci_expected_2, tolerance = 1e-7)
  expect_equal(unname(ci_res_1), unname(ci_res_2), tolerance = 1e-2)
})
