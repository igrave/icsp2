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
  expect_equal(result$coefficients, icr_result$coefficients, tolerance = 1e-7)

  expect_equal(result$intervals[[1]], icr_result$T_bull_Intervals)
})


test_that("PH model works for miceData", {
  data(miceData)

  result <- ic_sp_ph(Surv(l, u, type = 'interval2') ~ grp, data = miceData)

  icr_result <- icenReg::ic_sp(
    Surv(l, u, type = 'interval2') ~ grp,
    data = miceData
  )
  expect_equal(result$coefficients, icr_result$coefficients, tolerance = 1e-7)

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
  expect_equal(result$coefficients, icr_result$coefficients, tolerance = 1e-8)
  expect_equal(result$intervals[[1]], icr_result$T_bull_Intervals)
})


test_that("PH model works for with numerical derivative fallback", {
  set.seed(1281)
  n <- 4700
  sim_data <- simIC_weib(n = n, inspections = 5, inspectLength = .4)

  d1 <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data,
    control = ic_sp_control(derivMethod = 1)
  )

  expect_output(
    expect_warning(
      expect_warning(
        d121 <- ic_sp_ph(
          Surv(l, u, type = 'interval2') ~ x1 + x2,
          data = sim_data,
          control = ic_sp_control(derivMethod = c(12, 1))
        ),
        "observation log probability is <= 0",
        fixed = TRUE
      ),
      "Error encountered with derivative method",
      fixed = TRUE
    )
  )
  expect_equal(d1$coefficients, d121$coefficients, tolerance = 1e-6)
  expect_equal(d1$llk, d121$llk)
  expect_equal(d1$s, d121$s, tolerance = 1e-7)
  expect_equal(d1$scores, d121$scores)

  expect_error(
    expect_warning(
      expect_warning(
        ic_sp_ph(
          Surv(l, u, type = 'interval2') ~ x1 + x2,
          data = sim_data,
          control = ic_sp_control(derivMethod = 12)
        ),
        "observation log probability is <= 0",
        fixed = TRUE
      ),
      "Error encountered with derivative method",
      fixed = TRUE
    ),
    "Final log-likelihood is -Inf",
    fixed = TRUE
  )
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
