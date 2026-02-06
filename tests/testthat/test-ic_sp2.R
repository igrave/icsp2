test_that("PH model works for sim data", {
  sim_data <- simIC_weib(n = 500, inspections = 3, inspectLength = 1)

  result <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data,
    settings = ic_sp_settings(derivMethod = 1)
  )

  icr_result <- icenReg::ic_sp(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data
  )
  expect_equal(result$coefficients, icr_result$coefficients, tolerance = 1e-7)

  expect_equal(result$T_bull_Intervals[[1]], icr_result$T_bull_Intervals)
})


test_that("PH model works for miceData", {
  data(miceData)

  result <- ic_sp_ph(Surv(l, u, type = 'interval2') ~ grp, data = miceData)

  icr_result <- icenReg::ic_sp(
    Surv(l, u, type = 'interval2') ~ grp,
    data = miceData
  )
  expect_equal(result$coefficients, icr_result$coefficients, tolerance = 1e-8)

  expect_equal(result$T_bull_Intervals[[1]], icr_result$T_bull_Intervals)
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
    result$T_bull_Intervals[["ce"]],
    icr_result$fitList$ce$T_bull_Intervals
  )
  expect_equal(
    result$T_bull_Intervals[["ge"]],
    icr_result$fitList$ge$T_bull_Intervals
  )
})

test_that("PH model works for sim data with variance estimation", {
  set.seed(1951)
  sim_data <- simIC_weib(n = 500, inspections = 3, inspectLength = 1)

  result <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data,
    settings = ic_sp_settings(derivMethod = c(12, 1))
  )
  result_cov <- vcov(result, .2)
  dimnames(result_cov) <- list(
    names(result$coefficients),
    names(result$coefficients)
  )

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
    settings = ic_sp_settings(derivMethod = c(12, 1))
  )

  icr_result <- icenReg::ic_sp(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data,
    model = "po"
  )
  expect_equal(result$coefficients, icr_result$coefficients, tolerance = 1e-7)

  expect_equal(result$T_bull_Intervals[[1]], icr_result$T_bull_Intervals)
})


test_that("PO model works for sim data with variance estimation", {
  set.seed(1977)
  sim_data <- simIC_weib(n = 500, inspections = 3, inspectLength = 1)

  result <- ic_sp_po(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data,
    settings = ic_sp_settings(derivMethod = c(12, 1))
  )
  result_cov <- vcov(result, .2)
  dimnames(result_cov) <- list(
    names(result$coefficients),
    names(result$coefficients)
  )

  icr_result <- icenReg::ic_sp(
    Surv(l, u, type = 'interval2') ~ x1 + x2,
    data = sim_data,
    model = "po",
    bs_samples = 1000
  )
  icr_result_cov <- vcov(icr_result)
  expect_equal(result_cov, icr_result_cov, tolerance = 0.01)
  expect_equal(
    sqrt(diag(result_cov)),
    sqrt(diag(icr_result_cov)),
    tolerance = 0.05
  )
})
