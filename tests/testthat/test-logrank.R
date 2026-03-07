test_that("ic_logrank works with two groups", {
  data(miceData)

  # Simple two-group test
  result <- ic_logrank(Surv(l, u, type = "interval2") ~ grp, data = miceData)

  expect_s3_class(result, "ic_logrank")
  expect_true(!is.null(result$statistic))
  expect_true(!is.null(result$p.value))
  expect_equal(unname(result$parameter), 1) # 2 groups = 1 df
  expect_equal(length(result$score), 1)
  expect_equal(dim(result$information), c(1, 1))
})


test_that("ic_logrank works with multiple groups", {
  # Simulate data with 3 groups
  set.seed(123)
  n <- 150
  group <- factor(rep(c("A", "B", "C"), each = n / 3))

  # Generate interval-censored data
  true_time <- rexp(n, rate = 0.5)
  l <- pmax(0, true_time - runif(n, 0, 2))
  r <- true_time + runif(n, 0, 2)

  test_data <- data.frame(l = l, r = r, group = group)

  result <- ic_logrank(Surv(l, r, type = "interval2") ~ group, data = test_data)

  expect_s3_class(result, "ic_logrank")
  expect_equal(unname(result$parameter), 2) # 3 groups = 2 df
  expect_equal(length(result$score), 2)
  expect_equal(dim(result$information), c(2, 2))
  expect_equal(length(result$groups), 3)
})


test_that("ic_logrank works with stratification", {
  # Simulate stratified data
  set.seed(456)
  n <- 120
  group <- factor(rep(c("Control", "Treatment"), each = n / 2))
  site <- factor(rep(c("Site1", "Site2"), times = n / 2))

  true_time <- rexp(n, rate = 0.3 + 0.2 * (site == "Site2"))
  l <- pmax(0, true_time - runif(n, 0, 1.5))
  u <- true_time + runif(n, 0, 1.5)

  test_data <- data.frame(l = l, u = u, group = group, site = site)

  result <- ic_logrank(
    Surv(l, u, type = "interval2") ~ group + strata(site),
    data = test_data
  )

  expect_s3_class(result, "ic_logrank")
  expect_equal(unname(result$parameter), 1) # 2 groups = 1 df
  expect_true(!is.null(result$strata))
  expect_equal(length(result$strata), 2)
})


test_that("ic_logrank rejects invalid formulas", {
  data(miceData)

  # No group variable
  expect_error(
    ic_logrank(Surv(l, u, type = "interval2") ~ 1, data = miceData),
    "grouping variable"
  )

  # Wrong Surv type
  expect_error(
    ic_logrank(Surv(l) ~ grp, data = miceData),
    "interval2"
  )
})


test_that("ic_logrank handles factor conversion", {
  data(miceData)

  # Create numeric group variable
  test_data <- miceData
  test_data$grp_numeric <- as.numeric(test_data$grp)

  result <- ic_logrank(
    Surv(l, u, type = "interval2") ~ grp_numeric,
    data = test_data
  )

  expect_s3_class(result, "ic_logrank")
  expect_equal(unname(result$parameter), 1)
})


test_that("print.ic_logrank works", {
  data(miceData)
  result <- ic_logrank(Surv(l, u, type = "interval2") ~ grp, data = miceData)

  # Should not error
  expect_output(print(result), "Sun's log-rank test")
  expect_output(print(result), "Groups:")
  expect_output(print(result), "p-value")
})


test_that("summary.ic_logrank works", {
  data(miceData)
  result <- ic_logrank(Surv(l, u, type = "interval2") ~ grp, data = miceData)

  # Should not error
  expect_output(summary(result), "Sun's log-rank test")
  expect_output(summary(result), "Score vector")
  expect_output(summary(result), "Information matrix")
})


test_that("ic_logrank gives consistent results", {
  # Test with known data structure
  set.seed(789)
  n <- 100
  group <- factor(rep(c("A", "B"), each = n / 2))

  # Same distribution for both groups (H0 true)
  true_time <- rexp(n, rate = 0.5)
  l <- pmax(0, true_time - runif(n, 0, 1))
  r <- true_time + runif(n, 0, 1)

  test_data <- data.frame(l = l, r = r, group = group)

  result1 <- ic_logrank(
    Surv(l, r, type = "interval2") ~ group,
    data = test_data
  )
  result2 <- ic_logrank(
    Surv(l, r, type = "interval2") ~ group,
    data = test_data
  )

  # Results should be identical
  expect_equal(result1$statistic, result2$statistic)
  expect_equal(result1$p.value, result2$p.value)
})


test_that("ic_logrank test statistic is positive", {
  data(miceData)
  result <- ic_logrank(Surv(l, u, type = "interval2") ~ grp, data = miceData)

  expect_true(result$statistic >= 0)
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})


test_that("information matrix is positive definite", {
  data(miceData)
  result <- ic_logrank(Surv(l, u, type = "interval2") ~ grp, data = miceData)

  # Check eigenvalues are positive
  eigenvalues <- eigen(result$information, only.values = TRUE)$values
  expect_true(all(eigenvalues > 0))
})
