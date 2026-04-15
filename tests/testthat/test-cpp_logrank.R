test_that("cpp_logrank matches survdiff for 2 groups", {
  set.seed(123)
  n <- 200
  time <- rexp(n, rate = 1)
  event <- sample(0:1, n, replace = TRUE, prob = c(0.3, 0.7))
  group <- rep(1:2, each = n / 2)

  res <- cpp_logrank(time, event, group)
  sd_res <- survdiff(Surv(time, event) ~ group)

  expect_equal(res$observed, sd_res$obs, tolerance = 1e-10)
  expect_equal(res$expected, sd_res$exp, tolerance = 1e-10)
  expect_equal(res$variance, sd_res$var, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(res$n_by_group, sd_res$n, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("cpp_logrank matches survdiff for 3 groups", {
  set.seed(456)
  n <- 300
  time <- rexp(n, rate = rep(c(1, 1.5, 2), each = 100))
  event <- sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8))
  group <- rep(1:3, each = 100)

  res <- cpp_logrank(time, event, group)
  sd_res <- survdiff(Surv(time, event) ~ group)

  expect_equal(res$observed, sd_res$obs, tolerance = 1e-10)
  expect_equal(res$expected, sd_res$exp, tolerance = 1e-10)
  expect_equal(res$variance, sd_res$var, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(res$n_by_group, sd_res$n, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("cpp_logrank matches survdiff with heavy ties", {
  set.seed(789)
  n <- 150
  # Discrete times with many ties
  time <- sample(1:5, n, replace = TRUE)
  event <- sample(0:1, n, replace = TRUE, prob = c(0.4, 0.6))
  group <- sample(1:2, n, replace = TRUE)

  res <- cpp_logrank(time, event, group)
  sd_res <- survdiff(Surv(time, event) ~ group)

  expect_equal(res$observed, sd_res$obs, tolerance = 1e-10)
  expect_equal(res$expected, sd_res$exp, tolerance = 1e-10)
  expect_equal(res$variance, sd_res$var, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("cpp_logrank matches survdiff when all events observed", {
  set.seed(101)
  n <- 100
  time <- rexp(n)
  event <- rep(1L, n)
  group <- rep(1:2, each = 50)

  res <- cpp_logrank(time, event, group)
  sd_res <- survdiff(Surv(time, event) ~ group)

  expect_equal(res$observed, sd_res$obs, tolerance = 1e-10)
  expect_equal(res$expected, sd_res$exp, tolerance = 1e-10)
  expect_equal(res$variance, sd_res$var, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("cpp_logrank returns correct class and structure", {
  set.seed(42)
  time <- rexp(50)
  event <- sample(0:1, 50, replace = TRUE)
  group <- rep(1:2, each = 25)

  res <- cpp_logrank(time, event, group)
  expect_s3_class(res, "cpp_logrank")
  expect_named(res, c("observed", "expected", "variance", "n_by_group"))
  expect_length(res$observed, 2)
  expect_length(res$expected, 2)
  expect_equal(dim(res$variance), c(2, 2))
  expect_length(res$n_by_group, 2)
})

test_that("cpp_logrank errors on mismatched lengths", {
  expect_error(cpp_logrank(1:5, 0:1, 1:5), "same length")
  expect_error(cpp_logrank(1:5, rep(1, 5), 1:3), "same length")
})

test_that("cpp_logrank errors on invalid event values", {
  expect_error(cpp_logrank(1:5, c(0, 1, 2, 0, 1), rep(1L, 5)), "0 and 1")
})

test_that("cpp_logrank errors on non-consecutive groups", {
  expect_error(
    cpp_logrank(1:5, rep(1, 5), c(1L, 1L, 3L, 3L, 3L)),
    "consecutive integers"
  )
})

test_that("cpp_logrank errors on non-numeric time", {
  expect_error(
    cpp_logrank(letters[1:5], rep(1, 5), rep(1:2, length.out = 5)),
    "numeric"
  )
})
