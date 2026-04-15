set.seed(1995)
strata_a <- simIC_weib(
  n = 100,
  b1 = 0.1,
  b2 = -0.5,
  model = "ph",
  shape = 2,
  scale = 2,
  inspections = 2,
  inspectLength = 2.5,
  rndDigits = NULL,
  prob_cen = 1
)


strata_b <- simIC_weib(
  n = 100,
  b1 = 0.1,
  b2 = -0.5,
  model = "ph",
  shape = 1.5,
  scale = 2.5,
  inspections = 2,
  inspectLength = 2.5,
  rndDigits = NULL,
  prob_cen = 1
)

stratified_data <- rbind(
  data.frame(strata_a, group = "A"),
  data.frame(strata_b, group = "B")
)
stratified_data$arm <- factor(stratified_data$x2, labels = c("trtA", "trtB"))


test_that("ic test with sas type works with strata", {
  set.seed(123)
  result <- ic_logrank(
    Surv(l, u, type = "interval2") ~ arm + strata(group),
    data = stratified_data,
    type = "sas",
    n_samples = 1000
  )
  expect_s3_class(result, "ic_logrank")

  expected_logrank <- matrix(
    c(10.146555228, -10.146555228, 13.889014391, -13.889014391),
    nrow = 2,
    dimnames = list(c("trtA", "trtB"), c("A", "B"))
  )

  expect_equal(result$logrank, expected_logrank)

  expect_equal(result$logrank_overall, rowSums(expected_logrank))

  expected_p <- unname(pchisq(
    result$statistic,
    df = result$df,
    lower.tail = FALSE
  ))
  expect_equal(result$p.value, expected_p)

  expect_equal(result$var[1, 1], 20.4593484)
})


test_that("ic test with sas type works", {
  set.seed(1234)
  result <- ic_logrank(
    Surv(l, u, type = "interval2") ~ arm,
    data = stratified_data,
    type = "sas",
    n_samples = 1000
  )
  expect_s3_class(result, "ic_logrank")

  expected_logrank <- matrix(
    c(26.5067259, -26.5067259),
    nrow = 2,
    dimnames = list(c("trtA", "trtB"), "1")
  )

  expect_equal(result$logrank, expected_logrank)

  expect_equal(result$logrank_overall, rowSums(expected_logrank))

  expected_p <- unname(pchisq(
    result$statistic,
    df = result$df,
    lower.tail = FALSE
  ))

  expect_equal(result$p.value, expected_p)
  expect_equal(
    format.pval(result$p.value),
    format.pval(3.62680505591868e-08)
  )

  expect_equal(result$var[1, 1], 23.1582084)
})


test_that("ic test with hly type works", {
  set.seed(1234)
  result <- ic_logrank(
    Surv(l, u, type = "interval2") ~ arm,
    data = stratified_data,
    type = "hly",
    n_samples = 1000
  )
  expect_s3_class(result, "ic_logrank")
  expected_logrank <- matrix(
    c(26.4248745, -26.4248745),
    nrow = 2,
    dimnames = list(c("trtA", "trtB"), "1")
  )

  expect_equal(result$logrank, expected_logrank)

  expect_equal(result$logrank_overall, rowSums(expected_logrank))

  expected_p <- unname(pchisq(
    result$statistic,
    df = result$df,
    lower.tail = FALSE
  ))

  expect_equal(result$p.value, expected_p)
  expect_equal(
    format.pval(result$p.value),
    format.pval(3.99406586107093e-08)
  )

  expect_equal(result$var[1, 1], 23.1582084)
})

test_that("ic test with hly type works with strata", {
  set.seed(123)
  result <- ic_logrank(
    Surv(l, u, type = "interval2") ~ arm + strata(group),
    data = stratified_data,
    type = "hly",
    n_samples = 1000
  )
  expect_s3_class(result, "ic_logrank")
  expected_logrank <- matrix(
    c(10.1448717247994, -10.1448717247994, 13.9123283171416, -13.9123283171416),
    nrow = 2,
    dimnames = list(c("trtA", "trtB"), c("A", "B"))
  )

  expect_equal(result$logrank, expected_logrank)

  expect_equal(result$logrank_overall, rowSums(expected_logrank))

  expected_p <- unname(pchisq(
    result$statistic,
    df = result$df,
    lower.tail = FALSE
  ))
  expect_equal(result$p.value, expected_p)
  expect_equal(
    format.pval(result$p.value),
    format.pval(1.04555714969813e-07)
  )
  expect_equal(result$var[1, 1], 20.4593484)
})

test_that("ic_logrank fails with incorrect type", {
  expect_error(
    ic_logrank(
      Surv(l, u, type = "interval2") ~ arm,
      data = stratified_data,
      type = "invalid_type"
    ),
    "should be one of"
  )
})

test_that("ic_logrank fails with invalid n_samples", {
  expect_error(
    ic_logrank(
      Surv(l, u, type = "interval2") ~ arm,
      data = stratified_data,
      n_samples = -100
    ),
    "n_samples should be a positive integer"
  )
  expect_error(
    ic_logrank(
      Surv(l, u, type = "interval2") ~ arm,
      data = stratified_data,
      n_samples = 0
    ),
    "n_samples should be a positive integer"
  )
  expect_error(
    ic_logrank(
      Surv(l, u, type = "interval2") ~ arm,
      data = stratified_data,
      n_samples = NULL
    ),
    "n_samples should be a positive integer"
  )
  expect_error(
    ic_logrank(
      Surv(l, u, type = "interval2") ~ arm,
      data = stratified_data,
      n_samples = "a"
    ),
    "n_samples should be a positive integer"
  )
})

test_that("ic_logrank fails with missing strata variable", {
  expect_error(
    ic_logrank(
      Surv(l, u, type = "interval2") ~ arm + strata(nonexistent_group),
      data = stratified_data,
      type = "sas",
      n_samples = 1000
    ),
    "object 'nonexistent_group' not found"
  )
})

test_that("ic_logrank fails with missing group variable", {
  expect_error(
    ic_logrank(
      Surv(l, u, type = "interval2") ~ treatment,
      data = stratified_data,
      type = "sas",
      n_samples = 1000
    ),
    "object 'treatment' not found"
  )
})

test_that("ic_logrank fails with just intercept", {
  expect_error(
    ic_logrank(
      Surv(l, u, type = "interval2") ~ 1,
      data = stratified_data,
      type = "sas",
      n_samples = 1000
    ),
    "Formula should have exactly one grouping variable"
  )
})


test_that("ic_logrank fails with two group variables", {
  expect_error(
    ic_logrank(
      Surv(l, u, type = "interval2") ~ arm + group,
      data = stratified_data,
      type = "sas",
      n_samples = 1000
    ),
    "Formula should have exactly one grouping variable"
  )
})
