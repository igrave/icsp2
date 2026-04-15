test_that("check_weights catches invalid weights", {
  mf <- model.frame(
    x ~ y,
    data = data.frame(
      x = 1:10,
      y = 1:10,
      w = c(1, 2, -1, 4, 5, 6, 7, 8, 9, 10)
    ),
    weights = w
  )
  expect_error(
    check_weights(mf),
    "negative weights not allowed"
  )

  mf_na <- model.frame(
    x ~ y,
    data = data.frame(
      x = 1:10,
      y = 1:10,
      w = c(1, 2, 1, 4, 5, 6, 7, 8, NA, 10)
    ),
    na.action = na.pass,
    weights = w
  )
  expect_error(
    check_weights(mf_na),
    "NAs not allowed in weights"
  )
})

test_that("check_matrix catches rank deficient matrices", {
  mf <- model.frame(
    y ~ x1 + x2,
    data = data.frame(
      x1 = 1:10,
      x2 = 1:10,
      y = 1:10
    )
  )
  x_mat <- model.matrix(~ x1 + x2, data = mf)
  expect_error(
    check_matrix(x_mat),
    "Covariate matrix is rank deficient. Check covariates"
  )
})
