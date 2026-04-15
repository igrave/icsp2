test_that("summary works", {
  object <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ grp,
    data = miceData
  )

  result <- summary(object)
  expect_equal(result$coefficients, object$coefficients)
  expect_equal(result$var[[1]], 0.125079075)

  expected_coef_matrix <- matrix(
    c(0.678463889117337, 0.35366520141102, 1.9183789821856, 0.0550629741073061),
    nrow = 1,
    dimnames = list("grpge", c("Estimate", "Std.Error", "z.value", "p.value"))
  )

  expect_equal(result$coef_matrix, expected_coef_matrix, tolerance = 1e-7)
})

test_that("print.summary works", {
  object <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ grp,
    data = miceData
  )

  summary_object <- summary(object)
  expect_snapshot(print(summary_object))
})
