test_that("plot function works for stratified fit", {
  skip_if_not_installed("vdiffr")

  md <- miceData
  md$strata <- rep(c("A", "B"), length.out = nrow(md))
  fit <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ grp + strata(strata),
    data = md
  )

  vdiffr::expect_doppelganger(
    "stratified plot",
    plot(fit, col = 1:4, main = "Test Plot")
  )
})

test_that("plot function works for simple fit", {
  skip_if_not_installed("vdiffr")

  md <- miceData
  md$strata <- rep(c("A", "B"), length.out = nrow(md))
  fit <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ grp,
    data = md
  )

  vdiffr::expect_doppelganger(
    "simple plot",
    plot(fit, main = "Test Plot")
  )
})


test_that("lines function works for simple fit", {
  skip_if_not_installed("vdiffr")

  md <- miceData
  md$strata <- rep(c("A", "B"), length.out = nrow(md))
  fit <- ic_sp_ph(
    Surv(l, u, type = 'interval2') ~ grp,
    data = md
  )

  vdiffr::expect_doppelganger("lines plot", {
    plot(c(0, 1200), c(0, 2), type = "n", main = "Lines Plot")
    lines(fit)
  })
})
