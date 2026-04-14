test_that("simEventTime works for ph and weibull", {
  set.seed(1951)
  n <- 10
  x1 <- runif(n, -1, 1)
  x2 <- 1 - 2 * rbinom(n, 1, 0.5)
  b1 <- .3
  b2 <- -.3
  lin_pred <- x1 * b1 + x2 * b2

  set.seed(2001)
  sim_times_ph <- simEventTime(
    lin_pred,
    model = "ph",
    dist = qweibull,
    paramList = list(shape = 2, scale = 2)
  )
  expect_equal(
    sim_times_ph,
    c(
      1.05131276285532,
      1.40834876820508,
      2.17880592480362,
      2.4340307642724,
      2.83389253132937,
      2.05229572996039,
      1.45232583511554,
      2.47234605091589,
      1.32690276202114,
      1.98811584929
    )
  )

  set.seed(2001)
  sim_times_aft <- simEventTime(
    lin_pred,
    model = "aft",
    dist = qweibull,
    paramList = list(shape = 2, scale = 2)
  )
  expect_equal(
    sim_times_aft,
    c(
      2.38360092160946,
      1.94402494446242,
      1.27275832434661,
      1.01498187544161,
      0.762049433161491,
      1.27855851297724,
      1.78718498577176,
      0.968145936737901,
      2.06165894416794,
      1.42652724493612
    )
  )
})

test_that("simEventTime works for po and qexp", {
  set.seed(2001)
  lin_pred <- c(
    -0.000407785875722766,
    0.0050925994757563,
    0.248570929421112,
    -0.00522798215970399,
    -0.0668151799589395,
    -0.0681034242734313,
    -0.0946230549830943,
    -0.122658793022856,
    0.0233809104189276,
    0.129596768226475
  )
  sim_times_po <- simEventTime(
    lin_pred,
    model = "po",
    dist = qexp,
    paramList = list(rate = 2)
  )
  expect_equal(
    sim_times_po,
    c(
      0.138051788576318,
      0.250197568144732,
      0.860514144469606,
      0.734687869891502,
      0.910833623844616,
      0.47078227425592,
      0.222348792815792,
      0.631140533478268,
      0.229562533342731,
      0.607105708360802
    )
  )
})


test_that("simIC_weib works with rounding", {
  set.seed(2025)
  result <- simIC_weib(
    n = 5,
    b1 = .3,
    b2 = -.3,
    model = 'ph',
    shape = 2,
    scale = 2,
    inspections = 6,
    inspectLength = 1,
    rndDigits = 2
  )
  expect_equal(result$l, c(1.21, 0.37, 1.27, 2.89, 1.37))
  expect_equal(result$u, c(2.13, 1.18, 2.04, 3.88, 2.28))
})
