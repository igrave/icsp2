test_that("plot function creates expected output", {
  skip_if_not_installed("vdiffr")
  
  # Create test data
  test_data <- data.frame(
    x = 1:10,
    y = (1:10)^2
  )
  
  # Test the plot
  vdiffr::expect_doppelganger(
    "basic plot output",
    plot(test_data$x, test_data$y, main = "Test Plot")
  )
})