library(testthat)

test_that("pairwiseCorrelation calculates correct results", {
  # Test example data
  data <- rbind(c(2.5,4.2,1.6,6.2,3.1), c(5.2,3.5,8.2,7.4,5.3), c(0.3, 0.4, 1.4, 2.7, 3.1) )

  # Expected result
  expected_result <- c(-0.07970577, 0.35220623, 0.42153557)

  # Run the function
  result <- cormat(data)

  # Verify the result
  expect_equal(as.vector(result), expected_result)
})


test_that("pairwiseCorrelation matches base cor function", {
  # Test example data
  data <- rbind(c(1,1,3,4,NA),c(1,2,3,4,5) )

  # Expected result
  expected_result <- cor(data[1,],data[2,], use="pair")

  # Run the function and extract the first element
  result <- cormat(data)[1]

  # Verify the result
  expect_equal(result, expected_result)
})
