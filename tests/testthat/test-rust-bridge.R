# Tests for the extendr R-Rust bridge

test_that("rust_dot_product computes correct dot product", {
  # Basic case
  expect_equal(rust_dot_product(c(1, 2, 3), c(4, 5, 6)), 32)

  # Zero vectors
  expect_equal(rust_dot_product(c(0, 0, 0), c(1, 2, 3)), 0)

  # Single element
  expect_equal(rust_dot_product(c(5), c(3)), 15)

  # Negative values
  expect_equal(rust_dot_product(c(-1, 2, -3), c(4, -5, 6)), -32)

  # Floating point
  expect_equal(rust_dot_product(c(1.5, 2.5), c(3.0, 4.0)), 1.5 * 3.0 + 2.5 * 4.0)
})

test_that("rust_dot_product matches R's sum(x*y)", {
  set.seed(42)
  x <- rnorm(1000)
  y <- rnorm(1000)
  expect_equal(rust_dot_product(x, y), sum(x * y), tolerance = 1e-10)
})

test_that("rust_parallel_sum computes correct sum", {
  # Basic case
  expect_equal(rust_parallel_sum(c(1, 2, 3, 4, 5)), 15)

  # Single element
  expect_equal(rust_parallel_sum(c(42)), 42)

  # Negative values
  expect_equal(rust_parallel_sum(c(-1, -2, -3)), -6)

  # Mixed
  expect_equal(rust_parallel_sum(c(-1, 0, 1)), 0)
})

test_that("rust_parallel_sum matches R's sum()", {
  set.seed(42)
  x <- rnorm(10000)
  expect_equal(rust_parallel_sum(x), sum(x), tolerance = 1e-10)
})

test_that("rust_parallel_sum handles large vectors", {
  # Verify rayon parallel reduction works on larger data
  x <- rep(1.0, 100000)
  expect_equal(rust_parallel_sum(x), 100000)
})

test_that("Rust functions are accessible via the extendr bridge", {
  # Verify the functions exist and are callable
  expect_true(is.function(rust_dot_product))
  expect_true(is.function(rust_parallel_sum))
})
