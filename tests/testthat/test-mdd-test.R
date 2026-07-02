# Regression tests for the mdd_test wild bootstrap, in particular the
# parallel path: chunk sizes must sum to n_boot (no non-positive chunks),
# and `seed` must make parallel runs reproducible.

make_mdd_data <- function(n = 40) {
  set.seed(123)
  list(x = matrix(runif(n), ncol = 1), y = matrix(rnorm(n), ncol = 1))
}

test_that("mdd_test single-core is reproducible with seed", {
  d <- make_mdd_data()
  r1 <- mdd_test(d$x, d$y, n_boot = 50, seed = 1)
  r2 <- mdd_test(d$x, d$y, n_boot = 50, seed = 1)
  expect_identical(r1$wild_bootstrap_values, r2$wild_bootstrap_values)
  expect_identical(r1$p.value, r2$p.value)
  expect_length(r1$wild_bootstrap_values, 50)
})

test_that("mdd_test parallel path returns exactly n_boot values and is reproducible", {
  skip_on_os("windows") # mclapply forking is not available on Windows

  d <- make_mdd_data()
  r1 <- mdd_test(d$x, d$y, n_boot = 51, num_cores = 2, seed = 42)
  r2 <- mdd_test(d$x, d$y, n_boot = 51, num_cores = 2, seed = 42)

  expect_length(r1$wild_bootstrap_values, 51)
  expect_true(all(is.finite(r1$wild_bootstrap_values)))
  expect_identical(r1$wild_bootstrap_values, r2$wild_bootstrap_values)

  # The observed statistic must not depend on num_cores
  r0 <- mdd_test(d$x, d$y, n_boot = 51, num_cores = 1, seed = 42)
  expect_identical(r0$statistic, r1$statistic)
})

test_that("mdd_test handles num_cores exceeding n_boot", {
  skip_on_cran() # spawns more than 2 worker processes
  skip_on_os("windows")

  d <- make_mdd_data()
  # Pre-fix, ceiling-based chunking passed negative counts to the C++
  # bootstrap here and corrupted the result with a worker error.
  expect_no_warning({
    r <- mdd_test(d$x, d$y, n_boot = 2, num_cores = 4, seed = 1)
  })
  expect_length(r$wild_bootstrap_values, 2)
  expect_true(all(is.finite(r$wild_bootstrap_values)))
})
