# Regression tests for the C++ permutation kernels behind hsic_test /
# dcov_test, mmd_test / ed_test, and dhsic_test. The pinned values below were
# generated from the R-level permutation loops these kernels replaced and
# verified bit-identical; they must never drift.

test_that("hsic_test permutation distribution is reproducible and pinned", {
  set.seed(11)
  x <- matrix(rnorm(60), 60, 1)
  y <- x + matrix(rnorm(60), 60, 1)

  r <- hsic_test(x, y, type = "gaussian", n_perm = 50, seed = 3)
  expect_equal(r$statistic, 0.0211154912938505, tolerance = 1e-14)
  expect_identical(r$p.value, 0)
  expect_length(r$permutation_values, 50)
  expect_equal(r$permutation_values[1], 0.00293824854952654, tolerance = 1e-14)
  expect_equal(r$permutation_values[50], 0.00363583509606609, tolerance = 1e-14)

  # Same seed reproduces the full permutation vector exactly
  r2 <- hsic_test(x, y, type = "gaussian", n_perm = 50, seed = 3)
  expect_identical(r$permutation_values, r2$permutation_values)

  # num_cores no longer changes anything
  r4 <- hsic_test(x, y, type = "gaussian", n_perm = 50, seed = 3, num_cores = 4)
  expect_identical(r$permutation_values, r4$permutation_values)
})

test_that("mmd_test permutation distribution is reproducible and pinned", {
  set.seed(12)
  a <- matrix(rnorm(50 * 2), 50, 2)
  b <- matrix(rnorm(40 * 2, 0.5), 40, 2)

  r <- mmd_test(a, b, type = "euclidean", u_center = TRUE, n_perm = 50, seed = 4)
  expect_equal(r$statistic, 0.283998443457406, tolerance = 1e-14)
  expect_identical(r$p.value, 0)
  expect_equal(r$permutation_values[1], 0.00562114798914903, tolerance = 1e-14)
  expect_equal(r$permutation_values[50], -0.00751879917509735, tolerance = 1e-14)

  r2 <- mmd_test(a, b, type = "euclidean", u_center = TRUE, n_perm = 50, seed = 4)
  expect_identical(r$permutation_values, r2$permutation_values)
})

test_that("dhsic_test permutation distribution is reproducible and pinned", {
  set.seed(13)
  u <- runif(60, -pi, pi)
  xl <- list(matrix(sin(u)), matrix(cos(u)), matrix(sin(u) * cos(u)))

  r <- dhsic_test(xl, type = "gaussian", n_perm = 50, seed = 5)
  expect_equal(r$statistic, 0.0642048192002568, tolerance = 1e-14)
  expect_identical(r$p.value, 0)
  expect_equal(r$permutation_values[1], 0.00651284997468501, tolerance = 1e-14)
  expect_equal(r$permutation_values[50], 0.0073350473864085, tolerance = 1e-14)

  r2 <- dhsic_test(xl, type = "gaussian", n_perm = 50, seed = 5)
  expect_identical(r$permutation_values, r2$permutation_values)
})

test_that("hsic_test handles asymmetric user-supplied matrices (fallback path)", {
  set.seed(42)
  A <- matrix(abs(rnorm(30^2)), 30, 30); diag(A) <- 0  # deliberately asymmetric
  B <- matrix(abs(rnorm(30^2)), 30, 30); diag(B) <- 0

  r <- hsic_test(A, B, type = "euclidean", n_perm = 50, seed = 5, is_distance = TRUE)
  expect_true(is.numeric(r$p.value))
  expect_length(r$permutation_values, 50)

  r2 <- hsic_test(A, B, type = "euclidean", n_perm = 50, seed = 5, is_distance = TRUE)
  expect_identical(r$permutation_values, r2$permutation_values)
})

test_that("wrappers share their parent's permutation distribution exactly", {
  set.seed(15)
  x <- rnorm(40); y <- x^2 + rnorm(40)

  # dcov_test is hsic_test with type = "euclidean": identical everything
  r_wrap <- dcov_test(x, y, n_perm = 30, seed = 7)
  r_parent <- hsic_test(x, y, type = "euclidean", n_perm = 30, seed = 7)
  expect_identical(r_wrap$statistic, r_parent$statistic)
  expect_identical(r_wrap$permutation_values, r_parent$permutation_values)
  expect_identical(r_wrap$p.value, r_parent$p.value)

  # ed_test is mmd_test with type = "euclidean"
  a <- matrix(rnorm(30), 30, 1); b <- matrix(rnorm(25, 0.5), 25, 1)
  e_wrap <- ed_test(a, b, n_perm = 30, seed = 8)
  e_parent <- mmd_test(a, b, type = "euclidean", n_perm = 30, seed = 8)
  expect_identical(e_wrap$statistic, e_parent$statistic)
  expect_identical(e_wrap$permutation_values, e_parent$permutation_values)
})

test_that("n_perm is validated and edge sizes work", {
  set.seed(16)
  x <- rnorm(20); y <- rnorm(20)

  expect_error(hsic_test(x, y, n_perm = 0), "positive integer")
  expect_error(mmd_test(matrix(x), matrix(y), n_perm = -1), "positive integer")
  expect_error(dhsic_test(list(matrix(x), matrix(y)), n_perm = 0), "positive integer")

  # n_perm = 1 and small n run cleanly
  r <- hsic_test(x, y, n_perm = 1, seed = 1)
  expect_length(r$permutation_values, 1)
  r2 <- dcov_test(rnorm(10), rnorm(10), u_center = TRUE, n_perm = 20, seed = 2)
  expect_true(is.numeric(r2$p.value))
})

test_that("dhsic_test large-n automatic bandwidth uses the exact rebuild path", {
  skip_on_cran() # a few seconds of Gram-matrix rebuilds
  set.seed(17)
  n <- 1050 # above the bw_rcpp subsample limit
  xl <- list(matrix(rnorm(n)), matrix(rnorm(n)))
  r <- dhsic_test(xl, type = "gaussian", n_perm = 5, seed = 9)
  expect_length(r$permutation_values, 5)
  r2 <- dhsic_test(xl, type = "gaussian", n_perm = 5, seed = 9)
  expect_identical(r$permutation_values, r2$permutation_values)
})

test_that("observed statistics are unchanged by the permutation kernels", {
  # The observed statistic path (hsic/mmd/dhsic) is independent of the
  # permutation kernels; anchor it against direct statistic calls.
  set.seed(14)
  x <- matrix(rnorm(50), 50, 1)
  y <- x + matrix(rnorm(50), 50, 1)
  r <- hsic_test(x, y, type = "gaussian", n_perm = 10, seed = 1)
  expect_identical(r$statistic, hsic(x, y, type = "gaussian"))

  xl <- list(x, y)
  rd <- dhsic_test(xl, type = "gaussian", n_perm = 10, seed = 1)
  expect_identical(rd$statistic, dhsic(xl, type = "gaussian"))
})
