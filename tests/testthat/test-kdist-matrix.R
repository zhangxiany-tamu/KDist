# Smoke tests for the core distance/kernel constructors and base statistics

test_that("KDist_matrix produces correct dimensions", {
  set.seed(1)
  x <- matrix(rnorm(100), 20, 5)

  # Euclidean distance
  d_euc <- KDist_matrix(x, type = "euclidean")
  expect_equal(dim(d_euc), c(20, 20))
  expect_equal(diag(d_euc), rep(0, 20))
  expect_true(isSymmetric(d_euc))

  # Gaussian kernel
  k_gauss <- KDist_matrix(x, type = "gaussian")
  expect_equal(dim(k_gauss), c(20, 20))
  expect_equal(diag(k_gauss), rep(1, 20))
  expect_true(isSymmetric(k_gauss))

  # Laplacian kernel
  k_lap <- KDist_matrix(x, type = "laplacian")
  expect_equal(dim(k_lap), c(20, 20))
  expect_equal(diag(k_lap), rep(1, 20))
  expect_true(isSymmetric(k_lap))
})

test_that("bw_optim produces valid bandwidths", {
  set.seed(2)
  x <- matrix(rnorm(100), 20, 5)

  # Single bandwidth
  bw <- bw_optim(x)
  expect_length(bw, 1)
  expect_true(is.numeric(bw))
  expect_gt(bw, 0)

  # Grouped bandwidth
  group <- c(1, 1, 2, 2, 2)
  bw_group <- bw_optim(x, group = group)
  expect_length(bw_group, length(unique(group)))
  expect_true(all(bw_group > 0))
})

test_that("mmd from raw data matches mmd from a precomputed distance matrix", {
  set.seed(3)
  x <- matrix(rnorm(50), 10, 5)
  y <- matrix(rnorm(50, mean = 0.5), 10, 5)

  mmd_val <- mmd(x, y, type = "euclidean")
  expect_true(is.numeric(mmd_val))
  expect_length(mmd_val, 1)

  D <- KDist_matrix(rbind(x, y), type = "euclidean")
  mmd_val2 <- suppressMessages(mmd(D, type = "euclidean", n = 10, m = 10))
  expect_equal(mmd_val, mmd_val2, tolerance = 1e-10)
})

test_that("hsic produces valid values across input modes", {
  set.seed(4)
  x <- matrix(rnorm(50), 10, 5)
  y <- x + matrix(rnorm(50, sd = 0.1), 10, 5)

  hsic_val <- hsic(x, y)
  expect_true(is.numeric(hsic_val))
  expect_length(hsic_val, 1)
  expect_gt(hsic_val, 0)

  hsic_val_u <- hsic(x, y, u_center = TRUE)
  expect_true(is.numeric(hsic_val_u))

  Dx <- KDist_matrix(x)
  Dy <- KDist_matrix(y)
  hsic_val2 <- hsic(Dx, Dy, is_distance = TRUE)
  expect_true(is.numeric(hsic_val2))
})
