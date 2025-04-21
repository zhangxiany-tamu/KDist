library(testthat)
library(KDist)

context("Basic Distance and Kernel Tests")

# Test KDist_matrix function
test_that("KDist_matrix produces correct dimensions", {
  skip_if(!"KDistTest" %in% loadedNamespaces())

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

# Test bw_optim function
test_that("bw_optim produces valid bandwidths", {
  skip_if(!"KDistTest" %in% loadedNamespaces())

  x <- matrix(rnorm(100), 20, 5)

  # Single bandwidth
  bw <- bw_optim(x)
  expect_true(length(bw) == 1)
  expect_true(is.numeric(bw))
  expect_true(bw > 0)

  # Grouped bandwidth
  group <- c(1, 1, 2, 2, 2)
  bw_group <- bw_optim(x, group = group)
  expect_true(length(bw_group) == 2)
  expect_true(all(bw_group > 0))
})

# Test mmd function
test_that("mmd produces valid values", {
  skip_if(!"KDistTest" %in% loadedNamespaces())

  x <- matrix(rnorm(50), 10, 5)
  y <- matrix(rnorm(50, mean = 0.5), 10, 5)

  # Direct calculation
  mmd_val <- mmd(x, y, type = "euclidean")
  expect_true(is.numeric(mmd_val))
  expect_length(mmd_val, 1)

  # With distance matrix
  D <- KDist_matrix(rbind(x, y), type = "euclidean")
  mmd_val2 <- mmd(D, n = 10, m = 10, type = "euclidean")
  expect_true(is.numeric(mmd_val2))
  expect_length(mmd_val2, 1)

  # Values should be close
  expect_equal(mmd_val, mmd_val2, tolerance = 1e-10)
})

# Test hsic function
test_that("hsic produces valid values", {
  skip_if(!"KDistTest" %in% loadedNamespaces())

  x <- matrix(rnorm(50), 10, 5)
  y <- x + matrix(rnorm(50, sd = 0.1), 10, 5)  # Correlated data

  # Default calculation
  hsic_val <- hsic(x, y)
  expect_true(is.numeric(hsic_val))
  expect_length(hsic_val, 1)
  expect_true(hsic_val > 0)

  # With U-centering
  hsic_val_u <- hsic(x, y, u_center = TRUE)
  expect_true(is.numeric(hsic_val_u))

  # Distance matrices
  Dx <- KDist_matrix(x)
  Dy <- KDist_matrix(y)
  hsic_val2 <- hsic(Dx, Dy, is_distance = TRUE)
  expect_true(is.numeric(hsic_val2))
})

# Test matrix_v_center function
test_that("matrix_v_center produces correctly centered matrix", {
  skip_if(!"KDistTest" %in% loadedNamespaces())

  x <- matrix(1:16, 4, 4)
  centered <- matrix_v_center(x)

  # Sum of rows, columns and entire matrix should be approximately zero
  expect_equal(sum(rowSums(centered)), 0, tolerance = 1e-10)
  expect_equal(sum(colSums(centered)), 0, tolerance = 1e-10)
  expect_equal(sum(centered), 0, tolerance = 1e-10)
})
