# Regression tests: every documented `bandwidth` option of tcdt must be
# usable, and invalid inputs must error instead of silently falling back
# to the automatic bandwidth.

make_tcdt_data <- function(n = 50) {
  set.seed(123)
  df <- data.frame(X1 = rnorm(n), Y1 = rnorm(n), X2 = rnorm(n), Y2 = rnorm(n))
  df$Y1 <- df$X1 + rnorm(n)
  df$Y2 <- df$X2^2 + rnorm(n)
  df
}

test_that("tcdt works with all documented bandwidth options", {
  df <- make_tcdt_data()

  set.seed(1)
  r_null <- tcdt(Y1 | X1 ~ Y2 | X2, data = df, B = 20)
  expect_true(is.numeric(r_null$pvalue))

  set.seed(1)
  r_under <- tcdt(Y1 | X1 ~ Y2 | X2, data = df, B = 20,
                  bandwidth = "undersmoothing")
  expect_true(is.numeric(r_under$pvalue))

  # Numeric multiplier: bandwidths must scale exactly by the multiplier
  set.seed(1)
  r_mult <- tcdt(Y1 | X1 ~ Y2 | X2, data = df, B = 20, bandwidth = 2)
  expect_equal(r_mult$bandwidths, 2 * r_null$bandwidths)

  # Custom bandwidth function
  set.seed(1)
  r_fun <- tcdt(Y1 | X1 ~ Y2 | X2, data = df, B = 20,
                bandwidth = function(n, p, nu, A) 0.5 * A)
  expect_true(is.numeric(r_fun$pvalue))
  expect_true(all(r_fun$bandwidths > 0))

  # Fixed bandwidths as a list of two vectors: must actually be used
  set.seed(1)
  r_list <- tcdt(Y1 | X1 ~ Y2 | X2, data = df, B = 20,
                 bandwidth = list(0.3, 0.4))
  expect_equal(as.numeric(r_list$bandwidths), c(0.3, 0.4))
})

test_that("tcdt rejects invalid bandwidth inputs with clear errors", {
  df <- make_tcdt_data()

  expect_error(tcdt(Y1 | X1 ~ Y2 | X2, data = df, B = 10,
                    bandwidth = "undersmooth"), "Invalid 'bandwidth'")
  expect_error(tcdt(Y1 | X1 ~ Y2 | X2, data = df, B = 10,
                    bandwidth = -1), "positive")
  expect_error(tcdt(Y1 | X1 ~ Y2 | X2, data = df, B = 10,
                    bandwidth = c(1, 2)), "single positive value")
  expect_error(tcdt(Y1 | X1 ~ Y2 | X2, data = df, B = 10,
                    bandwidth = list(1, 2, 3)), "list")
  expect_error(tcdt(Y1 | X1 ~ Y2 | X2, data = df, B = 10,
                    bandwidth = TRUE), "Invalid 'bandwidth'")
})

test_that("tcdt validates per-column bandwidth lengths", {
  df <- make_tcdt_data()
  df$X1.2 <- rnorm(nrow(df))
  df$X2.2 <- rnorm(nrow(df))

  # Wrong length (3 bandwidths for 2 columns) must error, not corrupt results
  expect_error(
    tcdt(Y1 | X1 + X1.2 ~ Y2 | X2 + X2.2, data = df, B = 10,
         bandwidth = list(c(0.3, 0.3, 0.3), 0.4)),
    "length"
  )

  # Scalars recycle across columns
  set.seed(1)
  r <- tcdt(Y1 | X1 + X1.2 ~ Y2 | X2 + X2.2, data = df, B = 20,
            bandwidth = list(0.3, 0.4))
  expect_equal(dim(r$bandwidths), c(2L, 2L))
})

test_that("tcdt local testing works with a bandwidth multiplier", {
  df <- make_tcdt_data()
  set.seed(1)
  r <- tcdt(Y1 | X1 ~ Y2 | X2, data = df, B = 20, x0 = 0, bandwidth = 2)
  expect_true(is.numeric(r$pvalue))
  expect_identical(r$problem, "local")
})
