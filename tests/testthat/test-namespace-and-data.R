# Package-surface regression tests: the my_function template leftover must
# stay gone, and the internal asymptotic_stats reference table (R/sysdata.rda)
# must remain resolvable by the change-point functions after removing the
# duplicate copy from data/.

test_that("my_function is no longer exported", {
  expect_false("my_function" %in% getNamespaceExports("KDist"))
})

test_that("kcpd_single asymptotic calibration resolves the internal null table", {
  set.seed(1)
  p <- 5
  cp <- 15
  data <- rbind(matrix(rnorm(cp * p), cp, p),
                matrix(rnorm(cp * p, mean = 1), cp, p))

  r <- kcpd_single(data, type = "e-dist", method = "asymptotic")
  expect_true(is.numeric(r$pvalue))
  expect_gte(r$pvalue, 0)
  expect_lte(r$pvalue, 1)
  expect_true(is.numeric(r$locations))
})
