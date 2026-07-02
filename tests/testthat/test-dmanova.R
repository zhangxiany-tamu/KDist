# Regression tests: dmanova/dmanova2 set options(contrasts = ...) internally
# and must restore the user's global option on every exit path, including
# errors.

make_dmanova_data <- function(n = 40) {
  set.seed(1)
  df <- data.frame(g = factor(rep(1:2, each = n / 2)), z = rnorm(n))
  Y <- matrix(rnorm(n * 5), n, 5)
  Y[df$g == 2, ] <- Y[df$g == 2, ] + 0.5
  list(df = df, Y = Y)
}

test_that("dmanova restores the contrasts option on error and on success", {
  d <- make_dmanova_data()
  saved <- options()$contrasts

  expect_error(dmanova(d$Y ~ nonexistent_var, data = d$df))
  expect_identical(options()$contrasts, saved)

  r <- dmanova(d$Y ~ g | z, data = d$df, method = "asymptotic")
  expect_true(is.numeric(r$aov.tab[1, "Pr(>F)"]))
  expect_identical(options()$contrasts, saved)
})

test_that("dmanova2 restores the contrasts option on error and on success", {
  d <- make_dmanova_data()
  saved <- options()$contrasts

  expect_error(dmanova2(d$Y ~ nonexistent_var, data = d$df))
  expect_identical(options()$contrasts, saved)

  set.seed(2)
  r <- dmanova2(d$Y ~ g | z, data = d$df, n_perm = 49)
  expect_true(is.numeric(r$aov.tab[1, "Pr(>F)"]))
  expect_identical(options()$contrasts, saved)
})
