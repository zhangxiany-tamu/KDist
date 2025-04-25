#' Maximum Mean Discrepancy (MMD) Test
#'
#' Performs a permutation test based on the Maximum Mean Discrepancy to assess whether two samples come from
#' the same distribution. MMD is a powerful nonparametric two-sample test that can detect differences in
#' distributions.
#'
#' @param x First sample (matrix, data frame, or vector)
#' @param y Second sample (matrix, data frame, or vector)
#' @param type Type of kernel or distance to use (default: "gaussian").
#'   Options include "euclidean", "polynomial", "gaussian", "laplacian", "e-dist", "g-dist", or "l-dist".
#' @param bw Bandwidth parameter for kernel functions. If NULL, it will be automatically determined.
#' @param expo Exponent parameter for euclidean distance and polynomial kernel (default: 1).
#' @param scale_factor Scaling factor for automatic bandwidth calculation (default: 0.5).
#' @param group Optional vector specifying group membership for each column.
#'   Used for group-wise distance calculations in "e-dist", "g-dist", or "l-dist".
#' @param u_center Logical; use U-centering instead of V-centering (default: FALSE).
#' @param n_perm Number of permutations to use for the test (default: 1000).
#' @param seed Random seed for reproducibility (default: NULL).
#' @param num_cores Number of cores for parallel computing (default: 1).
#'
#' @return An object of class "mmd_test" containing:
#'   \item{statistic}{MMD test statistic value}
#'   \item{p.value}{Permutation-based p-value}
#'   \item{permutation_values}{Vector of MMD values from permutations}
#'
#' @details
#' The Maximum Mean Discrepancy (MMD) is a measure of the difference between two probability distributions
#' based on the mean embeddings of the distributions in a reproducing kernel Hilbert space (RKHS).
#'
#' For two samples \eqn{X = \{x_1, \ldots, x_n\}} and \eqn{Y = \{y_1, \ldots, y_m\}}, the MMD is defined as:
#'
#' For kernel-based types (gaussian, laplacian, polynomial):
#' \deqn{MMD^2(X, Y) = \frac{1}{n^2}\sum_{i,j=1}^{n}k(x_i, x_j) + \frac{1}{m^2}\sum_{i,j=1}^{m}k(y_i, y_j) - \frac{2}{nm}\sum_{i=1}^{n}\sum_{j=1}^{m}k(x_i, y_j)}
#'
#' For distance-based types (euclidean, e-dist, g-dist, l-dist):
#' \deqn{MMD^2(X, Y) = \frac{2}{nm}\sum_{i=1}^{n}\sum_{j=1}^{m}d(x_i, y_j) - \frac{1}{n^2}\sum_{i,j=1}^{n}d(x_i, x_j) - \frac{1}{m^2}\sum_{i,j=1}^{m}d(y_i, y_j)}
#'
#' where \eqn{k(x, y)} is a kernel function and \eqn{d(x, y)} is a distance function.
#'
#' The null hypothesis is that the two samples come from the same distribution. The alternative hypothesis
#' is that they come from different distributions. The p-value is calculated using a permutation test,
#' where the samples are randomly shuffled to create the null distribution of the test statistic.
#'
#' When u_center = TRUE, the U-statistic version of MMD is used, which provides an unbiased estimate.
#' When u_center = FALSE (default), the V-statistic version is used.
#'
#' @examples
#' # Example 1: Samples from the same distribution
#' set.seed(123)
#' x1 <- matrix(rnorm(100 * 3), ncol = 3)  # 100 observations in 3D from N(0,1)
#' y1 <- matrix(rnorm(80 * 3), ncol = 3)   # 80 observations in 3D from N(0,1)
#' test1 <- mmd_test(x1, y1, type = "gaussian", n_perm = 200)
#' print(test1)
#' plot(test1)
#'
#' # Example 2: Samples from different distributions
#' x2 <- matrix(rnorm(100 * 2), ncol = 2)  # From N(0,1)
#' y2 <- matrix(rnorm(100 * 2, mean = 1), ncol = 2)  # From N(1,1)
#' test2 <- mmd_test(x2, y2, type = "gaussian", n_perm = 200)
#' print(test2)
#'
#' # Example 3: Using a custom bandwidth
#' test3 <- mmd_test(x2, y2, type = "gaussian", bw = 0.5, n_perm = 200)
#' print(test3)
#'
#' # Example 4: Using U-centering (unbiased estimator)
#' test4 <- mmd_test(x2, y2, type = "gaussian", u_center = TRUE, n_perm = 200)
#' print(test4)
#'
#' # Example 5: Using vectors as input
#' x5 <- rnorm(50)
#' y5 <- rnorm(50, mean = 0.5)
#' test5 <- mmd_test(x5, y5, type = "euclidean", n_perm = 200)
#' print(test5)
#'
#' # Example 6: Using parallel computing (if multiple cores are available)
#' \dontrun{
#' test_parallel <- mmd_test(x2, y2, n_perm = 1000, num_cores = 4)
#' print(test_parallel)
#' }
#'
#' @references
#' Gretton, A., Borgwardt, K.M., Rasch, M.J., Schölkopf, B. and Smola, A. (2012).
#' A Kernel Two-Sample Test. Journal of Machine Learning Research, 13, pp.723-773.
#'
#' @seealso
#' \code{\link{mmd}} for calculating MMD without performing a hypothesis test
#' \code{\link{hsic_test}} for the HSIC independence test
#' \code{\link{dhsic_test}} for the d-variable HSIC test
#'
#' @export
mmd_test <- function(x, y, type = "gaussian", bw = NULL, expo = 1,
                     scale_factor = 0.5, group = NULL, u_center = FALSE,
                     n_perm = 1000, seed = NULL, num_cores = 1) {
  # Load parallel package
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("The 'parallel' package is required for this function. Please install it.")
  }

  # Set random seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Convert vectors to matrices if needed
  if (is.vector(x) && !is.matrix(x)) {
    x <- matrix(x, ncol = 1)
  }
  if (is.vector(y) && !is.matrix(y)) {
    y <- matrix(y, ncol = 1)
  }

  # Get dimensions
  n <- nrow(x)
  m <- nrow(y)
  total_n <- n + m

  # Ensure dimensions of x and y match in the second dimension (number of variables)
  if (ncol(x) != ncol(y)) {
    stop("Number of variables in x and y must be the same")
  }

  # Calculate distance matrix once
  D <- KDist_matrix(rbind(x, y), type = type, bw = bw, expo = expo, scale_factor = scale_factor, group = group)

  # Calculate observed test statistic
  observed_mmd <- suppressMessages(mmd(x = D, type = type, u_center = u_center, n = n, m = m))

  # Generate all permutation indices in advance for reproducibility
  all_perm_indices <- lapply(1:n_perm, function(i) {
    sample(total_n, total_n, replace = FALSE)
  })

  # Define the function to run for each permutation
  perm_function <- function(perm_indices) {
    # Reorder the distance matrix according to the permutation
    D_perm <- D[perm_indices, perm_indices]

    # Calculate MMD for the permuted data
    return(suppressMessages(mmd(x = D_perm, type = type, u_center = u_center, n = n, m = m)))
  }

  # Run permutations in parallel
  perm_results <- parallel::mclapply(all_perm_indices, perm_function, mc.cores = num_cores)
  perm_results <- unlist(perm_results)

  # Calculate p-value (proportion of permuted statistics >= observed)
  p.value <- mean(perm_results >= observed_mmd)

  # Create result list
  result <- list(
    statistic = observed_mmd,
    p.value = p.value,
    permutation_values = perm_results
  )

  class(result) <- "mmd_test"
  return(result)
}

#' Hilbert-Schmidt Independence Criterion (HSIC) Test
#'
#' Performs a permutation test based on the Hilbert-Schmidt Independence Criterion to assess
#' independence between two multivariate random variables. HSIC is a powerful nonparametric
#' measure of dependence that can detect both linear and non-linear associations.
#'
#' @param x First dataset (matrix, data frame, or vector) or pre-computed distance/kernel matrix
#' @param y Second dataset (matrix, data frame, or vector) or pre-computed distance/kernel matrix
#' @param type Type of kernel or distance to use (default: "gaussian").
#'   Options include "euclidean", "polynomial", "gaussian", "laplacian", "e-dist", "g-dist", or "l-dist".
#' @param bw Bandwidth parameter for kernel functions. If NULL, it will be automatically determined.
#' @param expo Exponent parameter for euclidean distance and polynomial kernel (default: 1).
#' @param scale_factor Scaling factor for automatic bandwidth calculation (default: 0.5).
#' @param group Optional list of length 2 specifying group membership for each column in x and y.
#'   The first element of the list should specify grouping for columns in x, and the
#'   second element should specify grouping for columns in y. Used for group-wise
#'   distance calculations in "e-dist", "g-dist", or "l-dist".
#' @param u_center Logical; use U-centering instead of V-centering (default: FALSE).
#'   U-centering provides an unbiased estimate of HSIC by excluding diagonal terms in the kernel matrices.
#' @param n_perm Number of permutations to use for the test (default: 1000).
#' @param seed Random seed for reproducibility (default: NULL).
#' @param num_cores Number of cores for parallel computing (default: 1).
#' @param is_distance Logical; whether input matrices x and y are already distance/kernel matrices (default: FALSE).
#'
#' @return An object of class "hsic_test" containing:
#'   \item{statistic}{HSIC test statistic value}
#'   \item{p.value}{Permutation-based p-value}
#'   \item{permutation_values}{Vector of HSIC values from permutations}
#'
#' @details
#' The Hilbert-Schmidt Independence Criterion (HSIC) measures the dependency between two random variables
#' by computing the squared Hilbert-Schmidt norm of the cross-covariance operator in reproducing kernel
#' Hilbert spaces (RKHS).
#'
#' For two samples \eqn{X = \{x_1, \ldots, x_n\}} and \eqn{Y = \{y_1, \ldots, y_n\}}, the HSIC is defined as:
#'
#' \deqn{HSIC(X, Y) = \frac{1}{n^2} \text{trace}(KHLH)}
#'
#' where \eqn{K} and \eqn{L} are kernel matrices for \eqn{X} and \eqn{Y} respectively, and \eqn{H} is the
#' centering matrix \eqn{H = I - \frac{1}{n}ee^T} with \eqn{e} being a vector of ones.
#'
#' For kernel-based types (gaussian, laplacian, polynomial), a higher value indicates stronger dependency.
#' For distance-based types (euclidean, e-dist, g-dist, l-dist), the interpretation may vary.
#'
#' The null hypothesis is that \eqn{X} and \eqn{Y} are independent. The alternative hypothesis
#' is that they are dependent. The p-value is calculated using a permutation test,
#' where one sample is randomly permuted while keeping the other fixed to create the null distribution.
#'
#' U-centering vs. V-centering:
#'
#' The HSIC can be estimated using either a V-statistic or a U-statistic:
#'
#' 1. V-statistic (u_center = FALSE, default):
#'    \deqn{HSIC_V(X, Y) = \frac{1}{n^2}\sum_{i,j=1}^{n} k(x_i,x_j)l(y_i,y_j) + \frac{1}{n^4}\sum_{i,j=1}^{n} k(x_i,x_j)\sum_{i,j=1}^{n} l(y_i,y_j) - \frac{2}{n^3}\sum_{i,j,q=1}^{n} k(x_i,x_j)l(y_i,y_q)}
#'    This is a biased estimator but often has better finite-sample properties.
#'
#' 2. U-statistic (u_center = TRUE):
#'    \deqn{HSIC_U(X, Y) = \frac{1}{(n)_2}\sum_{i \neq j} k(x_i,x_j)l(y_i,y_j) + \frac{1}{(n)_4}\sum_{i \neq j \neq q \neq r} k(x_i,x_j)l(y_q,y_r) -\frac{2}{(n)_3}\sum_{i \neq j \neq q} k(x_i,x_j)l(y_i,y_q)}
#'    where \eqn{(n)_k=n!/(n-k)!}. This removes the diagonal elements in the kernel matrices, providing an unbiased estimator of HSIC.
#'    U-centering is particularly important for small sample sizes, as it reduces estimation bias.
#'
#' When is_distance = TRUE, x and y are assumed to be pre-computed distance/kernel matrices. In this case,
#' the parameters bw, expo, scale_factor, and group are ignored.
#'
#' @examples
#' # Example 1: Independent variables
#' set.seed(123)
#' x1 <- matrix(rnorm(100), ncol = 1)
#' y1 <- matrix(rnorm(100), ncol = 1)
#' test1 <- hsic_test(x1, y1, type = "gaussian", n_perm = 200)
#' print(test1)
#' plot(test1)
#'
#' # Example 2: Linear relationship
#' x2 <- matrix(runif(100), ncol = 1)
#' y2 <- matrix(2*x2 + 0.1*rnorm(100), ncol = 1)  # Linear relationship with noise
#' test2 <- hsic_test(x2, y2, type = "gaussian", n_perm = 200)
#' print(test2)
#'
#' # Example 3: Non-linear relationship
#' x3 <- matrix(runif(100, -3, 3), ncol = 1)
#' y3 <- matrix(sin(x3) + 0.1*rnorm(100), ncol = 1)  # Sinusoidal relationship
#' test3 <- hsic_test(x3, y3, type = "gaussian", n_perm = 200)
#' print(test3)
#'
#' # Example 4: Using different kernel types
#' test_gaussian <- hsic_test(x3, y3, type = "gaussian", n_perm = 200)
#' test_laplacian <- hsic_test(x3, y3, type = "laplacian", n_perm = 200)
#' test_euclidean <- hsic_test(x3, y3, type = "euclidean", n_perm = 200)
#' print(c(gaussian = test_gaussian$p.value,
#'         laplacian = test_laplacian$p.value,
#'         euclidean = test_euclidean$p.value))
#'
#' # Example 5: Comparing V-statistic and U-statistic
#' # For a small sample size, the difference might be more noticeable
#' x_small <- matrix(rnorm(30), ncol = 1)
#' y_small <- matrix(0.5*x_small + rnorm(30, sd = 0.5), ncol = 1)
#'
#' test_v <- hsic_test(x_small, y_small, type = "gaussian", u_center = FALSE, n_perm = 200)
#' test_u <- hsic_test(x_small, y_small, type = "gaussian", u_center = TRUE, n_perm = 200)
#'
#' print("V-statistic (biased):")
#' print(test_v)
#' print("U-statistic (unbiased):")
#' print(test_u)
#'
#' # Example 6: Using multivariate data
#' x_multi <- matrix(rnorm(200), ncol = 2)
#' y_multi <- matrix(cbind(x_multi[,1] + rnorm(100, sd = 0.1),
#'                         x_multi[,2]^2 + rnorm(100, sd = 0.1)), ncol = 2)
#' test_multi <- hsic_test(x_multi, y_multi, type = "gaussian", n_perm = 200)
#' print(test_multi)
#'
#' # Example 7: Using pre-computed distance matrices
#' Dx <- KDist_matrix(x_multi, type = "gaussian")
#' Dy <- KDist_matrix(y_multi, type = "gaussian")
#' test_precomp <- hsic_test(Dx, Dy, type = "gaussian", is_distance = TRUE, n_perm = 200)
#' print(test_precomp)
#'
#' # Example 8: Using grouped variables with hsic_test
#' # Create sample data
#' set.seed(123)
#' n <- 100
#' x <- matrix(rnorm(n*4), ncol = 4)  # 4 variables in x
#' y <- matrix(rnorm(n*3), ncol = 3)  # 3 variables in y
#'
#' # Create a dependency between the variables in the first group
#' y[, 1:2] <- y[, 1:2] + 0.8 * x[, 1:2]
#'
#' # Define group structure for both datasets
#' x_groups <- c(1, 1, 2, 2)    # First 2 vars in group 1, last 2 in group 2
#' y_groups <- c(1, 1, 2)       # First 2 vars in group 1, last 1 in group 2
#'
#' # Combine into a list for group-specific analysis
#' group_list <- list(x_groups, y_groups)
#'
#' # Perform HSIC test with the group structure
#' test_grouped <- hsic_test(x, y, type = "e-dist", group = group_list, n_perm = 200)
#' print(test_grouped)
#'
#' # Compare with standard (non-grouped) HSIC test
#' test_standard <- hsic_test(x, y, type = "gaussian", n_perm = 200)
#' print(test_standard)
#'
#' # Example 9: Using parallel computing (if multiple cores are available)
#' \dontrun{
#' test_parallel <- hsic_test(x_multi, y_multi, n_perm = 1000, num_cores = 4)
#' print(test_parallel)
#' }
#'
#' @references
#' Gretton, A., Fukumizu, K., Teo, C. H., Song, L., Schölkopf, B., & Smola, A. J. (2007).
#' A kernel statistical test of independence. Advances in Neural Information Processing Systems, 20.
#'
#' @seealso
#' \code{\link{hsic}} for calculating HSIC without performing a hypothesis test
#' \code{\link{mmd_test}} for the MMD two-sample test
#' \code{\link{dhsic_test}} for the d-variable HSIC test
#' \code{\link{dcov_test}} for the distance covariance test (when type="euclidean")
#'
#' @export
hsic_test <- function(x, y, type = "gaussian", bw = NULL, expo = 1,
                      scale_factor = 0.5, group = NULL, u_center = FALSE,
                      n_perm = 1000, seed = NULL, num_cores = 1, is_distance = FALSE) {
  # Load parallel package
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("The 'parallel' package is required for this function. Please install it.")
  }

  # Set random seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!is_distance) {
    # Convert vectors to matrices if needed
    if (is.vector(x) && !is.matrix(x)) {
      x <- matrix(x, ncol = 1)
    }

    if (is.vector(y) && !is.matrix(y)) {
      y <- matrix(y, ncol = 1)
    }

    # Ensure dimensions of x and y match in the first dimension (number of observations)
    if (nrow(x) != nrow(y)) {
      stop("Number of observations in x and y must be the same")
    }

    # Handle group parameter appropriately for KDist_matrix
    group_x <- NULL
    group_y <- NULL

    # If group is a list of size 2, extract the appropriate group vectors
    if (is.list(group) && length(group) == 2) {
      group_x <- group[[1]]
      group_y <- group[[2]]
    } else {
      # If group is not a list or not of size 2, use the same group for both
      group_x <- group
      group_y <- group
    }

    # Calculate distance matrices with appropriate group parameters
    Dx <- KDist_matrix(x, type = type, bw = bw, expo = expo, scale_factor = scale_factor, group = group_x)
    Dy <- KDist_matrix(y, type = type, bw = bw, expo = expo, scale_factor = scale_factor, group = group_y)
  } else {
    # Use x and y directly as distance matrices
    Dx <- x
    Dy <- y

    # Check if matrices are square and of the same size
    if (!is.matrix(Dx) || !is.matrix(Dy) || nrow(Dx) != ncol(Dx) || nrow(Dy) != ncol(Dy)) {
      stop("When is_distance=TRUE, x and y must be square matrices")
    }

    if (nrow(Dx) != nrow(Dy)) {
      stop("Distance matrices x and y must have the same dimensions")
    }
  }

  # Get number of observations
  n <- nrow(Dx)

  # Calculate observed test statistic
  observed_hsic <- suppressMessages(
    hsic(Dx, Dy, type = type, u_center = u_center, is_distance = TRUE)
  )

  # Define the function to run for each permutation
  perm_function <- function(i) {
    # Generate a random permutation of indices
    perm_indices <- sample(n, n, replace = FALSE)

    # Reorder the second distance matrix according to the permutation
    Dy_perm <- Dy[perm_indices, perm_indices]

    # Calculate HSIC for the permuted data
    return(suppressMessages(
      hsic(Dx, Dy_perm, type = type, u_center = u_center, is_distance = TRUE)
    ))
  }

  # Run permutations in parallel
  perm_results <- parallel::mclapply(1:n_perm, function(i) perm_function(i), mc.cores = num_cores)
  perm_results <- unlist(perm_results)

  # Calculate p-value (proportion of permuted statistics >= observed)
  p.value <- mean(perm_results >= observed_hsic)

  # Create result list
  result <- list(
    statistic = observed_hsic,
    p.value = p.value,
    permutation_values = perm_results
  )

  class(result) <- "hsic_test"
  return(result)
}

#' Energy Distance Test
#'
#' Performs a permutation test for the equality of distributions based on energy distance.
#' This is a powerful nonparametric two-sample test that measures the difference between
#' probability distributions.
#'
#' @param x First sample (matrix, data frame, or vector) or full distance matrix if y is NULL
#' @param y Second sample (matrix, data frame, or vector)
#' @param a Exponent parameter for Euclidean distance (default: 1)
#' @param u_center Logical; use U-centering instead of V-centering (default: FALSE).
#'   U-centering provides an unbiased estimate by excluding diagonal terms in the distance matrices.
#' @param n_perm Number of permutations to use for the test (default: 1000)
#' @param seed Random seed for reproducibility (default: NULL)
#' @param num_cores Number of cores for parallel computing (default: 1)
#'
#' @return An object of class "ed_test" containing:
#'   \item{statistic}{Energy distance test statistic value}
#'   \item{p.value}{Permutation-based p-value}
#'   \item{permutation_values}{Vector of energy distance values from permutations}
#'
#' @details
#' The energy distance between two probability distributions F and G is defined as:
#'
#' \deqn{E(F, G) = 2E||X - Y||^a - E||X - X'||^a - E||Y - Y'||^a}
#'
#' where \eqn{X, X'} are independent random variables with distribution \eqn{F}, and \eqn{Y, Y'} are
#' independent random variables with distribution \eqn{G}, \eqn{||\cdot||} denotes Euclidean distance,
#' and \eqn{a} is the exponent parameter (usually set to 1).
#'
#' Energy distance has the property that \eqn{E(F, G) ≥ 0}, with equality if and only if \eqn{F = G}.
#' This makes it suitable for testing whether two samples come from the same distribution.
#'
#' The null hypothesis is that the two samples come from the same distribution. The alternative
#' hypothesis is that they come from different distributions. The p-value is calculated
#' using a permutation test, where the samples are randomly shuffled to create the null
#' distribution of the test statistic.
#'
#' U-centering vs. V-centering:
#'
#' The energy distance can be estimated using either a V-statistic or a U-statistic:
#'
#' 1. V-statistic (u_center = FALSE, default):
#'    \deqn{E_V(X, Y) = \frac{2}{nm}\sum_{i=1}^{n}\sum_{j=1}^{m}||x_i - y_j||^a - \frac{1}{n^2}\sum_{i=1}^{n}\sum_{j=1}^{n}||x_i - x_j||^a - \frac{1}{m^2}\sum_{i=1}^{m}\sum_{j=1}^{m}||y_i - y_j||^a}
#'    This is a biased estimator but often has better finite-sample properties.
#'
#' 2. U-statistic (u_center = TRUE):
#'    \deqn{E_U(X, Y) = \frac{2}{nm}\sum_{i=1}^{n}\sum_{j=1}^{m}||x_i - y_j||^a - \frac{1}{n(n-1)}\sum_{i \neq j}||x_i - x_j||^a - \frac{1}{m(m-1)}\sum_{i \neq j}||y_i - y_j||^a}
#'    This removes self-distances (i.e., diagonal elements in the distance matrices), providing an
#'    unbiased estimator of energy distance. U-centering is particularly important for small sample sizes,
#'    as it reduces estimation bias.
#'
#' Note: `ed_test` is implemented as a wrapper around `mmd_test` with type="euclidean", as energy
#' distance is equivalent to Maximum Mean Discrepancy with a specific kernel.
#'
#' @examples
#' # Example 1: Samples from the same distribution
#' set.seed(123)
#' x1 <- matrix(rnorm(100 * 3), ncol = 3)  # 100 observations in 3D from N(0,1)
#' y1 <- matrix(rnorm(80 * 3), ncol = 3)   # 80 observations in 3D from N(0,1)
#' test1 <- ed_test(x1, y1, n_perm = 200)
#' print(test1)
#' plot(test1)
#'
#' # Example 2: Samples from different distributions
#' x2 <- matrix(rnorm(100 * 2), ncol = 2)  # From N(0,1)
#' y2 <- matrix(rnorm(100 * 2, mean = 1), ncol = 2)  # From N(1,1)
#' test2 <- ed_test(x2, y2, n_perm = 200)
#' print(test2)
#'
#' # Example 3: Different exponent values
#' test_a1 <- ed_test(x2, y2, a = 1, n_perm = 200)  # Standard energy distance
#' test_a2 <- ed_test(x2, y2, a = 2, n_perm = 200)  # Squared energy distance
#' print(c(a1 = test_a1$p.value, a2 = test_a2$p.value))
#'
#' # Example 4: Comparing V-statistic and U-statistic
#' # For small sample sizes, the difference might be more noticeable
#' x_small <- matrix(rnorm(30, 0, 1), ncol = 1)
#' y_small <- matrix(rnorm(30, 0.5, 1), ncol = 1)
#'
#' test_v <- ed_test(x_small, y_small, u_center = FALSE, n_perm = 200)
#' test_u <- ed_test(x_small, y_small, u_center = TRUE, n_perm = 200)
#'
#' print("V-statistic (biased):")
#' print(test_v)
#' print("U-statistic (unbiased):")
#' print(test_u)
#'
#' # Example 5: Using vector inputs
#' x5 <- rnorm(50)
#' y5 <- rnorm(50, mean = 0.5)
#' test5 <- ed_test(x5, y5, n_perm = 200)
#' print(test5)
#'
#' # Example 6: Using parallel computing (if multiple cores are available)
#' \dontrun{
#' test_parallel <- ed_test(x2, y2, n_perm = 1000, num_cores = 4)
#' print(test_parallel)
#' }
#'
#' @references
#' Székely, G. J., & Rizzo, M. L. (2004). Testing for Equal Distributions in High Dimension.
#' InterStat, Nov. (5).
#'
#' Székely, G. J., & Rizzo, M. L. (2013). Energy statistics: A class of statistics based on
#' distances. Journal of Statistical Planning and Inference, 143(8), 1249-1272.
#'
#' Baringhaus, L., & Franz, C. (2004). On a new multivariate two-sample test.
#' Journal of Multivariate Analysis, 88(1), 190-206.
#'
#' @seealso
#' \code{\link{ed}} for calculating energy distance without performing a hypothesis test
#' \code{\link{mmd_test}} for the general Maximum Mean Discrepancy test
#' \code{\link{dcov_test}} for the distance covariance test of independence
#'
#' @export
ed_test <- function(x, y, a = 1, u_center = FALSE, n_perm = 1000,
                    seed = NULL, num_cores = 1) {
  # Call mmd_test with type="euclidean"
  result <- mmd_test(x = x, y = y, type = "euclidean", expo = a,
                     u_center = u_center, n_perm = n_perm,
                     seed = seed, num_cores = num_cores)

  # Update class to indicate this is an ED test
  class(result) <- c("ed_test", class(result))

  return(result)
}

#' Distance Covariance Test
#'
#' Performs a permutation test based on the distance covariance to assess independence between
#' two multivariate random variables. Distance covariance is a powerful nonparametric measure
#' of dependence that can detect both linear and non-linear associations.
#'
#' @param x First dataset (matrix, data frame, or vector) or pre-computed distance matrix
#' @param y Second dataset (matrix, data frame, or vector) or pre-computed distance matrix
#' @param a Exponent parameter for Euclidean distance (default: 1)
#' @param u_center Logical; use U-centering instead of V-centering (default: FALSE).
#'   U-centering provides an unbiased estimate by excluding diagonal terms in the distance matrices.
#' @param n_perm Number of permutations to use for the test (default: 1000)
#' @param seed Random seed for reproducibility (default: NULL)
#' @param num_cores Number of cores for parallel computing (default: 1)
#' @param is_distance Logical; whether input matrices x and y are already distance matrices (default: FALSE)
#'
#' @return An object of class "dcov_test" containing:
#'   \item{statistic}{Distance covariance test statistic value}
#'   \item{p.value}{Permutation-based p-value}
#'   \item{permutation_values}{Vector of distance covariance values from permutations}
#'
#' @details
#' The distance covariance between two random vectors \eqn{X} and \eqn{Y} is a measure of the
#' dependence between them. It is defined as:
#'
#' \deqn{dCov^2(X, Y) = E[||X - X'||^a \cdot ||Y - Y'||^a] + E[||X - X'||^a] \cdot E[||Y - Y'||^a] - 2E[||X - X'||^a \cdot ||Y - Y''||^a]}
#'
#' where \eqn{||.||} denotes Euclidean distance, \eqn{(X', Y')} and \eqn{(X'', Y'')} are
#' independent copies of \eqn{(X, Y)}, and \eqn{a} is the exponent parameter (usually set to 1).
#'
#' The distance covariance has the property that \eqn{dCov^2(X, Y) \geq 0}, with equality if and
#' only if \eqn{X} and \eqn{Y} are independent. This makes it suitable for testing independence
#' between random variables.
#'
#' The null hypothesis is that \eqn{X} and \eqn{Y} are independent. The alternative hypothesis
#' is that they are dependent. The p-value is calculated using a permutation test,
#' where one sample is randomly permuted while keeping the other fixed to create the null distribution.
#'
#' U-centering vs. V-centering:
#'
#' The distance covariance can be estimated using either a V-statistic or a U-statistic:
#'
#' 1. V-statistic (u_center = FALSE, default):
#'    \deqn{dCov_V^2(X, Y) = \frac{1}{n^2}\sum_{i,j=1}^{n} a_{ij}b_{ij}}
#'    where \eqn{a_{ij}} and \eqn{b_{ij}} are the \eqn{(i,j)}th element of the doubly-centered distance matrices.
#'
#' 2. U-statistic (u_center = TRUE):
#'    \deqn{dCov_U^2(X, Y) = \frac{1}{n(n-3)}\sum_{i \neq j} \tilde{a}_{ij}\tilde{b}_{ij}}
#'    where \eqn{\tilde{a}_{ij}} and \eqn{\tilde{b}_{ij}} are the \eqn{(i,j)}th element of the U-centered distance matrices.
#'    This removes diagonal elements in the distance matrices, providing an unbiased estimator
#'    of distance covariance. U-centering is particularly important for small sample sizes,
#'    as it reduces estimation bias.
#'
#' Note: `dcov_test` is implemented as a wrapper around `hsic_test` with type="euclidean", as distance
#' covariance is equivalent to HSIC with a specific kernel.
#'
#' When is_distance = TRUE, x and y are assumed to be pre-computed distance matrices. In this case,
#' the parameter a is ignored.
#'
#' @examples
#' # Example 1: Independent variables
#' set.seed(123)
#' x1 <- matrix(rnorm(100), ncol = 1)
#' y1 <- matrix(rnorm(100), ncol = 1)
#' test1 <- dcov_test(x1, y1, n_perm = 200)
#' print(test1)
#' plot(test1)
#'
#' # Example 2: Linear relationship
#' x2 <- matrix(runif(100), ncol = 1)
#' y2 <- matrix(2*x2 + 0.1*rnorm(100), ncol = 1)  # Linear relationship with noise
#' test2 <- dcov_test(x2, y2, n_perm = 200)
#' print(test2)
#'
#' # Example 3: Non-linear relationship
#' x3 <- matrix(runif(100, -3, 3), ncol = 1)
#' y3 <- matrix(sin(x3) + 0.1*rnorm(100), ncol = 1)  # Sinusoidal relationship
#' test3 <- dcov_test(x3, y3, n_perm = 200)
#' print(test3)
#'
#' # Example 4: Different exponent values
#' test_a1 <- dcov_test(x3, y3, a = 1, n_perm = 200)  # Standard distance covariance
#' test_a2 <- dcov_test(x3, y3, a = 2, n_perm = 200)  # Squared distance covariance
#' print(c(a1 = test_a1$p.value, a2 = test_a2$p.value))
#'
#' # Example 5: Comparing V-statistic and U-statistic
#' # For small sample sizes, the difference might be more noticeable
#' x_small <- matrix(rnorm(30), ncol = 1)
#' y_small <- matrix(0.5*x_small + rnorm(30, sd = 0.5), ncol = 1)
#'
#' test_v <- dcov_test(x_small, y_small, u_center = FALSE, n_perm = 200)
#' test_u <- dcov_test(x_small, y_small, u_center = TRUE, n_perm = 200)
#'
#' print("V-statistic (biased):")
#' print(test_v)
#' print("U-statistic (unbiased):")
#' print(test_u)
#'
#' # Example 6: Using multivariate data
#' x_multi <- matrix(rnorm(200), ncol = 2)
#' y_multi <- matrix(cbind(x_multi[,1] + rnorm(100, sd = 0.1),
#'                         x_multi[,2]^2 + rnorm(100, sd = 0.1)), ncol = 2)
#' test_multi <- dcov_test(x_multi, y_multi, n_perm = 200)
#' print(test_multi)
#'
#' # Example 7: Using pre-computed distance matrices
#' Dx <- as.matrix(dist(x_multi))  # Euclidean distance matrix
#' Dy <- as.matrix(dist(y_multi))  # Euclidean distance matrix
#' test_precomp <- dcov_test(Dx, Dy, is_distance = TRUE, n_perm = 200)
#' print(test_precomp)
#'
#' # Example 8: Using parallel computing (if multiple cores are available)
#' \dontrun{
#' test_parallel <- dcov_test(x_multi, y_multi, n_perm = 1000, num_cores = 4)
#' print(test_parallel)
#' }
#'
#' @references
#' Székely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing dependence by
#' correlation of distances. The Annals of Statistics, 35(6), 2769-2794.
#'
#' Székely, G. J., & Rizzo, M. L. (2009). Brownian distance covariance. The Annals of Applied
#' Statistics, 3(4), 1236-1265.
#'
#' @seealso
#' \code{\link{dcov}} for calculating distance covariance without performing a hypothesis test
#' \code{\link{dcor}} for calculating distance correlation
#' \code{\link{hsic_test}} for the general HSIC independence test
#' \code{\link{ed_test}} for the energy distance two-sample test
#'
#' @export
dcov_test <- function(x, y, a = 1, u_center = FALSE, n_perm = 1000,
                      seed = NULL, num_cores = 1, is_distance = FALSE) {
  # Call hsic_test with type="euclidean"
  result <- hsic_test(x = x, y = y, type = "euclidean", expo = a,
                      u_center = u_center, n_perm = n_perm,
                      seed = seed, num_cores = num_cores,
                      is_distance = is_distance)

  # Update class to indicate this is a dCov test
  class(result) <- c("dcov_test", class(result))

  return(result)
}

#' d-variable HSIC Independence Test (dHSIC)
#'
#' Performs a permutation test based on the d-variable Hilbert-Schmidt Independence Criterion
#' to assess mutual independence among multiple random variables. This implementation follows
#' the methodology described in Pfister et al. (2018).
#'
#' @param x List of matrices or vectors, where each element corresponds to a random variable
#' @param type Type of kernel or distance to use (default: "gaussian").
#'   Options include "euclidean", "polynomial", "gaussian", "laplacian", "e-dist", "g-dist", or "l-dist".
#' @param bw Bandwidth parameter for kernel functions. If NULL, it will be automatically determined.
#' @param expo Exponent parameter for Euclidean distance and polynomial kernel (default: 1).
#' @param scale_factor Scaling factor for automatic bandwidth calculation (default: 0.5).
#' @param group Optional list specifying group membership for each column in the corresponding input list elements.
#'   The length of the group list should match the length of the input list x, where each element
#'   contains the grouping information for the respective element in x. Used for group-wise
#'   distance calculations in "e-dist", "g-dist", or "l-dist".
#' @param n_perm Number of permutations to use for the test (default: 1000).
#' @param seed Random seed for reproducibility (default: NULL).
#' @param num_cores Number of cores for parallel computing (default: 1).
#'
#' @return An object of class "dhsic_test" containing:
#'   \item{statistic}{dHSIC test statistic value}
#'   \item{p.value}{Permutation-based p-value}
#'   \item{permutation_values}{Vector of dHSIC values from permutations}
#'
#' @details
#' The d-variable Hilbert-Schmidt Independence Criterion (dHSIC) extends the pairwise HSIC to
#' measure mutual independence among d random variables. While pairwise independence tests
#' can detect dependence between pairs of variables, they may miss higher-order dependencies
#' among three or more variables. The dHSIC addresses this limitation by providing a single
#' test statistic for mutual independence.
#'
#' For d random variables \eqn{X^{(1)}, X^{(2)}, \ldots, X^{(d)}}, the dHSIC is defined as:
#'
#' \deqn{dHSIC(X^{(1)}, \ldots, X^{(d)}) = ||\otimes_{j=1}^{d} \mu_{X^{(j)}} - \mu_{X^{(1)},\ldots,X^{(d)}}||^2_{\mathcal{H}}}
#'
#' where \eqn{\mu_{X^{(j)}}} is the mean embedding of the distribution of \eqn{X^{(j)}} in a
#' reproducing kernel Hilbert space (RKHS), \eqn{\mu_{X^{(1)},\ldots,X^{(d)}}} is the joint
#' embedding, and \eqn{\otimes} denotes the tensor product.
#'
#' The null hypothesis is that all random variables are mutually independent. The alternative
#' hypothesis is that at least one variable depends on the others in some way. The p-value is
#' calculated using a permutation test, where all variables except the first one are randomly
#' permuted to create the null distribution.
#'
#' The permutation scheme ensures that under the null distribution, all variables are independent,
#' while preserving the marginal distributions of each variable.
#'
#' @examples
#' # Example 1: Mutually independent variables
#' set.seed(123)
#' x1 <- matrix(rnorm(100), ncol = 1)
#' x2 <- matrix(rnorm(100), ncol = 1)
#' x3 <- matrix(rnorm(100), ncol = 1)
#' test1 <- dhsic_test(list(x1, x2, x3), type = "gaussian", n_perm = 200)
#' print(test1)
#' plot(test1)
#'
#' # Example 2: Variables with pairwise independence but mutual dependence
#' set.seed(456)
#' u <- runif(100, -pi, pi)
#' x1 <- matrix(sin(u), ncol = 1)
#' x2 <- matrix(cos(u), ncol = 1)
#' x3 <- matrix(sin(u) * cos(u), ncol = 1)
#' test2 <- dhsic_test(list(x1, x2, x3), type = "gaussian", n_perm = 200)
#' print(test2)
#'
#' # Example 3: Using different kernel types
#' test_gaussian <- dhsic_test(list(x1, x2, x3), type = "gaussian", n_perm = 200)
#' test_laplacian <- dhsic_test(list(x1, x2, x3), type = "laplacian", n_perm = 200)
#' test_euclidean <- dhsic_test(list(x1, x2, x3), type = "euclidean", n_perm = 200)
#' print(c(gaussian = test_gaussian$p.value,
#'         laplacian = test_laplacian$p.value,
#'         euclidean = test_euclidean$p.value))
#'
#' # Example 4: Using multivariate data
#' x1_multi <- matrix(rnorm(200), ncol = 2)
#' x2_multi <- matrix(rnorm(200), ncol = 2)
#' x3_multi <- matrix(rnorm(200), ncol = 2)
#' test_multi <- dhsic_test(list(x1_multi, x2_multi, x3_multi), n_perm = 200)
#' print(test_multi)
#'
#' # Example 5: Testing with correlated variables
#' x1 <- matrix(rnorm(100), ncol = 1)
#' x2 <- matrix(0.7 * x1 + 0.3 * rnorm(100), ncol = 1)  # Correlated with x1
#' x3 <- matrix(rnorm(100), ncol = 1)  # Independent of both
#' test_corr <- dhsic_test(list(x1, x2, x3), n_perm = 200)
#' print(test_corr)
#'
#' # Example 6: Using grouped variables with dhsic_test
#' # Create sample data with 3 datasets
#' set.seed(123)
#' n <- 100
#' x1 <- matrix(rnorm(n*3), ncol = 3)  # 3 variables in x1
#' x2 <- matrix(rnorm(n*4), ncol = 4)  # 4 variables in x2
#' x3 <- matrix(rnorm(n*2), ncol = 2)  # 2 variables in x3
#'
#' # Create dependencies between datasets (first groups are related)
#' z <- matrix(rnorm(n), ncol = 1)  # Shared latent factor
#' x1[, 1:2] <- x1[, 1:2] + z %*% matrix(1, 1, 2)
#' x2[, 1:2] <- x2[, 1:2] + z %*% matrix(1, 1, 2)
#' x3[, 1] <- x3[, 1] + z
#'
#' # Define group structure for each dataset
#' groups1 <- c(1, 1, 2)        # First 2 vars in group 1, last 1 in group 2
#' groups2 <- c(1, 1, 2, 2)     # First 2 vars in group 1, last 2 in group 2
#' groups3 <- c(1, 2)           # Each var in its own group
#'
#' # Combine into a list matching the input list structure
#' group_list <- list(groups1, groups2, groups3)
#'
#' # Perform dHSIC test with the group structure
#' test_grouped <- dhsic_test(list(x1, x2, x3), type = "e-dist",
#'                           group = group_list, n_perm = 200)
#' print(test_grouped)
#'
#' # Compare with standard (non-grouped) dHSIC test
#' test_standard <- dhsic_test(list(x1, x2, x3), type = "gaussian", n_perm = 200)
#' print(test_standard)
#'
#' # Example 7: Using parallel computing (if multiple cores are available)
#' \dontrun{
#' test_parallel <- dhsic_test(list(x1, x2, x3), n_perm = 1000, num_cores = 4)
#' print(test_parallel)
#' }
#'
#' @references
#' Pfister, N., Bühlmann, P., Schölkopf, B., & Peters, J. (2018).
#' Kernel-based tests for joint independence. Journal of the Royal Statistical Society:
#' Series B (Statistical Methodology), 80(1), 5-31.
#'
#' @seealso
#' \code{\link{dhsic}} for calculating dHSIC without performing a hypothesis test
#' \code{\link{hsic_test}} for the pairwise HSIC independence test
#' \code{\link{jhsic_test}} for an alternative joint independence test
#'
#' @export
dhsic_test <- function(x, type = "gaussian", bw = NULL, expo = 1,
                       scale_factor = 0.5, group = NULL,
                       n_perm = 1000, seed = NULL, num_cores = 1) {
  # Load parallel package
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("The 'parallel' package is required for this function. Please install it.")
  }

  # Set random seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Validate inputs
  if (!is.list(x)) {
    stop("Input 'x' must be a list of matrices or vectors")
  }

  if (length(x) < 2) {
    stop("At least two matrices are required for dHSIC")
  }

  # Convert all vectors to matrices
  x_matrices <- lapply(x, function(xi) {
    if(is.null(dim(xi)) || length(dim(xi)) == 1) {
      matrix(xi, ncol = 1)
    } else {
      as.matrix(xi)
    }
  })

  # Get dimensions
  d <- length(x_matrices)
  n <- nrow(x_matrices[[1]])

  # Ensure all matrices have the same number of rows
  for (i in 2:d) {
    if (nrow(x_matrices[[i]]) != n) {
      stop("All matrices must have the same number of rows")
    }
  }

  # Calculate observed test statistic
  observed_dhsic <- dhsic(x_matrices, type = type, bw = bw, expo = expo,
                         scale_factor = scale_factor, group = group)

  # Define the function to run for each permutation
  perm_function <- function(i) {
    # Create a permuted version of the list
    x_perm <- x_matrices

    # Permute all matrices except the first one
    for (j in 2:d) {
      perm_indices <- sample(n, n, replace = FALSE)
      x_perm[[j]] <- x_matrices[[j]][perm_indices, , drop = FALSE]
    }

    # Calculate dHSIC for the permuted data
    return(dhsic(x_perm, type = type, bw = bw, expo = expo,
                scale_factor = scale_factor, group = group))
  }

  # Run permutations in parallel
  perm_results <- parallel::mclapply(1:n_perm, function(i) perm_function(i), mc.cores = num_cores)
  perm_results <- unlist(perm_results)

  # Calculate p-value (proportion of permuted statistics >= observed)
  p.value <- mean(perm_results >= observed_dhsic)

  # Create result list
  result <- list(
    statistic = observed_dhsic,
    p.value = p.value,
    permutation_values = perm_results
  )

  class(result) <- "dhsic_test"
  return(result)
}

#' Joint HSIC and dCov Independence Test
#'
#' Performs a permutation test based on the joint Hilbert-Schmidt Independence Criterion and joint
#' covariance to assess mutual independence among multiple random variables. This implementation follows the methodology
#' described in Chakraborty & Zhang (2019), which extends traditional pairwise independence measures
#' to a joint measure capable of detecting complex dependence structures among multiple variables.
#'
#' @param x List of matrices or vectors, where each element corresponds to a random variable
#' @param cc Constant parameter (default: 1)
#' @param type Type of kernel or distance to use (default: "gaussian").
#'   Options include "euclidean", "polynomial", "gaussian", "laplacian", "e-dist", "g-dist", or "l-dist".
#' @param stat_type Type of statistic ("V", "U", "US", or "UR")
#'   \itemize{
#'     \item "V": V-statistic
#'     \item "U": U-statistic (partially unbiased)
#'     \item "US": Scale-free U-statistic
#'     \item "UR": Rank-based U-statistic (robust to monotonic transformations)
#'   }
#' @param bw Bandwidth parameter for kernel functions. If NULL, it will be automatically determined.
#' @param expo Exponent parameter for euclidean distance and polynomial kernel (default: 1).
#' @param scale_factor Scaling factor for automatic bandwidth calculation (default: 0.5).
#' @param group Optional list specifying group membership for each column in the corresponding input list elements.
#'   The length of the group list should match the length of the input list x, where each element
#'   contains the grouping information for the respective element in x. Used for group-wise
#'   distance calculations in "e-dist", "g-dist", or "l-dist".
#' @param n_perm Number of permutations to use for the test (default: 1000).
#' @param seed Random seed for reproducibility (default: NULL).
#' @param num_cores Number of cores for parallel computing (default: 1).
#'
#' @return An object of class "jhsic_test" containing:
#'   \item{statistic}{JHSIC test statistic value}
#'   \item{p.value}{Permutation-based p-value}
#'   \item{permutation_values}{Vector of JHSIC values from permutations}
#'   \item{stat_description}{Description of the statistic type used}
#'
#' @details
#' The Joint Hilbert-Schmidt Independence Criterion (JHSIC) extends the pairwise HSIC to detect complex
#' dependencies among multiple variables. When type = "euclidean", it becomes the joint distance covariance (JdCov)
#' proposed in Chakraborty & Zhang (2019), Unlike pairwise independence tests, JHSIC and JdCov can detect
#' higher-order dependencies that might be missed when only testing pairs of variables.
#'
#' The null hypothesis is that all random variables are mutually independent. The alternative hypothesis
#' is that there exists some form of dependence structure among the variables. The p-value is calculated
#' using a permutation test, where all variables except the first one are randomly permuted to create
#' the null distribution.
#'
#' Different statistic types offer various properties:
#'   \itemize{
#'     \item V-statistic: Standard, biased estimator
#'     \item U-statistic: Partially unbiased estimator that removes diagonal terms
#'     \item US-statistic: Scale-free U-statistic for normalized measure
#'     \item UR-statistic: Rank-based U-statistic, robust to monotonic transformations
#'   }
#'
#' @examples
#' # Example 1: Testing for independence among truly independent variables
#' set.seed(123)
#' x1 <- matrix(rnorm(100), ncol = 1)
#' x2 <- matrix(rnorm(100), ncol = 1)
#' x3 <- matrix(rnorm(100), ncol = 1)
#' test1 <- jhsic_test(list(x1, x2, x3), stat_type = "V", n_perm = 200)
#' print(test1)
#' plot(test1)
#'
#' # Example 2: Testing variables with pairwise independence but joint dependence
#' set.seed(456)
#' u <- runif(100, -pi, pi)
#' x1 <- matrix(sin(u), ncol = 1)
#' x2 <- matrix(cos(u), ncol = 1)
#' x3 <- matrix(sin(u) * cos(u), ncol = 1)
#' test2 <- jhsic_test(list(x1, x2, x3), stat_type = "V", n_perm = 200)
#' print(test2)
#'
#' # Example 3: Comparing different statistic types
#' test_v <- jhsic_test(list(x1, x2, x3), stat_type = "V", n_perm = 200)
#' test_u <- jhsic_test(list(x1, x2, x3), stat_type = "U", n_perm = 200)
#' test_us <- jhsic_test(list(x1, x2, x3), stat_type = "US", n_perm = 200)
#' test_ur <- jhsic_test(list(x1, x2, x3), stat_type = "UR", n_perm = 200)
#'
#' print(c(V = test_v$p.value, U = test_u$p.value,
#'         US = test_us$p.value, UR = test_ur$p.value))
#'
#' # Example 4: Using multivariate data
#' x_multi1 <- matrix(rnorm(200), ncol = 2)
#' x_multi2 <- matrix(rnorm(200), ncol = 2)
#' x_multi3 <- matrix(rnorm(200), ncol = 2)
#' test_multi <- jhsic_test(list(x_multi1, x_multi2, x_multi3), n_perm = 200)
#' print(test_multi)
#'
#' # Example 5: Testing with non-linear relationships
#' set.seed(789)
#' x <- matrix(runif(100, -3, 3), ncol = 1)
#' y <- matrix(x^2 + rnorm(100, sd = 0.5), ncol = 1)
#' z <- matrix(exp(x) + rnorm(100, sd = 0.5), ncol = 1)
#' test_nonlin <- jhsic_test(list(x, y, z), type = "gaussian", n_perm = 200)
#' print(test_nonlin)
#'
#' # Example 6: Using grouped variables with jhsic_test
#' # Create sample data with 3 datasets
#' set.seed(123)
#' n <- 100
#' x1 <- matrix(rnorm(n*3), ncol = 3)  # 3 variables in x1
#' x2 <- matrix(rnorm(n*4), ncol = 4)  # 4 variables in x2
#' x3 <- matrix(rnorm(n*2), ncol = 2)  # 2 variables in x3
#'
#' # Create complex dependency structure via a common factor
#' z <- runif(n, -pi, pi)
#' x1[, 1] <- 2*sin(z) + 0.1*rnorm(n)
#' x2[, 1] <- 2*cos(z) + 0.1*rnorm(n)
#' x3[, 1] <- sin(z)*cos(z) + 0.1*rnorm(n)
#'
#' # Define group structure for each dataset
#' groups1 <- c(1, 2, 2)        # First var in group 1, others in group 2
#' groups2 <- c(1, 2, 2, 2)     # First var in group 1, others in group 2
#' groups3 <- c(1, 2)           # Each var in its own group
#'
#' # Combine into a list matching the input list structure
#' group_list <- list(groups1, groups2, groups3)
#'
#' # Perform JHSIC test with the group structure
#' test_grouped <- jhsic_test(list(x1, x2, x3), type = "e-dist",
#'                           group = group_list, n_perm = 200)
#' print(test_grouped)
#'
#' # Compare with standard JHSIC test using V-statistic
#' test_standard <- jhsic_test(list(x1, x2, x3), type = "gaussian",
#'                            stat_type = "V", n_perm = 200)
#' print(test_standard)
#'
#' # Example 7: Using parallel computing for larger tests
#' \dontrun{
#' if(parallel::detectCores() > 1) {
#'   test_parallel <- jhsic_test(list(x_multi1, x_multi2, x_multi3),
#'                              n_perm = 1000, num_cores = 2)
#'   print(test_parallel)
#' }
#' }
#'
#' @references
#' Chakraborty, S., & Zhang, X. (2019). Distance Metrics for Measuring Joint Dependence with
#' Application to Causal Inference. \emph{Journal of the American Statistical Association}, 114, 1638-1650.
#'
#' @seealso
#' \code{\link{hsic_test}} for pairwise independence test
#' \code{\link{dhsic_test}} for an alternative joint independence test
#' \code{\link{mmd_test}} for two-sample testing
#'
#' @export
jhsic_test <- function(x, cc = 1, type = "gaussian", stat_type = "V", bw = NULL,
                       expo = 1, scale_factor = 0.5, group = NULL,
                       n_perm = 1000, seed = NULL, num_cores = 1) {
  # Load parallel package
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("The 'parallel' package is required for this function. Please install it.")
  }

  # Set random seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Validate inputs
  if (!is.list(x)) {
    stop("Input 'x' must be a list of matrices")
  }

  if (length(x) < 2) {
    stop("At least two matrices are required for jHSIC")
  }

  # Define valid statistic types
  valid_stats <- c("V", "U", "US", "UR")

  # Check if stat_type is valid
  if (!stat_type %in% valid_stats) {
    stop(sprintf("Invalid stat_type: '%s'. Must be one of: %s",
                 stat_type, paste(valid_stats, collapse = ", ")))
  }

  # Convert all vectors to matrices
  x_matrices <- lapply(x, function(xi) {
    if(is.null(dim(xi)) || length(dim(xi)) == 1) {
      matrix(xi, ncol = 1)
    } else {
      as.matrix(xi)
    }
  })

  # Get dimensions
  d <- length(x_matrices)
  n <- nrow(x_matrices[[1]])

  # Ensure all matrices have the same number of rows
  for (i in 2:d) {
    if (nrow(x_matrices[[i]]) != n) {
      stop("All matrices must have the same number of rows")
    }
  }

  # Calculate observed test statistic
  observed_jhsic <- jhsic(x_matrices, cc = cc, type = type, stat_type = stat_type,
                          bw = bw, expo = expo, scale_factor = scale_factor, group = group)

  # Define the function to run for each permutation
  perm_function <- function(i) {
    # Create a permuted version of the list
    x_perm <- x_matrices

    # Permute all matrices except the first one
    for (j in 2:d) {
      perm_indices <- sample(n, n, replace = FALSE)
      x_perm[[j]] <- x_matrices[[j]][perm_indices, , drop = FALSE]
    }

    # Calculate jHSIC for the permuted data
    return(jhsic(x_perm, cc = cc, type = type, stat_type = stat_type,
                bw = bw, expo = expo, scale_factor = scale_factor, group = group))
  }

  # Run permutations in parallel
  perm_results <- parallel::mclapply(1:n_perm, function(i) perm_function(i), mc.cores = num_cores)
  perm_results <- unlist(perm_results)

  # Calculate p-value (proportion of permuted statistics >= observed)
  p.value <- mean(perm_results >= observed_jhsic)

  # Get statistic description
  stat_description <- switch(stat_type,
                             "V" = "V-statistic",
                             "U" = "U-statistic",
                             "US" = "scale-free U-statistic",
                             "UR" = "rank-based U-statistic")

  # Create result list
  result <- list(
    statistic = observed_jhsic,
    p.value = p.value,
    permutation_values = perm_results,
    stat_description = stat_description
  )

  class(result) <- "jhsic_test"
  return(result)
}

#' High-Dimensional Two-Sample t-Test
#'
#' Performs a test for the difference between two high-dimensional distributions. This implementation
#' follows the methodology described in Chakraborty & Zhang (2021), which provides a framework
#' specifically designed for settings where the dimensionality is large compared to sample size.
#'
#' @param x First sample (matrix, data frame, or vector)
#' @param y Second sample (matrix, data frame, or vector)
#' @param type Type of kernel or distance to use (default: "e-dist").
#'   Options include "euclidean", "polynomial", "gaussian", "laplacian", "e-dist", "g-dist", or "l-dist".
#'   For high-dimensional data, "e-dist", "g-dist", and "l-dist" are specifically designed to provide
#'   better power.
#' @param bw Bandwidth parameter for kernel functions. If NULL, it will be automatically determined.
#' @param expo Exponent parameter for euclidean distance and polynomial kernel (default: 1).
#' @param scale_factor Scaling factor for automatic bandwidth calculation (default: 0.5).
#' @param group Optional vector specifying group membership for each column.
#'   Used for group-wise distance calculations in "e-dist", "g-dist", or "l-dist".
#'
#' @return A list with test statistic and p-value:
#'   \item{statistic}{The test statistic value}
#'   \item{p.value}{The p-value for the test}
#'
#' @details
#' The high-dimensional two-sample test is designed to compare two distributions when the number
#' of dimensions (p) is large relative to the sample size (n), potentially even p >> n. The test
#' is based on a modified Maximum Mean Discrepancy (MMD) statistic that has improved power in
#' high-dimensional settings.
#'
#' The test statistic incorporates bias correction terms and scaling factors to provide a
#' more stable and powerful test in high dimensions. The statistic follows a t-distribution
#' under the null hypothesis, enabling analytical p-value calculation without requiring
#' permutation procedures.
#'
#' When using "e-dist", "g-dist", or "l-dist" as the distance type, the function implements
#' specialized distance metrics designed for high-dimensional data:
#' \itemize{
#'   \item "e-dist": Euclidean-based aggregated distance
#'   \item "g-dist": Gaussian kernel-based aggregated distance
#'   \item "l-dist": Laplacian kernel-based aggregated distance
#' }
#'
#' These metrics address the concentration of distances phenomenon in high dimensions by
#' aggregating componentwise distances in a way that maintains discriminative power.
#'
#' @examples
#' # Example 1: High-dimensional data with true difference
#' set.seed(123)
#' p <- 500  # High dimensionality
#' n1 <- 30  # Small sample size
#' n2 <- 25  # Small sample size
#'
#' # Generate high-dimensional data with difference in first p/4 dimensions
#' x <- matrix(rnorm(n1 * p), nrow = n1, ncol = p)
#' y <- matrix(rnorm(n2 * p), nrow = n2, ncol = p)
#' y[, 1:(p/4)] <- y[, 1:(p/4)] + 0.7  # Add shift in first p/4 dimensions
#'
#' # Test using e-dist (recommended for high-dimensional data)
#' test1 <- hd_two_test(x, y, type = "e-dist")
#' print(test1)
#'
#' # Example 2: High-dimensional data with no difference (null hypothesis)
#' x2 <- matrix(rnorm(n1 * p), nrow = n1, ncol = p)
#' y2 <- matrix(rnorm(n2 * p), nrow = n2, ncol = p)
#'
#' test2 <- hd_two_test(x2, y2, type = "e-dist")
#' print(test2)
#'
#' # Example 3: Comparing different distance types for high-dimensional data
#' test_edist <- hd_two_test(x, y, type = "e-dist")
#' test_gdist <- hd_two_test(x, y, type = "g-dist")
#' test_ldist <- hd_two_test(x, y, type = "l-dist")
#'
#' cat("E-dist p-value:", test_edist$p.value, "\n")
#' cat("G-dist p-value:", test_gdist$p.value, "\n")
#' cat("L-dist p-value:", test_ldist$p.value, "\n")
#'
#' # Example 4: Using group structure for very high-dimensional data
#' p <- 1000
#' x3 <- matrix(rnorm(n1 * p), nrow = n1, ncol = p)
#' y3 <- matrix(rnorm(n2 * p), nrow = n2, ncol = p)
#' y3[, 1:(p/4)] <- y3[, 1:(p/4)] + 0.5  # Add shift in first p/4 dimensions
#'
#' # Create groups (e.g., 10 groups of 100 variables each)
#' groups <- rep(1:10, each = 100)
#'
#' # Test with group structure
#' test_grouped <- hd_two_test(x3, y3, type = "e-dist", group = groups)
#' print(test_grouped)
#'
#' @references
#' Chakraborty, S., & Zhang, X. (2021). A New Framework for Distance and Kernel-based Metrics in
#' High-dimension. \emph{Electronic Journal of Statistics}, 15, 5455-5522.
#'
#' @seealso
#' \code{\link{mmd_test}} for the standard MMD two-sample test
#' \code{\link{ed_test}} for the energy distance test
#' \code{\link{hd_dep_test}} for high-dimensional dependence testing
#'
#' @export
hd_two_test <- function(x, y, type = "e-dist", bw = NULL, expo = 1, scale_factor = 0.5, group = NULL) {
  # Convert vectors to single-column matrices if needed
  if (is.vector(x) && !is.matrix(x)) {
    x <- matrix(x, ncol = 1)
  }
  if (is.vector(y) && !is.matrix(y)) {
    y <- matrix(y, ncol = 1)
  }

  # Get dimensions
  n <- nrow(x)
  m <- nrow(y)

  # Calculate distance matrices
  D <- KDist_matrix(rbind(x, y), type = type, bw = bw, expo = expo, scale_factor = scale_factor, group = group)
  Dxx <- D[1:n, 1:n]
  Dyy <- D[(n+1):(n+m), (n+1):(n+m)]
  Dxy <- D[1:n, (n+1):(n+m)]

  # Calculate test statistics
  mmd_val <- suppressMessages(mmd(x = D, type = type, u_center = TRUE, n = n, m = m))
  hsic_val_x <- hsic(x = Dxx, y = Dxx, type = type, u_center = TRUE, is_distance = TRUE)
  hsic_val_y <- hsic(x = Dyy, y = Dyy, type = type, u_center = TRUE, is_distance = TRUE)
  chsic_val <- chsic(x = Dxy)

  # Calculate denominator scalars
  scalar_den <- 1/(n*m) + 1/(2*n*(n-1)) + 1/(2*m*(m-1))

  # Calculate the Snm statistic components
  n_factor <- n * (n - 3) / 2
  m_factor <- m * (m - 3) / 2

  # Numerator of Snm
  num_Snm <- 4 * (n - 1) * (m - 1) * chsic_val +
    4 * n_factor * hsic_val_x +
    4 * m_factor * hsic_val_y

  # Denominator of Snm
  den_Snm <- (n - 1) * (m - 1) + n_factor + m_factor

  # Compute Snm
  Snm <- num_Snm / den_Snm

  # Compute denominator and final statistic
  Den <- sqrt(scalar_den * Snm)

  # There seems to be a variable 'ed_all' used but not defined earlier
  # Assuming it should be mmd_val based on context
  stat <- mmd_val / Den

  # Calculate p-value
  pval <- 1 - pt(stat, df = den_Snm)

  return(list(statistic = stat, p.value = pval))
}

#' High-Dimensional Dependence Test
#'
#' Performs a test for dependence between two high-dimensional datasets. This implementation
#' follows the methodology described in Chakraborty & Zhang (2021), which provides a framework
#' specifically designed for settings where the dimensionality is large compared to sample size.
#'
#' @param x First dataset (matrix, data frame, or vector)
#' @param y Second dataset (matrix, data frame, or vector)
#' @param type Type of kernel or distance to use (default: "e-dist").
#'   Options include "euclidean", "polynomial", "gaussian", "laplacian", "e-dist", "g-dist", or "l-dist".
#'   For high-dimensional data, "e-dist", "g-dist", and "l-dist" are specifically designed to provide
#'   better power.
#' @param bw Bandwidth parameter for kernel functions. If NULL, it will be automatically determined.
#' @param expo Exponent parameter for euclidean distance and polynomial kernel (default: 1).
#' @param scale_factor Scaling factor for automatic bandwidth calculation (default: 0.5).
#' @param group Optional list specifying group membership for each column in the corresponding input list elements.
#'   The length of the group list should match the length of the input list x, where each element
#'   contains the grouping information for the respective element in x. Used for group-wise
#'   distance calculations in "e-dist", "g-dist", or "l-dist".
#' @param is_distance Logical; whether input matrices are already distance matrices (default: FALSE).
#'
#' @return A list containing:
#'   \item{statistic}{The test statistic value}
#'   \item{p.value}{The p-value for the test}
#'   \item{correlation}{The HSIC correlation value}
#'   \item{df}{Degrees of freedom for the test}
#'
#' @details
#' The high-dimensional dependence test is designed to detect dependence between two sets of variables
#' when the dimensionality is large relative to the sample size. The test is based on a modified
#' Hilbert-Schmidt Independence Criterion (HSIC) that has improved power in high-dimensional settings.
#'
#' The test uses a t-statistic based on the HSIC correlation, which follows a t-distribution under
#' the null hypothesis of independence. This allows analytical calculation of p-values without
#' requiring computationally intensive permutation procedures, making it particularly suitable for
#' high-dimensional applications.
#'
#' When using "e-dist", "g-dist", or "l-dist" as the distance type, the function implements
#' specialized distance metrics designed for high-dimensional data:
#' \itemize{
#'   \item "e-dist": Euclidean-based aggregated distance
#'   \item "g-dist": Gaussian kernel-based aggregated distance
#'   \item "l-dist": Laplacian kernel-based aggregated distance
#' }
#'
#' These metrics address the concentration of distances phenomenon in high dimensions by
#' aggregating componentwise distances in a way that maintains discriminative power for
#' dependence testing.
#'
#' @examples
#' # Example 1: High-dimensional data with true dependence
#' set.seed(123)
#' n <- 50      # Sample size
#' p_x <- 200   # Dimensions of first dataset
#' p_y <- 150   # Dimensions of second dataset
#'
#' # Generate first high-dimensional dataset
#' x <- matrix(rnorm(n * p_x), nrow = n, ncol = p_x)
#'
#' # Generate dependent second dataset
#' # Only the first 10 dimensions are related
#' y <- matrix(rnorm(n * p_y), nrow = n, ncol = p_y)
#' y[, 1:10] <- x[, 1:10] + 0.5 * matrix(rnorm(n * 10), nrow = n, ncol = 10)
#'
#' # Test using e-dist (recommended for high-dimensional data)
#' test1 <- hd_dep_test(x, y, type = "e-dist")
#' print(test1)
#'
#' # Example 2: High-dimensional data with independence (null hypothesis)
#' x2 <- matrix(rnorm(n * p_x), nrow = n, ncol = p_x)
#' y2 <- matrix(rnorm(n * p_y), nrow = n, ncol = p_y)
#'
#' test2 <- hd_dep_test(x2, y2, type = "e-dist")
#' print(test2)
#'
#' # Example 3: Comparing different distance types
#' test_edist <- hd_dep_test(x, y, type = "e-dist")
#' test_gdist <- hd_dep_test(x, y, type = "g-dist")
#' test_ldist <- hd_dep_test(x, y, type = "l-dist")
#'
#' cat("E-dist p-value:", test_edist$p.value,
#'     "correlation:", test_edist$correlation, "\n")
#' cat("G-dist p-value:", test_gdist$p.value,
#'     "correlation:", test_gdist$correlation, "\n")
#' cat("L-dist p-value:", test_ldist$p.value,
#'     "correlation:", test_ldist$correlation, "\n")
#'
#' # Example 4: Very high-dimensional scenario
#' n <- 30       # Smaller sample size
#' p_x <- 500    # Very high dimensions
#' p_y <- 400    # Very high dimensions
#'
#' # Create datasets with sparse dependence
#' x3 <- matrix(rnorm(n * p_x), nrow = n, ncol = p_x)
#' y3 <- matrix(rnorm(n * p_y), nrow = n, ncol = p_y)
#'
#' # 100 dimensions have dependencies
#' shared_dims <- sample(1:min(p_x, p_y), 100)
#' for (i in shared_dims) {
#'   y3[, i] <- x3[, i] + 0.3 * rnorm(n)
#' }
#'
#' test3 <- hd_dep_test(x3, y3, type = "e-dist")
#' print(test3)
#'
#' # Example 5: Using grouped variables with hd_dep_test
#' # Generate high-dimensional datasets with clustered variables
#' set.seed(123)
#' n <- 50      # Sample size
#' p1 <- 100    # Dimensions of x
#' p2 <- 80     # Dimensions of y
#'
#' # Create data matrices
#' x <- matrix(rnorm(n * p1), nrow = n, ncol = p1)
#' y <- matrix(rnorm(n * p2), nrow = n, ncol = p2)
#'
#' # Create dependencies between specific groups
#' # First 10 variables in x affect first 8 variables in y
#' y[, 1:8] <- y[, 1:8] + x[, 1:8]
#'
#' # Define group structure for x and y
#' # Group x variables into 5 groups of 20 each
#' x_groups <- rep(1:5, each = 20)
#' # Group y variables into 4 groups of 20 each
#' y_groups <- rep(1:4, each = 20)
#'
#' # Combine into a list
#' group_list <- list(x_groups, y_groups)
#'
#' # Perform high-dimensional dependence test with grouped variables
#' test_grouped <- hd_dep_test(x, y, type = "e-dist", group = group_list)
#' print(test_grouped)
#'
#' # Compare with standard test (without grouping)
#' test_standard <- hd_dep_test(x, y, type = "e-dist")
#' print(test_standard)
#'
#' # Example 6: Using pre-computed distance matrices
#' library(stats)
#' Dx <- KDist_matrix(x, type = "e-dist")
#' Dy <- KDist_matrix(y, type = "e-dist")
#'
#' test_dist <- hd_dep_test(Dx, Dy, is_distance = TRUE)
#' print(test_dist)
#'
#' @references
#' Chakraborty, S., & Zhang, X. (2021). A New Framework for Distance and Kernel-based Metrics in
#' High-dimension. \emph{Electronic Journal of Statistics}, 15, 5455-5522.
#'
#' @seealso
#' \code{\link{hsic_test}} for the standard HSIC independence test
#' \code{\link{dcov_test}} for the distance covariance test
#' \code{\link{hd_two_test}} for high-dimensional two-sample testing
#'
#' @export
hd_dep_test <- function(x, y, type = "e-dist", bw = NULL,
                        expo = 1, scale_factor = 0.5, group = NULL, is_distance = FALSE) {
  # Convert vectors to single-column matrices if needed
  if (is.vector(x) && !is.matrix(x)) {
    x <- matrix(x, ncol = 1)
  }
  if (is.vector(y) && !is.matrix(y)) {
    y <- matrix(y, ncol = 1)
  }

  # Check if dimensions are compatible
  n_x <- nrow(x)
  n_y <- nrow(y)

  if (n_x != n_y) {
    stop("The number of observations in x and y must be the same.")
  }

  # Calculate degrees of freedom
  df <- n_x * (n_x - 3) / 2

  # Calculate HSIC correlation with U-centering
  hsic_cor_val <- hsic_cor(x, y, type, bw, expo, scale_factor, group,
                           u_center = TRUE, is_distance)

  # Calculate test statistic
  t_statistic <- sqrt((df - 1) / (1 - hsic_cor_val^2)) * hsic_cor_val

  # Calculate p-value
  p_value <- 1 - pt(t_statistic, df = df - 1)

  # Return results
  return(list(statistic = t_statistic, p.value = p_value,
              correlation = hsic_cor_val, df = df - 1))
}
