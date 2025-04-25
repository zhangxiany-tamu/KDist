#' Joint Hilbert-Schmidt Independence Criterion
#'
#' Calculates the joint HSIC (JHSIC) for measuring mutual independence among multiple random variables.
#' JHSIC extends the pairwise HSIC to detect complex dependencies among three or more variables
#' that might be missed by testing only pairs of variables. Unlike pairwise tests that
#' can only detect dependence between pairs, JHSIC provides a single test statistic for
#' mutual independence.
#'
#' @param x List of matrices or data frames, where each element corresponds to a random variable
#' @param cc Constant parameter (default: 1)
#' @param type Type of kernel or distance to use (default: "gaussian").
#'   Options include "euclidean", "polynomial", "gaussian", "laplacian", "e-dist", "g-dist", or "l-dist".
#'   Note that when type = "euclidean", the function calculates the joint distance covariance (JdCov)
#'   as described in Chakraborty & Zhang (2019).
#' @param stat_type Type of statistic (default: "V")
#'   \itemize{
#'     \item "V": V-statistic (standard, biased estimator)
#'     \item "U": U-statistic (partially unbiased estimator that removes diagonal terms)
#'     \item "US": Scale-free U-statistic for normalized measure
#'     \item "UR": Rank-based U-statistic, robust to monotonic transformations
#'   }
#' @param bw Bandwidth parameter for kernel functions. If NULL, it will be automatically determined.
#' @param expo Exponent parameter for euclidean distance and polynomial kernel (default: 1).
#' @param scale_factor Scaling factor for automatic bandwidth calculation (default: 0.5).
#' @param group Optional list specifying group membership for each column in the corresponding input list elements.
#'   Used for group-wise distance calculations in "e-dist", "g-dist", or "l-dist".
#'
#' @return JHSIC value
#'
#' @details
#' The Joint Hilbert-Schmidt Independence Criterion (JHSIC) extends the pairwise HSIC to detect
#' complex dependencies among multiple variables. It is particularly useful for cases where
#' variables might have pairwise independence but joint dependence.
#'
#' When type = "euclidean", JHSIC becomes the joint distance covariance (JdCov) proposed in
#' Chakraborty & Zhang (2019). JdCov is a measure of joint dependence based on distances
#' that can detect complex dependency structures among multiple variables. When type = "polynomial",
#' "gaussian", or "laplacian", the kernel matrix will be transformed into a kernel-induced distance
#' matrix.
#'
#' Different statistic types offer various properties:
#'   \itemize{
#'     \item V-statistic: Standard, biased estimator
#'     \item U-statistic: Partially unbiased estimator that removes diagonal terms
#'     \item US-statistic: Scale-free U-statistic for normalized measure
#'     \item UR-statistic: Rank-based U-statistic, robust to monotonic transformations
#'   }
#'
#' When using "e-dist", "g-dist", or "l-dist" as the distance type, the function implements
#' specialized distance metrics designed for high-dimensional data:
#' \itemize{
#'   \item "e-dist": Euclidean-based aggregated distance
#'   \item "g-dist": Gaussian kernel-based aggregated distance
#'   \item "l-dist": Laplacian kernel-based aggregated distance
#' }
#'
#' @examples
#' # Example 1: Testing for independence among truly independent variables
#' set.seed(123)
#' x1 <- matrix(rnorm(100), ncol = 1)
#' x2 <- matrix(rnorm(100), ncol = 1)
#' x3 <- matrix(rnorm(100), ncol = 1)
#' jhsic_indep <- jhsic(list(x1, x2, x3), stat_type = "V")
#' print(jhsic_indep)
#'
#' # Example 2: Testing variables with pairwise independence but joint dependence
#' set.seed(456)
#' u <- runif(100, -pi, pi)
#' x1 <- matrix(sin(u), ncol = 1)
#' x2 <- matrix(cos(u), ncol = 1)
#' x3 <- matrix(sin(u) * cos(u), ncol = 1)
#' jhsic_dep <- jhsic(list(x1, x2, x3), stat_type = "V")
#' print(jhsic_dep)
#'
#' # Example 3: Comparing different statistic types
#' jhsic_v <- jhsic(list(x1, x2, x3), stat_type = "V")
#' jhsic_u <- jhsic(list(x1, x2, x3), stat_type = "U")
#' jhsic_us <- jhsic(list(x1, x2, x3), stat_type = "US")
#' jhsic_ur <- jhsic(list(x1, x2, x3), stat_type = "UR")
#'
#' print(c(V = jhsic_v, U = jhsic_u, US = jhsic_us, UR = jhsic_ur))
#'
#' # Example 4: Using different kernel types
#' jhsic_gaussian <- jhsic(list(x1, x2, x3), type = "gaussian", stat_type = "V")
#' jhsic_laplacian <- jhsic(list(x1, x2, x3), type = "laplacian", stat_type = "V")
#' # Using type = "euclidean" calculates the joint distance covariance (JdCov)
#' # as described in Chakraborty & Zhang (2019)
#' jhsic_euclidean <- jhsic(list(x1, x2, x3), type = "euclidean", stat_type = "V")
#'
#' print(c(gaussian = jhsic_gaussian,
#'         laplacian = jhsic_laplacian,
#'         euclidean = jhsic_euclidean))
#'
#' # Example 5: Using multivariate data
#' x_multi1 <- matrix(rnorm(200), ncol = 2)
#' x_multi2 <- matrix(rnorm(200), ncol = 2)
#' x_multi3 <- matrix(rnorm(200), ncol = 2)
#' jhsic_multi <- jhsic(list(x_multi1, x_multi2, x_multi3))
#' print(jhsic_multi)
#'
#' # Example 6: Testing with non-linear relationships
#' set.seed(789)
#' x <- matrix(runif(100, -3, 3), ncol = 1)
#' y <- matrix(x^2 + rnorm(100, sd = 0.5), ncol = 1)
#' z <- matrix(exp(x) + rnorm(100, sd = 0.5), ncol = 1)
#' jhsic_nonlin <- jhsic(list(x, y, z), type = "gaussian")
#' print(jhsic_nonlin)
#'
#' # Example 7: Using grouped variables with jhsic
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
#' # Calculate JHSIC with the group structure
#' jhsic_grouped <- jhsic(list(x1, x2, x3), type = "e-dist", group = group_list)
#' print(jhsic_grouped)
#'
#' @references
#' Chakraborty, S., & Zhang, X. (2019). Distance Metrics for Measuring Joint Dependence with
#' Application to Causal Inference. Journal of the American Statistical Association, 114, 1638-1650.
#'
#' @seealso
#' \code{\link{dhsic}} for an alternative joint independence test
#' \code{\link{hsic}} for pairwise independence test
#' \code{\link{jhsic_test}} for permutation-based testing using JHSIC
#'
#' @export
jhsic <- function(x, cc = 1, type = "gaussian", stat_type = "V", bw = NULL, expo = 1,
                  scale_factor = 0.5, group = NULL) {
  # Define valid statistic types
  valid_stats <- c("V", "U", "US", "UR")

  # Check if stat_type is valid
  if (!stat_type %in% valid_stats) {
    stop(sprintf("Invalid stat_type: '%s'. Must be one of: %s",
                 stat_type, paste(valid_stats, collapse = ", ")))
  }

  # Input validation
  if (!is.list(x)) {
    stop("Input 'x' must be a list of matrices")
  }

  if (length(x) < 2) {
    stop("At least two matrices are required in the list 'x'")
  }

  # Convert all elements to matrices if needed
  for (i in seq_along(x)) {
    if (is.vector(x[[i]])) {
      x[[i]] <- matrix(x[[i]], ncol = 1)
    } else if (is.data.frame(x[[i]])) {
      x[[i]] <- as.matrix(x[[i]])
    } else if (!is.matrix(x[[i]])) {
      stop("Each element in 'x' must be a vector, matrix or data frame")
    }
  }

  # Check that all matrices have the same number of rows
  n_rows <- nrow(x[[1]])
  for (i in 2:length(x)) {
    if (nrow(x[[i]]) != n_rows) {
      stop("All matrices in 'x' must have the same number of rows")
    }
  }

  # Use switch for cleaner code structure
  result <- switch(stat_type,
                   "V" = jhsic_v(x, cc = cc, type = type, bw_sexp = bw, expo = expo,
                                 scale_factor = scale_factor, group_ = group),
                   "U" = jhsic_u(x, cc = cc, type = type, bw_sexp = bw, expo = expo,
                                 scale_factor = scale_factor, group_ = group),
                   "US" = jhsic_u(x, cc = cc, type = type, bw_sexp = bw, expo = expo,
                                  scale = TRUE, scale_factor = scale_factor, group_ = group),
                   "UR" = jhsic_u(rank_list(x), cc = cc, type = type, bw_sexp = bw, expo = expo,
                                  scale_factor = scale_factor, group_ = group),
                   stop("Unhandled stat_type (should not happen due to earlier check)")
  )

  return(result)
}

#' d-variable Hilbert-Schmidt Independence Criterion (dHSIC)
#'
#' Calculates the dHSIC for measuring mutual independence among multiple random variables.
#' This implementation follows the methodology described in Pfister et al. (2018).
#' While pairwise independence tests can detect dependence between pairs of variables,
#' they may miss higher-order dependencies among three or more variables. The dHSIC
#' addresses this limitation by providing a single test statistic for mutual independence.
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
#'
#' @return dHSIC value
#'
#' @details
#' The d-variable Hilbert-Schmidt Independence Criterion (dHSIC) extends the pairwise HSIC to
#' measure mutual independence among d random variables. For d random variables
#' \eqn{X^{(1)}, X^{(2)}, \ldots, X^{(d)}}, the dHSIC is defined as:
#'
#' \deqn{dHSIC(X^{(1)}, \ldots, X^{(d)}) = ||\otimes_{j=1}^{d} \mu_{X^{(j)}} - \mu_{X^{(1)},\ldots,X^{(d)}}||^2_{\mathcal{H}}}
#'
#' where \eqn{\mu_{X^{(j)}}} is the mean embedding of the distribution of \eqn{X^{(j)}} in a
#' reproducing kernel Hilbert space (RKHS), \eqn{\mu_{X^{(1)},\ldots,X^{(d)}}} is the joint
#' embedding, and \eqn{\otimes} denotes the tensor product.
#'
#' dHSIC can detect higher-order dependencies among multiple variables that might be missed
#' when only testing pairs of variables.
#'
#' @examples
#' # Example 1: Variables with pairwise independence but mutual dependence
#' set.seed(456)
#' u <- runif(100, -pi, pi)
#' x1 <- matrix(sin(u), ncol = 1)
#' x2 <- matrix(cos(u), ncol = 1)
#' x3 <- matrix(sin(u) * cos(u), ncol = 1)
#' dhsic_val <- dhsic(list(x1, x2, x3), type = "gaussian")
#' print(dhsic_val)
#'
#' # Example 2: Using different kernel types
#' dhsic_gaussian <- dhsic(list(x1, x2, x3), type = "gaussian")
#' dhsic_laplacian <- dhsic(list(x1, x2, x3), type = "laplacian")
#' dhsic_euclidean <- dhsic(list(x1, x2, x3), type = "euclidean")
#' print(c(gaussian = dhsic_gaussian,
#'         laplacian = dhsic_laplacian,
#'         euclidean = dhsic_euclidean))
#'
#' # Example 3: Using multivariate data
#' x_multi1 <- matrix(rnorm(200), ncol = 2)  # 2D data
#' x_multi2 <- matrix(rnorm(200), ncol = 2)  # 2D data
#' # Create dependency: second dataset depends on first
#' x_multi2 <- x_multi2 + 0.7 * x_multi1
#' x_multi3 <- matrix(rnorm(200), ncol = 2)  # Independent 2D data
#' dhsic_multi <- dhsic(list(x_multi1, x_multi2, x_multi3))
#' print(dhsic_multi)
#'
#' @references
#' Pfister, N., Bühlmann, P., Schölkopf, B., & Peters, J. (2018).
#' Kernel-based tests for joint independence. Journal of the Royal Statistical Society:
#' Series B (Statistical Methodology), 80(1), 5-31.
#'
#' @seealso
#' \code{\link{jhsic}} for an alternative joint independence measure
#' \code{\link{hsic}} for pairwise independence test
#' \code{\link{dhsic_test}} for permutation-based testing using dHSIC
#'
#' @export
dhsic <- function(x, type = "gaussian", bw = NULL, expo = 1, scale_factor = 0.5, group = NULL) {
  # Validate inputs
  if (!is.list(x)) stop("'x' must be a list of matrices or vectors")

  d <- length(x)

  # Validate bandwidth parameter
  if (!is.null(bw) && length(bw) > 1 && length(bw) != d) {
    stop("'bw' must be NULL, a scalar, or a vector of length equal to the number of matrices")
  }

  # Convert inputs to matrices
  x_matrices <- lapply(x, function(xi) {
    if(is.null(dim(xi)) || length(dim(xi)) == 1) {
      matrix(xi, ncol = 1)
    } else {
      as.matrix(xi)
    }
  })

  # Check that all matrices have the same number of rows
  n_rows <- nrow(x_matrices[[1]])
  for (i in 2:d) {
    if (nrow(x_matrices[[i]]) != n_rows) {
      stop("All matrices in 'x' must have the same number of rows")
    }
  }

  # Call the C++ implementation
  dhsic_fast_rcpp(x_matrices, type, bw, expo, scale_factor, group)
}

#' A Parallel Implementation of HSIC and dCov for Mutual Independence Testing in High-dimension
#'
#' @param x Matrix with columns to test for mutual independence
#' @param type Type of kernel or distance
#' @param bw Bandwidth parameter
#' @param expo Exponent parameter
#' @param scale_factor Scaling factor for bandwidth
#' @param cores Number of cores for parallel processing
#'
#' @return A list containing test statistic and p-value
#'
#' @importFrom parallel mclapply detectCores
#' @importFrom stats pnorm
#' @keywords internal
.mhsic_parallel <- function(x, type = "gaussian", bw = NULL, expo = 1, scale_factor = 0.5, cores = NULL) {
  n <- dim(x)[1]
  p <- dim(x)[2]

  # Set number of cores
  if(is.null(cores)) cores <- detectCores() - 1

  # Step 1: Compute diagonal elements in parallel
  hsic_diag <- mclapply(1:p, function(i) {
    hsic(x[,i], x[,i], type = type, bw = bw, expo = expo, scale_factor = scale_factor, u_center = TRUE)
  }, mc.cores = cores)
  hsic_diag <- unlist(hsic_diag)

  # Step 2: Generate all pairs to compute
  pairs <- expand.grid(i = 2:p, j = 1:(p-1))
  pairs <- subset(pairs, j < i)

  # Step 3: Compute HSIC values for all pairs in parallel
  hsic_results <- mclapply(1:nrow(pairs), function(idx) {
    i <- pairs$i[idx]
    j <- pairs$j[idx]
    hsic_val <- hsic(x[,i], x[,j], type = type, bw = bw, expo = expo, scale_factor = scale_factor, u_center = TRUE)
    return(hsic_val)
  }, mc.cores = cores)

  # Combine results
  num <- sum(unlist(hsic_results))

  # Compute denominator
  den <- 0
  for (i in 2:p) {
    for (j in 1:(i-1)) {
      den <- den + hsic_diag[i] * hsic_diag[j]
    }
  }

  stat <- sqrt(n * (n - 1) / 2 / den) * num
  pval <- 1 - pnorm(stat)
  return(list(statistic = stat, p.value = pval))
}

#' HSIC and dCov for Mutual Independence Testing in High-dimension
#'
#' Performs a test for mutual independence among multiple variables in a high-dimensional setting.
#' This implementation is based on the methodology described in Yao, Zhang, & Shao (2018),
#' which was originally developed for distance covariance. The original method corresponds to
#' using this function with type = "euclidean", while other kernel options extend the approach
#' to different dependency measures.
#'
#' @param x Matrix where each column is a variable to test for mutual independence
#' @param type Type of kernel or distance to use (default: "gaussian").
#'   Options include "euclidean", "polynomial", "gaussian", "laplacian", "e-dist", "g-dist", or "l-dist".
#'   Note that type = "euclidean" corresponds to the original distance covariance method
#'   proposed in Yao, Zhang, & Shao (2018).
#' @param bw Bandwidth parameter for kernel functions. If NULL, it will be automatically determined.
#' @param expo Exponent parameter for Euclidean distance and polynomial kernel (default: 1).
#' @param scale_factor Scaling factor for automatic bandwidth calculation (default: 0.5).
#' @param num_cores Number of cores for parallel processing (default: NULL, which uses
#'   one less than the available number of cores).
#'
#' @return A list with test statistic and p-value:
#'   \item{statistic}{The test statistic value}
#'   \item{p.value}{The p-value for the test}
#'
#' @details
#' The mutual independence test is designed to detect dependencies among all
#' variables simultaneously in high-dimensional settings. The method is particularly effective
#' when the number of variables (p) is large, potentially even larger than the sample size (n).
#'
#' This implementation is based on the approach of Yao, Zhang, & Shao (2018), which was
#' specifically developed for distance covariance (corresponding to type = "euclidean").
#' The method uses pairwise dependency calculations to create a test statistic that follows
#' a standard normal distribution under the null hypothesis of mutual independence.
#' This approach avoids computationally intensive permutation procedures, making it
#' suitable for high-dimensional data.
#'
#' When type is set to other values, the function extends the original methodology to use
#' different kernel or distance measures, which may be more appropriate for specific types
#' of dependencies.
#'
#' When the number of variables is large (p > 60), the function automatically uses
#' parallel computing to speed up the calculation.
#'
#' @examples
#' # Example 1: High-dimensional data with mutual independence
#' set.seed(123)
#' n <- 100  # Sample size
#' p <- 200  # High dimensionality
#' # Generate independent high-dimensional data
#' x_indep <- matrix(rnorm(n * p), ncol = p)
#' # Using the original distance covariance method (Yao et al., 2018)
#' test_indep <- mhsic(x_indep, type = "euclidean")
#' print(test_indep)
#'
#' # Example 2: High-dimensional data with dependence
#' # Create dependency structure: first 20 variables are dependent
#' x_dep <- matrix(rnorm(n * p), ncol = p)
#' # Add a common factor to first 20 variables
#' common_factor <- rnorm(n)
#' for(i in 1:20) {
#'   x_dep[, i] <- x_dep[, i] + 0.7 * common_factor
#' }
#' test_dep <- mhsic(x_dep, type = "euclidean")
#' print(test_dep)
#'
#' # Example 3: Using different kernels with high-dimensional data
#' p_moderate <- 300  # Moderately high dimensionality
#' x_mod <- matrix(rnorm(n * p_moderate), ncol = p_moderate)
#' # Make first 20 variables dependent
#' common_factor <- rnorm(n)
#' for(i in 1:20) {
#'   x_mod[, i] <- x_mod[, i] + 0.8 * common_factor
#' }
#' # Original distance covariance method (Yao et al., 2018)
#' test_dCov <- mhsic(x_mod, type = "euclidean", num_cores = 2)
#' # Extensions using other kernels
#' test_gaussian <- mhsic(x_mod, type = "gaussian", num_cores = 2)
#' test_laplacian <- mhsic(x_mod, type = "laplacian", num_cores = 2)
#' cat("Distance covariance p-value:", test_dCov$p.value, "\n")
#' cat("Gaussian kernel p-value:", test_gaussian$p.value, "\n")
#' cat("Laplacian kernel p-value:", test_laplacian$p.value, "\n")
#'
#' @references
#' Yao, S., Zhang, X., & Shao, X. (2018). Testing Mutual Independence in High Dimension
#' via Distance Covariance. Journal of the Royal Statistical Society: Series B, 80, 455-480.
#'
#' @seealso
#' \code{\link{hsic}} for pairwise independence test
#' \code{\link{jhsic}} for joint independence test
#' \code{\link{dhsic}} for d-variable HSIC test
#'
#' @export
mhsic <- function(x, type = "gaussian", bw = NULL, expo = 1, scale_factor = 0.5, num_cores = 1) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame")
  }

  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  p <- dim(x)[2]

  # Use parallel processing only for larger matrices
  if (p > 60) {
    return(.mhsic_parallel(x, type, bw, expo, scale_factor, num_cores))
  } else {
    return(mhsic_cpp(x, type, bw, expo, scale_factor))
  }
}
