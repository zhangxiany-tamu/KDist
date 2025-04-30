#' Calculate Kernel or Distance Matrix
#'
#' This function calculates various kernel or distance matrices for input data.
#'
#' @param data Matrix or data frame. Each row represents an observation.
#' @param type Character string specifying the type of distance/kernel.
#'   Options include "euclidean", "polynomial", "gaussian", "laplacian",
#'   "e-dist", "g-dist", or "l-dist".
#' @param bw Bandwidth parameter for gaussian and laplacian kernels.
#'   If NULL, it will be automatically determined.
#' @param expo Exponent parameter for euclidean distance and polynomial kernel.
#' @param scale_factor Scaling factor for automatic bandwidth calculation.
#' @param group Optional vector specifying group membership for each column.
#'   Used for group-wise distance calculations in "e-dist", "g-dist", or "l-dist".
#'
#' @return A distance or kernel matrix.
#'
#' @examples
#' # Create sample data
#' set.seed(123)
#' x <- matrix(rnorm(50), ncol = 5)
#'
#' # Calculate Euclidean distance matrix
#' D_euclidean <- KDist_matrix(x, type = "euclidean")
#'
#' # Calculate Gaussian kernel matrix with automatic bandwidth
#' D_gaussian <- KDist_matrix(x, type = "gaussian")
#'
#' # Calculate Laplacian kernel with custom bandwidth
#' D_laplacian <- KDist_matrix(x, type = "laplacian", bw = 0.5)
#'
#' # Calculate polynomial kernel with degree = 2
#' D_poly <- KDist_matrix(x, type = "polynomial", expo = 2)
#'
#' # Using group structure - columns 1-2 in group 1, columns 3-5 in group 2
#' groups <- c(1, 1, 2, 2, 2)
#' D_grouped <- KDist_matrix(x, type = "e-dist", group = groups)
#'
#' # Working with a data frame
#' df <- data.frame(a = rnorm(10), b = rnorm(10), c = rnorm(10))
#' D_df <- KDist_matrix(df, type = "gaussian", scale_factor = 0.7)
#'
#' # Working with a vector
#' vec <- rnorm(15)
#' D_vec <- KDist_matrix(vec)
#' @export
KDist_matrix <- function(data, type = "euclidean", bw = NULL, expo = 1, scale_factor = 0.5, group = NULL) {
  # Ensure data is a matrix
  if (is.vector(data)) {
    data <- as.matrix(data)
    p <- 1
  } else if (is.data.frame(data)) {
    data <- as.matrix(data)
    p <- ncol(data)
  } else if (is.matrix(data)) {
    p <- ncol(data)
  } else {
    stop("Input 'data' must be a vector, matrix or data frame.")
  }

  # Check the format of group
  if (!is.null(group)) {
    if (length(group) != p) {
      stop("Length of group vector must match number of columns in data")
    }
  }

  if (type == "euclidean") {
    # Use Rcpp function for Euclidean distance
    return(euclidean_dist_rcpp(data, expo))
  }

  if (type == "polynomial") {
    # Use Rcpp function for polynomial kernel
    return(polynomial_kernel_rcpp(data, expo))
  }

  if (type %in% c("gaussian", "laplacian")) {
    # Determine bandwidth parameter if not provided
    if (is.null(bw)) {
      bw <- bw_optim(data, scale_factor = scale_factor)
    } else if (length(bw) > 1) {
      bw <- bw[1]
    }

    # Apply kernel transformation using Rcpp functions
    if (type == "laplacian") {
      return(laplacian_kernel_rcpp(data, bw))
    } else { # type == "gaussian"
      return(gaussian_kernel_rcpp(data, bw))
    }
  }

  if (type == "e-dist") {
    return(e_dist_rcpp(data, group))
  }

  if (type %in% c("g-dist", "l-dist")) {
    if (is.null(bw)) {
      # Automatically calculate bandwidth if not provided
      bw <- bw_optim(
        data,
        scale_factor = scale_factor,
        group = if(is.null(group)) 1:p else group
      )
    } else if (!is.vector(bw) || length(bw) != p) {
      stop("For induced kernel distances, bw must be a vector with length equal to the number of columns")
    }
    return(kernel_induced_distance_matrix(data, bw, type, group))
  }

  # Strict error for unknown distance types
  stop(paste("Unknown distance type:", type))
}

#' Calculate optimal bandwidth
#'
#' @param x Matrix or vector of data
#' @param sample_limit Maximum number of samples to use for calculation
#' @param scale_factor Scaling factor for the bandwidth
#' @param group Optional grouping of variables
#'
#' @return Bandwidth parameter value or vector
#'
#' @examples
#' # Working with a data frame (must have all numeric columns)
#' df <- data.frame(a = rnorm(20), b = rnorm(20), c = rnorm(20))
#' bw_df <- bw_optim(df)
#'
#' # The function will throw an error if non-numeric columns are present
#' \dontrun{
#' df_with_factors <- data.frame(a = rnorm(20), b = as.factor(1:20))
#' bw_error <- bw_optim(df_with_factors) # This will error
#' }
#' @export
bw_optim <- function(x, sample_limit = 1000, scale_factor = 0.5, group = NULL) {
  # Ensure x is a numeric matrix
  if (is.vector(x)) {
    if (!is.numeric(x)) {
      stop("Input vector must be numeric")
    }
    x <- matrix(x, ncol = 1)
  } else if (is.data.frame(x)) {
    # Check if all columns are numeric
    if (!all(sapply(x, is.numeric))) {
      stop("All data frame columns must be numeric")
    }
    x <- as.matrix(x)
  } else if (!is.matrix(x)) {
    stop("Input must be a numeric vector, matrix, or data.frame with numeric columns")
  }

  # Additional check for numerical matrix
  if (!is.numeric(x)) {
    stop("Input matrix must contain numeric values only")
  }

  # If no group is provided, treat the entire matrix as one group
  if (is.null(group)) {
    return(bw_rcpp(x, sample_limit, scale_factor))
  }

  # Validate group
  if (length(group) != ncol(x)) {
    stop("Length of group vector must match number of columns in x")
  }

  # Find maximum group and all used groups
  max_group <- max(group)
  used_groups <- sort(unique(group))

  # Create result vector of length max_group, all filled with NAs initially
  result <- rep(NA_real_, max_group)

  # Calculate bandwidth for each used group
  for (g in used_groups) {
    # Extract columns for this group
    group_cols <- which(group == g)
    group_data <- x[, group_cols, drop = FALSE]

    # Calculate bandwidth for this group
    result[g] <- bw_rcpp(group_data, sample_limit, scale_factor)
  }

  return(result)
}

#' Maximum Mean Discrepancy (MMD)
#'
#' Calculates the Maximum Mean Discrepancy between two samples, which is a measure of the
#' difference between distributions. This implementation is based on the method described in
#' Gretton et al. (2012).
#'
#' @param x First sample or full distance matrix if y is NULL
#' @param y Second sample (optional)
#' @param type Type of kernel or distance
#' @param bw Bandwidth parameter
#' @param expo Exponent parameter
#' @param scale_factor Scaling factor for bandwidth
#' @param group Optional grouping of variables
#' @param u_center Logical; use U-centering instead of V-centering
#' @param n Size of first sample (required if x is a distance matrix)
#' @param m Size of second sample (required if x is a distance matrix)
#'
#' @return MMD value
#' @examples
#' # Create two sample datasets
#' set.seed(123)
#' x <- matrix(rnorm(100), nrow = 20, ncol = 5)  # First sample
#' y <- matrix(rnorm(100, mean = 0.5), nrow = 20, ncol = 5)  # Second sample with different mean
#'
#' # Calculate MMD with Euclidean distance (default)
#' mmd_euclidean <- mmd(x, y)
#' print(mmd_euclidean)
#'
#' # Calculate MMD with Gaussian kernel
#' mmd_gaussian <- mmd(x, y, type = "gaussian")
#' print(mmd_gaussian)
#'
#' # Calculate MMD with Laplacian kernel
#' mmd_laplacian <- mmd(x, y, type = "laplacian")
#' print(mmd_laplacian)
#'
#' # Calculate MMD with polynomial kernel
#' mmd_polynomial <- mmd(x, y, type = "polynomial", expo = 2)
#' print(mmd_polynomial)
#'
#' # Using custom bandwidth parameter
#' mmd_custom_bw <- mmd(x, y, type = "gaussian", bw = 0.7)
#' print(mmd_custom_bw)
#'
#' # Compare V-centering and U-centering
#' mmd_v_center <- mmd(x, y, u_center = FALSE)  # Default
#' mmd_u_center <- mmd(x, y, u_center = TRUE)
#' print(c("V-centering:", mmd_v_center, "U-centering:", mmd_u_center))
#'
#' # Using vector inputs
#' x_vec <- rnorm(30)
#' y_vec <- rnorm(30, mean = 0.5)
#' mmd_vectors <- mmd(x_vec, y_vec)
#' print(mmd_vectors)
#'
#' # Using the group parameter for grouped variables
#' # Suppose the 5 columns are in 3 groups
#' group_vec <- c(1, 1, 2, 2, 3)
#' mmd_grouped <- mmd(x, y, type = "e-dist", group = group_vec)
#' print(mmd_grouped)
#'
#' # IMPORTANT: When providing a pre-computed distance matrix, make sure to specify
#' # the correct 'type' parameter that matches how the matrix was calculated
#'
#' # Example with a pre-computed Euclidean distance matrix
#' combined_data <- rbind(x, y)
#' D_euclidean <- KDist_matrix(combined_data, type = "euclidean")
#' mmd_euclidean <- mmd(D_euclidean, n = nrow(x), m = nrow(y), type = "euclidean")
#'
#' # Example with a pre-computed Gaussian kernel matrix
#' D_gaussian <- KDist_matrix(combined_data, type = "gaussian")
#' mmd_gaussian <- mmd(D_gaussian, n = nrow(x), m = nrow(y), type = "gaussian")
#'
#' # Using the wrong type parameter will give incorrect results!
#' # For example, this would be INCORRECT:
#' \dontrun{
#' # D_gaussian is a Gaussian kernel matrix but type="euclidean" is specified
#' mmd_wrong <- mmd(D_gaussian, n = nrow(x), m = nrow(y), type = "euclidean")
#' }
#'
#' @references
#' Gretton, A., Borgwardt, K. M., Rasch, M. J., Schölkopf, B., & Smola, A. (2012).
#' A kernel two-sample test. \emph{Journal of Machine Learning Research, 13}, 723-773.
#'
#' @export
mmd <- function(x, y = NULL, type = "euclidean", bw = NULL, expo = 1, scale_factor = 0.5, group = NULL, u_center = FALSE, n = NULL, m = NULL) {
  # Validate the type parameter first
  valid_types <- c("euclidean", "polynomial", "gaussian", "laplacian", "e-dist", "g-dist", "l-dist")
  if (!type %in% valid_types) {
    stop(paste("Invalid type parameter:", type,
               "Valid options are:", paste(valid_types, collapse = ", ")))
  }

  if (is.null(y)) {
    # If y is NULL, x is treated as the full distance matrix D
    D <- x

    # Check if the input is a valid distance matrix (should be square)
    if (!is.matrix(D) || nrow(D) != ncol(D)) {
      stop("When y is NULL, x must be a square distance matrix")
    }

    # Check if n and m are provided
    if (is.null(n) || is.null(m)) {
      stop("When y is NULL, n and m must be provided to specify the dimensions of the samples")
    }

    # Check if n and m are valid
    if (n <= 0 || m <= 0 || !is.numeric(n) || !is.numeric(m)) {
      stop("n and m must be positive numeric values")
    }

    # Check if n+m equals the size of D
    total_size <- nrow(D)
    if (n + m != total_size) {
      stop(paste("The sum of n and m must equal the size of the distance matrix (", total_size, ")", sep = ""))
    }

    # Add a reminder about specifying the correct type
    message("When providing a pre-computed distance matrix, make sure the 'type' parameter matches how the matrix was computed.")

    # Extract submatrices
    Dxx <- D[1:n, 1:n]
    Dyy <- D[(n+1):(n+m), (n+1):(n+m)]
    Dxy <- D[1:n, (n+1):(n+m)]
  } else {
    if (is.vector(x) && !is.matrix(x)) {
      x <- matrix(x, ncol = 1)
    }
    if (is.vector(y) && !is.matrix(y)) {
      y <- matrix(y, ncol = 1)
    }
    n <- nrow(x)
    m <- nrow(y)
    # Calculate the full kernel distance matrices
    D <- KDist_matrix(rbind(x,y), type = type, bw = bw, expo = expo, scale_factor = scale_factor, group = group)
    Dxx <- D[1:n, 1:n]
    Dyy <- D[(n+1):(n+m), (n+1):(n+m)]
    Dxy <- D[1:n, (n+1):(n+m)]
  }

  if(!u_center) {
    # Calculate Generalized Energy Distance based on distance type
    if (type %in% c("euclidean", "e-dist", "l-dist", "g-dist")) {
      # For distance-based types, MMD = 2*Dxy - Dxx - Dyy
      return(2 * mean(Dxy) - mean(Dxx) - mean(Dyy))
    } else {
      # For kernel-based types, MMD = Dxx + Dyy - 2*Dxy
      return(mean(Dxx) + mean(Dyy) - 2 * mean(Dxy))
    }
  } else {
    if (type %in% c("euclidean", "e-dist", "l-dist", "g-dist")) {
      # For distance-based types, MMD = 2*Dxy - Dxx - Dyy
      return(2 * mean(Dxy) - sum(Dxx)/n/(n-1) - sum(Dyy)/m/(m-1))
    } else {
      # For kernel-based types, MMD = Dxx + Dyy - 2*Dxy
      diag(Dxx) <- 0
      diag(Dyy) <- 0
      return(sum(Dxx)/n/(n-1) + sum(Dyy)/m/(m-1) - 2 * mean(Dxy))
    }
  }
}

#' Hilbert-Schmidt Independence Criterion (HSIC)
#'
#' Calculates the HSIC between two datasets or matrices, which measures the
#' dependence between random variables. This implementation is based on the method
#' described in Gretton et al. (2007).
#'
#' @param x First dataset or distance matrix
#' @param y Second dataset or distance matrix
#' @param type Type of kernel or distance
#' @param bw Bandwidth parameter
#' @param expo Exponent parameter
#' @param scale_factor Scaling factor for bandwidth
#' @param group Optional list of length 2 specifying group membership for each column in x and y.
#'   The first element of the list should specify grouping for columns in x, and the
#'   second element should specify grouping for columns in y. Used for group-wise
#'   distance calculations in "e-dist", "g-dist", or "l-dist".
#' @param u_center Logical; use U-centering instead of V-centering
#' @param is_distance Logical; whether input matrices are already distance matrices
#'
#' @return HSIC value
#'
#' @examples
#' # Create example datasets
#' set.seed(123)
#'
#' # Example 1: Independent variables
#' x1 <- matrix(rnorm(100), ncol = 1)
#' y1 <- matrix(rnorm(100), ncol = 1)
#' hsic_indep <- hsic(x1, y1, type = "gaussian")
#' cat("HSIC for independent variables:", hsic_indep, "\n")
#'
#' # Example 2: Linear relationship
#' x2 <- matrix(1:100, ncol = 1)
#' y2 <- matrix(2*x2 + rnorm(100, sd = 3), ncol = 1)
#' hsic_linear <- hsic(x2, y2, type = "gaussian")
#' cat("HSIC for linear relationship:", hsic_linear, "\n")
#'
#' # Example 3: Non-linear relationship
#' x3 <- matrix(seq(-3, 3, length.out = 100), ncol = 1)
#' y3 <- matrix(sin(x3) + rnorm(100, sd = 0.2), ncol = 1)
#' hsic_nonlinear <- hsic(x3, y3, type = "gaussian")
#' cat("HSIC for non-linear relationship:", hsic_nonlinear, "\n")
#'
#' # Example 4: Using different kernels/distances
#' hsic_gaussian <- hsic(x3, y3, type = "gaussian")
#' hsic_laplacian <- hsic(x3, y3, type = "laplacian")
#' hsic_euclidean <- hsic(x3, y3, type = "euclidean")
#' cat("HSIC values with different kernels:",
#'     "\n  Gaussian:", hsic_gaussian,
#'     "\n  Laplacian:", hsic_laplacian,
#'     "\n  Euclidean:", hsic_euclidean, "\n")
#'
#' # Example 5: Using custom bandwidth
#' hsic_bw_small <- hsic(x3, y3, type = "gaussian", bw = 0.1)
#' hsic_bw_large <- hsic(x3, y3, type = "gaussian", bw = 2.0)
#' cat("HSIC with different bandwidths:",
#'     "\n  Small bandwidth (0.1):", hsic_bw_small,
#'     "\n  Large bandwidth (2.0):", hsic_bw_large, "\n")
#'
#' # Example 6: Using U-centering (unbiased estimator)
#' hsic_v <- hsic(x3, y3, type = "gaussian", u_center = FALSE)
#' hsic_u <- hsic(x3, y3, type = "gaussian", u_center = TRUE)
#' cat("HSIC with different centering:",
#'     "\n  V-centering (biased):", hsic_v,
#'     "\n  U-centering (unbiased):", hsic_u, "\n")
#'
#' # Example 7: Using pre-computed distance matrices
#' Dx <- KDist_matrix(x3, type = "gaussian")
#' Dy <- KDist_matrix(y3, type = "gaussian")
#' hsic_precomp <- hsic(Dx, Dy, type = "gaussian", is_distance = TRUE)
#' cat("HSIC with pre-computed matrices:", hsic_precomp, "\n")
#'
#' # Example 8: Using grouped variables with hsic
#' # Create sample data
#' x <- matrix(rnorm(100*5), ncol = 5)  # 5 variables in x
#' y <- matrix(rnorm(100*4), ncol = 4)  # 4 variables in y
#'
#' # Define group structure: x has 2 groups, y has 2 groups
#' x_groups <- c(1, 1, 2, 2, 2)  # First 2 vars in group 1, last 3 in group 2
#' y_groups <- c(1, 1, 2, 2)     # First 2 vars in group 1, last 2 in group 2
#'
#' # Combine into a list as required by hsic
#' group_list <- list(x_groups, y_groups)
#'
#' # Calculate HSIC with the group structure
#' hsic_result <- hsic(x, y, type = "e-dist", group = group_list)
#' print(hsic_result)
#'
#' @references
#' Gretton, A., Fukumizu, K., Teo, C., Song, L., Schölkopf, B., & Smola, A. (2007).
#' A kernel statistical test of independence. \emph{Advances in neural information processing systems, 20}.
#'
#' @export
hsic <- function(x, y, type = "gaussian", bw = NULL, expo = 1, scale_factor = 0.5, group = NULL, u_center = FALSE, is_distance = FALSE) {
  # Convert vectors to single-column matrices if needed
  if (is.vector(x) && !is.matrix(x)) {
    x <- matrix(x, ncol = 1)
  }

  if (is.vector(y) && !is.matrix(y)) {
    y <- matrix(y, ncol = 1)
  }

  # Call the C++ function with matrices
  hsic_cpp(x, y, type, bw, expo, scale_factor, group, u_center, is_distance)
}

#' Cross-Sample Hilbert-Schmidt Independence Criterion
#'
#' Calculates the cross HSIC between two samples. This is an extension of the notion of
#' cross distance covariance (cdCov) proposed in Chakraborty & Zhang (2021). This quantity
#' appears in the asymptotic variance of MMD or (generalized) energy distance in the
#' high-dimensional setting.
#'
#' @param x First sample or cross-distance matrix
#' @param y Second sample (optional)
#' @param type Type of kernel or distance
#' @param bw Bandwidth parameter
#' @param expo Exponent parameter
#' @param scale_factor Scaling factor for bandwidth
#' @param group Optional grouping of variables
#'
#' @return Cross HSIC value
#'
#' @references
#' Chakraborty, S., & Zhang, X. (2021). A New Framework for Distance and
#' Kernel-based Metrics in High-dimension. \emph{Electronic Journal of Statistics, 15},
#' 5455-5522.
#'
#' @examples
#' # Example 1: Basic usage with randomly generated data
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol = 1)
#' y <- matrix(rnorm(100), ncol = 1)
#' chsic_value <- chsic(x, y)
#' print(chsic_value)
#'
#' # Example 2: Using dependent data
#' x2 <- matrix(1:100, ncol = 1)
#' y2 <- matrix(1:100 + rnorm(100, sd = 0.5), ncol = 1)
#' chsic_value2 <- chsic(x2, y2)
#' print(chsic_value2)
#'
#' # Example 3: Using different distance/kernel types
#' chsic_euclidean <- chsic(x2, y2, type = "euclidean")
#' chsic_gaussian <- chsic(x2, y2, type = "gaussian")
#' chsic_laplacian <- chsic(x2, y2, type = "laplacian")
#' print(c(euclidean = chsic_euclidean,
#'        gaussian = chsic_gaussian,
#'        laplacian = chsic_laplacian))
#'
#' # Example 4: Using a pre-computed cross-distance matrix
#' # First combine the samples
#' combined_x <- rbind(x2, y2)
#' # Calculate the full distance matrix
#' D <- KDist_matrix(combined_x, type = "euclidean")
#' # Extract the cross-distance matrix (upper-right quadrant)
#' n <- nrow(x2)
#' m <- nrow(y2)
#' D_cross <- D[1:n, (n+1):(n+m)]
#' # Calculate chsic using the cross-distance matrix
#' chsic_precomputed <- chsic(D_cross)
#' print(chsic_precomputed)
#'
#' # Example 5: Using chsic with vector inputs
#' x_vec <- rnorm(100)
#' y_vec <- rnorm(100)
#' chsic_vec <- chsic(x_vec, y_vec)
#' print(chsic_vec)
#'
#' # Example 6: Using multivariate data
#' x_multi <- matrix(rnorm(200), ncol = 2)
#' y_multi <- matrix(rnorm(200), ncol = 2)
#' chsic_multi <- chsic(x_multi, y_multi)
#' print(chsic_multi)
#'
#' @export
chsic <- function(x, y = NULL, type = "gaussian", bw = NULL, expo = 1, scale_factor = 0.5, group = NULL) {
  if (is.null(y)) {
    # If y is NULL, x is treated as the distance matrix D
    D <- x

    # Check if the input is a valid matrix
    if (!is.matrix(D)) {
      stop("When y is NULL, x must be a matrix")
    }

    # Get dimensions
    n <- nrow(D)
    m <- ncol(D)
  } else {
    # If both x and y are provided, calculate the cross-distance matrix
    # Convert vectors to single-column matrices if needed
    if (is.vector(x) && !is.matrix(x)) {
      x <- matrix(x, ncol = 1)
    }

    if (is.vector(y) && !is.matrix(y)) {
      y <- matrix(y, ncol = 1)
    }
    n <- nrow(x)
    m <- nrow(y)

    # Calculate the full distance matrix
    A <- KDist_matrix(rbind(x, y), type = type, bw = bw, expo = expo, scale_factor = scale_factor, group = group)

    # Extract the cross-distance matrix
    D <- A[1:n, (n+1):(n+m)]
  }

  # Use the optimized C++ function for calculation
  return(chsic_cpp(D))
}

#' HSIC-based Correlation Coefficient
#'
#' Calculates a correlation coefficient based on the Hilbert-Schmidt Independence Criterion,
#' providing a normalized measure of dependence between variables. This is proposed in
#' Gretton et al. (2007) as the kernel version of the distance correlation later developed
#' by Székely et al. (2007).
#'
#' @param x First dataset or distance matrix
#' @param y Second dataset or distance matrix
#' @param type Type of kernel or distance
#' @param bw Bandwidth parameter
#' @param expo Exponent parameter
#' @param scale_factor Scaling factor for bandwidth
#' @param group Optional list of length 2 specifying group membership for each column in x and y.
#'   The first element of the list should specify grouping for columns in x, and the
#'   second element should specify grouping for columns in y. Used for group-wise
#'   distance calculations in "e-dist", "g-dist", or "l-dist".
#' @param u_center Logical; use U-centering instead of V-centering
#' @param is_distance Logical; whether input matrices are already distance matrices
#'
#' @return HSIC correlation coefficient
#'
#' @references
#' Gretton, A., Fukumizu, K., Teo, C., Song, L., Schölkopf, B., & Smola, A. (2007).
#' A kernel statistical test of independence. \emph{Advances in neural information processing systems, 20}.
#'
#' Székely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing dependence by
#' correlation of distances. \emph{The Annals of Statistics, 35}(6), 2769-2794.
#'
#' @examples
#' # Example 1: Independent variables
#' set.seed(123)
#' x1 <- matrix(rnorm(100), ncol = 1)
#' y1 <- matrix(rnorm(100), ncol = 1)
#' hsic_corr1 <- hsic_cor(x1, y1, type = "gaussian")
#' # Should be close to 0
#' print(hsic_corr1)
#'
#' # Example 2: Linear relationship
#' x2 <- matrix(1:100, ncol = 1)
#' y2 <- matrix(2*x2 + rnorm(100, sd = 5), ncol = 1)
#' hsic_corr2 <- hsic_cor(x2, y2, type = "gaussian")
#' # Should be close to 1
#' print(hsic_corr2)
#'
#' # Example 3: Non-linear relationship
#' x3 <- matrix(seq(-3, 3, length.out = 100), ncol = 1)
#' y3 <- matrix(x3^2 + rnorm(100, sd = 0.5), ncol = 1)
#' hsic_corr3 <- hsic_cor(x3, y3, type = "gaussian")
#' # Should detect the non-linear relationship
#' print(hsic_corr3)
#'
#' # Compare with Pearson correlation
#' cor_pearson3 <- cor(x3, y3)
#' # Pearson correlation may miss the non-linear relationship
#' print(paste("HSIC correlation:", hsic_corr3, "Pearson correlation:", cor_pearson3))
#'
#' # Example 4: Using different kernel types
#' hsic_gaussian <- hsic_cor(x3, y3, type = "gaussian")
#' hsic_laplacian <- hsic_cor(x3, y3, type = "laplacian")
#' hsic_euclidean <- hsic_cor(x3, y3, type = "euclidean")
#' print(c(gaussian = hsic_gaussian, laplacian = hsic_laplacian, euclidean = hsic_euclidean))
#'
#' # Example 5: Using pre-computed distance matrices
#' Dx <- KDist_matrix(x3, type = "gaussian")
#' Dy <- KDist_matrix(y3, type = "gaussian")
#' hsic_precomp <- hsic_cor(Dx, Dy, type = "gaussian", is_distance = TRUE)
#' print(hsic_precomp)
#'
#' # Example 6: Perfect correlation (identical variables)
#' hsic_perfect <- hsic_cor(x2, x2, type = "gaussian")
#' # Should be exactly 1
#' print(hsic_perfect)
#'
#' # Example 7: Sinusoidal relationship
#' t <- seq(0, 4*pi, length.out = 100)
#' x7 <- matrix(sin(t), ncol = 1)
#' y7 <- matrix(cos(t), ncol = 1)
#' hsic_sin <- hsic_cor(x7, y7, type = "gaussian")
#' cor_pearson7 <- cor(x7, y7)
#' # HSIC should detect this non-linear relationship
#' print(paste("HSIC correlation:", hsic_sin, "Pearson correlation:", cor_pearson7))
#'
#' # Example 8: Using grouped variables with hsic_cor
#' # Create sample data
#' x <- matrix(rnorm(100*4), ncol = 4)  # 4 variables in x
#' y <- matrix(rnorm(100*3), ncol = 3)  # 3 variables in y
#'
#' # Create some dependency between the datasets
#' y[, 1:2] <- y[, 1:2] + 0.8 * x[, 1:2]
#'
#' # Define group structure for both datasets
#' x_groups <- c(1, 1, 2, 2)    # First 2 vars in group 1, last 2 in group 2
#' y_groups <- c(1, 1, 2)       # First 2 vars in group 1, last 1 in group 2
#'
#' # Combine into a list as required for grouped analysis
#' group_list <- list(x_groups, y_groups)
#'
#' # Calculate HSIC correlation with the group structure
#' cor_grouped <- hsic_cor(x, y, type = "e-dist", group = group_list)
#' print(paste("HSIC correlation with grouped variables:", cor_grouped))
#'
#' # Compare with standard (non-grouped) HSIC correlation
#' cor_standard <- hsic_cor(x, y, type = "gaussian")
#' print(paste("Standard HSIC correlation:", cor_standard))
#'
#' @export
hsic_cor <- function(x, y, type = "gaussian", bw = NULL, expo = 1, scale_factor = 0.5, group = NULL, u_center = FALSE, is_distance = FALSE) {
  # Check if x and y are the same object (by memory address)
  same_data <- identical(x, y)

  # If x and y are identical, return 1.0 immediately
  if (same_data) {
    return(1.0)  # Perfect correlation when x and y are identical
  }

  # Convert vectors to single-column matrices if needed
  if (is.vector(x) && !is.matrix(x)) {
    x <- matrix(x, ncol = 1)
  }

  if (is.vector(y) && !is.matrix(y)) {
    y <- matrix(y, ncol = 1)
  }

  if (is_distance) {
    # Use x and y directly as distance matrices
    Dxx <- x
    Dyy <- y
  } else {
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

    # Compute distance matrices with appropriate group parameter
    Dxx <- KDist_matrix(x, type, bw, expo, scale_factor, group_x)
    Dyy <- KDist_matrix(y, type, bw, expo, scale_factor, group_y)
  }

  # Calculate HSIC values
  hsic_val <- hsic_cpp(Dxx, Dyy, type, u_center = u_center, is_distance = TRUE)
  hsic_val_x <- hsic_cpp(Dxx, Dxx, type, u_center = u_center, is_distance = TRUE)
  hsic_val_y <- hsic_cpp(Dyy, Dyy, type, u_center = u_center, is_distance = TRUE)

  # Check for numerical stability
  denominator <- sqrt(hsic_val_x * hsic_val_y)
  if (denominator < .Machine$double.eps) {
    warning("Denominator close to zero, returning NA")
    return(NA)
  }

  # Calculate correlation and ensure it's in [-1, 1] range
  cor_val <- hsic_val / denominator
  return(max(min(cor_val, 1), -1))
}

#' Energy Distance Between Two Samples
#'
#' Calculates the energy distance between two multivariate samples, which is a measure of the
#' difference between distributions.
#'
#' @param x First sample or full distance matrix if y is NULL
#' @param y Second sample (optional)
#' @param a Exponent parameter for Euclidean distance, default is 1
#' @param u_center Logical; use U-centering instead of V-centering
#' @param n Size of first sample (required if x is a distance matrix)
#' @param m Size of second sample (required if x is a distance matrix)
#'
#' @details
#' The energy distance between two samples is defined as:
#'
#' \deqn{E(X, Y) = 2E[||X - Y||^a] - E[||X - X'||^a] - E[||Y - Y'||^a]}
#'
#' where \eqn{||\cdot||} denotes Euclidean distance, and \eqn{(X', Y')} is an independent copy of \eqn{(X, Y)}.
#'
#' This implementation is based on the method proposed by Székely and Rizzo (2004).
#' Energy distance is a measure of statistical distance between probability
#' distributions, and has the property that \eqn{E(F, G)\geq 0} with equality if and
#' only if \eqn{F = G} (where \eqn{F} and \eqn{G} are the distributions from which \eqn{X} and \eqn{Y} are sampled).
#'
#' @return The energy distance value between the two samples
#'
#' @references
#' Székely, G. J., & Rizzo, M. L. (2004). Testing for Equal Distributions in High Dimension.
#' InterStat, Nov. (5).
#'
#' Székely, G. J., & Rizzo, M. L. (2013). Energy statistics: A class of statistics based on
#' distances. Journal of Statistical Planning and Inference, 143(8), 1249-1272.
#'
#' @examples
#' # Example 1: Energy distance between samples from the same distribution
#' set.seed(123)
#' x <- matrix(rnorm(100 * 3), ncol = 3)  # 100 observations in 3D from N(0,1)
#' y <- matrix(rnorm(80 * 3), ncol = 3)   # 80 observations in 3D from N(0,1)
#' ed_same <- ed(x, y)
#' print(paste("Energy distance (same distribution):", round(ed_same, 4)))
#'
#' # Example 2: Energy distance between samples from different distributions
#' set.seed(456)
#' x <- matrix(rnorm(100 * 2), ncol = 2)  # 100 observations in 2D from N(0,1)
#' y <- matrix(rnorm(100 * 2, mean = 1), ncol = 2)  # 100 observations in 2D from N(1,1)
#' ed_diff <- ed(x, y)
#' print(paste("Energy distance (different means):", round(ed_diff, 4)))
#'
#' # Example 3: Using different exponent values
#' ed_a1 <- ed(x, y, a = 1)  # Default
#' ed_a2 <- ed(x, y, a = 2)  # Squared Euclidean distance
#' print(c("a=1" = ed_a1, "a=2" = ed_a2))
#'
#' # Example 4: Using U-centering instead of V-centering
#' ed_v <- ed(x, y, u_center = FALSE)  # Default, V-statistic
#' ed_u <- ed(x, y, u_center = TRUE)   # U-statistic (unbiased)
#' print(c("V-statistic" = ed_v, "U-statistic" = ed_u))
#'
#' # Example 5: Using vector inputs
#' x_vec <- rnorm(50)
#' y_vec <- rnorm(50, mean = 0.5)
#' ed_vec <- ed(x_vec, y_vec)
#' print(paste("Energy distance with vectors:", round(ed_vec, 4)))
#'
#' # Example 6: Using a pre-computed distance matrix
#' combined <- rbind(x, y)
#' D <- as.matrix(dist(combined)^1)  # ^1 for a=1
#' n_x <- nrow(x)
#' m_y <- nrow(y)
#' ed_precomp <- ed(D, n = n_x, m = m_y)
#' print(paste("Using pre-computed matrix:", round(ed_precomp, 4)))
#'
#' @seealso \code{\link{mmd}} for Maximum Mean Discrepancy, \code{\link{dcov}} for distance covariance
#' @export
ed <- function(x, y = NULL, a = 1, u_center = FALSE, n = NULL, m = NULL) {
  # Handle the pre-computed distance matrix case directly
  if (is.null(y) && is.matrix(x)) {
    # Check if n and m are provided
    if (is.null(n) || is.null(m)) {
      stop("When y is NULL, n and m must be provided to specify the dimensions of the samples")
    }

    # Validate inputs
    total_size <- nrow(x)
    if (n + m != total_size) {
      stop(paste("The sum of n and m must equal the size of the distance matrix (",
                 total_size, ")", sep = ""))
    }

    # Extract submatrices
    Dxx <- x[1:n, 1:n]
    Dyy <- x[(n+1):(n+m), (n+1):(n+m)]
    Dxy <- x[1:n, (n+1):(n+m)]

    # Calculate energy distance
    if (!u_center) {
      # V-statistic version
      return(2 * mean(Dxy) - mean(Dxx) - mean(Dyy))
    } else {
      # U-statistic version
      return(2 * mean(Dxy) - sum(Dxx)/n/(n-1) - sum(Dyy)/m/(m-1))
    }
  } else {
    # Normal case - just call mmd
    return(mmd(x = x, y = y, type = "euclidean", expo = a,
               u_center = u_center, n = n, m = m))
  }
}

#' Distance Covariance Between Two Samples
#'
#' Calculates the squared distance covariance between two multivariate samples, which
#' measures the dependence between random vectors.
#'
#' @param x First dataset or distance matrix
#' @param y Second dataset or distance matrix
#' @param a Exponent parameter for Euclidean distance, default is 1
#' @param u_center Logical; use U-centering instead of V-centering (default: FALSE)
#' @param is_distance Logical; whether input matrices are already distance matrices (default: FALSE)
#'
#' @details
#' The squared distance covariance between two samples is defined as:
#'
#' \deqn{dcov^2(X, Y) = E[||X - X'||^a · ||Y - Y'||^a] + E[||X - X'||^a] · E[||Y - Y'||^a]
#'                - 2E[||X - X'||^a · ||Y - Y''||^a]}
#'
#' where \eqn{||\cdot||} denotes Euclidean distance, \eqn{(X', Y')} and \eqn{(X'', Y'')} are independent copies
#' of \eqn{(X, Y)}, and a is typically set to 1.
#'
#' Distance covariance has several important properties:
#' \itemize{
#'   \item \eqn{dcov^2(X, Y) \geq 0} for all distributions of \eqn{X} and \eqn{Y}
#'   \item \eqn{dcov^2(X, Y) = 0} if and only if \eqn{X} and \eqn{Y} are independent
#'   \item It can detect both linear and non-linear dependence
#'   \item It is applicable to vectors of arbitrary dimensions
#' }
#'
#' When u_center = FALSE (default), the function computes the V-statistic version.
#' When u_center = TRUE, it computes the U-statistic version that's unbiased.
#'
#' @return The squared distance covariance value between the two samples
#'
#' @references
#' Székely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing dependence by
#' correlation of distances. *The Annals of Statistics*, 35(6), 2769-2794.
#'
#' @examples
#' # Example 1: Distance covariance between independent samples
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol = 1)
#' y <- matrix(rnorm(100), ncol = 1)
#' dcov_indep <- dcov(x, y)
#' print(paste("Distance covariance (independent):", round(dcov_indep, 6)))
#'
#' # Example 2: Distance covariance with linear relationship
#' x <- matrix(runif(100), ncol = 1)
#' y <- matrix(x + 0.1*rnorm(100), ncol = 1)  # Linear relationship with noise
#' dcov_linear <- dcov(x, y)
#' print(paste("Distance covariance (linear):", round(dcov_linear, 6)))
#'
#' # Example 3: Distance covariance with non-linear relationship
#' x <- matrix(runif(100, -3, 3), ncol = 1)
#' y <- matrix(x^2 + 0.1*rnorm(100), ncol = 1)  # Quadratic relationship
#' dcov_nonlinear <- dcov(x, y)
#' print(paste("Distance covariance (non-linear):", round(dcov_nonlinear, 6)))
#'
#' # Example 4: Comparing V-statistic and U-statistic estimators
#' dcov_v <- dcov(x, y, u_center = FALSE)  # Default, V-statistic
#' dcov_u <- dcov(x, y, u_center = TRUE)   # U-statistic (unbiased)
#' print(c("V-statistic" = dcov_v, "U-statistic" = dcov_u))
#'
#' # Example 5: Using multivariate data
#' x_multi <- matrix(rnorm(200), ncol = 2)
#' y_multi <- cbind(x_multi[,1] + rnorm(100, sd = 0.1),
#'                  x_multi[,2]^2 + rnorm(100, sd = 0.1))
#' dcov_multi <- dcov(x_multi, y_multi)
#' print(paste("Multivariate distance covariance:", round(dcov_multi, 6)))
#'
#' # Example 6: Using pre-computed distance matrices
#' Dx <- as.matrix(dist(x_multi))
#' Dy <- as.matrix(dist(y_multi))
#' dcov_precomp <- dcov(Dx, Dy, is_distance = TRUE)
#' print(paste("Using pre-computed matrices:", round(dcov_precomp, 6)))
#'
#' @seealso
#' \code{\link{dcor}} for distance correlation,
#' \code{\link{hsic}} for Hilbert-Schmidt Independence Criterion,
#' \code{\link{ed}} for energy distance
#'
#' @export
dcov <- function(x, y, a = 1, u_center = FALSE, is_distance = FALSE) {
  # Call hsic with type="euclidean"
  return(hsic(x = x, y = y, type = "euclidean", expo = a,
              u_center = u_center, is_distance = is_distance))
}

#' Distance Correlation Between Two Samples
#'
#' Calculates the squared distance correlation between two multivariate samples, which
#' provides a scaled measure of dependence between random vectors.
#'
#' @param x First dataset or distance matrix
#' @param y Second dataset or distance matrix
#' @param a Exponent parameter for Euclidean distance, default is 1
#' @param u_center Logical; use U-centering instead of V-centering (default: FALSE)
#' @param is_distance Logical; whether input matrices are already distance matrices (default: FALSE)
#'
#' @details
#' The squared distance correlation between two samples is defined as:
#'
#' \deqn{dcor^2(X, Y) = dcov^2(X, Y) / \sqrt{dcov^2(X, X) dcov^2(Y, Y)}}
#'
#' where \eqn{dcov^2} is the squared distance covariance. Distance correlation normalizes
#' the distance covariance to produce a measure between 0 and 1.
#'
#' Distance correlation has these important properties:
#' \itemize{
#'   \item \eqn{0 \leq dcor^2(X, Y) \leq 1}
#'   \item \eqn{dcor^2(X, Y) = 0} if and only if \eqn{X} and \eqn{Y} are independent
#'   \item \eqn{dcor^2(X, Y) = 1} only when \eqn{X} and \eqn{Y} are related by a linear transformation
#'   \item It can detect both linear and non-linear dependence
#'   \item It is applicable to vectors of arbitrary and different dimensions
#' }
#'
#' When u_center = FALSE (default), the function computes the V-statistic version.
#' When u_center = TRUE, it computes the U-statistic version that's unbiased.
#'
#' @return The squared distance correlation value between the two samples (between 0 and 1)
#'
#' @references
#' Székely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing dependence by
#' correlation of distances. *The Annals of Statistics*, 35(6), 2769-2794.
#'
#' @examples
#' # Example 1: Distance correlation between independent samples
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol = 1)
#' y <- matrix(rnorm(100), ncol = 1)
#' dcor_indep <- dcor(x, y)
#' print(paste("Distance correlation (independent):", round(dcor_indep, 6)))
#'
#' # Example 2: Distance correlation with linear relationship
#' x <- matrix(runif(100), ncol = 1)
#' y <- matrix(2*x + 0.1*rnorm(100), ncol = 1)  # Linear relationship with noise
#' dcor_linear <- dcor(x, y)
#' print(paste("Distance correlation (linear):", round(dcor_linear, 6)))
#'
#' # Example 3: Distance correlation with non-linear relationship
#' x <- matrix(runif(100, -3, 3), ncol = 1)
#' y <- matrix(sin(x) + 0.1*rnorm(100), ncol = 1)  # Sinusoidal relationship
#' dcor_nonlinear <- dcor(x, y)
#' print(paste("Distance correlation (non-linear):", round(dcor_nonlinear, 6)))
#'
#' # Compare with Pearson correlation
#' pearson_cor <- cor(x, y)[1,1]
#' print(paste("Distance correlation:", round(dcor_nonlinear, 6),
#'             "Pearson correlation:", round(pearson_cor, 6)))
#'
#' # Example 4: Comparing V-statistic and U-statistic estimators
#' dcor_v <- dcor(x, y, u_center = FALSE)  # Default, V-statistic
#' dcor_u <- dcor(x, y, u_center = TRUE)   # U-statistic (unbiased)
#' print(c("V-statistic" = dcor_v, "U-statistic" = dcor_u))
#'
#' # Example 5: Using multivariate data with different dimensions
#' x_multi <- matrix(rnorm(100*3), ncol = 3)  # 3-dimensional
#' y_multi <- matrix(rnorm(100*2), ncol = 2)  # 2-dimensional
#' # Make them dependent
#' y_multi[,1] <- x_multi[,1] + rnorm(100, sd = 0.1)
#' y_multi[,2] <- x_multi[,2]^2 + rnorm(100, sd = 0.1)
#' dcor_multi <- dcor(x_multi, y_multi)
#' print(paste("Multivariate distance correlation:", round(dcor_multi, 6)))
#'
#' # Example 6: Using pre-computed distance matrices
#' Dx <- as.matrix(dist(x_multi))
#' Dy <- as.matrix(dist(y_multi))
#' dcor_precomp <- dcor(Dx, Dy, is_distance = TRUE)
#' print(paste("Using pre-computed matrices:", round(dcor_precomp, 6)))
#'
#' # Example 7: Perfect linear relationship
#' x_perfect <- matrix(1:100, ncol = 1)
#' y_perfect <- matrix(3*x_perfect + 5, ncol = 1)  # Perfect linear relationship
#' dcor_perfect <- dcor(x_perfect, y_perfect)
#' print(paste("Perfect linear relationship:", round(dcor_perfect, 6)))
#'
#' @seealso
#' \code{\link{dcov}} for distance covariance,
#' \code{\link{hsic_cor}} for kernel-based correlation,
#' \code{\link{cor}} for Pearson correlation
#'
#' @export
dcor <- function(x, y, a = 1, u_center = FALSE, is_distance = FALSE) {
  # Call hsic_cor with type="euclidean"
  return(hsic_cor(x = x, y = y, type = "euclidean", expo = a,
                  u_center = u_center, is_distance = is_distance))
}

#' Martingale Difference Divergence (MDD)
#'
#' Calculates the martingale difference divergence between a predictor matrix
#' and response variable or matrix. MDD measures the conditional mean dependence
#' and can be used for testing conditional mean independence, nonlinear feature screening,
#' and dimension reduction for stationary multivariate time series.
#'
#' @param x The predictor dataset or distance matrix
#' @param y The response variable (vector) or matrix
#' @param type Type of kernel or distance for the predictor (default: "euclidean")
#' @param bw Bandwidth parameter (default: NULL for automatic selection)
#' @param expo Exponent parameter (default: 1)
#' @param scale_factor Scaling factor for bandwidth (default: 0.5)
#' @param group_x Optional grouping of predictor variables
#' @param u_center Logical; use U-centering instead of V-centering (default: FALSE)
#' @param is_distance Logical; whether the x input is already a distance matrix (default: FALSE)
#'
#' @details
#' Martingale Difference Divergence (MDD) is a measure of conditional mean dependence.
#' MDD = 0 if and only if E(Y|X) = E(Y) almost surely, which characterizes conditional
#' mean independence. This makes it particularly useful for:
#'
#' 1. Testing conditional mean independence
#' 2. High-dimensional variable screening
#' 3. Feature selection based on conditional mean dependence
#' 4. Dimension reduction for stationary multivariate time series
#'
#' For vector y, MDD returns a scalar measuring the conditional mean dependence.
#' For matrix y, MDD returns a p×p martingale difference divergence matrix.
#'
#' When u_center = TRUE, the function uses an unbiased estimator based on U-statistics,
#' which can improve performance in smaller sample sizes.
#'
#' @return
#' If y is a vector, returns a scalar measuring conditional mean dependence.
#' If y is a matrix, returns a p×p martingale difference divergence matrix.
#'
#' @examples
#' # Example 1: MDD with vector response
#' set.seed(123)
#' x <- matrix(rnorm(100*3), ncol = 3)
#' y <- x[,1]^2 + x[,2] + rnorm(100, sd = 0.5)
#' mdd_val <- mdd(x, y)
#' print(mdd_val)
#'
#' # Example 2: MDD with matrix response
#' y_mat <- cbind(x[,1]^2, x[,3] + rnorm(100, sd = 0.5))
#' mdd_mat <- mdd(x, y_mat)
#' print(mdd_mat)
#'
#' # Example 3: Using U-centering for improved performance with smaller samples
#' mdd_u <- mdd(x, y, u_center = TRUE)
#' print(mdd_u)
#'
#' # Example 4: Variable screening example
#' p <- 10
#' n <- 100
#' x_large <- matrix(rnorm(n*p), ncol = p)
#' y <- x_large[,1]^2 + 0.5*x_large[,3] + rnorm(n, sd = 0.5)
#' # Calculate MDD for each predictor
#' mdd_values <- numeric(p)
#' for(j in 1:p) {
#'   mdd_values[j] <- mdd(x_large[,j,drop=FALSE], y)
#' }
#' # Identify most important predictors
#' important_vars <- order(mdd_values, decreasing = TRUE)[1:3]
#' print(important_vars)  # Should identify variables 1 and 3 as important
#'
#' @references
#' Shao, X., & Zhang, J. (2014). Martingale difference correlation and its use in
#' high-dimensional variable screening. Journal of the American Statistical
#' Association, 109(507), 1302-1318.
#'
#' Zhang, X., Yao, S., & Shao, X. (2018). Conditional mean and quantile dependence
#' testing in high dimension. Journal of the American Statistical Association,
#' 113(524), 1763-1776.
#'
#' Lee, C. E., & Shao, X. (2018). Martingale difference divergence matrix and its
#' application to dimension reduction for stationary multivariate time series.
#' Journal of the American Statistical Association, 113(521), 216-229.
#'
#' @seealso
#' \code{\link{dcov}} for distance covariance,
#' \code{\link{hsic}} for Hilbert-Schmidt Independence Criterion
#'
#' @export
mdd <- function(x, y, type = "euclidean", bw = NULL, expo = 1, scale_factor = 0.5,
                group_x = NULL, u_center = FALSE, is_distance = FALSE) {
  # Handle vectors
  if (is.vector(x) && !is.matrix(x)) {
    x <- as.matrix(x)
  }

  # Call the C++ function
  mdd_cpp(x, y, type, bw, expo, scale_factor, group_x, u_center, is_distance)
}

