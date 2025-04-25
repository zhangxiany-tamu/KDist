#' Retrieve a statistic from the internal asymptotic_stats list
#'
#' This function retrieves a statistic from the internal dataset `asymptotic_stats`
#' by a given key.
#'
#' @param key A character string specifying the name of the statistic to retrieve.
#'
#' @return The corresponding statistic (could be a number, list, etc.).
#'
#' @examples
#' \dontrun{
#' my_function("mean_stat")
#' }
#'
#' @export
my_function <- function(key) {
  if (!key %in% names(asymptotic_stats)) {
    stop("Invalid key.")
  }
  return(asymptotic_stats[[key]])
}

#' Single Change-Point Detection Using Kernel and Distance-Based Metrics
#'
#' @description
#' Detects a single change point in sequential data using kernel-based generalized
#' homogeneity metrics. The method is based on the high-dimensional t-test framework proposed in
#' Chakraborty & Zhang (2021) and implemented for change-point detection as described in
#' Chakraborty & Zhang (2021, arXiv:2105.08976).
#'
#' @details
#' This function implements a kernel and distance-based approach for detecting a single structural change in multivariate
#' time series data, building on the high-dimensional t-test framework for distance and kernel-based metrics.
#' The method works by:
#'
#' 1. Computing a distance/kernel matrix based on the input data
#' 2. Calculating test statistics based on the high-dimensional t-test framework using MMD (Maximum Mean Discrepancy)
#'    and HSIC (Hilbert-Schmidt Independence Criterion) measures
#' 3. Identifying the location with the maximum statistic value as the estimated change point
#' 4. Computing a p-value through permutation testing or asymptotic approximation
#'
#' The approach is especially effective for high-dimensional data where traditional methods might fail,
#' as it operates in kernel space rather than directly on the raw features, allowing for detection of changes in higher
#' order moments even when the number of variables exceeds the sample size.
#'
#' @param data A matrix or data frame with rows representing time points (observations) and
#'        columns representing variables or features.
#' @param type Type of distance or kernel to use (default: "e-dist"). Options include:
#'        \itemize{
#'          \item "euclidean": Euclidean distance
#'          \item "gaussian": Gaussian kernel
#'          \item "laplacian": Laplacian kernel
#'          \item "polynomial": Polynomial kernel
#'          \item "e-dist": Euclidean-based aggregated distance (for high-dimensional data)
#'          \item "g-dist": Gaussian kernel-based aggregated distance
#'          \item "l-dist": Laplacian kernel-based aggregated distance
#'        }
#' @param bw Bandwidth parameter for kernel calculations. If NULL (default), it will be automatically determined.
#' @param expo Exponent parameter for distance calculation (default: 1).
#' @param scale_factor Scaling factor for automatic bandwidth calculation (default: 0.5).
#' @param group Optional vector specifying group membership for each variable/column.
#'        Used for group-wise distance calculations in "e-dist", "g-dist", or "l-dist".
#' @param B Number of permutations for p-value calculation (default: 199).
#' @param alpha Significance level for determining whether a change point is significant (default: 0.05).
#' @param seeds Random seed for reproducibility. If NULL, no seed is set.
#' @param method Testing method for p-value calculation: "permutation" or "asymptotic" (default: "permutation").
#' @param num_cores Number of cores for parallel processing of permutations (default: 1).
#'
#' @return A list containing:
#'   \item{locations}{The estimated change point location (index)}
#'   \item{pvalue}{P-value for the detected change point}
#'   \item{statistic}{The test statistic value at the detected change point}
#'   \item{cluster}{A vector of cluster assignments (1 for observations before the change point,
#'                  2 for observations after). If no significant change point is found (p-value > alpha),
#'                  all observations are assigned to cluster 1.}
#'
#' @examples
#' # Example 1: Change point in univariate time series
#' set.seed(123)
#' # Generate data with a change in mean at t=50
#' n <- 100
#' x1 <- rnorm(50, mean = 0, sd = 1)
#' x2 <- rnorm(50, mean = 2, sd = 1)
#' x <- c(x1, x2)
#' data <- matrix(x, ncol = 1)
#'
#' # Detect change point
#' result <- kcpd_single(data, method = "permutation", B = 99)
#' print(result$locations)  # Should be close to 50
#' print(result$pvalue)     # Should be < 0.05
#'
#' # Example 2: Change point in high-dimensional time series
#' set.seed(456)
#' # Generate high-dimensional data (p > n) with a change at t=60
#' n <- 120
#' p <- 150  # Number of variables exceeds sample size
#'
#' # First segment: independent variables
#' X1 <- matrix(rnorm(60 * p), ncol = p)
#'
#' # Second segment: shifted mean
#' X2 <- matrix(rnorm(60 * p, mean = 0.5), ncol = p)
#'
#' # Combined data
#' X <- rbind(X1, X2)
#'
#' # Detect change point using e-dist (efficient for high-dimensional data)
#' result <- kcpd_single(X, type = "e-dist", method = "permutation", B = 99)
#' print(result$locations)  # Should be close to 60
#' print(result$pvalue)     # Should be < 0.05
#' print(table(result$cluster))  # Distribution of cluster assignments
#'
#' @references
#' Chakraborty, S., & Zhang, X. (2021). A New Framework for Distance and Kernel-based Metrics
#' in High-dimension. \emph{Electronic Journal of Statistics}, 15, 5455-5522.
#'
#' Chakraborty, S., & Zhang, X. (2021). High-dimensional Change-point Detection Using
#' Generalized Homogeneity Metrics. arXiv:2105.08976.
#'
#' @seealso
#' \code{\link{kcpd_sbs}} for detecting multiple change points using seeded binary segmentation
#' \code{\link{KDist_matrix}} for details on distance/kernel matrix calculation
#' \code{\link{chsic}} for computing cross HSIC
#'
#' @importFrom parallel mclapply detectCores
#'
#' @export
kcpd_single <- function(data, type = "e-dist", bw = NULL, expo = 1, scale_factor = 0.5, group = NULL,
                          B = 199, alpha = 0.05, seeds = NULL,
                          method = "permutation", num_cores = 1) {
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Input data must be a matrix or data frame")
  }
  n <- nrow(data)

  # Set number of cores
  mc.cores <- detectCores() - 1

  # Calculate the distance matrix for the original data
  full_dist <- KDist_matrix(data = data, type = type, bw = bw, expo = expo, scale_factor = scale_factor, group = group)

  # Apply stat_recur function to get test statistic and change point
  result <- .stat_recur(full_dist, type = type)

  # Extract values from the unnamed list returned by stat_recur
  val_test_stat <- result[[1]]  # First element is the statistic value
  cp <- result[[2]]             # Second element is the change point

  # Calculate p-value based on specified method
  if (method == "permutation") {
    # Permutation test - use parallel processing
    if (!is.null(seeds)) {
      set.seed(seeds)
    }

    # Process permutations in parallel
    stat_vec <- unlist(mclapply(1:B, function(b) {
      # Generate permuted indices
      perm_indices <- sample(1:n)

      # Compute distance matrix for permuted data
      perm_dist <- full_dist[perm_indices, perm_indices]

      # Calculate test statistic for permuted data
      perm_result <- .stat_recur(perm_dist, type = type)

      # Return the test statistic (first element of the result list)
      return(perm_result[[1]])
    }, mc.cores = num_cores))

    # Calculate p-value
    pval <- (1 + sum(stat_vec > val_test_stat)) / (1 + B)

  } else if (method == "asymptotic") {
    # If method is asymptotic, use the package's asymptotic_stats
    stat_val <- asymptotic_stats  # Use the globally loaded data

    # Calculate p-value using asymptotic approximation
    pval <- (1 + sum(stat_val > val_test_stat)) / (1 + length(stat_val))

  } else {
    stop("Invalid method. Choose either 'permutation' or 'asymptotic'")
  }

  # Create and return results
  return(list(
    locations = cp,
    pvalue = pval,
    statistic = val_test_stat,
    cluster = if(pval < alpha) c(rep(1, cp), rep(2, n-cp)) else rep(1, n)
  ))
}

#' Seeded Binary Segmentation for Multiple Change Point Detection
#'
#' @description
#' Detects multiple change points in time series or sequential data using the Seeded Binary Segmentation (SeedBS)
#' method combined with kernel-based metrics. This approach is particularly effective for high-dimensional data
#' and can detect various types of changes in distribution beyond just mean shifts.
#'
#' @details
#' This function implements a multiple change point detection approach that combines:
#' 1. The Seeded Binary Segmentation (SeedBS) methodology from Kovács et al. (2023), which uses
#'    a deterministic set of intervals at multiple scales to search for change points
#' 2. The kernel-based testing framework from Chakraborty & Zhang (2021), which enables
#'    detection of general distributional changes even in high-dimensional settings
#'
#' The algorithm works by:
#' 1. Generating a structured set of intervals using the `get_seeded_intervals` function
#' 2. Applying the `kcpd_single` detection method to each interval independently
#' 3. Selecting significant change points in order of p-value (smallest first)
#' 4. Removing overlapping intervals to avoid redundant change points
#' 5. Returning the final set of change points with their p-values and segment assignments
#'
#' This method is more computationally efficient than exhaustive search approaches while
#' maintaining statistical power and accuracy, especially in high-dimensional settings.
#'
#' @param data A matrix or data frame with rows representing time points (observations) and
#'        columns representing variables or features.
#' @param type Type of distance or kernel to use (default: "e-dist"). Options include:
#'        \itemize{
#'          \item "euclidean": Euclidean distance
#'          \item "gaussian": Gaussian kernel
#'          \item "laplacian": Laplacian kernel
#'          \item "polynomial": Polynomial kernel
#'          \item "e-dist": Euclidean-based aggregated distance (for high-dimensional data)
#'          \item "g-dist": Gaussian kernel-based aggregated distance
#'          \item "l-dist": Laplacian kernel-based aggregated distance
#'        }
#' @param bw Bandwidth parameter for kernel calculations. If NULL (default), it will be automatically determined.
#' @param expo Exponent parameter for distance calculation (default: 1).
#' @param scale_factor Scaling factor for automatic bandwidth calculation (default: 0.5).
#' @param group Optional vector specifying group membership for each variable/column.
#'        Used for group-wise distance calculations in "e-dist", "g-dist", or "l-dist".
#' @param B Number of permutations for p-value calculation (default: 199).
#' @param alpha Significance level for determining whether a change point is significant (default: 0.05).
#' @param seeds Random seed for reproducibility. If NULL, no seed is set.
#' @param method Testing method for p-value calculation: "permutation" or "asymptotic" (default: "asymptotic").
#'        The asymptotic method is typically faster and recommended for larger datasets.
#' @param decay Decay factor for interval generation (default: sqrt(2)).
#'        Controls the density of intervals; smaller values create more intervals.
#' @param unique_int Logical; whether to use only unique intervals (default: FALSE).
#' @param bound Minimum segment size to consider (default: 20).
#'        Intervals shorter than this threshold are discarded. It should be greater than 5.
#' @param num_cores Number of cores for parallel processing (default: 1).
#'        If NULL, will use one less than the available number of cores.
#'
#' @return A list containing:
#'   \item{locations}{A vector of detected change point locations (indices)}
#'   \item{pvalues}{A vector of p-values corresponding to each detected change point}
#'   \item{cluster}{A vector of cluster assignments for each observation, where each segment
#'                  between change points is assigned a unique integer}
#'
#' @examples
#' # Example 1: High-dimensional data with change in distribution
#' set.seed(123)
#' n <- 100  # Total sample size
#' p <- 200  # Number of dimensions (high-dimensional setting)
#' cp <- 50  # Change point location
#'
#' # Generate high-dimensional data
#' # Before change: standard normal distribution
#' X1 <- matrix(rnorm(cp * p), nrow = cp, ncol = p)
#'
#' # After change: centered standard exponential distribution
#' X2 <- matrix(rexp(n = (n-cp) * p) - 1, nrow = (n-cp), ncol = p)
#'
#' # Combine data
#' X <- rbind(X1, X2)
#'
#' # Detect change points
#' result <- kcpd_sbs(X, type = "e-dist", method = "permutation")
#'
#' print(result$locations)  # Should be close to 50
#' print(result$pvalues)    # Should be < 0.05
#' print(table(result$cluster))  # Should show 2 segments
#'
#' # Visualize first few dimensions (for demonstration only)
#' matplot(1:n, X[, 1:5], type = "l", col = 1:5, lty = 1,
#'         xlab = "Time", ylab = "Value", main = "First 5 Dimensions")
#' if(length(result$locations) > 0) {
#'   abline(v = result$locations, col = "red", lwd = 2, lty = 2)
#' }
#'
#' # Example 2: Multiple change points in high-dimensional data
#' set.seed(456)
#' # Generate high-dimensional data with changes at t=40 and t=80
#' n <- 120
#' p <- 50  # High-dimensional
#'
#' # Three segments with different distributions
#' X1 <- matrix(rnorm(40 * p), ncol = p)
#' X2 <- matrix(rnorm(40 * p, mean = 0.5), ncol = p)
#'
#' # Third segment has different covariance structure
#' Sigma <- matrix(0.5, p, p)
#' diag(Sigma) <- 1
#' X3 <- matrix(rnorm(40 * p), ncol = p)
#' X3 <- X3 %*% chol(Sigma)  # Apply covariance structure
#'
#' # Combined data
#' X <- rbind(X1, X2, X3)
#'
#' # Detect multiple change points using e-dist (efficient for high-dimensional data)
#' result <- kcpd_sbs(X, type = "e-dist", method = "asymptotic", decay = 1.5)
#' print(result$locations)  # Should be close to 40 and 80
#' print(result$pvalues)
#' print(table(result$cluster))
#'
#' # Simply visualize the results with base R
#' mean_vals <- rowMeans(X)  # Just for visualization
#' plot(1:n, mean_vals, pch = 16, col = result$cluster,
#'      xlab = "Time", ylab = "Mean Value",
#'      main = "Multiple Change Point Detection")
#' # Add vertical lines at change points
#' if(length(result$locations) > 0) {
#'   abline(v = result$locations, lty = 2, col = "red")
#' }
#'
#' @references
#' Kovács, S., Bühlmann, P., Li, H., & Munk, A. (2023). Seeded binary segmentation: a
#' general methodology for fast and optimal changepoint detection. *Biometrika*, *110*(1), 249-256.
#'
#' Chakraborty, S., & Zhang, X. (2021). A New Framework for Distance and Kernel-based Metrics
#' in High-dimension. \emph{Electronic Journal of Statistics}, 15, 5455-5522.
#'
#' Chakraborty, S., & Zhang, X. (2021). High-dimensional Change-point Detection Using
#' Generalized Homogeneity Metrics. arXiv:2105.08976.
#'
#' @seealso
#' \code{\link{kcpd_single}} for single change point detection
#' \code{\link{get_seeded_intervals}} for generating the seeded intervals
#' \code{\link{adjustedRandIndex}} for evaluating change point detection results
#' \code{\link{KDist_matrix}} for details on distance/kernel matrix calculation
#'
#' @importFrom parallel mclapply detectCores
#' @export
kcpd_sbs <- function(data, type = "e-dist", bw = NULL, expo = 1, scale_factor = 0.5,
                     group = NULL, B = 199, alpha = 0.05, seeds = NULL,
                     method = "asymptotic", decay = sqrt(2), unique_int = TRUE, bound = 20, num_cores = 1) {

  # Check if bound is large enough
  if (bound <= 5) {
    stop("Parameter 'bound' must be greater than 5 to ensure numerical stability in change point detection.")
  }

  # Input validation
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Input data must be a matrix or data frame")
  }

  n <- nrow(data)

  # If method is asymptotic, use the package's asymptotic_stats
  if (method == "asymptotic") {
    stat_val <- asymptotic_stats  # Use the globally loaded data
  }

  # Get seeded intervals
  intervals <- get_seeded_intervals(n, decay = decay, unique_int = unique_int, bound = bound)

  # Initialize results matrix for all intervals
  interval_results <- matrix(NA, nrow = nrow(intervals), ncol = 4)
  colnames(interval_results) <- c("start", "end", "cp", "pvalue")

  # Process intervals in parallel using mclapply
  if(is.null(num_cores)) num_cores <- detectCores() - 1

  interval_results <- mclapply(1:nrow(intervals), function(i) {
    start_idx <- intervals[i, 1]
    end_idx <- intervals[i, 2]

    # Check if interval size meets minimum requirement
    if (end_idx - start_idx + 1 < bound) {
      return(c(start_idx, end_idx, NA, 1))
    }

    # Extract segment data
    segment_data <- data[start_idx:end_idx, , drop = FALSE]

    # Apply kcpd_single to the segment
    result <- tryCatch({
      kcpd_single(
        data = segment_data,
        type = type,
        bw = bw,
        expo = expo,
        scale_factor = scale_factor,
        group = group,
        B = B,
        alpha = alpha,
        seeds = seeds,
        method = method,
        num_cores = num_cores
      )
    }, error = function(e) {
      list(pvalue = 1, locations = NA, statistic = 0)
    })

    # Compute global change point index
    global_cp <- ifelse(is.na(result$locations), NA, start_idx + result$locations - 1)

    # Return results for this interval
    return(c(start_idx, end_idx, global_cp, result$pvalue))
  }, mc.cores = num_cores)

  # Convert list to matrix
  interval_results <- do.call(rbind, interval_results)
  colnames(interval_results) <- c("start", "end", "cp", "pvalue")

  # Initialize results vectors
  change_points <- numeric(0)
  p_values <- numeric(0)

  # Find all significant intervals
  significant <- interval_results[interval_results[, "pvalue"] < alpha & !is.na(interval_results[, "cp"]), , drop = FALSE]

  # If no significant change points, return empty result
  if (nrow(significant) == 0) {
    return(list(
      locations = numeric(0),
      pvalues = numeric(0),
      cluster = rep(1, n)
    ))
  }

  # Process significant change points in order of significance (lowest p-value first)
  while (nrow(significant) > 0) {
    # Get interval with smallest p-value
    best_idx <- which.min(significant[, "pvalue"])
    best_cp <- significant[best_idx, "cp"]
    best_pvalue <- significant[best_idx, "pvalue"]

    # Add to results
    change_points <- c(change_points, best_cp)
    p_values <- c(p_values, best_pvalue)

    # Remove all intervals containing this change point
    to_remove <- apply(significant, 1, function(interval) {
      best_cp >= interval[1] && best_cp <= interval[2]
    })

    significant <- significant[!to_remove, , drop = FALSE]
  }

  # Sort change points
  if (length(change_points) > 0) {
    # Sort the change points
    sorted_indices <- order(change_points)
    change_points <- change_points[sorted_indices]
    p_values <- p_values[sorted_indices]

    # Create cluster labels
    clusters <- numeric(n)
    prev_cp <- 0

    for (i in 1:(length(change_points) + 1)) {
      current_cp <- if (i <= length(change_points)) change_points[i] else n
      clusters[(prev_cp + 1):current_cp] <- i
      prev_cp <- current_cp
    }
  } else {
    # No change points found
    clusters <- rep(1, n)
  }

  # Return results
  return(list(
    locations = change_points,
    pvalues = p_values,
    cluster = clusters
  ))
}

#' Generate Seeded Intervals for Binary Segmentation
#'
#' @description
#' Creates a structured set of intervals for use in seeded binary segmentation (SeedBS)
#' for change point detection, based on the methodology described in Kovács et al. (2023).
#'
#' @details
#' This function implements the interval generation scheme for Seeded Binary Segmentation (SeedBS),
#' a general methodology for fast and optimal changepoint detection. The approach generates a
#' carefully chosen set of intervals at different scales, which are used to search for change points.
#'
#' The intervals are created using a dyadic structure controlled by the decay parameter, starting
#' with the full interval \[1,n\] and then generating intervals of decreasing lengths. This creates
#' a multi-scale representation that enables efficient change point detection while maintaining
#' statistical power.
#'
#' SeedBS offers computational advantages over traditional approaches like Wild Binary Segmentation
#' (WBS) by using a deterministic set of intervals rather than random ones, while still maintaining
#' optimal statistical properties.
#'
#' @param n Length of time series (number of observations)
#' @param decay Decay factor that controls the interval generation (default: sqrt(2)).
#'        Smaller values create more intervals and finer granularity.
#' @param unique_int Logical; whether to return only unique intervals (default: TRUE).
#'        When TRUE, duplicate intervals are removed.
#' @param bound Minimum interval size to consider (default: 2).
#'        Intervals shorter than this threshold are discarded.
#'
#' @return A matrix with two columns representing the start and end points of each interval.
#'         Each row defines one interval \[start, end\] to be considered in the change point detection.
#'
#' @examples
#' # Example 1: Generate intervals for a time series of length 100
#' intervals <- get_seeded_intervals(100)
#' head(intervals)  # View the first few intervals
#' dim(intervals)   # See how many intervals were generated
#'
#' # Example 2: Generate intervals with a smaller decay factor (more intervals)
#' intervals_fine <- get_seeded_intervals(100, decay = 1.2)
#' dim(intervals_fine)  # Should generate more intervals than the default
#'
#' # Example 3: Generate only unique intervals with minimum size 10
#' intervals_unique <- get_seeded_intervals(200, unique_int = TRUE, bound = 10)
#' head(intervals_unique)
#'
#' # Example 4: Visualize the interval structure
#' intervals_small <- get_seeded_intervals(50, decay = 1.5)
#' plot(1, 1, type = "n", xlim = c(1, 50), ylim = c(1, nrow(intervals_small)),
#'      xlab = "Time points", ylab = "Interval index", main = "Seeded Intervals")
#' for (i in 1:nrow(intervals_small)) {
#'   lines(c(intervals_small[i, 1], intervals_small[i, 2]), c(i, i), lwd = 2)
#' }
#'
#' @references
#' Kovács, S., Bühlmann, P., Li, H., & Munk, A. (2023). Seeded binary segmentation: a
#' general methodology for fast and optimal changepoint detection. *Biometrika*, *110*(1), 249-256.
#'
#' @seealso
#' \code{\link{kcpd_sbs}} for detecting multiple change points using these seeded intervals
#'
#' @export
get_seeded_intervals <- function(n, decay = sqrt(2), unique_int = TRUE, bound = 2) {
  # Convert n to integer
  n <- as.integer(n)

  # Calculate the maximum depth based on decay rate
  depth <- log(n, base = decay)
  depth <- ceiling(depth)

  # Initialize the boundary matrix with the full interval
  boundary_mtx <- matrix(NA, ncol = 2)
  colnames(boundary_mtx) <- c("start", "end")
  boundary_mtx[1, ] <- c(1, n)

  # Create intervals at each depth level
  for (i in 2:depth) {
    int_length <- n * (1 / decay)^(i - 1)

    # Skip this depth level if intervals would be smaller than bound
    if (int_length < bound) {
      break
    }

    n_int <- ceiling(round(n / int_length, 14)) * 2 - 1

    boundary_mtx <- rbind(
      boundary_mtx,
      cbind(
        floor(seq(1, n - int_length + 1, length.out = (n_int))),
        ceiling(seq(int_length, n, length.out = (n_int)))
      )
    )
  }

  # Filter out intervals shorter than bound
  # Calculate interval lengths
  interval_lengths <- boundary_mtx[, "end"] - boundary_mtx[, "start"] + 1

  # Keep only intervals with length >= bound
  boundary_mtx <- boundary_mtx[interval_lengths >= bound, ]

  # Return unique intervals if required
  if (unique_int) {
    return(unique(boundary_mtx))
  }

  boundary_mtx
}

#' Calculate Statistics Recursively for Change Point Detection
#'
#' @description
#' Internal function that calculates test statistics for various potential change point
#' locations. This implements the kernel-based test statistic calculation as described in
#' Chakraborty & Zhang (2021).
#'
#' @details
#' This function computes the test statistics for all possible change point locations by:
#' 1. Calculating MMD (Maximum Mean Discrepancy) values
#' 2. Computing HSIC (Hilbert-Schmidt Independence Criterion) components
#' 3. Combining these to form the final test statistic
#'
#' The implementation uses efficient recursion strategies and C++ implementations for
#' computationally intensive parts.
#'
#' @param full_dist Full distance/kernel matrix computed from the input data
#' @param type Type of distance or kernel ("euclidean", "gaussian", etc.)
#'
#' @return A list containing:
#'   \item{statistic}{Maximum test statistic value}
#'   \item{change_point}{Estimated change point location}
#'   \item{all_statistics}{Vector of test statistics for all potential change points}
#'
#' @keywords internal
.stat_recur <- function(full_dist, type) {
  # Get the dimensions of the distance matrix
  n <- nrow(full_dist)

  # Create index vectors for partitioning
  # Ensuring we have at least 4 elements on each side
  n1 <- 4:(n-4)
  n2 <- n - n1

  # Calculate all necessary statistics
  mmd_all <- .mmd_recur(full_dist, type = type)
  chsic_all <- chsic_recur_cpp(full_dist)
  hsic_x_all <- hsic_recur_cpp(full_dist)

  # Calculate hsic_y by reversing the matrix and then reversing the result
  reversed_dist <- full_dist[n:1, n:1]
  hsic_y_all <- hsic_recur_cpp(reversed_dist)
  hsic_y_all <- hsic_y_all[length(hsic_y_all):1]

  # Calculate denominator scalars
  scalar_den <- 1/(n1*n2) + 1/(2*n1*(n1-1)) + 1/(2*n2*(n2-1))

  # Calculate the Snm statistic components
  n1_factor <- n1 * (n1 - 3) / 2
  n2_factor <- n2 * (n2 - 3) / 2

  # Numerator of Snm
  num_Snm <- 4 * (n1 - 1) * (n2 - 1) * chsic_all +
    4 * n1_factor * hsic_x_all +
    4 * n2_factor * hsic_y_all

  # Denominator of Snm
  den_Snm <- (n1 - 1) * (n2 - 1) + n1_factor + n2_factor

  # Compute Snm
  Snm <- num_Snm / den_Snm

  # Compute denominator and final statistic
  Den <- sqrt(scalar_den * Snm)
  stat <- mmd_all * n1 * n2/n^2 / Den

  # Find the index with maximum statistic value
  max_idx <- which.max(stat)

  # Calculate the estimated change point
  cp_est <- 3 + max_idx
  stat_val <- stat[max_idx]

  # Return results as a list
  return(list(
    statistic = stat_val,
    change_point = cp_est,
    all_statistics = stat
  ))
}

#' Calculate MMD Recursively for Change Point Detection
#'
#' @description
#' Internal function that calculates Maximum Mean Discrepancy (MMD) values recursively
#' for all potential change point locations.
#'
#' @details
#' This function computes MMD values by calculating partial averages of kernel/distance
#' matrix elements using efficient C++ implementations. The computation differs based
#' on whether the input is a distance-based or kernel-based metric.
#'
#' @param full_dist Full distance/kernel matrix computed from the input data
#' @param type Type of distance or kernel ("euclidean", "gaussian", etc.)
#'
#' @return A vector of MMD values for all potential change points
#'
#' @keywords internal
.mmd_recur <- function(full_dist, type) {
  n <- nrow(full_dist)
  # Calculate the partial sums using C++ implementations
  Axx <- calculate_partial_aves_cpp(full_dist)
  Ayy <- calculate_partial_aves_cpp(full_dist[n:1,n:1])
  Ayy <- Ayy[length(Ayy):1]
  Axy <- calculate_partial_aves_cross_cpp(full_dist)

  # Calculate Generalized Energy Distance based on distance type
  if (type %in% c("euclidean", "e-dist", "l-dist", "g-dist")) {
    # For distance-based types, MMD = 2*ave_Axy - ave_Axx - ave_Ayy
    ed <- 2 * Axy - Axx - Ayy
  } else {
    # For kernel-based types, MMD = ave_Axx + ave_Ayy - 2*ave_Axy
    ed <- Axx + Ayy - 2 * Axy
  }
  return(ed)
}

#' Adjusted Rand Index for Comparing Clusterings
#'
#' @description
#' Calculates the Adjusted Rand Index (ARI) between two clusterings or partitions of the same data.
#' The ARI is a measure of the similarity between two data clusterings, adjusted for chance.
#'
#' @details
#' The Adjusted Rand Index is a modified version of the Rand Index that is adjusted for the
#' chance groupings of elements. It has an expected value of 0 for random partitions and a
#' maximum value of 1 for identical clusterings.
#'
#' The formula for the ARI is:
#'
#' \deqn{ARI = \frac{\sum_{ij} \binom{n_{ij}}{2} - \big[\sum_i \binom{a_i}{2} \sum_j \binom{b_j}{2} \big] / \binom{n}{2}}
#' {\frac{1}{2}\big[\sum_i \binom{a_i}{2} + \sum_j \binom{b_j}{2}\big] - \big[\sum_i \binom{a_i}{2} \sum_j \binom{b_j}{2}\big] / \binom{n}{2}}}
#'
#' where:
#' - \eqn{n_{ij}} is the number of objects in both class \eqn{i} of clustering \eqn{x} and class \eqn{j} of clustering \eqn{y}
#' - \eqn{a_i} is the number of objects in class \eqn{i} of clustering \eqn{x}
#' - \eqn{b_j} is the number of objects in class \eqn{j} of clustering \eqn{y}
#' - \eqn{n} is the total number of objects
#'
#' The ARI is especially useful for evaluating change point detection algorithms by comparing the
#' detected segmentation with a ground truth segmentation.
#'
#' @param x First clustering (vector of cluster assignments)
#' @param y Second clustering (vector of cluster assignments)
#'
#' @return A numeric value between -1 and 1:
#'   \itemize{
#'     \item 1: Perfect agreement between the two clusterings
#'     \item ~0: Agreement equivalent to random chance
#'     \item Negative values: Agreement less than random chance
#'   }
#'
#' @examples
#' # Example 1: Identical clusterings
#' x <- c(1, 1, 1, 2, 2, 3, 3, 3)
#' y <- c(1, 1, 1, 2, 2, 3, 3, 3)
#' adjustedRandIndex(x, y)  # Should return 1
#'
#' # Example 2: Different label names but same structure
#' x <- c(1, 1, 1, 2, 2, 3, 3, 3)
#' y <- c(5, 5, 5, 9, 9, 2, 2, 2)
#' adjustedRandIndex(x, y)  # Should return 1
#'
#' # Example 3: Similar but not identical clusterings
#' x <- c(1, 1, 1, 2, 2, 3, 3, 3)
#' y <- c(1, 1, 2, 2, 2, 3, 3, 3)
#' adjustedRandIndex(x, y)  # Should return a value between 0 and 1
#'
#' # Example 4: Very different clusterings
#' x <- c(1, 1, 1, 1, 2, 2, 2, 2)
#' y <- c(1, 2, 1, 2, 1, 2, 1, 2)
#' adjustedRandIndex(x, y)  # Should return a value close to 0
#'
#' # Example 5: Comparing change point detection results with ground truth
#' true_segments <- c(rep(1, 50), rep(2, 30), rep(3, 40))
#' detected_segments <- c(rep(1, 45), rep(2, 35), rep(3, 40))
#' ari <- adjustedRandIndex(true_segments, detected_segments)
#' print(paste("ARI between true and detected segments:", ari))
#'
#' @references
#' Hubert, L., & Arabie, P. (1985). Comparing partitions. Journal of Classification, 2(1), 193-218.
#'
#' @seealso
#' \code{\link{kcpd_single}} and \code{\link{kcpd_sbs}} for change point detection methods
#' that can be evaluated using the Adjusted Rand Index
#'
#' @export
adjustedRandIndex <- function(x, y)
{
  x <- as.vector(x)
  y <- as.vector(y)
  if (length(x) != length(y))
    stop("arguments must be vectors of the same length")
  tab <- table(x, y)
  if (all(dim(tab) == c(1, 1)))
    return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}
