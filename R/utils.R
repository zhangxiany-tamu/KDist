#' V-center a Matrix
#'
#' Performs V-centering on a matrix by subtracting row means, column means,
#' and adding the overall mean. V-centering is a technique used in various
#' statistical methods including distance covariance, HSIC, and kernel-based
#' independence tests.
#'
#' @param A A numeric matrix to be V-centered
#'
#' @return A V-centered matrix of the same dimensions as the input
#'
#' @details
#' For a matrix A, the V-centered matrix is calculated as:
#'
#' V(A) = A - row_means - col_means + total_mean
#'
#' where row_means is a matrix where each element in row i is the mean of row i in A,
#' col_means is a matrix where each element in column j is the mean of column j in A,
#' and total_mean is the overall mean of all elements in A.
#'
#' This function uses an efficient Rcpp implementation for improved performance,
#' which is particularly beneficial for large matrices.
#'
#' @examples
#' # Create a sample matrix
#' set.seed(123)
#' A <- matrix(rnorm(100), nrow = 10)
#'
#' # Apply V-centering
#' A_centered <- v_center(A)
#'
#' # Verify properties of V-centered matrix
#' # Row sums should be approximately zero
#' row_sums <- rowSums(A_centered)
#' print(round(row_sums, 10))
#'
#' # Column sums should be approximately zero
#' col_sums <- colSums(A_centered)
#' print(round(col_sums, 10))
#'
#' # Total sum should be approximately zero
#' total_sum <- sum(A_centered)
#' print(round(total_sum, 10))
#'
#' @seealso
#' \code{\link{hsic}} for Hilbert-Schmidt Independence Criterion
#' \code{\link{dcov}} for distance covariance
#'
#' @export
v_center <- function(A) {
  # Check if input is a matrix
  if (!is.matrix(A)) {
    stop("Input must be a matrix")
  }

  # Check if matrix is numeric
  if (!is.numeric(A)) {
    stop("Input matrix must contain numeric values")
  }

  # Call the Rcpp function
  return(matrix_v_center(A))
}

#' U-center a Symmetric Matrix
#'
#' Performs U-centering on a symmetric matrix, which is a technique that provides
#' an unbiased estimator by removing diagonal elements and applying appropriate
#' centering. U-centering is commonly used in statistical tests based on distance
#' or kernel matrices, such as HSIC and distance covariance.
#'
#' @param A A symmetric numeric matrix to be U-centered
#'
#' @return A U-centered matrix of the same dimensions as the input, with zeros on the diagonal
#'
#' @details
#' For a symmetric matrix A, the U-centered matrix is calculated as:
#'
#' U(A) = A + (total sum)/((n-1)(n-2)) - (row sum)/(n-2) - (col sum)/(n-2)
#'
#' where row sum and column sum are identical since A is symmetric, total sum is the
#' sum of all elements in A, and n is the size of the matrix. The diagonal elements
#' are set to zero.
#'
#' U-centering provides an unbiased estimator for kernel and distance-based measures,
#' which is particularly important for small sample sizes. This implementation uses Rcpp
#' for improved computational efficiency with large matrices.
#'
#' @examples
#' # Create sample data
#' set.seed(123)
#' X <- matrix(rnorm(100), nrow = 20)  # 20 observations in 5D
#'
#' # Generate a Gaussian kernel matrix using KDist_matrix
#' K <- KDist_matrix(X)
#'
#' # Apply U-centering to the kernel matrix
#' K_u_centered <- u_center(K)
#'
#' # Verify properties of U-centered matrix
#' # Diagonal should be zero
#' diag_elements <- diag(K_u_centered)
#' print(all(diag_elements == 0))
#'
#' # Row and column sums (off-diagonal elements) should be approximately zero
#' row_sums <- rowSums(K_u_centered)
#' print(all(abs(row_sums) < 1e-10))
#'
#' @seealso
#' \code{\link{v_center}} for V-centering
#' \code{\link{hsic}} for Hilbert-Schmidt Independence Criterion
#' \code{\link{dcov}} for distance covariance
#'
#' @export
u_center <- function(A) {
  # Check if input is a matrix
  if (!is.matrix(A)) {
    stop("Input must be a matrix")
  }

  # Check if matrix is numeric
  if (!is.numeric(A)) {
    stop("Input matrix must contain numeric values")
  }

  # Check if matrix is square
  if (nrow(A) != ncol(A)) {
    stop("Input matrix must be square (symmetric)")
  }

  # Check for symmetry (with some tolerance for floating-point operations)
  if (!isSymmetric(A, tol = 1e-8)) {
    warning("Input matrix should be symmetric; results may be unexpected")
  }

  # Call the Rcpp function
  return(matrix_u_center(A))
}

#' Print a summary of a mmd_test object
#'
#' @param x A mmd_test object
#' @param ... Additional arguments to be passed to methods
#'
#' @return Invisibly returns the object
#' @export
print.mmd_test <- function(x, ...) {
  cat("Maximum Mean Discrepancy Test\n\n")
  cat("Statistic:", format(x$statistic, digits = 4), "\n")
  cat("P-value:", format(x$p.value, digits = 4), "\n")
  cat("\n")
  cat("Based on", length(x$permutation_values), "permutations\n")
  invisible(x)
}

#' Print a summary of a hsic_test object
#'
#' @param x A hsic_test object
#' @param ... Additional arguments to be passed to methods
#'
#' @return Invisibly returns the object
#' @export
print.hsic_test <- function(x, ...) {
  cat("Hilbert-Schmidt Independence Criterion Test\n\n")
  cat("Statistic:", format(x$statistic, digits = 4), "\n")
  cat("P-value:", format(x$p.value, digits = 4), "\n")
  cat("\n")
  cat("Based on", length(x$permutation_values), "permutations\n")
  invisible(x)
}

#' Print a summary of a dhsic_test object
#'
#' @param x A dhsic_test object
#' @param ... Additional arguments to be passed to methods
#'
#' @return Invisibly returns the object
#' @export
print.dhsic_test <- function(x, ...) {
  cat("Distance Hilbert-Schmidt Independence Criterion Test\n\n")
  cat("Statistic:", format(x$statistic, digits = 4), "\n")
  cat("P-value:", format(x$p.value, digits = 4), "\n")
  cat("\n")
  cat("Based on", length(x$permutation_values), "permutations\n")
  invisible(x)
}

#' Print a summary of a jhsic_test object
#'
#' @param x A jhsic_test object
#' @param ... Additional arguments to be passed to methods
#'
#' @return Invisibly returns the object
#' @export
print.jhsic_test <- function(x, ...) {
  cat("Joint Hilbert-Schmidt Independence Criterion Test\n\n")
  cat("Statistic type:", x$stat_description, "\n")
  cat("Statistic:", format(x$statistic, digits = 4), "\n")
  cat("P-value:", format(x$p.value, digits = 4), "\n")
  cat("\n")
  cat("Based on", length(x$permutation_values), "permutations\n")
  invisible(x)
}

#' Plot a histogram of permutation values for independence tests
#'
#' @param x A test object (mmd_test, hsic_test, dhsic_test, or jhsic_test)
#' @param ... Additional arguments to be passed to methods
#'
#' @return Invisibly returns the plot object
#' @export
plot.mmd_test <- function(x, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  hist(x$permutation_values,
       main = "Histogram of Permutation Values",
       xlab = "MMD Statistic",
       freq = FALSE,
       ...)

  # Draw a vertical line for the observed statistic
  abline(v = x$statistic, col = "red", lwd = 2)

  # Add a legend
  legend("topright",
         legend = c(paste("Observed =", format(x$statistic, digits = 4)),
                    paste("P-value =", format(x$p.value, digits = 4))),
         bty = "n")

  invisible(x)
}

#' Plot method for hsic_test
#'
#' @param x A hsic_test object
#' @param ... Additional arguments to be passed to methods
#'
#' @return Invisibly returns the plot object
#' @export
plot.hsic_test <- function(x, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  hist(x$permutation_values,
       main = "Histogram of Permutation Values",
       xlab = "HSIC Statistic",
       freq = FALSE,
       ...)

  # Draw a vertical line for the observed statistic
  abline(v = x$statistic, col = "red", lwd = 2)

  # Add a legend
  legend("topright",
         legend = c(paste("Observed =", format(x$statistic, digits = 4)),
                    paste("P-value =", format(x$p.value, digits = 4))),
         bty = "n")

  invisible(x)
}

#' Plot method for dhsic_test
#'
#' @param x A dhsic_test object
#' @param ... Additional arguments to be passed to methods
#'
#' @return Invisibly returns the plot object
#' @export
plot.dhsic_test <- function(x, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  hist(x$permutation_values,
       main = "Histogram of Permutation Values",
       xlab = "dHSIC Statistic",
       freq = FALSE,
       ...)

  # Draw a vertical line for the observed statistic
  abline(v = x$statistic, col = "red", lwd = 2)

  # Add a legend
  legend("topright",
         legend = c(paste("Observed =", format(x$statistic, digits = 4)),
                    paste("P-value =", format(x$p.value, digits = 4))),
         bty = "n")

  invisible(x)
}

#' Plot method for jhsic_test
#'
#' @param x A jhsic_test object
#' @param ... Additional arguments to be passed to methods
#'
#' @return Invisibly returns the plot object
#' @export
plot.jhsic_test <- function(x, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  hist(x$permutation_values,
       main = "Histogram of Permutation Values",
       xlab = "JHSIC Statistic",
       freq = FALSE,
       ...)

  # Draw a vertical line for the observed statistic
  abline(v = x$statistic, col = "red", lwd = 2)

  # Add a legend
  legend("topright",
         legend = c(paste("Test type:", x$stat_description),
                    paste("Observed =", format(x$statistic, digits = 4)),
                    paste("P-value =", format(x$p.value, digits = 4))),
         bty = "n")

  invisible(x)
}

#' Print a summary of a ed_test object
#'
#' @param x An ed_test object
#' @param ... Additional arguments to be passed to methods
#'
#' @return Invisibly returns the object
#' @export
print.ed_test <- function(x, ...) {
  cat("Energy Distance Test\n\n")
  cat("Statistic:", format(x$statistic, digits = 4), "\n")
  cat("P-value:", format(x$p.value, digits = 4), "\n")
  cat("\n")
  cat("Based on", length(x$permutation_values), "permutations\n")
  invisible(x)
}

#' Print a summary of a dcov_test object
#'
#' @param x A dcov_test object
#' @param ... Additional arguments to be passed to methods
#'
#' @return Invisibly returns the object
#' @export
print.dcov_test <- function(x, ...) {
  cat("Distance Covariance Test\n\n")
  cat("Statistic:", format(x$statistic, digits = 4), "\n")
  cat("P-value:", format(x$p.value, digits = 4), "\n")
  cat("\n")
  cat("Based on", length(x$permutation_values), "permutations\n")
  invisible(x)
}

#' Plot a histogram of permutation values for a ed_test
#'
#' @param x An ed_test object
#' @param ... Additional arguments to be passed to methods
#'
#' @return Invisibly returns the plot object
#' @export
plot.ed_test <- function(x, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  hist(x$permutation_values,
       main = "Histogram of Permutation Values",
       xlab = "Energy Distance Statistic",
       freq = FALSE,
       ...)

  # Draw a vertical line for the observed statistic
  abline(v = x$statistic, col = "red", lwd = 2)

  # Add a legend
  legend("topright",
         legend = c(paste("Observed =", format(x$statistic, digits = 4)),
                    paste("P-value =", format(x$p.value, digits = 4))),
         bty = "n")

  invisible(x)
}

#' Plot a histogram of permutation values for a dcov_test
#'
#' @param x A dcov_test object
#' @param ... Additional arguments to be passed to methods
#'
#' @return Invisibly returns the plot object
#' @export
plot.dcov_test <- function(x, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  hist(x$permutation_values,
       main = "Histogram of Permutation Values",
       xlab = "Distance Covariance Statistic",
       freq = FALSE,
       ...)

  # Draw a vertical line for the observed statistic
  abline(v = x$statistic, col = "red", lwd = 2)

  # Add a legend
  legend("topright",
         legend = c(paste("Observed =", format(x$statistic, digits = 4)),
                    paste("P-value =", format(x$p.value, digits = 4))),
         bty = "n")

  invisible(x)
}
