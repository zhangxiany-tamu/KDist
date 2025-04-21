#' Joint Hilbert-Schmidt Independence Criterion
#'
#' Calculates the joint HSIC for multiple datasets.
#'
#' @param x List of matrices or data frames
#' @param cc Constant parameter, default is 1
#' @param type Type of kernel or distance
#' @param stat_type Type of statistic ("V", "U", "US", or "UR")
#' @param bw Bandwidth parameter
#' @param expo Exponent parameter
#' @param scale_factor Scaling factor for bandwidth
#' @param group Optional grouping of variables
#'
#' @return Joint HSIC value
#'
#' @export
jhsic <- function(x, cc = 1, type = "euclidean", stat_type = "V", bw = NULL, expo = 1,
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
#' Calculate the dHSIC for multiple datasets.
#'
#' @param x List of matrices or data frames
#' @param type Type of kernel or distance
#' @param bw Bandwidth parameter
#' @param expo Exponent parameter
#' @param scale_factor Scaling factor for bandwidth
#' @param group Optional grouping of variables
#'
#' @return dHSIC value
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

#' A Parallel Implementation of HSIC for Mutual Independence Testing in High-dimension
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
#' @export
mhsic_parallel <- function(x, type = "euclidean", bw = NULL, expo = 1, scale_factor = 0.5, cores = NULL) {
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

#' HSIC for Mutual Independence Testing in High-dimension
#'
#' @param x Matrix where each column is a variable to test for mutual independence
#' @param type Type of kernel or distance
#' @param bw Bandwidth parameter
#' @param expo Exponent parameter
#' @param scale_factor Scaling factor for bandwidth
#' @param cores Number of cores for parallel processing
#'
#' @return A list with test statistic and p-value
#'
#' @export
mhsic <- function(x, type = "euclidean", bw = NULL, expo = 1, scale_factor = 0.5, cores = NULL) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame")
  }

  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  p <- dim(x)[2]

  # Use parallel processing only for larger matrices
  if (p > 60) {
    return(mhsic_parallel(x, type, bw, expo, scale_factor, cores))
  } else {
    return(mhsic_cpp(x, type, bw, expo, scale_factor))
  }
}
