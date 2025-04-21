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
#' @export
bw_optim <- function(x, sample_limit = 1000, scale_factor = 0.5, group = NULL) {
  # Input validation
  if (!is.numeric(x)) {
    stop("Input must be numeric")
  }
  
  # Ensure x is a matrix
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
  } else if (is.data.frame(x)) {
    x <- as.matrix(x)
  } else if (!is.matrix(x)) {
    stop("Input must be a vector, matrix, or data.frame")
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

#' Calculate Maximum Mean Discrepancy (MMD)
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
#'
#' @export
mmd <- function(x, y = NULL, type = "euclidean", bw = NULL, expo = 1, scale_factor = 0.5, group = NULL, u_center = FALSE, n = NULL, m = NULL) {
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
#' Calculates the HSIC between two datasets or matrices.
#'
#' @param x First dataset or distance matrix
#' @param y Second dataset or distance matrix
#' @param type Type of kernel or distance
#' @param bw Bandwidth parameter
#' @param expo Exponent parameter
#' @param scale_factor Scaling factor for bandwidth
#' @param group Optional grouping of variables
#' @param u_center Logical; use U-centering instead of V-centering
#' @param is_distance Logical; whether input matrices are already distance matrices
#'
#' @return HSIC value
#'
#' @export
hsic <- function(x, y, type = "euclidean", bw = NULL, expo = 1, scale_factor = 0.5, group = NULL, u_center = FALSE, is_distance = FALSE) {
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
#' Calculates the cross HSIC between two samples.
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
#' @export
chsic <- function(x, y = NULL, type = "euclidean", bw = NULL, expo = 1, scale_factor = 0.5, group = NULL) {
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
#' @param x First dataset or distance matrix
#' @param y Second dataset or distance matrix 
#' @param type Type of kernel or distance
#' @param bw Bandwidth parameter
#' @param expo Exponent parameter
#' @param scale_factor Scaling factor for bandwidth
#' @param group Optional grouping of variables
#' @param u_center Logical; use U-centering instead of V-centering
#' @param is_distance Logical; whether input matrices are already distance matrices
#'
#' @return HSIC correlation coefficient between -1 and 1
#'
#' @export
hsic_cor <- function(x, y, type = "euclidean", bw = NULL, expo = 1, scale_factor = 0.5, group = NULL, u_center = FALSE, is_distance = FALSE) {
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
    # Compute distance matrices
    Dxx <- KDist_matrix(x, type, bw, expo, scale_factor, group)
    Dyy <- KDist_matrix(y, type, bw, expo, scale_factor, group)
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
