#' Maximum Mean Discrepancy (MMD) Test
#'
#' Performs a permutation test based on the Maximum Mean Discrepancy.
#'
#' @param x First sample
#' @param y Second sample
#' @param type Type of kernel or distance
#' @param bw Bandwidth parameter
#' @param expo Exponent parameter
#' @param scale_factor Scaling factor for bandwidth
#' @param group Optional grouping of variables
#' @param u_center Logical; use U-centering instead of V-centering
#' @param n_perm Number of permutations
#' @param seed Random seed for reproducibility
#' @param cores Number of cores for parallel computing
#'
#' @return A list with test statistic, p-value, and permutation values
#'
#' @export
mmd_test <- function(x, y, type = "euclidean", bw = NULL, expo = 1, 
                     scale_factor = 0.5, group = NULL, u_center = FALSE, 
                     n_perm = 1000, seed = NULL, cores = parallel::detectCores() - 1) {
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
  observed_mmd <- mmd(x = D, u_center = u_center, n = n, m = m)
  
  # Generate all permutation indices in advance for reproducibility
  all_perm_indices <- lapply(1:n_perm, function(i) {
    sample(total_n, total_n, replace = FALSE)
  })
  
  # Define the function to run for each permutation
  perm_function <- function(perm_indices) {
    # Reorder the distance matrix according to the permutation
    D_perm <- D[perm_indices, perm_indices]
    
    # Calculate MMD for the permuted data
    return(mmd(x = D_perm, u_center = u_center, n = n, m = m))
  }
  
  # Run permutations in parallel
  perm_results <- parallel::mclapply(all_perm_indices, perm_function, mc.cores = cores)
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

#' HSIC Independence Test
#'
#' Performs a permutation test based on the Hilbert-Schmidt Independence Criterion.
#'
#' @param x First dataset
#' @param y Second dataset
#' @param type Type of kernel or distance
#' @param bw Bandwidth parameter
#' @param expo Exponent parameter
#' @param scale_factor Scaling factor for bandwidth
#' @param group Optional grouping of variables
#' @param u_center Logical; use U-centering instead of V-centering
#' @param n_perm Number of permutations
#' @param seed Random seed for reproducibility
#' @param cores Number of cores for parallel computing
#'
#' @return A list with test statistic, p-value, and permutation values
#'
#' @export
hsic_test <- function(x, y, type = "euclidean", bw = NULL, expo = 1, 
                      scale_factor = 0.5, group = NULL, u_center = FALSE, 
                      n_perm = 1000, seed = NULL, cores = parallel::detectCores() - 1) {
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
  
  # Ensure dimensions of x and y match in the first dimension (number of observations)
  if (nrow(x) != nrow(y)) {
    stop("Number of observations in x and y must be the same")
  }
  
  # Calculate distance matrices once
  Dx <- KDist_matrix(x, type = type, bw = bw, expo = expo, scale_factor = scale_factor, group = group)
  Dy <- KDist_matrix(y, type = type, bw = bw, expo = expo, scale_factor = scale_factor, group = group)
  
  # Calculate observed test statistic
  observed_hsic <- hsic(Dx, Dy, u_center = u_center, is_distance = TRUE)
  
  # Define the function to run for each permutation
  perm_function <- function(i) {
    # Generate a random permutation of indices
    perm_indices <- sample(n, n, replace = FALSE)
    
    # Reorder the second distance matrix according to the permutation
    Dy_perm <- Dy[perm_indices, perm_indices]
    
    # Calculate HSIC for the permuted data
    return(hsic(Dx, Dy_perm, u_center = u_center, is_distance = TRUE))
  }
  
  # Run permutations in parallel
  perm_results <- parallel::mclapply(1:n_perm, function(i) perm_function(i), mc.cores = cores)
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

#' d-variable HSIC Independence Test (dHSIC)
#'
#' Performs a permutation test based on the d-variable Hilbert-Schmidt Independence Criterion with multiple datasets.
#'
#' @param x List of matrices or vectors
#' @param type Type of kernel or distance
#' @param bw Bandwidth parameter
#' @param expo Exponent parameter
#' @param scale_factor Scaling factor for bandwidth
#' @param group Optional grouping of variables
#' @param n_perm Number of permutations
#' @param seed Random seed for reproducibility
#' @param cores Number of cores for parallel computing
#'
#' @return A list with test statistic, p-value, and permutation values
#'
#' @export
dhsic_test <- function(x, type = "gaussian", bw = NULL, expo = 1, 
                       scale_factor = 0.5, group = NULL, 
                       n_perm = 1000, seed = NULL, cores = parallel::detectCores() - 1) {
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
  perm_results <- parallel::mclapply(1:n_perm, function(i) perm_function(i), mc.cores = cores)
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

#' Joint HSIC Independence Test
#'
#' Performs a permutation test based on the joint Hilbert-Schmidt Independence Criterion.
#'
#' @param x List of matrices or vectors
#' @param cc Constant parameter, default is 1
#' @param type Type of kernel or distance
#' @param stat_type Type of statistic ("V", "U", "US", or "UR")
#' @param bw Bandwidth parameter
#' @param expo Exponent parameter
#' @param scale_factor Scaling factor for bandwidth
#' @param group Optional grouping of variables
#' @param n_perm Number of permutations
#' @param seed Random seed for reproducibility
#' @param cores Number of cores for parallel computing
#'
#' @return A list with test statistic, p-value, and permutation values
#'
#' @export
jhsic_test <- function(x, cc = 1, type = "euclidean", stat_type = "V", bw = NULL, 
                       expo = 1, scale_factor = 0.5, group = NULL, 
                       n_perm = 1000, seed = NULL, cores = parallel::detectCores() - 1) {
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
  perm_results <- parallel::mclapply(1:n_perm, function(i) perm_function(i), mc.cores = cores)
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

#' High-Dimensional Two-Sample Test
#'
#' Performs a test for the difference between two high-dimensional distributions.
#'
#' @param x First sample
#' @param y Second sample
#' @param type Type of kernel or distance
#' @param bw Bandwidth parameter
#' @param expo Exponent parameter
#' @param scale_factor Scaling factor for bandwidth
#' @param group Optional grouping of variables
#'
#' @return A list with test statistic and p-value
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
  mmd_val <- mmd(x = D, type = type, u_center = TRUE, n = n, m = m)
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
#' Performs a test for dependence between two high-dimensional datasets.
#'
#' @param x First dataset
#' @param y Second dataset
#' @param type Type of kernel or distance
#' @param bw Bandwidth parameter
#' @param expo Exponent parameter
#' @param scale_factor Scaling factor for bandwidth
#' @param group Optional grouping of variables
#' @param is_distance Logical; whether input matrices are already distance matrices
#'
#' @return A list with test statistic, p-value, correlation value, and degrees of freedom
#'
#' @export
hd_dep_test <- function(x, y, type = "euclidean", bw = NULL, 
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