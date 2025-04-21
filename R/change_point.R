#' Single Change-Point Detection
#'
#' Detects a single change point in time series or sequential data.
#'
#' @param data A matrix or data frame with rows representing time points
#' @param type Type of distance to use (default: "e-dist")
#' @param bw Bandwidth parameter
#' @param expo Exponent parameter
#' @param scale_factor Scaling factor for bandwidth
#' @param group Optional grouping of variables
#' @param B Number of permutations
#' @param alpha Significance level
#' @param seeds Random seed
#' @param method Testing method: "permutation" or "asymptotic"
#' @param stat_val Reference statistics for asymptotic method
#'
#' @return A list with change point location, p-value, statistic, and clustering
#'
#' @export
kcpd_single <- function(data, type = "e-dist", bw = NULL, expo = 1, scale_factor = 0.5, group = NULL,  
                          B = 199, alpha = 0.05, seeds = NULL, 
                          method = "permutation", stat_val = NULL) {
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Input data must be a matrix or data frame")
  }
  n <- nrow(data)
  
  # Set number of cores 
  mc.cores <- detectCores() - 1
  
  # Calculate the distance matrix for the original data
  full_dist <- KDist_matrix(data = data, type = type, bw = bw, expo = expo, scale_factor = scale_factor, group = group)
    
  # Apply stat_recur function to get test statistic and change point
  result <- stat_recur(full_dist, type = type)
  
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
      perm_result <- stat_recur(perm_dist, type = type)
      
      # Return the test statistic (first element of the result list)
      return(perm_result[[1]])
    }, mc.cores = mc.cores))
    
    # Calculate p-value
    pval <- (1 + sum(stat_vec > val_test_stat)) / (1 + B)
    
  } else if (method == "asymptotic") {
    # Check if stat_val is provided
    if (is.null(stat_val) || !is.numeric(stat_val)) {
      stop("When method='asymptotic', you must provide the stat_val parameter as a numeric vector")
    }
    
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

#' Seeded Binary Segmentation for Multiple Change Points
#'
#' Detects multiple change points using seeded binary segmentation.
#'
#' @param data A matrix or data frame with rows representing time points
#' @param type Type of distance to use (default: "e-dist")
#' @param bw Bandwidth parameter
#' @param expo Exponent parameter
#' @param scale_factor Scaling factor for bandwidth
#' @param group Optional grouping of variables
#' @param B Number of permutations
#' @param alpha Significance level
#' @param seeds Random seed
#' @param method Testing method: "permutation" or "asymptotic"
#' @param stat_val Reference statistics for asymptotic method
#' @param decay Decay factor for interval selection
#' @param unique_int Whether to use unique intervals
#' @param bound Minimum segment size
#' @param num_cores Number of cores for parallel processing
#'
#' @return A list with change point locations, p-values, and clustering
#'
#' @export
kcpd_sbs <- function(data, type = "e-dist", bw = NULL, expo = 1, scale_factor = 0.5, 
                     group = NULL, B = 199, alpha = 0.05, seeds = NULL,
                     method = "asymptotic", stat_val = NULL, 
                     decay = sqrt(2), unique_int = FALSE, bound = 20, num_cores = NULL) {
  
  # Input validation
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Input data must be a matrix or data frame")
  }
  
  n <- nrow(data)
  
  # Get seeded intervals
  intervals <- get_seeded_intervals(n, decay = decay, unique.int = unique_int, bound = bound)
  
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
        stat_val = stat_val
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

#' Generate Seeded Intervals
#'
#' Creates seeded intervals for binary segmentation.
#'
#' @param n Length of time series
#' @param decay Decay factor for interval generation
#' @param unique.int Whether to return unique intervals
#' @param bound Minimum interval size
#'
#' @return A matrix of interval start and end points
#'
#' @export
get_seeded_intervals <- function(n, decay = sqrt(2), unique.int = FALSE, bound = 2) {
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
  if (unique.int) {
    return(unique(boundary_mtx))
  }
  
  boundary_mtx
}

#' Calculate Statistics Recursively
#'
#' Calculates test statistics for change point detection recursively.
#' This is an internal function used by kcpd_single.
#' 
#' @param full_dist Full distance matrix
#' @param type Type of distance
#'
#' @return A list with test statistic, change point, and all statistics
stat_recur <- function(full_dist, type) {
  # Get the dimensions of the distance matrix
  n <- nrow(full_dist)
  
  # Create index vectors for partitioning
  # Ensuring we have at least 4 elements on each side
  n1 <- 4:(n-4)
  n2 <- n - n1
  
  # Calculate all necessary statistics
  mmd_all <- mmd_recur(full_dist, type = type)
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

#' Calculate MMD Recursively
#'
#' Calculates MMD values recursively for various split points.
#'
#' @param full_dist Full distance matrix
#' @param type Type of distance
#'
#' @return A vector of MMD values
mmd_recur <- function(full_dist, type) {
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

#' Adjusted Rand Index
#'
#' Calculates the Adjusted Rand Index between two clusterings.
#'
#' @param x First clustering
#' @param y Second clustering
#'
#' @return Adjusted Rand Index value
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