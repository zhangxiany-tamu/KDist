#' Kernel Conditional Independence Test
#'
#' Tests for conditional independence between two variables given a third variable.
#'
#' @param data_x First dataset
#' @param data_y Second dataset
#' @param data_z Conditioning dataset
#' @param type_x Type of kernel for first dataset
#' @param type_y Type of kernel for second dataset
#' @param type_z Type of kernel for conditioning dataset
#' @param bw_x Bandwidth parameter for first dataset
#' @param bw_y Bandwidth parameter for second dataset
#' @param bw_z Bandwidth parameter for conditioning dataset
#' @param expo_x Exponent parameter for first dataset
#' @param expo_y Exponent parameter for second dataset
#' @param expo_z Exponent parameter for conditioning dataset
#' @param scale_factor Scaling factor for bandwidth
#' @param group_x Grouping for first dataset variables
#' @param group_y Grouping for second dataset variables
#' @param group_z Grouping for conditioning dataset variables
#' @param method Resampling method: "bootstrap" or "knn"
#' @param B Number of bootstrap replicates
#' @param knn_k Number of nearest neighbors for KNN method
#' @param knn_weighted Whether to use weighted KNN
#' @param num_cores Number of cores for parallel computing
#'
#' @return A list with p-value, test statistic, and other results
#'
#' @export
kcid <- function(data_x, data_y, data_z, type_x = "gaussian", type_y = "gaussian", type_z = "gaussian",
                 bw_x = NULL, bw_y = NULL, bw_z = NULL, expo_x = 1, expo_y = 1, expo_z = 1, scale_factor = 0.5, 
                 group_x = NULL, group_y = NULL, group_z = NULL,
                 method = c("bootstrap", "knn"), B = 299, 
                 knn_k = 5, knn_weighted = FALSE,
                 num_cores = parallel::detectCores() - 1) {
  
  # Match method argument
  method <- match.arg(method)
  
  # Load required packages for parallel processing
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("The 'parallel' package is required for this function. Please install it.")
  }
  
  # Regularization parameters
  epsilon_x <- 1e-3
  epsilon_y <- 1e-3
  
  # Standardize data
  data_x <- scale(data_x)
  data_x[is.na(data_x)] <- 0
  
  data_y <- scale(data_y)
  data_y[is.na(data_y)] <- 0
  
  data_z <- scale(data_z)
  data_z[is.na(data_z)] <- 0
  
  # Concatenate X and Z
  data_x_with_z <- cbind(data_x, 0.5 * data_z)
  
  # Compute kernel matrices
  Kx <- KDist_matrix(data_x_with_z, type = type_x, bw = bw_x, expo = expo_x, scale_factor = scale_factor, group = group_x)
  Ky <- KDist_matrix(data_y, type = type_y, bw = bw_y, expo = expo_y, scale_factor = scale_factor, group = group_y)
  Kz <- KDist_matrix(data_z, type = type_z, bw = bw_z, expo = expo_z, scale_factor = scale_factor, group = group_z)
  
  # Center kernel matrices
  Kx <- matrix_v_center(Kx)
  Ky <- matrix_v_center(Ky)
  Kz <- matrix_v_center(Kz)
  
  # Center kernel matrices for regression
  n <- nrow(data_x)
  I_n <- diag(n)
  
  # Using ridge regression for centering
  Rzx <- epsilon_x * solve(Kz + epsilon_x * I_n)
  KxR <- Rzx %*% Kx %*% Rzx
  
  if (epsilon_x != epsilon_y) {
    Rzy <- epsilon_y * solve(Kz + epsilon_y * I_n)
    KyR <- Rzy %*% Ky %*% Rzy
  } else {
    # Use the same residual matrix if epsilon values are the same
    KyR <- Rzx %*% Ky %*% Rzx
  }
  
  # Calculate test statistic
  V_stat <- sum(KxR * KyR)
  
  # Setup for different resampling methods
  if (method == "bootstrap") {
    # Create kernel matrix for Z for local smoothing
    G_b <- kernel_for_smooth("gaussian", nu = 2)
    
    # Generate bandwidth for local bootstrap
    A <- pmin(
      apply(X = data_z, MARGIN = 2, FUN = IQR) / 1.34,
      apply(X = data_z, MARGIN = 2, FUN = sd)
    )
    h_b <- n^(-1 / (ncol(data_z)/2 + 4)) * A
    
    # Create local smoothing kernel matrix
    Gz <- matrix(1, nrow = n, ncol = n)
    for (j in 1:ncol(data_z)) {
      Gz <- Gz * G_b(outer(data_z[, j], data_z[, j], "-"), h_b[j])
    }
    
    # Normalize weights
    w <- sweep(Gz, 2, colSums(Gz), "/", check.margin = FALSE)
  }
  
  # Define resampling function for both methods
  resampling_iteration <- function(b) {
    # Generate resampled indices
    if (method == "bootstrap") {
      # Local bootstrap sampling
      idx_b <- numeric(n)
      for (i in 1:n) {
        idx_b[i] <- sample(1:n, 1, prob = w[, i])
      }
    } else if (method == "knn") {
      # KNN conditional sampling
      idx_b <- knn_conditional_sampling(data_z, k = knn_k, weighted = knn_weighted)
    }
    
    # Get bootstrap samples
    x_b <- data_x[idx_b, , drop = FALSE]
    y_b <- data_y
    z_b <- data_z
    
    # Compute kernel matrices for resampled data
    x_with_z_b <- cbind(x_b, 0.5 * z_b)
    Kx_b <- KDist_matrix(x_with_z_b, type = type_x, bw = bw_x, expo = expo_x, scale_factor = scale_factor, group = group_x)
    Ky_b <- KDist_matrix(y_b, type = type_y, bw = bw_y, expo = expo_y, scale_factor = scale_factor, group = group_y)
    Kz_b <- KDist_matrix(z_b, type = type_z, bw = bw_z, expo = expo_z, scale_factor = scale_factor, group = group_z)
    
    # Center resampled kernel matrices
    Kx_b <- matrix_v_center(Kx_b)
    Ky_b <- matrix_v_center(Ky_b)
    Kz_b <- matrix_v_center(Kz_b)
    
    # Center resampled kernel matrices for regression
    Rzx_b <- epsilon_x * solve(Kz_b + epsilon_x * I_n)
    KxR_b <- Rzx_b %*% Kx_b %*% Rzx_b
    
    if (epsilon_x != epsilon_y) {
      Rzy_b <- epsilon_y * solve(Kz_b + epsilon_y * I_n)
      KyR_b <- Rzy_b %*% Ky_b %*% Rzy_b
    } else {
      KyR_b <- Rzx_b %*% Ky_b %*% Rzx_b
    }
    
    # Calculate resampled test statistic
    return(sum(KxR_b * KyR_b))
  }
  
  # Parallel computation of resampling statistics
  V_stat_b <- unlist(mclapply(1:B, resampling_iteration, mc.cores = num_cores))
  pvalue <- (1 + sum(V_stat_b > V_stat)) / (1 + B)
  
  return(list(
    pvalue = pvalue,
    statistic = V_stat,
    resampling_stats = V_stat_b,
    method = method,
    B = B,
    method_params = if(method == "knn") list(k = knn_k, weighted = knn_weighted) else NULL
  ))
}

#' Two-Sample Conditional Distribution Test
#'
#' Performs a test for differences in conditional distributions.
#'
#' @param X1 First conditioning dataset
#' @param X2 Second conditioning dataset
#' @param Y1 First response dataset
#' @param Y2 Second response dataset
#' @param x0 Local point for local testing (NULL for global testing)
#' @param B Number of bootstrap replications
#' @param h Bandwidth for local smoothing of X
#' @param adj_bw Adjustment factor for automatic bandwidth
#' @param const_factor Constant factor for automatic bandwidth
#' @param h_boot Bandwidth for bootstrap
#' @param adj_bw_boot Adjustment factor for bootstrap bandwidth
#' @param const_factor_boot Constant factor for bootstrap bandwidth
#' @param scale_factor Scaling factor for bandwidth
#' @param kern_smooth Kernel type for smoothing
#' @param nu Order of smoothing kernel
#' @param stat Test statistic type
#' @param kern_mmd Kernel type for MMD
#' @param sampling_method Sampling method
#' @param knn_k Number of nearest neighbors for KNN
#' @param knn_weighted Whether to use weighted KNN
#' @param num_cores Number of cores for parallel computing
#'
#' @return A list with test results
#'
#' @export
tcdt <- function(X1, X2, Y1, Y2, x0 = NULL, B = 299,
                 h = NULL, adj_bw = 0.1, const_factor = 1,
                 h_boot = NULL, adj_bw_boot = 0, const_factor_boot = 1, scale_factor = 1,
                 kern_smooth = c("gaussian", "epanechnikov", "uniform"), nu = 2,
                 stat = c("cmmd", "cged"), kern_mmd = c("gaussian", "laplacian"),
                 sampling_method = c("bootstrap", "knn"), knn_k = 5, knn_weighted = FALSE,
                 num_cores = parallel::detectCores() - 1) {
  
  # Load required packages
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("The 'parallel' package is required for this function. Please install it.")
  }
  
  # Match arguments
  stat <- match.arg(stat)
  kern_mmd <- match.arg(kern_mmd)
  kern_smooth <- match.arg(kern_smooth)
  sampling_method <- match.arg(sampling_method)
  
  # Coerce and check
  for (name_var in c("X1", "X2", "Y1", "Y2")) {
    obj <- get(name_var)
    if (inherits(obj, "data.frame")) {
      assign(name_var, as.matrix(obj))
    } else if (inherits(obj, "numeric")) {
      assign(name_var, matrix(obj, nrow = length(obj)))
    } else if (!inherits(obj, "matrix")) {
      stop(paste(name_var, "should be a vector, matrix or data frame."))
    }
  }
  if (ncol(X1) != ncol(X2)) {
    stop("`X1` and `X2` should have the same number of columns.")
  }
  if (ncol(Y1) != ncol(Y2)) {
    stop("`Y1` and `Y2` should have the same number of columns.")
  }
  
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  p <- ncol(X1)
  q <- ncol(Y1)
  
  n <- n1 + n2
  idx1 <- (1:n1)
  idx2 <- ((n1+1):n)
  
  if (is.null(x0)) {
    problem <- "global"
  } else {
    problem <- "local"
    x0 <- as.numeric(x0)
    if (length(x0) != p) {
      stop("The dimension of `x0` should be matched with the dimension of `X1` and `X2`.")
    }
  }
  
  list_return <- list(
    problem = problem,
    B = B,
    kern_smooth = kern_smooth,
    stat = stat,
    sampling_method = sampling_method
  )
  if (stat == "cmmd") {
    list_return$kern_mmd <- kern_mmd
  }
  
  # Kernel matrix for Y
  Y_pool <- rbind(Y1, Y2)
  if (stat == "cged") {
    dY <- KDist_matrix(Y_pool, type = "euclidean")
  } else if (stat == "cmmd" && kern_mmd == "gaussian") {
    dY <- KDist_matrix(Y_pool, type = "gaussian", scale_factor = scale_factor)
  } else if (stat == "cmmd" && kern_mmd == "laplacian") {
    dY <- KDist_matrix(Y_pool, type = "laplace", scale_factor = scale_factor)
  }
  kY11 <- dY[idx1, idx1]
  kY22 <- dY[idx2, idx2]
  kY12 <- dY[idx1, idx2]
  
  if (stat == "cged") {
    # Conditional energy distance
    kY11 <- -kY11
    kY22 <- -kY22
    kY12 <- -kY12
  } 
  
  # Local smoothing kernel and bandwidths
  G <- kernel_for_smooth(kern_smooth, nu)
  nu <- attr(G, "nu")
  
  A1 <- pmin(
    apply(X = X1, MARGIN = 2, FUN = IQR) / 1.34,
    apply(X = X1, MARGIN = 2, FUN = sd)
  )
  A2 <- pmin(
    apply(X = X2, MARGIN = 2, FUN = IQR) / 1.34,
    apply(X = X2, MARGIN = 2, FUN = sd)
  )
  
  if (is.null(h) || h == "undersmoothing") {
    if (problem == "global") {
      h1 <- const_factor * n1^(-1 / (p / 2 + nu - adj_bw)) * A1
      h2 <- const_factor * n2^(-1 / (p / 2 + nu - adj_bw)) * A2
    } else if (problem == "local") {
      h1 <- const_factor * n1^(-1 / (p + nu - adj_bw)) * A1
      h2 <- const_factor * n2^(-1 / (p + nu - adj_bw)) * A2
    }
  } else if (is.function(h)) {
    h1 <- h(n1, p, nu, A1)
    h2 <- h(n2, p, nu, A2)
  } else {
    # Using provided bandwidth
    if (is.list(h) && length(h) == 2) {
      h1 <- h[[1]]
      h2 <- h[[2]]
    } else {
      stop("If h is provided, it should be a list of two vectors with appropriate lengths.")
    }
  }
  list_return$h <- cbind(h1, h2)
  
  # Test statistic
  if (problem == "global") {
    G1_11 <- matrix(1, nrow = n1, ncol = n1)
    G1_12 <- matrix(1, nrow = n1, ncol = n2)
    G2_21 <- matrix(1, nrow = n2, ncol = n1)
    G2_22 <- matrix(1, nrow = n2, ncol = n2)
    for (j in 1:p) {
      G1_11 <- G1_11 * G(outer(X1[, j], X1[, j], "-"), h1[j])
      G1_12 <- G1_12 * G(outer(X1[, j], X2[, j], "-"), h1[j])
      G2_21 <- G2_21 * G(outer(X2[, j], X1[, j], "-"), h2[j])
      G2_22 <- G2_22 * G(outer(X2[, j], X2[, j], "-"), h2[j])
    }
    
    S1_1 <- colSums(G2_21)
    S1_2 <- colSums(G1_12)
    S2_1 <- S1_1^2 - colSums(G2_21^2) # n1-vector
    S2_2 <- S1_2^2 - colSums(G1_12^2) # n2-vector
    
    Tn <- teststatg(n1, n2, kY11, kY22, kY12, G1_11, G1_12, G2_22, G2_21, S1_1, S1_2, S2_1, S2_2)
  } else if (problem == "local") {
    G1X <- rep(1, n1)
    G2X <- rep(1, n2)
    for (j in 1:p) {
      G1X <- G1X * G(drop(outer(X1[, j], x0[j], "-")), h1[j])
      G2X <- G2X * G(drop(outer(X2[, j], x0[j], "-")), h2[j])
    }
    
    S1_1 <- sum(G1X)
    S1_2 <- sum(G2X)
    S2_1 <- S1_1^2 - sum(G1X^2) # scalar
    S2_2 <- S1_2^2 - sum(G2X^2) # scalar
    
    Tn <- teststatl(n1, n2, kY11, kY22, kY12, G1X, G2X, S1_1, S1_2, S2_1, S2_2)
  }
  list_return$Tn <- Tn
  
  # Bootstrap / KNN sampling for p-value calculation
  if (B > 0) {
    # Prepare for sampling methods
    if (sampling_method == "bootstrap") {
      # Kernel for local bootstrap
      G_b <- kernel_for_smooth("gaussian", nu = 2)
      
      # Kernel matrix used for local bootstrap
      if (is.null(h_boot) || h_boot == "thumb") {
        h1_b <- const_factor_boot * n1^(-1 / (p + 4 + adj_bw_boot)) * A1
        h2_b <- const_factor_boot * n2^(-1 / (p + 4 + adj_bw_boot)) * A2
        
      } else if (h_boot == "undersmoothing") {
        h1_b <- const_factor_boot * n1^(-1 / (p / 2 + 2)) * A1
        h2_b <- const_factor_boot * n2^(-1 / (p / 2 + 2)) * A2
        
      } else if (is.function(h_boot)) {
        h1_b <- h_boot(n1, p, nu, A1)
        h2_b <- h_boot(n2, p, nu, A2)
      } else {
        # Using provided bootstrap bandwidth
        if (is.list(h_boot) && length(h_boot) == 2) {
          h1_b <- h_boot[[1]]
          h2_b <- h_boot[[2]]
        } else {
          stop("If h_boot is provided, it should be a list of two vectors with appropriate lengths.")
        }
      }
      list_return$h_b <- cbind(h1_b, h2_b)
      
      G1_11_b <- matrix(1, nrow = n1, ncol = n1)
      G1_12_b <- matrix(1, nrow = n1, ncol = n2)
      G2_21_b <- matrix(1, nrow = n2, ncol = n1)
      G2_22_b <- matrix(1, nrow = n2, ncol = n2)
      for (j in 1:p) {
        G1_11_b <- G1_11_b * G_b(outer(X1[, j], X1[, j], "-"), h1_b[j])
        G1_12_b <- G1_12_b * G_b(outer(X1[, j], X2[, j], "-"), h1_b[j])
        G2_21_b <- G2_21_b * G_b(outer(X2[, j], X1[, j], "-"), h2_b[j])
        G2_22_b <- G2_22_b * G_b(outer(X2[, j], X2[, j], "-"), h2_b[j])
      }
      Gp <- rbind(cbind(G1_11_b, G1_12_b), cbind(G2_21_b, G2_22_b))
      
      w <- sweep(Gp, 2, colSums(Gp), "/", check.margin = FALSE)
    } else if (sampling_method == "knn") {
      # For KNN, we need to pool the X data
      X_pool <- rbind(X1, X2)
      list_return$knn_params <- list(k = knn_k, weighted = knn_weighted)
    }
    
    # Define bootstrap function to be used with mclapply
    bootstrap_iteration <- function(b) {
      # Choose sampling method
      if (sampling_method == "bootstrap") {
        # Local bootstrap sampling
        idx_b <- numeric(n)
        for (i in 1:n) {
          idx_b[i] <- sample(1:n, 1, prob = w[, i])
        }
      } else if (sampling_method == "knn") {
        # KNN sampling
        idx_b <- knn_conditional_sampling(X_pool, k = knn_k, weighted = knn_weighted)
      }
      
      # Extract indices for each sample
      idx1_b <- idx_b[idx1]
      idx2_b <- idx_b[idx2]
      
      # Calculate kernel matrices for resampled data
      if (stat == "cged") {
        kYb11 <- -dY[idx1_b, idx1_b]
        kYb22 <- -dY[idx2_b, idx2_b]
        kYb12 <- -dY[idx1_b, idx2_b]
      } else if (stat == "cmmd") {
        kYb11 <- dY[idx1_b, idx1_b]
        kYb22 <- dY[idx2_b, idx2_b]
        kYb12 <- dY[idx1_b, idx2_b]
      }
      
      # Calculate test statistic for bootstrap sample
      if (problem == "global") {
        return(teststatg(n1, n2, kYb11, kYb22, kYb12, G1_11, G1_12, G2_22, G2_21, S1_1, S1_2, S2_1, S2_2))
      } else if (problem == "local") {
        return(teststatl(n1, n2, kYb11, kYb22, kYb12, G1X, G2X, S1_1, S1_2, S2_1, S2_2))
      }
    }
    
    # Parallel computation of bootstrap statistics
    Tn_b <- unlist(parallel::mclapply(1:B, bootstrap_iteration, mc.cores = num_cores))
    pvalue <- (1 + sum(Tn_b > Tn)) / (1 + B)
    
    list_return$Tn_b <- Tn_b
    list_return$pvalue <- pvalue
  }
  
  return(list_return)
}

#' K-Nearest Neighbor Conditional Sampling
#'
#' Performs KNN-based conditional sampling.
#'
#' @param z_data Conditioning dataset
#' @param k Number of nearest neighbors
#' @param weighted Whether to use weighted sampling
#'
#' @return A vector of sampled indices
#'
#' @export
knn_conditional_sampling <- function(z_data, k = 5, weighted = TRUE) {
  n <- nrow(z_data)
  idx <- numeric(n)
  
  # Calculate all pairwise distances at once
  # dist() is more efficient than calculating manually
  dist_matrix <- as.matrix(dist(z_data)^2)
  
  # Set self-distances to infinity
  diag(dist_matrix) <- Inf
  
  for (i in 1:n) {
    # Get k nearest neighbors
    nn_indices <- order(dist_matrix[i, ])[1:k]
    
    if (weighted) {
      # Get distances to k nearest neighbors
      nn_dists <- dist_matrix[i, nn_indices]
      
      # Convert distances to weights (smaller distance = larger weight)
      weights <- 1/(nn_dists + 1e-10)
      weights <- weights / sum(weights)
      
      # Sample with probability proportional to weights
      idx[i] <- sample(nn_indices, 1, prob = weights)
    } else {
      # Uniform sampling from neighbors
      idx[i] <- sample(nn_indices, 1)
    }
  }
  return(idx)
}

#' Kernel Function for Local Smoothing
#'
#' Creates a kernel function for local smoothing with the specified properties.
#'
#' @param type Type of kernel: "gaussian", "epanechnikov", or "uniform"
#' @param nu Order of the kernel
#' @param param Additional parameter (for uniform kernel)
#'
#' @return A kernel function
#'
#' @export
kernel_for_smooth <- function(type = c("gaussian", "epanechnikov", "uniform"), 
                              nu = 2, param = NULL) {
  type <- match.arg(type)
  
  if (type == "gaussian") {
    if (nu == 2) {
      kernel <- function(x, h) {
        exp(-(x / h)^2 / 2) / sqrt(2 * pi) / h
      }
    } else if (nu == 4) {
      kernel <- function(x, h) {
        (3 - (x/h)^2) / 2 * exp(-(x / h)^2 / 2) / sqrt(2 * pi) / h
      }
    } else if (nu == 6) {
      kernel <- function(x, h) {
        (15 - 10 * (x/h)^2 + (x/h)^4) / 8 * exp(-(x / h)^2 / 2) / sqrt(2 * pi) / h
      }
    } else {
      stop("For gaussian smoothing kernel, only support `nu` = 2, 4, or 6 currently.")
    }
  } else if (type == "epanechnikov") {
    if (nu == 2) {
      kernel <- function(x, h) {
        ifelse(
          abs(x / h) <= 1,
          (1 - (x / h)^2) * 3 / 4 / h,
          0
        )
      }
    } else if (nu == 4) {
      kernel <- function(x, h) {
        ifelse(
          abs(x / h) <= 1,
          (1 - (x / h)^2)^2 * 15 / 16 / h,
          0
        )
      } 
    } else if (nu == 6) {
      kernel <- function(x, h) {
        ifelse(
          abs(x / h) <= 1,
          (1 - (x / h)^2)^3 * 35 / 32 / h,
          0
        )
      } 
    } else {
      stop("For epanechnikov smoothing kernel, only support `nu` = 2, 4, or 6 currently.")
    }
  } else if (type == "uniform") {
    if (nu == 2) {
      if (is.null(param)) {
        boundary <- 1
      } else {
        boundary <- param
      }
      kernel <- function(x, h) {
        ifelse(
          (x / h) >= -boundary & (x / h) <= boundary,
          1 / (2 * boundary * h),
          0
        )
      }
    } else {
      stop("For uniform smoothing kernel, only support `nu` = 2 currently.")
    }
  }
  
  attr(kernel, "nu") <- nu
  return(kernel)
}