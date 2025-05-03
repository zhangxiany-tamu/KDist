#' Calculate Kernel or Distance Matrix for Functional Data
#'
#' This function calculates distance or kernel matrices for functional data objects.
#'
#' @param fd_obj A functional data object from the fda package, or a list of individual fd objects
#' @param type Character string specifying the type of distance calculation.
#'   Options include "L2", "L2-deriv", "L2-deriv-weighted"
#' @param output Character string specifying the type of output: "distance", "gaussian", or "laplacian"
#' @param deriv_order Order of derivative to use when type involves derivatives (default: 1)
#' @param weight_deriv Weight for derivative component in weighted distances (default: 0.5)
#' @param bw Bandwidth parameter for Gaussian and Laplacian kernels.
#'   If NULL, it will be automatically determined using median heuristic.
#' @param scale_factor Scaling factor for automatic bandwidth calculation (default: 0.5)
#'
#' @return A distance or kernel matrix depending on the output parameter
#'
#' @examples
#' library(fda)
#'
#' # Create sample data with 50 curves
#' t_points <- seq(0, 1, length.out = 101)
#' n_curves <- 50
#' curves <- matrix(0, n_curves, length(t_points))
#'
#' # Generate 50 different curves
#' for (i in 1:n_curves) {
#'   # Create curves with different frequencies and phases
#'   freq <- 1 + (i %% 5)
#'   phase <- (i %% 10) * pi/10
#'   amplitude <- 0.8 + 0.4 * (i %% 3) / 3
#'   curves[i,] <- amplitude * sin(2*pi*freq*t_points + phase)
#' }
#'
#' # Create functional data object
#' basis <- create.bspline.basis(c(0, 1), 23)
#' fd_obj <- Data2fd(t_points, t(curves), basis)
#'
#' # Calculate different types of matrices
#'
#' # 1. Standard L2 distance matrix
#' dist_mat_L2 <- KDist_fda_matrix(fd_obj, type = "L2", output = "distance")
#'
#' # 2. L2 distance of first derivatives
#' dist_mat_deriv <- KDist_fda_matrix(fd_obj, type = "L2-deriv",
#'                                    deriv_order = 1, output = "distance")
#'
#' # 3. Gaussian kernel based on weighted combination of function and derivatives
#' kernel_mat <- KDist_fda_matrix(fd_obj, type = "L2-deriv-weighted",
#'                               output = "gaussian", weight_deriv = 0.3)
#'
#' # 4. Laplacian kernel with custom bandwidth
#' laplacian_mat <- KDist_fda_matrix(fd_obj, type = "L2",
#'                                  output = "laplacian", bw = 0.5)
#'
#' # Visualize the kernel matrices
#' par(mfrow=c(1,2))
#' image(kernel_mat, main="Gaussian Kernel", xlab="Curve Index", ylab="Curve Index")
#' image(laplacian_mat, main="Laplacian Kernel", xlab="Curve Index", ylab="Curve Index")
#'
#' @export
KDist_fda_matrix <- function(fd_obj, type = "L2", output = "distance", deriv_order = 1,
                             weight_deriv = 0.5, bw = NULL, scale_factor = 0.5) {
  # Check if fda package is available
  if (!requireNamespace("fda", quietly = TRUE)) {
    stop("Package 'fda' needed for this function. Please install it.")
  }

  # Validate output parameter
  if (!output %in% c("distance", "gaussian", "laplacian")) {
    stop("Parameter 'output' must be one of: 'distance', 'gaussian', or 'laplacian'")
  }

  # Convert list of individual fd objects to a single fd object if needed
  if (is.list(fd_obj) && !inherits(fd_obj, "fd")) {
    # Check if all elements are fd objects
    if (!all(sapply(fd_obj, inherits, "fd"))) {
      stop("All elements in the list must be functional data objects")
    }

    # Get the basis from the first fd object
    basis <- fd_obj[[1]]$basis

    # Extract coefficients from each fd object
    coefs <- sapply(fd_obj, function(fd) fd$coefs)

    # Create a single fd object with all curves
    fd_obj <- fda::fd(coefs, basis)
  }

  # Now we have a single fd object, possibly with multiple curves
  n <- ncol(fd_obj$coefs)  # Number of curves

  # Initialize distance matrix
  dist_matrix <- matrix(0, n, n)

  # Calculate distance matrix based on the selected type
  if (type == "L2") {
    # Calculate L2 distances
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        # Create fd objects for individual curves
        fd_i <- fda::fd(fd_obj$coefs[,i], fd_obj$basis)
        fd_j <- fda::fd(fd_obj$coefs[,j], fd_obj$basis)

        # Calculate the difference
        fd_diff <- fd_i - fd_j

        # Calculate L2 norm using inprod
        dist_matrix[i,j] <- sqrt(fda::inprod(fd_diff, fd_diff))
        dist_matrix[j,i] <- dist_matrix[i,j]  # Symmetry
      }
    }
  } else if (type == "L2-deriv") {
    # Calculate the derivative fd object
    fd_deriv <- fda::deriv.fd(fd_obj, deriv_order)

    # Calculate derivative-based distances
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        # Create fd objects for individual derivative curves
        fd_deriv_i <- fda::fd(fd_deriv$coefs[,i], fd_deriv$basis)
        fd_deriv_j <- fda::fd(fd_deriv$coefs[,j], fd_deriv$basis)

        # Calculate the difference
        fd_deriv_diff <- fd_deriv_i - fd_deriv_j

        # Calculate L2 norm of derivatives
        dist_matrix[i,j] <- sqrt(fda::inprod(fd_deriv_diff, fd_deriv_diff))
        dist_matrix[j,i] <- dist_matrix[i,j]  # Symmetry
      }
    }
  } else if (type == "L2-deriv-weighted") {
    # Calculate the derivative fd object
    fd_deriv <- fda::deriv.fd(fd_obj, deriv_order)

    # Calculate weighted distances
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        # Create fd objects for individual curves
        fd_i <- fda::fd(fd_obj$coefs[,i], fd_obj$basis)
        fd_j <- fda::fd(fd_obj$coefs[,j], fd_obj$basis)

        # Create fd objects for individual derivative curves
        fd_deriv_i <- fda::fd(fd_deriv$coefs[,i], fd_deriv$basis)
        fd_deriv_j <- fda::fd(fd_deriv$coefs[,j], fd_deriv$basis)

        # Calculate the differences
        fd_diff <- fd_i - fd_j
        fd_deriv_diff <- fd_deriv_i - fd_deriv_j

        # Calculate weighted L2 norm
        dist_func <- fda::inprod(fd_diff, fd_diff)
        dist_deriv <- fda::inprod(fd_deriv_diff, fd_deriv_diff)

        dist_matrix[i,j] <- sqrt(
          (1 - weight_deriv) * dist_func + weight_deriv * dist_deriv
        )
        dist_matrix[j,i] <- dist_matrix[i,j]  # Symmetry
      }
    }
  } else {
    stop("Unsupported distance type: ", type)
  }

  # Apply kernel transformation if requested
  if (output %in% c("gaussian", "laplacian")) {
    # Determine bandwidth parameter if not provided
    if (is.null(bw)) {
      # Extract upper triangle of the distance matrix (excluding diagonal)
      dist_upper <- dist_matrix[upper.tri(dist_matrix)]

      # Apply median heuristic
      if (length(dist_upper) > 0) {
        bw <- scale_factor * median(dist_upper)
      } else {
        # Default if there's only one curve
        bw <- 1.0
      }
    }

    # Apply kernel transformation
    if (output == "gaussian") {
      # Gaussian kernel: exp(-dist^2 / bw^2)
      kernel_matrix <- exp(-(dist_matrix^2) / (bw^2))
      return(kernel_matrix)
    } else if (output == "laplacian") {
      # Laplacian kernel: exp(-dist / bw)
      kernel_matrix <- exp(-dist_matrix / bw)
      return(kernel_matrix)
    }
  }

  # Return the distance matrix if output = "distance"
  return(dist_matrix)
}
