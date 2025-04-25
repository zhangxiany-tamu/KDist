#' Kernel Conditional Independence Test
#'
#' Performs a kernel-based conditional independence test to assess whether two variables
#' are conditionally independent given a third variable. This implementation follows the
#' methodology described in Zhang et al. (2012).
#'
#' @param formula A formula of the form x ~ y | z, where x is the first variable,
#'   y is the second variable, and z is the conditioning variable or variables.
#'   Multiple variables can be included using the standard formula notation,
#'   e.g., x ~ y | z1 + z2 to condition on both z1 and z2.
#' @param data A data frame containing the variables in the formula.
#' @param type A character vector of length 3 specifying the kernel types for x, y, and z.
#'   Options include "gaussian", "laplacian", "polynomial", "euclidean", "e-dist", "g-dist", or "l-dist".
#'   Default is c("gaussian", "gaussian", "gaussian").
#' @param bw A numeric vector of length 3 specifying bandwidth parameters for x, y, and z.
#'   If NULL (default), bandwidths will be automatically determined.
#' @param expo A numeric vector of length 3 specifying exponent parameters for x, y, and z.
#'   Default is c(1, 1, 1).
#' @param scale_factor Scaling factor for automatic bandwidth calculation (default: 0.5).
#' @param groups A list of length 3, containing optional group membership vectors
#'   for columns in x, y, and z. Used for group-wise distance calculations.
#' @param method Resampling method for obtaining the null distribution (default: "bootstrap").
#'   Options are "bootstrap" or "knn".
#' @param B Number of bootstrap or resampling iterations (default: 299).
#' @param knn_k Number of nearest neighbors to use when method = "knn" (default: 5).
#' @param knn_weighted Logical; whether to use weighted KNN when method = "knn" (default: FALSE).
#' @param num_cores Number of cores for parallel processing (default: 1).
#'
#' @return A list containing:
#'   \item{pvalue}{P-value for the conditional independence test}
#'   \item{statistic}{Test statistic value}
#'   \item{resampling_stats}{Vector of test statistics from resampling}
#'   \item{method}{The resampling method used}
#'   \item{B}{Number of resampling iterations}
#'   \item{method_params}{Parameters used for the resampling method (when applicable)}
#'
#' @details
#' The kernel conditional independence test (KCID) evaluates whether X and Y are conditionally
#' independent given Z, denoted as X ⊥⊥ Y | Z. The test uses the conditional cross-covariance
#' operator in reproducing kernel Hilbert spaces (RKHS) as a measure of conditional dependence.
#'
#' The formula interface uses the '|' operator to specify conditioning variables. For example,
#' the formula 'x ~ y | z' tests whether x and y are conditionally independent given z.
#'
#' @examples
#' # Example 1: Conditional Independence
#' set.seed(123)
#' n <- 200
#' df <- data.frame(
#'   z = runif(n),
#'   stringsAsFactors = FALSE
#' )
#' df$x <- df$z + 0.5*rnorm(n)  # X depends on Z
#' df$y <- df$z + 0.5*rnorm(n)  # Y depends on Z, not on X given Z
#'
#' # Test conditional independence with formula interface
#' result1 <- kcid(x ~ y | z, data = df, method = "bootstrap", B = 100)
#' print(paste("Example 1 p-value:", result1$pvalue))
#' # The p-value should be high, indicating conditional independence
#'
#' # Example 2: Conditional Dependence
#' set.seed(456)
#' n <- 200
#' df2 <- data.frame(
#'   z = runif(n),
#'   x = rnorm(n),
#'   stringsAsFactors = FALSE
#' )
#' # Y depends on both X and Z - conditional dependence
#' df2$y <- 0.7*df2$x + 0.7*df2$z + 0.3*rnorm(n)
#'
#' # Test conditional dependence
#' result2 <- kcid(x ~ y | z, data = df2, method = "bootstrap", B = 100)
#' print(paste("Example 2 p-value:", result2$pvalue))
#' # The p-value should be low, indicating conditional dependence
#'
#' # Example 3: Nonlinear Conditional Dependence
#' set.seed(789)
#' n <- 200
#' df3 <- data.frame(
#'   z = runif(n, -1, 1),
#'   x = rnorm(n),
#'   stringsAsFactors = FALSE
#' )
#' # Y depends on X^2 and Z - nonlinear conditional dependence
#' df3$y <- df3$x^2 + df3$z + 0.3*rnorm(n)
#'
#' # Test nonlinear conditional dependence
#' result3 <- kcid(x ~ y | z, data = df3, method = "bootstrap", B = 100)
#' print(paste("Example 3 p-value:", result3$pvalue))
#' # The p-value should be low, indicating conditional dependence
#'
#' # Example 4: Conditioning on multiple variables
#' set.seed(101)
#' n <- 200
#' df4 <- data.frame(
#'   z1 = runif(n),
#'   z2 = runif(n),
#'   stringsAsFactors = FALSE
#' )
#' df4$x <- df4$z1 + 0.5*rnorm(n)
#' df4$y <- df4$z1 + df4$z2 + 0.3*rnorm(n)
#'
#' # Test conditioning on multiple variables
#' result4 <- kcid(x ~ y | z1 + z2, data = df4, method = "knn", knn_k = 10, B = 100)
#' print(paste("Example 4 p-value:", result4$pvalue))
#'
#' @references
#' Zhang, K., Peters, J., Janzing, D., & Schölkopf, B. (2012). Kernel-based conditional
#' independence test and application in causal discovery. \emph{arXiv preprint arXiv:1202.3775}.
#'
#' @export
kcid <- function(formula, data, type = c("gaussian", "gaussian", "gaussian"),
                 bw = NULL, expo = c(1, 1, 1), scale_factor = 0.5,
                 groups = NULL, method = c("bootstrap", "knn"), B = 299,
                 knn_k = 5, knn_weighted = FALSE, num_cores = 1) {

  # Parse the formula
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
  }

  if (length(formula) != 3) {
    stop("Formula must be of the form 'X ~ Y | Z'")
  }

  # Extract the parts of the formula
  formula_terms <- terms(formula)

  # Handle the conditional part
  rhs <- formula[[3]]
  if (!inherits(rhs, "call") || rhs[[1]] != as.name("|")) {
    stop("Formula must contain the '|' operator to specify conditional variables")
  }

  x_var <- deparse(formula[[2]])
  y_var <- deparse(rhs[[2]])
  z_var <- deparse(rhs[[3]])

  # Function to extract column patterns from dataframe
  get_columns <- function(var_name, df) {
    # Check if var_name contains multiple variables (with + operator)
    if (grepl("\\+", var_name)) {
      # Split by + and trim whitespace
      var_names <- strsplit(var_name, "\\+")[[1]]
      var_names <- trimws(var_names)

      # Extract each variable and combine them
      result <- NULL
      for (v in var_names) {
        if (is.null(result)) {
          result <- get_columns(v, df)
        } else {
          result <- cbind(result, get_columns(v, df))
        }
      }
      return(result)
    }

    # Original code for single variable case
    pattern <- paste0("^", var_name, "\\.[0-9]+$")
    col_names <- grep(pattern, names(df), value = TRUE)

    if (length(col_names) == 0) {
      if (var_name %in% names(df)) {
        return(df[, var_name, drop = FALSE])
      } else {
        stop(paste("Variable", var_name, "not found in data"))
      }
    } else {
      col_names <- col_names[order(as.numeric(sub(paste0(var_name, "\\."), "", col_names)))]
      return(as.matrix(df[, col_names, drop = FALSE]))
    }
  }

  # Extract the data matrices
  data_x <- get_columns(x_var, data)
  data_y <- get_columns(y_var, data)
  data_z <- get_columns(z_var, data)

  # Match method argument
  method <- match.arg(method, c("bootstrap", "knn"))

  # Process type, bw, and expo parameters
  if (length(type) == 1) {
    type <- rep(type, 3)
  } else if (length(type) != 3) {
    stop("'type' must be a character vector of length 1 or 3")
  }
  type_x <- type[1]
  type_y <- type[2]
  type_z <- type[3]

  if (is.null(bw)) {
    bw_x <- NULL
    bw_y <- NULL
    bw_z <- NULL
  } else {
    if (length(bw) == 1) {
      bw <- rep(bw, 3)
    } else if (length(bw) != 3) {
      stop("'bw' must be NULL or a numeric vector of length 1 or 3")
    }
    bw_x <- bw[1]
    bw_y <- bw[2]
    bw_z <- bw[3]
  }

  if (length(expo) == 1) {
    expo <- rep(expo, 3)
  } else if (length(expo) != 3) {
    stop("'expo' must be a numeric vector of length 1 or 3")
  }
  expo_x <- expo[1]
  expo_y <- expo[2]
  expo_z <- expo[3]

  # Process groups parameter
  if (is.null(groups)) {
    group_x <- NULL
    group_y <- NULL
    group_z <- NULL
    group_xz <- NULL
  } else {
    if (!is.list(groups) || length(groups) != 3) {
      stop("'groups' must be NULL or a list of length 3")
    }
    group_x <- groups[[1]]
    group_y <- groups[[2]]
    group_z <- groups[[3]]
    group_xz <- c(group_x, group_z)

    # Check that group lengths match variable dimensions
    if (!is.null(group_x) && length(group_x) != ncol(data_x)) {
      stop(paste("Length of group_x must match number of columns in", x_var,
                 "(", ncol(data_x), "columns )"))
    }
    if (!is.null(group_y) && length(group_y) != ncol(data_y)) {
      stop(paste("Length of group_y must match number of columns in", y_var,
                 "(", ncol(data_y), "columns )"))
    }
    if (!is.null(group_z) && length(group_z) != ncol(data_z)) {
      stop(paste("Length of group_z must match number of columns in", z_var,
                 "(", ncol(data_z), "columns )"))
    }
  }

  # Call the implementation function
  result <- .kcid_impl(data_x, data_y, data_z,
                       type_x = type_x, type_y = type_y, type_z = type_z,
                       bw_x = bw_x, bw_y = bw_y, bw_z = bw_z,
                       expo_x = expo_x, expo_y = expo_y, expo_z = expo_z,
                       scale_factor = scale_factor,
                       group_xz = group_xz, group_y = group_y, group_z = group_z,
                       method = method, B = B,
                       knn_k = knn_k, knn_weighted = knn_weighted,
                       num_cores = num_cores)

  # Add formula and variable info to the result
  result$formula <- formula
  result$variables <- list(
    x = list(name = x_var, cols = ncol(data_x)),
    y = list(name = y_var, cols = ncol(data_y)),
    z = list(name = z_var, cols = ncol(data_z))
  )
  result$call <- match.call()

  return(result)
}

# Internal implementation function - not exported
.kcid_impl <- function(data_x, data_y, data_z, type_x = "gaussian", type_y = "gaussian", type_z = "gaussian",
                       bw_x = NULL, bw_y = NULL, bw_z = NULL, expo_x = 1, expo_y = 1, expo_z = 1, scale_factor = 0.5,
                       group_xz = NULL, group_y = NULL, group_z = NULL,
                       method = c("bootstrap", "knn"), B = 299,
                       knn_k = 5, knn_weighted = FALSE,
                       num_cores = 1) {

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
  Kx <- KDist_matrix(data_x_with_z, type = type_x, bw = bw_x, expo = expo_x, scale_factor = scale_factor, group = group_xz)
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
    Kx_b <- KDist_matrix(x_with_z_b, type = type_x, bw = bw_x, expo = expo_x, scale_factor = scale_factor, group = group_xz)
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
#' @description
#' Performs a test for differences in conditional distributions using kernel-based methods.
#' This function implements both global and local two-sample conditional distribution testing
#' as described in Yan, Li, & Zhang (2024).
#'
#' @details
#' The TCDT tests whether two conditional distributions are equal, either globally:
#'   \deqn{H_0: P(Y_1|X_1) = P(Y_2|X_2)}
#'   for all \eqn{X}, or locally at a specific point \eqn{x_0}:
#'   \deqn{H_0: P(Y_1|X_1=x_0) = P(Y_2|X_2=x_0)}
#'
#' The formula interface uses the syntax `Y1 | X1 ~ Y2 | X2`, where:
#' - Y1 and Y2 are the response variables
#' - X1 and X2 are the conditioning variables
#' - The `|` operator separates response and conditioning variables in each sample
#' - The `~` operator separates the two samples
#'
#' @param formula A formula of the form `Y1 | X1 ~ Y2 | X2` or `Y1 | X1 ~ Y2 | X2 | local(x0)`.
#'        For the local testing case, `x0` specifies the point at which to test.
#'        Multi-dimensional variables can be included using the `+` operator, e.g.,
#'        `Y1.1 + Y1.2 | X1.1 + X1.2 ~ Y2.1 + Y2.2 | X2.1 + X2.2`.
#' @param data A data frame containing all variables specified in the formula.
#' @param x0 Local point for local testing (default: NULL for global testing). Can be specified
#'        either in the formula using `local(x0)` or as a separate parameter.
#' @param stat Test statistic: "cmmd" for conditional maximum mean discrepancy or
#'        "cged" for conditional generalized energy distance (default: "cmmd").
#' @param kernel Type of kernel for Y. Options: "gaussian" or "laplacian" for CMMD,
#'        automatically set to "euclidean" for CGED (default: "gaussian").
#' @param smooth_kernel Type of kernel for X smoothing: "gaussian", "epanechnikov",
#'        or "uniform" (default: "gaussian").
#' @param sampling_method Method for resampling: "bootstrap" or "knn" (default: "bootstrap").
#' @param B Number of bootstrap or resampling iterations (default: 299).
#' @param bandwidth Controls the bandwidth for X smoothing:
#'        - NULL (default): automatic selection
#'        - "undersmoothing": use undersmoothing
#'        - A numeric value: use this value as a multiplier for automatic bandwidth
#'        - A function: specify custom bandwidth function
#' @param bandwidth_adjust Adjustment factor for automatic bandwidth (default: 0.1).
#' @param knn_k Number of nearest neighbors when sampling_method="knn" (default: 5).
#' @param num_cores Number of cores for parallel computation (default: 1).
#'
#' @return A list containing:
#'   \item{problem}{Type of problem: "global" or "local"}
#'   \item{stat}{Test statistic value}
#'   \item{pvalue}{P-value for the test}
#'   \item{bandwidths}{Bandwidths used for smoothing}
#'   \item{method}{Testing method details}
#'   \item{formula}{The original formula}
#'   \item{call}{The function call}
#'
#' @examples
#' # Example 1: Global test with formula interface
#' set.seed(123)
#' n <- 100
#' df <- data.frame(
#'   X1 = rnorm(n),
#'   Y1 = rnorm(n),
#'   X2 = rnorm(n),
#'   Y2 = rnorm(n)
#' )
#'
#' # Make Y1|X1 and Y2|X2 follow different distributions
#' df$Y1 <- df$X1 + rnorm(n)
#' df$Y2 <- df$X2^2 + rnorm(n)
#'
#' # Test the conditional distributions
#' result <- tcdt(Y1 | X1 ~ Y2 | X2, data = df, B = 100)
#' print(result)
#'
#' # Example 2: Local test at a specific point
#' result2 <- tcdt(Y1 | X1 ~ Y2 | X2 | local(0), data = df, B = 100)
#' # OR use the x0 parameter
#' result2b <- tcdt(Y1 | X1 ~ Y2 | X2, data = df, x0 = 0, B = 100)
#'
#' # Example 3: Multi-dimensional variables
#' df$X1.2 <- rnorm(n)
#' df$X2.2 <- rnorm(n)
#' result3 <- tcdt(Y1 | X1 + X1.2 ~ Y2 | X2 + X2.2, data = df, B = 100)
#'
#' @references
#' Yan, J., Li, Z., & Zhang, X. (2024). Distance and Kernel-Based Measures for Global and Local
#' Two-Sample Conditional Distribution Testing. *arXiv preprint arXiv:2210.08149*.
#'
#' @seealso
#' \code{\link{kcid}} for kernel conditional independence testing
#' \code{\link{knn_conditional_sampling}} for the KNN sampling method used in resampling
#' \code{\link{kernel_for_smooth}} for the kernel functions used in local smoothing
#'
#' @export
tcdt <- function(formula, data, x0 = NULL, stat = c("cmmd", "cged"),
                 kernel = c("gaussian", "laplacian"),
                 smooth_kernel = c("gaussian", "epanechnikov", "uniform"),
                 sampling_method = c("bootstrap", "knn"),
                 B = 299, bandwidth = NULL, bandwidth_adjust = 0.1,
                 knn_k = 5, num_cores = 1) {

  # Match arguments
  stat <- match.arg(stat)
  kernel <- match.arg(kernel)
  smooth_kernel <- match.arg(smooth_kernel)
  sampling_method <- match.arg(sampling_method)

  # Parse the formula
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
  }

  # Process the formula to extract Y1|X1 ~ Y2|X2
  formula_str <- deparse(formula)

  # Check for local testing in formula
  local_pattern <- "local\\((.+?)\\)"
  if (grepl(local_pattern, formula_str)) {
    # Extract x0 from formula
    x0_match <- regexpr(local_pattern, formula_str)
    if (x0_match > 0) {
      x0_str <- regmatches(formula_str, regexpr(local_pattern, formula_str))
      x0_str <- gsub("local\\((.+?)\\)", "\\1", x0_str)

      # Parse x0
      x0 <- eval(parse(text = x0_str), envir = data)

      # Remove local part from formula
      formula_str <- gsub("\\s*\\|\\s*local\\(.+?\\)", "", formula_str)
      formula <- as.formula(formula_str)
    }
  }

  # Split formula at the tilde
  formula_parts <- strsplit(formula_str, "~")[[1]]
  if (length(formula_parts) != 2) {
    stop("Formula must be of the form 'Y1 | X1 ~ Y2 | X2'")
  }

  # Split each part at the pipe
  lhs <- trimws(formula_parts[1])
  rhs <- trimws(formula_parts[2])

  lhs_parts <- strsplit(lhs, "\\|")[[1]]
  rhs_parts <- strsplit(rhs, "\\|")[[1]]

  if (length(lhs_parts) != 2 || length(rhs_parts) < 1) {
    stop("Formula must be of the form 'Y1 | X1 ~ Y2 | X2'")
  }

  # Extract variable expressions
  y1_expr <- trimws(lhs_parts[1])
  x1_expr <- trimws(lhs_parts[2])

  y2_expr <- trimws(rhs_parts[1])
  if (length(rhs_parts) >= 2) {
    x2_expr <- trimws(rhs_parts[2])
  } else {
    # If X2 is not specified, assume it's the same as X1
    x2_expr <- x1_expr
  }

  # Helper function to extract variables
  extract_data <- function(expr, data) {
    if (grepl("\\+", expr)) {
      # Multiple variables
      var_names <- strsplit(expr, "\\+")[[1]]
      var_names <- trimws(var_names)
      return(as.matrix(data[, var_names, drop = FALSE]))
    } else {
      # Single variable
      return(as.matrix(data[, expr, drop = FALSE]))
    }
  }

  # Extract data
  Y1 <- extract_data(y1_expr, data)
  X1 <- extract_data(x1_expr, data)
  Y2 <- extract_data(y2_expr, data)
  X2 <- extract_data(x2_expr, data)

  # Set up parameters for .tcdt_impl
  if (is.null(bandwidth)) {
    h <- NULL
  } else if (is.numeric(bandwidth)) {
    h <- bandwidth
  } else if (bandwidth == "undersmoothing") {
    h <- "undersmoothing"
  } else {
    h <- NULL
  }

  # Check if x0 is provided and has correct dimensions
  if (!is.null(x0)) {
    if (is.numeric(x0)) {
      if (length(x0) == 1 && ncol(X1) > 1) {
        # Scalar x0 but multi-dimensional X1, replicate x0
        x0 <- rep(x0, ncol(X1))
      } else if (length(x0) != ncol(X1)) {
        stop("Dimension of x0 does not match dimension of X1")
      }
    } else {
      stop("x0 must be numeric")
    }
  }

  # Map kernel parameters
  kern_mmd <- ifelse(stat == "cmmd", kernel, "euclidean")

  # Call the implementation function with the extracted data
  result <- .tcdt_impl(
    X1 = X1, X2 = X2, Y1 = Y1, Y2 = Y2,
    x0 = x0, B = B,
    h = h, adj_bw = bandwidth_adjust,
    kern_smooth = smooth_kernel,
    stat = stat, kern_mmd = kern_mmd,
    sampling_method = sampling_method,
    knn_k = knn_k,
    num_cores = num_cores
  )

  # Add formula and call information
  result$formula <- formula
  result$call <- match.call()

  # Clean up the results to make a more concise output
  clean_result <- list(
    problem = result$problem,
    stat = result$Tn,
    pvalue = result$pvalue,
    bandwidths = result$h,
    method = list(
      stat_type = result$stat,
      kernel = ifelse(stat == "cmmd", result$kern_mmd, "euclidean"),
      smooth_kernel = result$kern_smooth,
      sampling = result$sampling_method,
      B = result$B
    ),
    formula = formula,
    call = match.call()
  )

  class(clean_result) <- c("tcdt", class(result))

  return(clean_result)
}

#' Implementation of Two-Sample Conditional Distribution Test
#'
#' @description
#' Internal implementation function for the two-sample conditional distribution test.
#' This function should not be called directly by users unless they need
#' fine-grained control over all parameters.
#'
#' @param X1 First conditioning dataset (matrix, data frame, or vector)
#' @param X2 Second conditioning dataset with the same number of columns as X1
#' @param Y1 First response dataset (matrix, data frame, or vector)
#' @param Y2 Second response dataset with the same number of columns as Y1
#' @param x0 Local point for local testing (default: NULL for global testing)
#' @param B Number of bootstrap or resampling iterations (default: 299)
#' @param h Bandwidth for local smoothing of X (default: NULL for automatic selection)
#' @param adj_bw Adjustment factor for automatic bandwidth (default: 0.1)
#' @param const_factor Constant scaling factor for automatic bandwidth (default: 1)
#' @param h_boot Bandwidth for bootstrap resampling (default: NULL)
#' @param adj_bw_boot Adjustment factor for bootstrap bandwidth (default: 0)
#' @param const_factor_boot Constant scaling factor for bootstrap bandwidth (default: 1)
#' @param scale_factor Scaling factor for kernel bandwidth in MMD calculation (default: 1)
#' @param kern_smooth Kernel type for smoothing (default: "gaussian")
#' @param nu Order of smoothing kernel (default: 2)
#' @param stat Test statistic type (default: "cmmd")
#' @param kern_mmd Kernel type for MMD when stat="cmmd" (default: "gaussian")
#' @param sampling_method Resampling method (default: "bootstrap")
#' @param knn_k Number of nearest neighbors for KNN (default: 5)
#' @param knn_weighted Whether to use weighted KNN (default: FALSE)
#' @param num_cores Number of cores for parallel computation (default: 1)
#'
#' @return A list containing test results similar to the original tcdt function
#'
#' @keywords internal
.tcdt_impl <- function(X1, X2, Y1, Y2, x0 = NULL, B = 299,
                       h = NULL, adj_bw = 0.1, const_factor = 1,
                       h_boot = NULL, adj_bw_boot = 0, const_factor_boot = 1, scale_factor = 1,
                       kern_smooth = "gaussian", nu = 2,
                       stat = "cmmd", kern_mmd = "gaussian",
                       sampling_method = "bootstrap", knn_k = 5, knn_weighted = FALSE,
                       num_cores = 1) {

  # Load required packages
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("The 'parallel' package is required for this function. Please install it.")
  }

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
    dY <- KDist_matrix(Y_pool, type = "laplacian", scale_factor = scale_factor)
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
#' @description
#' Performs KNN-based conditional sampling from a dataset. This function samples from
#' the conditional distribution at each data point by randomly selecting from its K
#' nearest neighbors, which serves as a nonparametric approximation to sampling from
#' the true conditional distribution.
#'
#' @details
#' For each observation in the input dataset, the function:
#' 1. Finds its K nearest neighbors based on Euclidean distance
#' 2. Samples one of these neighbors, either uniformly or with weights inversely
#'    proportional to their distances
#'
#' This approach is particularly useful for resampling in conditional independence
#' and distribution tests, as it preserves the conditional structure in the data.
#'
#' @param z_data Conditioning dataset (matrix or data frame). Each row represents an observation.
#' @param k Number of nearest neighbors to consider (default: 5).
#' @param weighted Logical; whether to use weighted sampling based on distances (default: TRUE).
#'   When TRUE, closer neighbors have higher probability of being selected.
#'
#' @return A vector of sampled indices of the same length as the number of rows in z_data.
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
#' @description
#' Creates a kernel function for local smoothing with the specified properties.
#' The returned function can be used to apply various types of kernel smoothing
#' in nonparametric conditional testing procedures.
#'
#' @details
#' This function returns a kernel function with the specified type and order.
#' The returned kernel is properly normalized and has appropriate support.
#' Higher order kernels (higher values of nu) provide better asymptotic
#' theoretical properties but may be less stable in small samples.
#'
#' The returned function takes two arguments:
#' - x: The input values (typically differences between data points)
#' - h: The bandwidth parameter controlling the smoothing amount
#'
#' @param type Type of kernel. Options:
#'   - "gaussian": Gaussian kernel (normal density)
#'   - "epanechnikov": Epanechnikov kernel (quadratic function with compact support)
#'   - "uniform": Uniform kernel (rectangular function)
#' @param nu Order of the kernel. Higher values provide higher-order approximations (default: 2).
#'   Available values depend on the kernel type:
#'   - Gaussian: 2, 4, or 6
#'   - Epanechnikov: 2, 4, or 6
#'   - Uniform: only 2
#' @param param Additional parameter. For uniform kernel, specifies the boundary (default: NULL,
#'   which uses boundary = 1 for uniform kernel).
#'
#' @return A function that computes the kernel values with the given specifications.
#'   The returned function has an additional attribute "nu" storing the order.
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
