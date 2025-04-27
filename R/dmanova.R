#' Distance-Based Multivariate Analysis of Variance (D-MANOVA)
#'
#' @description
#' Performs distance-based multivariate analysis of variance (D-MANOVA) for testing
#' the dependence of a multivariate response Y on predictors X after adjusting for
#' covariates Z.
#'
#' @details
#' D-MANOVA transforms the multivariate response into a kernel or distance matrix and
#' then tests whether this distance structure depends on a set of predictors after
#' accounting for covariates. The approach is based on partitioning the total sum of
#' squares of the centered distance matrix.
#'
#' The formula interface can be specified in two ways:
#' \itemize{
#'   \item \code{Y ~ X | Z} - Test the dependence of Y on X after adjusting for Z.
#'   \item \code{Y ~ X} - Test the dependence of Y on X with only an intercept term as covariate.
#' }
#'
#' By default, Z will contain an intercept. If Z is not specified in the formula (i.e., using
#' \code{Y ~ X} format), then Z will be set to an intercept-only model.
#'
#' If the input Y is already a distance matrix, it will be used directly, bypassing the distance calculation step.
#'
#' When \code{method = "asymptotic"}, the function uses an asymptotic approximation to
#' calculate the p-value. When \code{method = "permutation"}, a permutation procedure
#' is used to obtain the p-value, which is more robust but computationally more intensive.
#'
#' @param formula A formula of the form Y ~ X | Z or Y ~ X, where Y is the response, X contains
#'        the predictors of interest, and Z contains the covariates to adjust for (if specified).
#' @param data An optional data frame containing the variables in the formula.
#' @param type Type of distance or kernel to use for the response (default: "euclidean").
#'        Options include "euclidean", "polynomial", "gaussian", "laplacian", "e-dist",
#'        "g-dist", or "l-dist".
#' @param bw Bandwidth parameter for kernel distances. If NULL, it will be automatically determined.
#' @param expo Exponent parameter for Euclidean distance and polynomial kernel (default: 1).
#' @param scale_factor Scaling factor for automatic bandwidth calculation (default: 0.5).
#' @param group Optional vector specifying group membership for each variable in the response.
#'        The length should match the number of columns in the response matrix.
#' @param method Testing method (default: "asymptotic"). Options are:
#'        \itemize{
#'          \item \code{"asymptotic"}: Uses an asymptotic approximation for the p-value
#'          \item \code{"permutation"}: Uses permutation testing for the p-value
#'        }
#' @param n_perm Number of permutations to use when \code{method = "permutation"} (default: 999).
#' @param parallel Logical; whether to use parallel processing for permutation testing (default: FALSE).
#' @param num_cores Number of cores for parallel processing. If NULL, uses all available cores minus one.
#' @param contrasts Contrasts for factor variables (default: NULL, which uses contr.sum for unordered factors
#'        and contr.poly for ordered factors).
#' @param returnG Logical; whether to return the centered distance matrix G (default: FALSE).
#' @param is_distance Logical; whether the input Y is already a distance matrix (default: FALSE).
#'
#' @return A list containing:
#'   \item{aov.tab}{ANOVA table with F-statistic, R-squared values, and p-value}
#'   \item{df}{Degrees of freedom parameter for the asymptotic approximation (when applicable)}
#'   \item{call}{The function call}
#'   \item{G}{The centered distance matrix (if \code{returnG = TRUE})}
#'   \item{perm.F}{Vector of F-statistics from permutations (when \code{method = "permutation"})}
#'   \item{n_perm}{Number of permutations used (when \code{method = "permutation"})}
#'
#' @examples
#' # Example 1: Basic usage with multivariate response
#' set.seed(123)
#' n <- 100
#' X1 <- factor(sample(1:3, n, replace = TRUE))
#' X2 <- runif(n)
#' Z1 <- runif(n)
#' Z2 <- factor(sample(1:2, n, replace = TRUE))
#'
#' # Create multivariate response
#' Y <- matrix(rnorm(n*5), ncol = 5)
#' Y[,1] <- Y[,1] + as.numeric(X1) + Z1
#' Y[,2] <- Y[,2] + X2 + as.numeric(Z2)
#'
#' # Combine data
#' dataset <- data.frame(X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2)
#'
#' # Test dependence of Y on X1 and X2 after adjusting for Z1 and Z2
#' result1 <- dmanova(Y ~ X1 + X2 | Z1 + Z2, data = dataset)
#' print(result1$aov.tab)
#'
#' # Example 2: Using a simplified formula without covariates
#' # (This will use an intercept-only model for Z)
#' result2 <- dmanova(Y ~ X1 + X2, data = dataset)
#' print(result2$aov.tab)
#'
#' # Example 3: Using a pre-computed distance matrix
#' # Calculate a distance matrix first
#' D <- as.matrix(dist(Y))
#' # Use the distance matrix directly
#' result3 <- dmanova(D ~ X1 + X2 | Z1 + Z2, data = dataset, is_distance = TRUE)
#' print(result3$aov.tab)
#'
#' @references
#' Chen, J., & Zhang, X. (2022). D-MANOVA: Fast Distance-based Multivariate Analysis
#' of Variance for Large-scale Microbiome Association Studies. \emph{Bioinformatics}, 38, 286-288.
#'
#' Legendre, P., & Anderson, M. J. (1999). Distance-based redundancy analysis:
#' Testing multispecies responses in multifactorial ecological experiments.
#' \emph{Ecological Monographs}, 69(1), 1-24.
#'
#' McArdle, B. H., & Anderson, M. J. (2001). Fitting multivariate models to community data:
#' A comment on distance-based redundancy analysis. \emph{Ecology}, 82(1), 290-297.
#'
#' @export
dmanova <- function(formula, data = NULL, type = "euclidean", bw = NULL, expo = 1,
                    scale_factor = 0.5, group = NULL, method = c("asymptotic", "permutation"),
                    n_perm = 999, parallel = FALSE, num_cores = NULL, contrasts = NULL,
                    returnG = FALSE, is_distance = FALSE) {
  # Match the method argument
  method <- match.arg(method)

  # Default contrasts
  if(is.null(contrasts)) {
    contrasts <- list(unordered = "contr.sum", ordered = "contr.poly")
  }

  # Set contrasts
  op.c <- options()$contrasts
  options(contrasts = c(contrasts$unordered, contrasts$ordered))

  # Parse the formula
  formula_str <- as.character(formula)

  # Check if the formula contains the | operator
  has_covariates <- length(formula_str) >= 3 && grepl("\\|", formula_str[3])

  # Extract parts of the formula
  lhs <- formula[[2]]  # Y (response)

  if (has_covariates) {
    # Split the right side at the | symbol
    rhs_parts <- strsplit(formula_str[3], "\\|")[[1]]
    X_part <- trimws(rhs_parts[1])
    Z_part <- trimws(rhs_parts[2])
  } else {
    # If no | symbol, all the right side is X, and Z is just an intercept
    X_part <- formula_str[3]
    Z_part <- "1"  # Intercept only
  }

  # Evaluate response
  Y <- eval(lhs, data, parent.frame())

  # Create formulas with intercept
  Z_formula <- as.formula(paste("~", Z_part))
  full_formula <- as.formula(paste("~", paste(c(X_part, Z_part), collapse = "+")))

  # Create model frames and matrices
  Z_frame <- model.frame(Z_formula, data, drop.unused.levels = TRUE)
  Z <- model.matrix(Z_formula, Z_frame)

  full_frame <- model.frame(full_formula, data, drop.unused.levels = TRUE)
  XZ <- model.matrix(full_formula, full_frame)

  # Check for issues
  if (ncol(XZ) <= ncol(Z)) {
    warning("The number of parameters in X+Z is not greater than Z alone. This suggests X variables may already be included in Z.")
  }

  n <- nrow(XZ)

  # Handle distance matrix input
  if (is_distance) {
    # If Y is already a distance matrix
    if (!is.matrix(Y)) {
      stop("When is_distance=TRUE, Y must be a matrix")
    }

    if (nrow(Y) != n || ncol(Y) != n) {
      stop("Distance matrix Y must be square with dimensions matching the number of observations")
    }

    # Use Y directly as the distance matrix D
    D <- -Y
  } else {
    # Check dimensions of Y if it's not a distance matrix
    if (is.null(nrow(Y))) {
      # If Y is a vector, convert to matrix
      Y <- as.matrix(Y, ncol = 1)
    } else if (is.data.frame(Y)) {
      Y <- as.matrix(Y)
    }

    if (nrow(Y) != n) {
      stop("The number of observations in the response must match the predictors")
    }

    # Check group parameter if provided
    if (!is.null(group) && length(group) != ncol(Y)) {
      stop("Length of group vector must match number of columns in response")
    }

    # Calculate the distance matrix for Y
    D <- KDist_matrix(Y, type = type, bw = bw, expo = expo, scale_factor = scale_factor, group = group)

    if(type %in% c("euclidean", "e-dist", "l-dist", "g-dist")) {
      D <- -D
    }
  }

  # Center the distance matrix
  G <- matrix_v_center(D)

  # Calculate projection matrices
  XZi <- solve(t(XZ) %*% XZ)
  Zi <- solve(t(Z) %*% Z)
  HZ <- Z %*% Zi %*% t(Z)
  HXZ <- XZ %*% XZi %*% t(XZ)
  HX <- HXZ - HZ
  HIXZ <- diag(n) - HXZ
  HIX <- diag(n) - HX

  # Calculate degrees of freedom
  df1 <- ncol(XZ) - ncol(Z)
  df2 <- n - ncol(XZ)

  # Calculate sums of squares
  MSS <- sum(G * HX)
  RSS <- sum(G * HIXZ)
  TSS <- MSS + RSS

  # Calculate F-statistic
  f_stat <- (MSS / df1) / (RSS / df2)

  # Calculate p-value
  if(method == "asymptotic") {
    # Use asymptotic method
    G_tilde <- HIXZ %*% G %*% HIXZ
    dG <- diag(G_tilde)
    mu1 <- sum(dG) / df2
    mu2 <- (sum(G_tilde^2) - sum(dG^2)) / (df2^2 + sum(HIXZ^4) - 2 * sum(diag(HIXZ^2)))
    K <- mu1^2 / mu2
    pvalue <- pchisq(f_stat * K * df1, df = K * df1, lower.tail = FALSE)
  } else {
    # Permutation method
    K <- NA  # K is not applicable for permutation method

    # Determine number of cores
    if (parallel) {
      if (is.null(num_cores)) {
        num_cores <- parallel::detectCores() - 1
        if (is.na(num_cores) || num_cores < 1) num_cores <- 1
      }
    } else {
      num_cores <- 1  # sequential computation
    }

    # Create permutation indices in advance
    perm_indices <- replicate(n_perm, sample(1:n, n, replace = FALSE), simplify = FALSE)

    # Define the function to calculate F-statistic for permuted data
    calc_perm_f <- function(perm_idx) {
      # Permute the rows of G
      Gp <- G[perm_idx, perm_idx]
      # Calculate F-statistic on permuted data
      (sum(Gp * HX) / df1) / (sum(Gp * HIXZ) / df2)
    }

    # Use mclapply for parallel processing or lapply for sequential
    if (!requireNamespace("parallel", quietly = TRUE) || num_cores == 1) {
      f_perm_list <- lapply(perm_indices, calc_perm_f)
      f_perm <- unlist(f_perm_list)
    } else {
      f_perm_list <- parallel::mclapply(perm_indices, calc_perm_f, mc.cores = num_cores)
      f_perm <- unlist(f_perm_list)
    }

    # Calculate p-value as proportion of permuted F-statistics >= observed F-statistic
    pvalue <- (sum(f_perm >= f_stat) + 1) / (n_perm + 1)
  }

  # Create output table
  SumsOfSqs <- c(MSS, RSS, TSS)
  tab <- data.frame(
    Df = c(df1, df2, n - ncol(Z)),
    SumsOfSqs = SumsOfSqs,
    MeanSqs = c(MSS/df1, RSS/df2, NA),
    F.Model = c(f_stat, NA, NA),
    R2 = c(MSS/TSS, RSS/TSS, NA),
    "Pr(>F)" = c(pvalue, NA, NA)
  )

  rownames(tab) <- c("Model(Adjusted)", "Residuals", "Total")
  colnames(tab)[ncol(tab)] <- "Pr(>F)"
  class(tab) <- c("anova", class(tab))
  attr(tab, "heading") <- c("Distance-based MANOVA testing the dependence of response on predictors after adjusting for covariates\n")

  # Reset contrasts option
  options(contrasts = op.c)

  # Add information about the model for diagnostics
  model_info <- list(
    X_dims = ncol(XZ) - ncol(Z),
    Z_dims = ncol(Z),
    is_distance = is_distance,
    XZ_formula = deparse(full_formula),
    Z_formula = deparse(Z_formula)
  )

  # Build output list
  if (returnG) {
    out <- list(aov.tab = tab, df = K, G = G, model_info = model_info, call = match.call())
  } else {
    out <- list(aov.tab = tab, df = K, model_info = model_info, call = match.call())
  }

  # If using permutation method, add permutation results
  if (method == "permutation") {
    out$perm.F <- f_perm
    out$n_perm <- n_perm
  }

  return(out)
}

#' Fully Distance-Based Multivariate Analysis of Variance (D-MANOVA2)
#'
#' @description
#' Performs a fully distance-based multivariate analysis of variance for testing
#' the dependence of a multivariate response Y on predictors X after adjusting for
#' covariates Z. Unlike the standard dmanova, this function transforms all variables
#' (Y, X, and Z) into distance matrices and uses a ridge regression approach.
#'
#' @details
#' This implementation uses a fully distance-based approach where all variables
#' (response, predictors, and covariates) are transformed into distance or kernel matrices.
#' The test then assesses whether the distance structure of Y depends on that of X
#' after accounting for the distance structure of Z.
#'
#' The formula interface can be specified in two ways:
#' \itemize{
#'   \item \code{Y ~ X | Z} - Test the dependence of Y on X after adjusting for Z.
#'   \item \code{Y ~ X} - Test the dependence of Y on X with only an intercept term as covariate.
#' }
#'
#' By default, Z will contain an intercept. If Z is not specified in the formula (i.e., using
#' \code{Y ~ X} format), then Z will be set to an intercept-only model.
#'
#' If the input Y is already a distance matrix, it will be used directly, bypassing the distance calculation step.
#'
#' This function differs from the standard dmanova in that it uses ridge regression on
#' the distance matrices rather than projection matrices on the original variables, and
#' allows for separate distance/kernel configurations for response and covariates.
#'
#' @param formula A formula of the form Y ~ X | Z or Y ~ X, where Y is the response, X contains
#'        the predictors of interest, and Z contains the covariates to adjust for (if specified).
#' @param data An optional data frame containing the variables in the formula.
#' @param response_params A list of parameters for calculating the distance/kernel matrix of the response:
#'        \itemize{
#'          \item \code{type}: Type of distance or kernel to use (default: "euclidean").
#'                Options include "euclidean", "polynomial", "gaussian", "laplacian", "e-dist",
#'                "g-dist", or "l-dist".
#'          \item \code{bw}: Bandwidth parameter for kernel distances. If NULL, it will be automatically determined.
#'          \item \code{expo}: Exponent parameter for Euclidean distance and polynomial kernel (default: 2).
#'          \item \code{group}: Optional vector specifying group membership for each variable in the response.
#'                The length should match the number of columns in the response matrix.
#'        }
#' @param covariate_params A list of parameters for calculating the distance/kernel matrices of predictors and covariates:
#'        \itemize{
#'          \item \code{type}: Type of distance or kernel to use (default: "euclidean").
#'                Options include "euclidean", "polynomial", "gaussian", "laplacian", "e-dist",
#'                "g-dist", or "l-dist".
#'          \item \code{bw}: Bandwidth parameter for kernel distances. If NULL, it will be automatically determined.
#'          \item \code{expo}: Exponent parameter for Euclidean distance and polynomial kernel (default: 2).
#'        }
#'        Note: For model matrices (Z and XZ), no grouping structure is applied regardless of
#'        what parameters are provided.
#' @param scale_factor Scaling factor for automatic bandwidth calculation (default: 0.5).
#' @param epsilon Regularization parameter for ridge regression (default: 1e-3).
#' @param n_perm Number of permutations to use for the test (default: 999).
#' @param parallel Logical; whether to use parallel processing for permutation testing (default: FALSE).
#' @param num_cores Number of cores for parallel processing. If NULL, uses all available cores minus one.
#' @param contrasts Contrasts for factor variables (default: NULL, which uses contr.sum for unordered factors
#'        and contr.poly for ordered factors).
#' @param returnD Logical; whether to return the centered distance matrices (default: FALSE).
#' @param is_distance Logical; whether the input Y is already a distance matrix (default: FALSE).
#'
#' @return A list containing:
#'   \item{aov.tab}{ANOVA table with F-statistic, R-squared values, and p-value}
#'   \item{call}{The function call}
#'   \item{D_matrices}{The centered distance matrices (if \code{returnD = TRUE})}
#'   \item{perm.F}{Vector of F-statistics from permutations}
#'   \item{n_perm}{Number of permutations used}
#'   \item{model_info}{Information about the model fitted}
#'
#' @examples
#' # Example 1: Basic usage with multivariate response
#' set.seed(123)
#' n <- 100
#' X1 <- factor(sample(1:3, n, replace = TRUE))
#' X2 <- runif(n)
#' Z1 <- runif(n)
#' Z2 <- factor(sample(1:2, n, replace = TRUE))
#'
#' # Create multivariate response
#' Y <- matrix(rnorm(n*5), ncol = 5)
#' Y[,1] <- Y[,1] + as.numeric(X1) + Z1
#' Y[,2] <- Y[,2] + X2 + as.numeric(Z2)
#'
#' # Combine data
#' df <- data.frame(X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2)
#'
#' # Test dependence of Y on X1 and X2 after adjusting for Z1 and Z2
#' # Using default parameter lists
#' result1 <- dmanova2(Y ~ X1 + X2 | Z1 + Z2, data = df)
#' print(result1$aov.tab)
#'
#' # Example 2: Using a simplified formula without covariates
#' # (This will use an intercept-only model for Z)
#' result2 <- dmanova2(Y ~ X1 + X2, data = df)
#' print(result2$aov.tab)
#'
#' # Example 3: Using different parameters for response and covariates
#' result3 <- dmanova2(Y ~ X1 + X2 | Z1 + Z2, data = df,
#'                    response_params = list(type = "gaussian", expo = 1),
#'                    covariate_params = list(type = "euclidean", expo = 2))
#' print(result3$aov.tab)
#'
#' # Example 4: Using a pre-computed distance matrix
#' # Calculate a distance matrix first
#' D <- as.matrix(dist(Y))
#' # Use the distance matrix directly
#' result4 <- dmanova2(D ~ X1 + X2 | Z1 + Z2, data = df, is_distance = TRUE)
#' print(result4$aov.tab)
#'
#' # Example 5: Using group structure for multivariate response
#' # Suppose the 5 columns of Y belong to two groups
#' group_y <- c(1, 1, 1, 2, 2)
#' result5 <- dmanova2(Y ~ X1 + X2 | Z1 + Z2, data = df,
#'                    response_params = list(type = "e-dist", group = group_y))
#' print(result5$aov.tab)
#'
#' # Example 6: Using parallel processing for permutation tests
#' result6 <- dmanova2(Y ~ X1 + X2 | Z1 + Z2, data = df,
#'                    parallel = TRUE, n_perm = 499, num_cores = 2)
#' print(result6$aov.tab)
#'
#' @seealso
#' \code{\link{dmanova}} for the standard distance-based MANOVA approach.
#' \code{\link{KDist_matrix}} for details on distance/kernel matrix calculation.
#'
#' @export
dmanova2 <- function(formula, data = NULL,
                     response_params = list(type = "euclidean", bw = NULL, expo = 2, group = NULL),
                     covariate_params = list(type = "euclidean", bw = NULL, expo = 2),
                     scale_factor = 0.5, epsilon = 1e-3,
                     n_perm = 999, parallel = FALSE, num_cores = NULL, contrasts = NULL,
                     returnD = FALSE, is_distance = FALSE) {

  # Extract parameters for response
  type_rep <- response_params$type
  bw_rep <- response_params$bw
  expo_rep <- response_params$expo
  group_rep <- response_params$group

  # Extract parameters for covariates
  type_cov <- covariate_params$type
  bw_cov <- covariate_params$bw
  expo_cov <- covariate_params$expo

  # Default contrasts
  if(is.null(contrasts)) {
    contrasts <- list(unordered = "contr.sum", ordered = "contr.poly")
  }

  # Set contrasts
  op.c <- options()$contrasts
  options(contrasts = c(contrasts$unordered, contrasts$ordered))

  # Parse the formula
  formula_str <- as.character(formula)

  # Check if the formula contains the | operator
  has_covariates <- length(formula_str) >= 3 && grepl("\\|", formula_str[3])

  # Extract parts of the formula
  lhs <- formula[[2]]  # Y (response)

  if (has_covariates) {
    # Split the right side at the | symbol
    rhs_parts <- strsplit(formula_str[3], "\\|")[[1]]
    X_part <- trimws(rhs_parts[1])
    Z_part <- trimws(rhs_parts[2])
  } else {
    # If no | symbol, all the right side is X, and Z is just an intercept
    X_part <- formula_str[3]
    Z_part <- "1"  # Intercept only
  }

  # Evaluate response
  Y <- eval(lhs, data, parent.frame())

  # Create formulas with intercept
  Z_formula <- as.formula(paste("~", Z_part))
  X_formula <- as.formula(paste("~", X_part))
  full_formula <- as.formula(paste("~", paste(c(X_part, Z_part), collapse = "+")))

  # Create model frames and matrices
  Z_frame <- model.frame(Z_formula, data, drop.unused.levels = TRUE)
  Z <- model.matrix(Z_formula, Z_frame)

  X_frame <- model.frame(X_formula, data, drop.unused.levels = TRUE)
  X <- model.matrix(X_formula, X_frame)

  full_frame <- model.frame(full_formula, data, drop.unused.levels = TRUE)
  XZ <- model.matrix(full_formula, full_frame)

  # Check for issues
  if (ncol(XZ) <= ncol(Z)) {
    warning("The number of parameters in X+Z is not greater than Z alone. This suggests X variables may already be included in Z.")
  }

  n <- nrow(XZ)

  # Handle distance matrix input for Y
  if (is_distance) {
    # If Y is already a distance matrix
    if (!is.matrix(Y)) {
      stop("When is_distance=TRUE, Y must be a matrix")
    }

    if (nrow(Y) != n || ncol(Y) != n) {
      stop("Distance matrix Y must be square with dimensions matching the number of observations")
    }

    # Use Y directly as the distance matrix DY
    DY <- -Y
  } else {
    # Check dimensions of Y if it's not a distance matrix
    if (is.null(nrow(Y))) {
      # If Y is a vector, convert to matrix
      Y <- as.matrix(Y, ncol = 1)
    } else if (is.data.frame(Y)) {
      Y <- as.matrix(Y)
    }

    if (nrow(Y) != n) {
      stop("The number of observations in the response must match the predictors")
    }

    # Check group parameter if provided
    if (!is.null(group_rep) && length(group_rep) != ncol(Y)) {
      stop("Length of group_rep vector must match number of columns in response")
    }

    # Calculate the distance matrix for Y
    DY <- KDist_matrix(Y, type = type_rep, bw = bw_rep, expo = expo_rep,
                       scale_factor = scale_factor, group = group_rep)

    if(type_rep %in% c("euclidean", "e-dist", "l-dist", "g-dist")) {
      DY <- -DY
    }
  }

  # Calculate the distance matrices for Z and XZ
  DZ <- KDist_matrix(Z, type = type_cov, bw = bw_cov, expo = expo_cov,
                     scale_factor = scale_factor, group = NULL)
  DXZ <- KDist_matrix(XZ, type = type_cov, bw = bw_cov, expo = expo_cov,
                      scale_factor = scale_factor, group = NULL)

  if(type_cov %in% c("euclidean", "e-dist", "l-dist", "g-dist")) {
    DZ <- -DZ
    DXZ <- -DXZ
  }

  # Center the distance matrices
  I_n <- diag(n)
  DY <- matrix_v_center(DY)
  DZ <- matrix_v_center(DZ)
  DXZ <- matrix_v_center(DXZ)

  # Calculate the regularized matrices
  RXZ <- epsilon * solve(DXZ + epsilon * I_n)
  RZ <- epsilon * solve(DZ + epsilon * I_n)

  # Calculate the kernel matrices
  KXZ <- RXZ %*% DY %*% RXZ
  KZ <- RZ %*% DY %*% RZ

  # Calculate degrees of freedom
  df1 <- ncol(XZ) - ncol(Z)
  df2 <- n - ncol(XZ)

  # Calculate sums of squares
  MSS <- sum(diag(KZ - KXZ))
  RSS <- sum(diag(KXZ))
  TSS <- sum(diag(KZ))

  # Calculate F-statistic
  f_stat <- (MSS / df1) / (RSS / df2)

  # Calculate permutation p-value
  # Determine number of cores
  if (parallel) {
    if (is.null(num_cores)) {
      num_cores <- parallel::detectCores() - 1
      if (is.na(num_cores) || num_cores < 1) num_cores <- 1
    }
  } else {
    num_cores <- 1  # sequential computation
  }

  # Create permutation indices in advance
  perm_indices <- replicate(n_perm, sample(1:n, n, replace = FALSE), simplify = FALSE)

  # Define the function to calculate F-statistic for permuted data
  calc_perm_f <- function(perm_idx) {
    # Permute the rows of DY
    DYp <- DY[perm_idx, perm_idx]
    KXZp <- RXZ %*% DYp %*% RXZ
    KZp <- RZ %*% DYp %*% RZ
    # Calculate F-statistic on permuted data
    (sum(diag(KZp - KXZp)) / df1) / (sum(diag(KXZp)) / df2)
  }

  # Use mclapply for parallel processing or lapply for sequential
  if (!requireNamespace("parallel", quietly = TRUE) || num_cores == 1) {
    f_perm_list <- lapply(perm_indices, calc_perm_f)
    f_perm <- unlist(f_perm_list)
  } else {
    f_perm_list <- parallel::mclapply(perm_indices, calc_perm_f, mc.cores = num_cores)
    f_perm <- unlist(f_perm_list)
  }

  # Calculate p-value as proportion of permuted F-statistics >= observed F-statistic
  pvalue <- (sum(f_perm >= f_stat) + 1) / (n_perm + 1)

  # Create output table
  SumsOfSqs <- c(MSS, RSS, TSS)
  tab <- data.frame(
    Df = c(df1, df2, n - ncol(Z)),
    SumsOfSqs = SumsOfSqs,
    MeanSqs = c(MSS/df1, RSS/df2, NA),
    F.Model = c(f_stat, NA, NA),
    R2 = c(MSS/TSS, RSS/TSS, NA),
    "Pr(>F)" = c(pvalue, NA, NA)
  )

  rownames(tab) <- c("Model(Adjusted)", "Residuals", "Total")
  colnames(tab)[ncol(tab)] <- "Pr(>F)"
  class(tab) <- c("anova", class(tab))
  attr(tab, "heading") <- c("Fully Distance-based MANOVA testing the dependence of response on predictors after adjusting for covariates\n")

  # Reset contrasts option
  options(contrasts = op.c)

  # Add information about the model for diagnostics
  model_info <- list(
    X_dims = ncol(XZ) - ncol(Z),
    Z_dims = ncol(Z),
    is_distance = is_distance,
    XZ_formula = deparse(full_formula),
    Z_formula = deparse(Z_formula)
  )

  # Build output list
  if (returnD) {
    out <- list(
      aov.tab = tab,
      D_matrices = list(DY = DY, DZ = DZ, DXZ = DXZ),
      model_info = model_info,
      perm.F = f_perm,
      n_perm = n_perm,
      call = match.call()
    )
  } else {
    out <- list(
      aov.tab = tab,
      model_info = model_info,
      perm.F = f_perm,
      n_perm = n_perm,
      call = match.call()
    )
  }

  return(out)
}
