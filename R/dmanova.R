#' Kernel and Distance-Based MANOVA
#'
#' Performs kernel/distance-based multivariate analysis of variance.
#' Only the response variable will be transformed into a kernel or distance matrix.
#'
#' @param formula A formula with the response on the left and predictors on the right
#' @param data An optional data frame containing the variables in the formula
#' @param type Type of distance to use (default: "euclidean")
#' @param bw Bandwidth parameter (for kernel distances)
#' @param expo Exponent parameter
#' @param scale_factor Scaling factor for bandwidth
#' @param group Optional grouping of variables
#' @param method Testing method: "asymptotic" or "permutation"
#' @param n_perm Number of permutations when method="permutation"
#' @param parallel Logical; whether to use parallel processing
#' @param n_cores Number of cores for parallel processing
#' @param contr.unordered Contrast method for unordered factors
#' @param contr.ordered Contrast method for ordered factors
#' @param returnG Logical; whether to return the G matrix
#'
#' @return A list containing ANOVA table and other test results
#'
#' @export
dmanova <- function(formula, data = NULL, type = "euclidean", bw = NULL, expo = 1, 
                    scale_factor = 0.5, group = NULL, method = c("asymptotic", "permutation"), 
                    n_perm = 999, parallel = FALSE, n_cores = NULL, contr.unordered = "contr.sum", 
                    contr.ordered = "contr.poly", returnG = FALSE) {
  # Match the method argument
  method <- match.arg(method)
  
  # Parse the formula
  Terms <- terms(formula, data = data)
  lhs <- formula[[2]]
  lhs <- eval(lhs, data, parent.frame())
  
  formula[[2]] <- NULL
  rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE)
  op.c <- options()$contrasts
  options(contrasts = c(contr.unordered, contr.ordered))
  rhs <- model.matrix(formula, rhs.frame)
  options(contrasts = op.c)
  grps <- attr(rhs, "assign")
  
  # The first group includes the intercept, nterms count the intercept
  u.grps <- unique(grps)
  nterms <- length(u.grps)
  
  Z <- rhs[, grps %in% u.grps[1:(nterms-1)], drop=FALSE]
  XZ <- rhs[, grps %in% u.grps[1:nterms], drop=FALSE]
  
  n <- nrow(XZ)
  D <- KDist_matrix(lhs, type = type, bw = bw, expo = expo, scale_factor = scale_factor, group = group)
  
  if(type %in% c("euclidean", "e-dist", "l-dist", "g-dist")) D <- -D
  G <- matrix_v_center(D)
  XZi <- solve(t(XZ) %*% XZ)
  Zi <- solve(t(Z) %*% (Z))
  HZ <- Z %*% Zi %*% t(Z)
  HXZ <- XZ %*% XZi %*% t(XZ)
  HX <- HXZ - HZ
  HIXZ <- diag(n) - HXZ
  HIX <- diag(n) - HX
  
  df1 <- ncol(XZ) - ncol(Z)
  df2 <- n - ncol(XZ)
  
  MSS <- sum(G * HX)
  RSS <- sum(G * HIXZ)
  TSS <- sum(diag(G))
  
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
      if (is.null(n_cores)) {
        n_cores <- parallel::detectCores() - 1
        if (is.na(n_cores) || n_cores < 1) n_cores <- 1
      }
    } else {
      n_cores <- 1  # sequential computation
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
    if (!requireNamespace("parallel", quietly = TRUE) || n_cores == 1) {
      f_perm_list <- lapply(perm_indices, calc_perm_f)
      f_perm <- unlist(f_perm_list)
    } else {
      f_perm_list <- parallel::mclapply(perm_indices, calc_perm_f, mc.cores = n_cores)
      f_perm <- unlist(f_perm_list)
    }
    
    # Calculate p-value as proportion of permuted F-statistics >= observed F-statistic
    pvalue <- (sum(f_perm >= f_stat) + 1) / (n_perm + 1)
  }
  
  # Create output table
  SumsOfSqs <- c(MSS, RSS, TSS)
  tab <- data.frame(
    Df = c(df1, df2, n - 1), 
    SumsOfSqs = SumsOfSqs, 
    MeanSqs = c(MSS/df1, RSS/df2, NA), 
    F.Model = c(f_stat, NA, NA), 
    R2 = c(MSS/TSS, RSS/TSS, NA), 
    "Pr(>F)" = c(pvalue, NA, NA)
  )
  
  rownames(tab) <- c("Model(Adjusted)", "Residuals", "Total")
  colnames(tab)[ncol(tab)] <- "Pr(>F)"
  class(tab) <- c("anova", class(tab))
  attr(tab, "heading") <- c("F stat and P value of the last term is adjusted by preceding terms!\n")
  
  # Build output list
  if (returnG) {
    out <- list(aov.tab = tab, df = K, G = G, call = match.call())
  } else {
    out <- list(aov.tab = tab, df = K, call = match.call())
  }
  
  # If using permutation method, add permutation results
  if (method == "permutation") {
    out$perm.F <- f_perm
    out$n_perm <- n_perm
  }
  
  return(out)
}

#' Kernel and Distance-Based MANOVA
#'
#' An alternative implementation of kernel and distance-based multivariate analysis of variance.
#' All the variables will be transformed into kernel or distance matrices.
#'
#' @param formula A formula with the response on the left and predictors on the right
#' @param data An optional data frame containing the variables in the formula
#' @param type_rep Type of distance for response
#' @param type_cov Type of distance for covariates
#' @param expo_rep Exponent parameter for response
#' @param expo_cov Exponent parameter for covariates
#' @param n_perm Number of permutations
#' @param parallel Logical; whether to use parallel processing
#' @param n_cores Number of cores for parallel processing
#' @param contr.unordered Contrast method for unordered factors
#' @param contr.ordered Contrast method for ordered factors
#'
#' @return A list containing test statistics and p-value
#'
#' @export
dmanova2 <- function(formula, data = NULL, type_rep = "euclidean", type_cov = "euclidean",  
                    expo_rep = 1, expo_cov = 1, n_perm = 299, parallel = FALSE, 
                    n_cores = NULL, contr.unordered = "contr.sum", 
                    contr.ordered = "contr.poly") {
  # Parse the formula
  Terms <- terms(formula, data = data)
  lhs <- formula[[2]]
  lhs <- eval(lhs, data, parent.frame())
  
  formula[[2]] <- NULL
  rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE)
  op.c <- options()$contrasts
  options(contrasts = c(contr.unordered, contr.ordered))
  rhs <- model.matrix(formula, rhs.frame)
  options(contrasts = op.c)
  grps <- attr(rhs, "assign")
  
  # The first group includes the intercept, nterms count the intercept
  u.grps <- unique(grps)
  nterms <- length(u.grps)
  
  Z <- rhs[, grps %in% u.grps[1:(nterms-1)], drop=FALSE]
  XZ <- rhs[, grps %in% u.grps[1:nterms], drop=FALSE]
  
  n <- nrow(XZ)
  DY <- KDist_matrix(lhs, type = type_rep, bw = NULL, expo = expo_rep, scale_factor = 0.5, group = NULL)
  DZ <- KDist_matrix(Z, type = type_cov, bw = NULL, expo = expo_cov, scale_factor = 0.5, group = NULL)
  DXZ <- KDist_matrix(XZ, type = type_cov, bw = NULL, expo = expo_cov, scale_factor = 0.5, group = NULL)
  
  if(type_rep %in% c("euclidean", "e-dist", "l-dist", "g-dist")) {
    DY <- -DY
  }
  if(type_cov %in% c("euclidean", "e-dist", "l-dist", "g-dist")) {
    DZ <- -DZ
    DXZ <- -DXZ
  }
  
  epsilon <- 1e-3
  I_n <- diag(n)
  DY <- matrix_v_center(DY)
  DZ <- matrix_v_center(DZ)
  DXZ <- matrix_v_center(DXZ)
  
  RXZ <- epsilon*solve(DXZ + epsilon * I_n)
  RZ <- epsilon*solve(DZ + epsilon * I_n)
  KXZ <- RXZ %*% DY %*% RXZ
  KZ <- RZ %*% DY %*% RZ
  
  df1 <- ncol(XZ) - ncol(Z)
  df2 <- n - ncol(XZ)
  f_stat <- (sum(diag(KZ - KXZ)) / df1) / (sum(diag(KXZ)) / df2)
  
  # Calculate p-value
    # Determine number of cores
    if (parallel) {
      if (is.null(n_cores)) {
        n_cores <- parallel::detectCores() - 1
        if (is.na(n_cores) || n_cores < 1) n_cores <- 1
      }
    } else {
      n_cores <- 1  # sequential computation
    }
    
    # Create permutation indices in advance
    perm_indices <- replicate(n_perm, sample(1:n, n, replace = FALSE), simplify = FALSE)
    
    # Define the function to calculate F-statistic for permuted data
    calc_perm_f <- function(perm_idx) {
      # Permute the rows of G
      DYp <- DY[perm_idx, perm_idx]
      KXZp <- RXZ %*% DYp %*% RXZ
      KZp <- RZ %*% DYp %*% RZ
      # Calculate F-statistic on permuted data
      (sum(diag(KZp - KXZp)) / df1) / (sum(diag(KXZp)) / df2)
    }
    
    # Use mclapply for parallel processing or lapply for sequential
    if (!requireNamespace("parallel", quietly = TRUE) || n_cores == 1) {
      f_perm_list <- lapply(perm_indices, calc_perm_f)
      f_perm <- unlist(f_perm_list)
    } else {
      f_perm_list <- parallel::mclapply(perm_indices, calc_perm_f, mc.cores = n_cores)
      f_perm <- unlist(f_perm_list)
    }
    
    # Calculate p-value as proportion of permuted F-statistics >= observed F-statistic
    pvalue <- (sum(f_perm >= f_stat) + 1) / (n_perm + 1)
    
    # Build output list
    out <- list(statistics = f_stat, pvalue = pvalue)
    return(out)
}