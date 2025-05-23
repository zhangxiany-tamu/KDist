% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/change_point.R
\name{kcpd_sbs}
\alias{kcpd_sbs}
\title{Seeded Binary Segmentation for Multiple Change Point Detection}
\usage{
kcpd_sbs(
  data,
  type = "e-dist",
  bw = NULL,
  expo = 1,
  scale_factor = 0.5,
  group = NULL,
  B = 199,
  alpha = 0.05,
  seeds = NULL,
  method = "asymptotic",
  decay = sqrt(2),
  unique_int = TRUE,
  bound = 20,
  num_cores = 1
)
}
\arguments{
\item{data}{A matrix or data frame with rows representing time points (observations) and
columns representing variables or features.}

\item{type}{Type of distance or kernel to use (default: "e-dist"). Options include:
\itemize{
\item "euclidean": Euclidean distance
\item "gaussian": Gaussian kernel
\item "laplacian": Laplacian kernel
\item "polynomial": Polynomial kernel
\item "e-dist": Euclidean-based aggregated distance (for high-dimensional data)
\item "g-dist": Gaussian kernel-based aggregated distance
\item "l-dist": Laplacian kernel-based aggregated distance
}}

\item{bw}{Bandwidth parameter for kernel calculations. If NULL (default), it will be automatically determined.}

\item{expo}{Exponent parameter for distance calculation (default: 1).}

\item{scale_factor}{Scaling factor for automatic bandwidth calculation (default: 0.5).}

\item{group}{Optional vector specifying group membership for each variable/column.
Used for group-wise distance calculations in "e-dist", "g-dist", or "l-dist".}

\item{B}{Number of permutations for p-value calculation (default: 199).}

\item{alpha}{Significance level for determining whether a change point is significant (default: 0.05).}

\item{seeds}{Random seed for reproducibility. If NULL, no seed is set.}

\item{method}{Testing method for p-value calculation: "permutation" or "asymptotic" (default: "asymptotic").
The asymptotic method is typically faster and recommended for larger datasets.}

\item{decay}{Decay factor for interval generation (default: sqrt(2)).
Controls the density of intervals; smaller values create more intervals.}

\item{unique_int}{Logical; whether to use only unique intervals (default: FALSE).}

\item{bound}{Minimum segment size to consider (default: 20).
Intervals shorter than this threshold are discarded. It should be greater than 5.}

\item{num_cores}{Number of cores for parallel processing (default: 1).
If NULL, will use one less than the available number of cores.}
}
\value{
A list containing:
\item{locations}{A vector of detected change point locations (indices)}
\item{pvalues}{A vector of p-values corresponding to each detected change point}
\item{cluster}{A vector of cluster assignments for each observation, where each segment
between change points is assigned a unique integer}
}
\description{
Detects multiple change points in time series or sequential data using the Seeded Binary Segmentation (SeedBS)
method combined with kernel-based metrics. This approach is particularly effective for high-dimensional data
and can detect various types of changes in distribution beyond just mean shifts.
}
\details{
This function implements a multiple change point detection approach that combines:
\enumerate{
\item The Seeded Binary Segmentation (SeedBS) methodology from Kovács et al. (2023), which uses
a deterministic set of intervals at multiple scales to search for change points
\item The kernel-based testing framework from Chakraborty & Zhang (2021), which enables
detection of general distributional changes even in high-dimensional settings
}

The algorithm works by:
\enumerate{
\item Generating a structured set of intervals using the \code{get_seeded_intervals} function
\item Applying the \code{kcpd_single} detection method to each interval independently
\item Selecting significant change points in order of p-value (smallest first)
\item Removing overlapping intervals to avoid redundant change points
\item Returning the final set of change points with their p-values and segment assignments
}

This method is more computationally efficient than exhaustive search approaches while
maintaining statistical power and accuracy, especially in high-dimensional settings.

This version includes Bonferroni correction to control the family-wise error rate
when testing multiple intervals.
}
\examples{
# Example 1: High-dimensional data with change in distribution
set.seed(123)
n <- 100  # Total sample size
p <- 200  # Number of dimensions (high-dimensional setting)
cp <- 50  # Change point location

# Generate high-dimensional data
# Before change: standard normal distribution
X1 <- matrix(rnorm(cp * p), nrow = cp, ncol = p)

# After change: centered standard exponential distribution
X2 <- matrix(rexp(n = (n-cp) * p) - 1, nrow = (n-cp), ncol = p)

# Combine data
X <- rbind(X1, X2)

# Detect change points
result <- kcpd_sbs(X, type = "e-dist")

print(result$locations)  # Should be close to 50
print(result$pvalues)    # Should be < 0.05
print(table(result$cluster))  # Should show 2 segments

# Visualize first few dimensions (for demonstration only)
matplot(1:n, X[, 1:5], type = "l", col = 1:5, lty = 1,
        xlab = "Time", ylab = "Value", main = "First 5 Dimensions")
if(length(result$locations) > 0) {
  abline(v = result$locations, col = "red", lwd = 2, lty = 2)
}

# Example 2: Multiple change points in high-dimensional data
set.seed(456)
# Generate high-dimensional data with changes at t=40 and t=80
n <- 120
p <- 50  # High-dimensional

# Three segments with different distributions
X1 <- matrix(rnorm(40 * p), ncol = p)
X2 <- matrix(rnorm(40 * p, mean = 0.5), ncol = p)

# Third segment has different covariance structure
Sigma <- matrix(0.5, p, p)
diag(Sigma) <- 1
X3 <- matrix(rnorm(40 * p), ncol = p)
X3 <- X3 \%*\% chol(Sigma)  # Apply covariance structure

# Combined data
X <- rbind(X1, X2, X3)

# Detect multiple change points using e-dist (efficient for high-dimensional data)
result <- kcpd_sbs(X, type = "e-dist", decay = 1.5)
print(result$locations)  # Should be close to 40 and 80
print(result$pvalues)
print(table(result$cluster))

# Simply visualize the results with base R
mean_vals <- rowMeans(X)  # Just for visualization
plot(1:n, mean_vals, pch = 16, col = result$cluster,
     xlab = "Time", ylab = "Mean Value",
     main = "Multiple Change Point Detection")
# Add vertical lines at change points
if(length(result$locations) > 0) {
  abline(v = result$locations, lty = 2, col = "red")
}

}
