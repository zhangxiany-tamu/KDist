% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hypothesis_tests.R
\name{mdd_test}
\alias{mdd_test}
\title{Martingale Difference Divergence (MDD) Test}
\usage{
mdd_test(
  x,
  y,
  type = "euclidean",
  bw = NULL,
  expo = 1,
  scale_factor = 0.5,
  group_x = NULL,
  u_center = FALSE,
  n_boot = 1000,
  seed = NULL,
  num_cores = 1,
  is_distance = FALSE,
  boot_type = c("normal", "mammen", "rademacher")
)
}
\arguments{
\item{x}{Predictor dataset (matrix, data frame, or vector) or pre-computed distance matrix}

\item{y}{Response variable (vector)}

\item{type}{Type of kernel or distance to use for the predictor (default: "euclidean").
Options include "euclidean", "polynomial", "gaussian", "laplacian", "e-dist",
"g-dist", or "l-dist".}

\item{bw}{Bandwidth parameter for kernel distances. If NULL, it will be automatically determined.}

\item{expo}{Exponent parameter for euclidean distance and polynomial kernel (default: 1).}

\item{scale_factor}{Scaling factor for automatic bandwidth calculation (default: 0.5).}

\item{group_x}{Optional vector specifying group membership for each column in predictor data.
Used for group-wise distance calculations in "e-dist", "g-dist", or "l-dist".}

\item{u_center}{Logical; use U-centering instead of V-centering (default: FALSE).
U-centering provides an unbiased estimate by excluding diagonal terms.}

\item{n_boot}{Number of bootstrap iterations to use for the test (default: 1000).}

\item{seed}{Random seed for reproducibility (default: NULL).}

\item{num_cores}{Number of cores for parallel computing (default: 1).}

\item{is_distance}{Logical; whether the input x is already a distance matrix (default: FALSE).}

\item{boot_type}{Type of wild bootstrap distribution (default: "normal").
Options are "normal" (standard normal), "mammen" (Mammen's two-point distribution),
or "rademacher" (Rademacher distribution).}
}
\value{
An object of class "mdd_test" containing:
\item{statistic}{MDD test statistic value}
\item{p.value}{Wild bootstrap p-value}
\item{wild_bootstrap_values}{Vector of bootstrap statistic values}
\item{boot_type}{The bootstrap distribution type used}
\item{n}{Sample size}
\item{n_boot}{Number of bootstrap iterations}
}
\description{
Performs a wild bootstrap test for conditional mean independence based on the
Martingale Difference Divergence. This test assesses whether a response variable
is conditionally mean independent of a predictor, which is a weaker condition than
full conditional independence.
}
\details{
The MDD test checks whether E(Y|X) = E(Y), which is the conditional mean independence
condition. Unlike tests for full conditional independence, the MDD focuses specifically
on the conditional mean. The test is based on a distance/kernel representation that
transforms the original problem into one of checking the covariance between certain
distance matrices.

The wild bootstrap is used to calculate the p-value, with three distribution options:
\enumerate{
\item Normal: Standard normal distribution
\item Mammen: Mammen's two-point distribution, which has specific statistical properties
\item Rademacher: Takes values -1 or 1 with equal probability
}

When u_center = TRUE, the function uses an unbiased estimator based on U-statistics
that removes diagonal terms from the distance matrices, which can improve performance
for small sample sizes.
}
\examples{
# Example 1: Independent variables
set.seed(123)
x1 <- matrix(rnorm(100), ncol = 1)
y1 <- matrix(rnorm(100), ncol = 1)
test1 <- mdd_test(x1, y1, boot_type = "normal", n_boot = 200)
print(test1)
plot(test1)

# Example 2: Conditional mean dependence
x2 <- matrix(runif(100), ncol = 1)
y2 <- matrix(x2^2 + 0.5*rnorm(100), ncol = 1)  # Y depends on X^2
test2 <- mdd_test(x2, y2, boot_type = "mammen", n_boot = 200)
print(test2)

# Example 3: Parallel computation with multiple cores
\dontrun{
test3 <- mdd_test(x2, y2, boot_type = "rademacher", n_boot = 500, num_cores = 2)
print(test3)
}

}
\references{
Lee, C. E., Zhang, X., & Shao, X. (2020). Testing conditional mean independence for
functional data. \emph{Biometrika}, 107(2), 331-346.
}
\seealso{
\code{\link{mdd}} for calculating MDD without performing a test
\code{\link{hsic_test}} for general independence testing
}
