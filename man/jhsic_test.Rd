% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hypothesis_tests.R
\name{jhsic_test}
\alias{jhsic_test}
\title{Joint HSIC and dCov Independence Test}
\usage{
jhsic_test(
  x,
  cc = 1,
  type = "gaussian",
  stat_type = "V",
  bw = NULL,
  expo = 1,
  scale_factor = 0.5,
  group = NULL,
  n_perm = 1000,
  seed = NULL,
  num_cores = 1
)
}
\arguments{
\item{x}{List of matrices or vectors, where each element corresponds to a random variable}

\item{cc}{Constant parameter (default: 1)}

\item{type}{Type of kernel or distance to use (default: "gaussian").
Options include "euclidean", "polynomial", "gaussian", "laplacian", "e-dist", "g-dist", or "l-dist".}

\item{stat_type}{Type of statistic ("V", "U", "US", or "UR")
\itemize{
\item "V": V-statistic
\item "U": U-statistic (partially unbiased)
\item "US": Scale-free U-statistic
\item "UR": Rank-based U-statistic (robust to monotonic transformations)
}}

\item{bw}{Bandwidth parameter for kernel functions. If NULL, it will be automatically determined.}

\item{expo}{Exponent parameter for euclidean distance and polynomial kernel (default: 1).}

\item{scale_factor}{Scaling factor for automatic bandwidth calculation (default: 0.5).}

\item{group}{Optional list specifying group membership for each column in the corresponding input list elements.
The length of the group list should match the length of the input list x, where each element
contains the grouping information for the respective element in x. Used for group-wise
distance calculations in "e-dist", "g-dist", or "l-dist".}

\item{n_perm}{Number of permutations to use for the test (default: 1000).}

\item{seed}{Random seed for reproducibility (default: NULL).}

\item{num_cores}{Number of cores for parallel computing (default: 1).}
}
\value{
An object of class "jhsic_test" containing:
\item{statistic}{JHSIC test statistic value}
\item{p.value}{Permutation-based p-value}
\item{permutation_values}{Vector of JHSIC values from permutations}
\item{stat_description}{Description of the statistic type used}
}
\description{
Performs a permutation test based on the joint Hilbert-Schmidt Independence Criterion and joint
covariance to assess mutual independence among multiple random variables. This implementation follows the methodology
described in Chakraborty & Zhang (2019), which extends traditional pairwise independence measures
to a joint measure capable of detecting complex dependence structures among multiple variables.
}
\details{
The Joint Hilbert-Schmidt Independence Criterion (JHSIC) extends the pairwise HSIC to detect complex
dependencies among multiple variables. When type = "euclidean", it becomes the joint distance covariance (JdCov)
proposed in Chakraborty & Zhang (2019), Unlike pairwise independence tests, JHSIC and JdCov can detect
higher-order dependencies that might be missed when only testing pairs of variables.

The null hypothesis is that all random variables are mutually independent. The alternative hypothesis
is that there exists some form of dependence structure among the variables. The p-value is calculated
using a permutation test, where all variables except the first one are randomly permuted to create
the null distribution.

Different statistic types offer various properties:
\itemize{
\item V-statistic: Standard, biased estimator
\item U-statistic: Partially unbiased estimator that removes diagonal terms
\item US-statistic: Scale-free U-statistic for normalized measure
\item UR-statistic: Rank-based U-statistic, robust to monotonic transformations
}
}
\examples{
# Example 1: Testing for independence among truly independent variables
set.seed(123)
x1 <- matrix(rnorm(100), ncol = 1)
x2 <- matrix(rnorm(100), ncol = 1)
x3 <- matrix(rnorm(100), ncol = 1)
test1 <- jhsic_test(list(x1, x2, x3), stat_type = "V", n_perm = 200)
print(test1)
plot(test1)

# Example 2: Testing variables with pairwise independence but joint dependence
set.seed(456)
u <- runif(100, -pi, pi)
x1 <- matrix(sin(u), ncol = 1)
x2 <- matrix(cos(u), ncol = 1)
x3 <- matrix(sin(u) * cos(u), ncol = 1)
test2 <- jhsic_test(list(x1, x2, x3), stat_type = "V", n_perm = 200)
print(test2)

# Example 3: Comparing different statistic types
test_v <- jhsic_test(list(x1, x2, x3), stat_type = "V", n_perm = 200)
test_u <- jhsic_test(list(x1, x2, x3), stat_type = "U", n_perm = 200)
test_us <- jhsic_test(list(x1, x2, x3), stat_type = "US", n_perm = 200)
test_ur <- jhsic_test(list(x1, x2, x3), stat_type = "UR", n_perm = 200)

print(c(V = test_v$p.value, U = test_u$p.value,
        US = test_us$p.value, UR = test_ur$p.value))

# Example 4: Using multivariate data
x_multi1 <- matrix(rnorm(200), ncol = 2)
x_multi2 <- matrix(rnorm(200), ncol = 2)
x_multi3 <- matrix(rnorm(200), ncol = 2)
test_multi <- jhsic_test(list(x_multi1, x_multi2, x_multi3), n_perm = 200)
print(test_multi)

# Example 5: Testing with non-linear relationships
set.seed(789)
x <- matrix(runif(100, -3, 3), ncol = 1)
y <- matrix(x^2 + rnorm(100, sd = 0.5), ncol = 1)
z <- matrix(exp(x) + rnorm(100, sd = 0.5), ncol = 1)
test_nonlin <- jhsic_test(list(x, y, z), type = "gaussian", n_perm = 200)
print(test_nonlin)

# Example 6: Using grouped variables with jhsic_test
# Create sample data with 3 datasets
set.seed(123)
n <- 100
x1 <- matrix(rnorm(n*3), ncol = 3)  # 3 variables in x1
x2 <- matrix(rnorm(n*4), ncol = 4)  # 4 variables in x2
x3 <- matrix(rnorm(n*2), ncol = 2)  # 2 variables in x3

# Create complex dependency structure via a common factor
z <- runif(n, -pi, pi)
x1[, 1] <- 2*sin(z) + 0.1*rnorm(n)
x2[, 1] <- 2*cos(z) + 0.1*rnorm(n)
x3[, 1] <- sin(z)*cos(z) + 0.1*rnorm(n)

# Define group structure for each dataset
groups1 <- c(1, 2, 2)        # First var in group 1, others in group 2
groups2 <- c(1, 2, 2, 2)     # First var in group 1, others in group 2
groups3 <- c(1, 2)           # Each var in its own group

# Combine into a list matching the input list structure
group_list <- list(groups1, groups2, groups3)

# Perform JHSIC test with the group structure
test_grouped <- jhsic_test(list(x1, x2, x3), type = "e-dist",
                          group = group_list, n_perm = 200)
print(test_grouped)

# Compare with standard JHSIC test using V-statistic
test_standard <- jhsic_test(list(x1, x2, x3), type = "gaussian",
                           stat_type = "V", n_perm = 200)
print(test_standard)

# Example 7: Using parallel computing for larger tests
\dontrun{
if(parallel::detectCores() > 1) {
  test_parallel <- jhsic_test(list(x_multi1, x_multi2, x_multi3),
                             n_perm = 1000, num_cores = 2)
  print(test_parallel)
}
}

}
\references{
Chakraborty, S., & Zhang, X. (2019). Distance Metrics for Measuring Joint Dependence with
Application to Causal Inference. \emph{Journal of the American Statistical Association}, 114, 1638-1650.
}
\seealso{
\code{\link{hsic_test}} for pairwise independence test
\code{\link{dhsic_test}} for an alternative joint independence test
\code{\link{mmd_test}} for two-sample testing
}
