% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KDist.R
\name{dcor}
\alias{dcor}
\title{Distance Correlation Between Two Samples}
\usage{
dcor(x, y, a = 1, u_center = FALSE, is_distance = FALSE)
}
\arguments{
\item{x}{First dataset or distance matrix}

\item{y}{Second dataset or distance matrix}

\item{a}{Exponent parameter for Euclidean distance, default is 1}

\item{u_center}{Logical; use U-centering instead of V-centering (default: FALSE)}

\item{is_distance}{Logical; whether input matrices are already distance matrices (default: FALSE)}
}
\value{
The squared distance correlation value between the two samples (between 0 and 1)
}
\description{
Calculates the squared distance correlation between two multivariate samples, which
provides a scaled measure of dependence between random vectors.
}
\details{
The squared distance correlation between two samples is defined as:

\deqn{dcor^2(X, Y) = dcov^2(X, Y) / \sqrt{dcov^2(X, X) dcov^2(Y, Y)}}

where \eqn{dcov^2} is the squared distance covariance. Distance correlation normalizes
the distance covariance to produce a measure between 0 and 1.

Distance correlation has these important properties:
\itemize{
\item \eqn{0 \leq dcor^2(X, Y) \leq 1}
\item \eqn{dcor^2(X, Y) = 0} if and only if \eqn{X} and \eqn{Y} are independent
\item \eqn{dcor^2(X, Y) = 1} only when \eqn{X} and \eqn{Y} are related by a linear transformation
\item It can detect both linear and non-linear dependence
\item It is applicable to vectors of arbitrary and different dimensions
}

When u_center = FALSE (default), the function computes the V-statistic version.
When u_center = TRUE, it computes the U-statistic version that's unbiased.
}
\examples{
# Example 1: Distance correlation between independent samples
set.seed(123)
x <- matrix(rnorm(100), ncol = 1)
y <- matrix(rnorm(100), ncol = 1)
dcor_indep <- dcor(x, y)
print(paste("Distance correlation (independent):", round(dcor_indep, 6)))

# Example 2: Distance correlation with linear relationship
x <- matrix(runif(100), ncol = 1)
y <- matrix(2*x + 0.1*rnorm(100), ncol = 1)  # Linear relationship with noise
dcor_linear <- dcor(x, y)
print(paste("Distance correlation (linear):", round(dcor_linear, 6)))

# Example 3: Distance correlation with non-linear relationship
x <- matrix(runif(100, -3, 3), ncol = 1)
y <- matrix(sin(x) + 0.1*rnorm(100), ncol = 1)  # Sinusoidal relationship
dcor_nonlinear <- dcor(x, y)
print(paste("Distance correlation (non-linear):", round(dcor_nonlinear, 6)))

# Compare with Pearson correlation
pearson_cor <- cor(x, y)[1,1]
print(paste("Distance correlation:", round(dcor_nonlinear, 6),
            "Pearson correlation:", round(pearson_cor, 6)))

# Example 4: Comparing V-statistic and U-statistic estimators
dcor_v <- dcor(x, y, u_center = FALSE)  # Default, V-statistic
dcor_u <- dcor(x, y, u_center = TRUE)   # U-statistic (unbiased)
print(c("V-statistic" = dcor_v, "U-statistic" = dcor_u))

# Example 5: Using multivariate data with different dimensions
x_multi <- matrix(rnorm(100*3), ncol = 3)  # 3-dimensional
y_multi <- matrix(rnorm(100*2), ncol = 2)  # 2-dimensional
# Make them dependent
y_multi[,1] <- x_multi[,1] + rnorm(100, sd = 0.1)
y_multi[,2] <- x_multi[,2]^2 + rnorm(100, sd = 0.1)
dcor_multi <- dcor(x_multi, y_multi)
print(paste("Multivariate distance correlation:", round(dcor_multi, 6)))

# Example 6: Using pre-computed distance matrices
Dx <- as.matrix(dist(x_multi))
Dy <- as.matrix(dist(y_multi))
dcor_precomp <- dcor(Dx, Dy, is_distance = TRUE)
print(paste("Using pre-computed matrices:", round(dcor_precomp, 6)))

# Example 7: Perfect linear relationship
x_perfect <- matrix(1:100, ncol = 1)
y_perfect <- matrix(3*x_perfect + 5, ncol = 1)  # Perfect linear relationship
dcor_perfect <- dcor(x_perfect, y_perfect)
print(paste("Perfect linear relationship:", round(dcor_perfect, 6)))

}
\references{
Székely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing dependence by
correlation of distances. \emph{The Annals of Statistics}, 35(6), 2769-2794.
}
\seealso{
\code{\link{dcov}} for distance covariance,
\code{\link{hsic_cor}} for kernel-based correlation,
\code{\link{cor}} for Pearson correlation
}
