% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KDist.R
\name{hsic_cor}
\alias{hsic_cor}
\title{HSIC-based Correlation Coefficient}
\usage{
hsic_cor(
  x,
  y,
  type = "gaussian",
  bw = NULL,
  expo = 1,
  scale_factor = 0.5,
  group = NULL,
  u_center = FALSE,
  is_distance = FALSE
)
}
\arguments{
\item{x}{First dataset or distance matrix}

\item{y}{Second dataset or distance matrix}

\item{type}{Type of kernel or distance}

\item{bw}{Bandwidth parameter}

\item{expo}{Exponent parameter}

\item{scale_factor}{Scaling factor for bandwidth}

\item{group}{Optional list of length 2 specifying group membership for each column in x and y.
The first element of the list should specify grouping for columns in x, and the
second element should specify grouping for columns in y. Used for group-wise
distance calculations in "e-dist", "g-dist", or "l-dist".}

\item{u_center}{Logical; use U-centering instead of V-centering}

\item{is_distance}{Logical; whether input matrices are already distance matrices}
}
\value{
HSIC correlation coefficient
}
\description{
Calculates a correlation coefficient based on the Hilbert-Schmidt Independence Criterion,
providing a normalized measure of dependence between variables. This is proposed in
Gretton et al. (2007) as the kernel version of the distance correlation later developed
by Székely et al. (2007).
}
\examples{
# Example 1: Independent variables
set.seed(123)
x1 <- matrix(rnorm(100), ncol = 1)
y1 <- matrix(rnorm(100), ncol = 1)
hsic_corr1 <- hsic_cor(x1, y1, type = "gaussian")
# Should be close to 0
print(hsic_corr1)

# Example 2: Linear relationship
x2 <- matrix(1:100, ncol = 1)
y2 <- matrix(2*x2 + rnorm(100, sd = 5), ncol = 1)
hsic_corr2 <- hsic_cor(x2, y2, type = "gaussian")
# Should be close to 1
print(hsic_corr2)

# Example 3: Non-linear relationship
x3 <- matrix(seq(-3, 3, length.out = 100), ncol = 1)
y3 <- matrix(x3^2 + rnorm(100, sd = 0.5), ncol = 1)
hsic_corr3 <- hsic_cor(x3, y3, type = "gaussian")
# Should detect the non-linear relationship
print(hsic_corr3)

# Compare with Pearson correlation
cor_pearson3 <- cor(x3, y3)
# Pearson correlation may miss the non-linear relationship
print(paste("HSIC correlation:", hsic_corr3, "Pearson correlation:", cor_pearson3))

# Example 4: Using different kernel types
hsic_gaussian <- hsic_cor(x3, y3, type = "gaussian")
hsic_laplacian <- hsic_cor(x3, y3, type = "laplacian")
hsic_euclidean <- hsic_cor(x3, y3, type = "euclidean")
print(c(gaussian = hsic_gaussian, laplacian = hsic_laplacian, euclidean = hsic_euclidean))

# Example 5: Using pre-computed distance matrices
Dx <- KDist_matrix(x3, type = "gaussian")
Dy <- KDist_matrix(y3, type = "gaussian")
hsic_precomp <- hsic_cor(Dx, Dy, type = "gaussian", is_distance = TRUE)
print(hsic_precomp)

# Example 6: Perfect correlation (identical variables)
hsic_perfect <- hsic_cor(x2, x2, type = "gaussian")
# Should be exactly 1
print(hsic_perfect)

# Example 7: Sinusoidal relationship
t <- seq(0, 4*pi, length.out = 100)
x7 <- matrix(sin(t), ncol = 1)
y7 <- matrix(cos(t), ncol = 1)
hsic_sin <- hsic_cor(x7, y7, type = "gaussian")
cor_pearson7 <- cor(x7, y7)
# HSIC should detect this non-linear relationship
print(paste("HSIC correlation:", hsic_sin, "Pearson correlation:", cor_pearson7))

# Example 8: Using grouped variables with hsic_cor
# Create sample data
x <- matrix(rnorm(100*4), ncol = 4)  # 4 variables in x
y <- matrix(rnorm(100*3), ncol = 3)  # 3 variables in y

# Create some dependency between the datasets
y[, 1:2] <- y[, 1:2] + 0.8 * x[, 1:2]

# Define group structure for both datasets
x_groups <- c(1, 1, 2, 2)    # First 2 vars in group 1, last 2 in group 2
y_groups <- c(1, 1, 2)       # First 2 vars in group 1, last 1 in group 2

# Combine into a list as required for grouped analysis
group_list <- list(x_groups, y_groups)

# Calculate HSIC correlation with the group structure
cor_grouped <- hsic_cor(x, y, type = "e-dist", group = group_list)
print(paste("HSIC correlation with grouped variables:", cor_grouped))

# Compare with standard (non-grouped) HSIC correlation
cor_standard <- hsic_cor(x, y, type = "gaussian")
print(paste("Standard HSIC correlation:", cor_standard))

}
\references{
Gretton, A., Fukumizu, K., Teo, C., Song, L., Schölkopf, B., & Smola, A. (2007).
A kernel statistical test of independence. \emph{Advances in neural information processing systems, 20}.

Székely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing dependence by
correlation of distances. \emph{The Annals of Statistics, 35}(6), 2769-2794.
}
