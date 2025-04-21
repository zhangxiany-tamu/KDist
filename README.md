# KDist: A Collection of Kernel and Distance Methods for Statistical Inference

[![R-CMD-check](https://github.com/zhangxiany-tamu/KDist/workflows/R-CMD-check/badge.svg)](https://github.com/zhangxiany-tamu/KDist/actions)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/KDist)](https://cran.r-project.org/package=KDist)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

KDist is a comprehensive R package that implements a wide range of kernel and distance-based methods for statistical inference. These non-parametric methods are particularly useful for:

- Two-sample testing
- Testing independence between random vectors
- Multivariate hypothesis testing
- High-dimensional data analysis
- Detecting structural changes in sequential data
- Conditional independence and distribution testing

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("zhangxiany-tamu/KDist")
```

## Features

KDist offers a versatile toolkit with the following key components:

### Distance and Kernel Metrics

- Multiple distance options: Euclidean, Kernel-induced distances
- Various kernel options: Gaussian, Laplacian, Polynomial
- Automatic bandwidth selection
- Group-based distance calculations

### Testing Methods

- **Two-Sample Testing**: Test for differences between multivariate distributions using Maximum Mean Discrepancy (MMD) or Energy Distance (ED)
- **Independence Testing**: Hilbert-Schmidt Independence Criterion (HSIC), Distance Covariance (dCov), d-variable HSIC (dHSIC), Joint HSIC (jHSIC)
- **High-Dimensional Testing**: Methods developed for high-dimensional data settings
- **Multivariate Analysis**: Kernel and distance-based MANOVA
- **Conditional Testing**: Tests for conditional independence and conditional distribution differences
- **Change-Point Detection**: Methods for single and multiple change-point detection

### Implementation Features

- Efficient C++ implementations for core functions (via Rcpp)
- Parallel computing support to speed up permutation tests

## Example Usage

### Two-Sample Testing

Test whether two samples come from the same distribution:

```r
library(KDist)

# Generate two samples
set.seed(123)
n <- 100
p <- 5
x <- matrix(rnorm(n*p), n, p)
y <- matrix(rnorm(n*p, mean = 0.5), n, p)

# MMD test with Gaussian kernel
result <- mmd_test(x, y, type = "gaussian", n_perm = 199)
print(result)
plot(result)
```

### Independence Testing

Test independence between two random variables:

```r
# Generate dependent data
set.seed(456)
n <- 100
x <- matrix(rnorm(n*2), n, 2)
y <- x + matrix(rnorm(n*2, sd = 0.5), n, 2)

# HSIC test
test_result <- hsic_test(x, y, type = "gaussian", n_perm = 199)
print(test_result)
```

### Testing Mutual Independence Among Multiple Variables

```r
# Generate data with correlated variables
set.seed(789)
n <- 100
p <- 10
X <- matrix(rnorm(n*p), n, p)
X[,1:3] <- X[,1:3] + matrix(rnorm(n*3, sd=0.2), n, 3)  # Create dependency

# Test for mutual independence
result <- mhsic(X, type = "gaussian")
print(result)
```

### Distance-based MANOVA

```r
# Create example data
set.seed(123)
n <- 100
p <- 5
group <- rep(1:3, each=n/3)
x <- matrix(0, n, p)
x[group == 1, ] <- matrix(rnorm(n/3*p), n/3, p)
x[group == 2, ] <- matrix(rnorm(n/3*p, mean=0.5), n/3, p)
x[group == 3, ] <- matrix(rnorm(n/3*p, mean=1.0), n/3, p)
data <- data.frame(y=x, group=as.factor(group))

# Distance-based MANOVA
result <- dmanova(y ~ group, data=data, method="permutation", n_perm=99)
print(result$aov.tab)
```

### Change-point Detection

```r
# Generate data with change point
set.seed(123)
n <- 200
cp <- 100  # True change point
data1 <- matrix(rnorm(cp*5), cp, 5)
data2 <- matrix(rnorm(cp*5, mean=0.8), cp, 5) 
data <- rbind(data1, data2)

# Detect single change point
result <- kcpd_single(data, type="e-dist", method="permutation", B=99)
print(result)

# Multiple change-point detection
multi_result <- kcpd_sbs(data, type="e-dist", method="permutation", B=99)
print(multi_result)
```

### Conditional Independence Testing

```r
# Generate data
set.seed(123)
n <- 100
Z <- matrix(rnorm(n*2), n, 2)
X <- Z + matrix(rnorm(n*2, sd=0.5), n, 2)
Y <- Z + matrix(rnorm(n*2, sd=0.5), n, 2)

# Test conditional independence X ⊥ Y | Z
result <- kcid(X, Y, Z, type_x="gaussian", type_y="gaussian", type_z="gaussian", method="knn")
print(result)
```

## Available Functions

### Core Functions

- `KDist_matrix()`: Compute distance or kernel matrices
- `bw_optim()`: Automatic bandwidth selection for kernel methods

### Two-Sample Testing

- `mmd()`: Calculate Maximum Mean Discrepancy (MMD)
- `mmd_test()`: Permutation test based on MMD
- `hd_two_test()`: High-dimensional two-sample test

### Independence Testing

- `hsic()`: Calculate Hilbert-Schmidt Independence Criterion (HSIC) or Distance Covariance (dCov)
- `hsic_test()`: Permutation test for independence using HSIC/dCov
- `hsic_cor()`: HSIC-based correlation coefficient or squared distance correlation
- `dhsic()`: d-variable HSIC/dCov for multiple datasets
- `dhsic_test()`: Multi-variable independence test
- `jhsic()`: Joint independence measure
- `jhsic_test()`: Joint independence test
- `mhsic()`: Mutual independence test for high-dimensional data

### Multivariate Analysis

- `dmanova()`: Distance-based MANOVA
- `dmanova2()`: Alternative implementation of distance-based MANOVA

### Conditional Testing

- `kcid()`: Kernel conditional independence test
- `tcdt()`: Two-sample conditional distribution test
- `knn_conditional_sampling()`: KNN-based conditional sampling

### Change-point Detection

- `kcpd_single()`: Single change-point detection
- `kcpd_sbs()`: Multiple change-point detection with seeded binary segmentation
- `get_seeded_intervals()`: Generate intervals for binary segmentation

## How to Cite

If you use KDist in your research, please consider citing it as:

```
Zhang, X. (2025). KDist: A Collection of Kernel and Distance Methods for Statistical Inference. R package version 0.1.0.
```

## License

This package is available under the [MIT License](LICENSE).

## References

- Chen, J., & Zhang, X. (2022). D-MANOVA: Fast Distance-based Multivariate Analysis of Variance for Large-scale Microbiome Association Studies. *Bioinformatics, 38*, 286-288.
- Chakraborty, S., & Zhang, X. (2021). A New Framework for Distance and Kernel-based Metrics in High-dimension. *Electronic Journal of Statistics, 15*, 5455-5522.
- Chakraborty, S., & Zhang, X. (2021). High-dimensional Change-point Detection Using Generalized Homogeneity Metrics. arXiv:2105.08976.
- Zhu, C., Zhang, X., Yao, S., & Shao, X. (2020). Distance-based and RKHS-based Dependence Metrics in High-dimension. *Annals of Statistics, 48*, 3366-3394.
- Chakraborty, S., & Zhang, X. (2019). Distance Metrics for Measuring Joint Dependence with Application to Causal Inference. *Journal of the American Statistical Association, 114*, 1638-1650.
- Yao, S., Zhang, X., & Shao, X. (2018). Testing Mutual Independence in High Dimension via Distance Covariance. *Journal of the Royal Statistical Society Series B, 80*, 455-480.
- Gretton, A., Borgwardt, K. M., Rasch, M. J., Schölkopf, B., & Smola, A. (2012). A kernel two-sample test. *Journal of Machine Learning Research, 13*, 723-773.
- Székely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing dependence by correlation of distances. *The Annals of Statistics, 35*(6), 2769-2794.
- Pfister, N., Bühlmann, P., Schölkopf, B., & Peters, J. (2018). Kernel-based tests for joint independence. *Journal of the Royal Statistical Society: Series B, 80*(1), 5-31.

