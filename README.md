# KDistTest: Kernel-Based Distance Methods for Statistical Testing

## Overview

KDistTest is a comprehensive R package for kernel-based distance methods in statistical testing. It implements a wide range of non-parametric tests using distance and kernel-based approaches, including:

- Generalized Energy Distance
- Maximum Mean Discrepancy (MMD)
- Hilbert-Schmidt Independence Criterion (HSIC)
- Joint HSIC
- Distance-based MANOVA
- Conditional Independence Testing
- Change-point Detection

These methods are particularly useful for high-dimensional data analysis, testing independence between random variables, multivariate hypothesis testing, and detecting structural changes in sequential data.

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("YourUsername/KDistTest")
```

## Key Features

- **Flexible Distance Metrics**: Various distance and kernel options including Euclidean, Gaussian, Laplacian, and more.
- **High-dimensional Testing**: Methods optimized for high-dimensional settings.
- **Parallel Computing**: Support for parallel processing to speed up permutation tests.
- **Change-point Detection**: Binary segmentation methods for multiple change-point detection.
- **Conditional Testing**: Tests for conditional independence and conditional distribution differences.

## Main Functions

### Distance and Kernel Functions

- `KDist_matrix()`: Compute distance or kernel matrices with various metrics
- `bw_optim()`: Automatic bandwidth selection for kernel methods
- `matrix_v_center()`: Double-centering for distance/kernel matrices

### Two-Sample Testing

- `mmd_test()`: Test for differences between two multivariate distributions
- `hd_two_test()`: High-dimensional two-sample test

### Independence Testing

- `hsic_test()`: Test for independence between two variables
- `dhsic_test()`: Multi-variable independence test with distance-based HSIC
- `jhsic_test()`: Joint independence test with various statistics
- `mhsic()`: Mutual independence test for many variables
- `hd_dep_test()`: High-dimensional dependence test

### Multivariate Analysis

- `dmanova()`: Distance-based MANOVA with permutation option
- `dmanova2()`: Alternative implementation of distance-based MANOVA

### Conditional Testing

- `kcid()`: Kernel conditional independence test
- `tcdt()`: Two-sample conditional distribution test
- `knn_conditional_sampling()`: K-nearest neighbor conditional sampling

### Change-point Detection

- `kcpd_single()`: Single change-point detection
- `kcpd_sbs()`: Multiple change-point detection with seeded binary segmentation

## Examples

### Two-Sample Testing

```r
# Generate two samples
set.seed(123)
n <- 100
p <- 10
x <- matrix(rnorm(n*p), n, p)
y <- matrix(rnorm(n*p, mean = 0.5), n, p)

# MMD test
result <- mmd_test(x, y, type = "gaussian", n_perm = 199)
print(result)
plot(result)
```

### Independence Testing

```r
# Generate dependent data
set.seed(456)
n <- 200
x <- matrix(rnorm(n*2), n, 2)
y <- x + matrix(rnorm(n*2, sd = 0.5), n, 2)

# HSIC test
result <- hsic_test(x, y, type = "gaussian", n_perm = 199)
print(result)
```

### Distance-based MANOVA

```r
# Create example data
set.seed(789)
n <- 100
p <- 5
group <- rep(1:2, each = n/2)
x <- matrix(0, n, p)
x[group == 1, ] <- matrix(rnorm(n/2*p), n/2, p)
x[group == 2, ] <- matrix(rnorm(n/2*p, mean = 0.8), n/2, p)
data <- data.frame(y = x, group = as.factor(group))

# Distance-based MANOVA
out <- dmanova(y ~ group, data = data, method = "permutation", n_perm = 99)
print(out$aov.tab)
```

### Change-point Detection

```r
# Generate data with change point
set.seed(123)
n <- 200
cp <- 100  # True change point
data1 <- matrix(rnorm(cp*5), cp, 5)
data2 <- matrix(rnorm(cp*5, mean = 0.8), cp, 5) 
data <- rbind(data1, data2)

# Detect single change point
result <- kcpd_single(data, type = "e-dist", method = "permutation", B = 99)
print(result)
```

## References

- Shubhadeep, C., & Zhang, X. (2021). *A general framework for high-dimensional two-sample and independence tests.*
- Szekely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). *Measuring and testing dependence by correlation of distances.*
- Gretton, A., Bousquet, O., Smola, A., & Schölkopf, B. (2005). *Measuring statistical dependence with Hilbert-Schmidt norms.*

## License

This package is available under the MIT License.
