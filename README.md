# KDist: Kernel and Distance Methods for Statistical Inference

<div align="center">

<img src="images/kdist-logo.png" alt="KDist Logo" width="200"/>

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/KDist)](https://cran.r-project.org/package=KDist)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

</div>

## 📋 Overview

KDist provides a comprehensive collection of kernel and distance-based methods for nonparametric statistical inference. These powerful methods excel in scenarios where traditional parametric approaches may fail, particularly with complex, high-dimensional data.

### 🔍 Key Applications

- **Two-sample testing**: Detect differences between multivariate distributions
- **Independence testing**: Measure dependence between random vectors
- **Multivariate hypothesis testing**: Perform nonparametric MANOVA
- **High-dimensional analysis**: Analyze data where dimensions exceed sample size
- **Change-point detection**: Identify structural changes in sequential data
- **Conditional testing**: Test conditional independence and distribution differences

## 🚀 Installation

```r
# Development version from GitHub
# install.packages("devtools")
devtools::install_github("zhangxiany-tamu/KDist")
```

## ✨ Features

### 🧮 Distance and Kernel Options

- **Distance metrics**: Euclidean, kernel-induced distances
- **Kernel functions**: Gaussian, Laplacian, Polynomial
- **Optimization features**:
  - Automatic bandwidth selection
  - Group-based distance calculations
  - Efficient C++ implementations
  - Parallel computing support

### 🔬 Testing Methods

| Category | Methods |
|----------|---------|
| **Two-Sample Testing** | Maximum Mean Discrepancy (MMD), Energy Distance (ED) |
| **Independence Testing** | Hilbert-Schmidt Independence Criterion (HSIC), Distance Covariance (dCov) |
| **Multi-Variable Testing** | d-variable HSIC (dHSIC), Joint HSIC (jHSIC), Joint dCov (jdCov) |
| **High-Dimensional Testing** | Specialized methods for p >> n settings |
| **Multivariate Analysis** | Kernel and distance-based MANOVA |
| **Conditional Testing** | Conditional independence and distribution tests |
| **Change-Point Analysis** | Single and multiple change-point detection |

## 🏎️ Performance Comparison

KDist offers comparable or superior performance in both accuracy and computational efficiency compared to other popular packages:

### Accuracy Comparison

```r
library(energy)
library(dHSIC)
library(KDist)

# Compare distance covariance implementation
set.seed(123)
n <- 1000
x <- rnorm(n)
y <- x + rnorm(n)

# Standard distance covariance
energy::dcov(x,y)^2  # Result: ~0.265
KDist::dcov(x,y)     # Result: ~0.265

# Distance correlation
energy::dcor(x,y)^2  # Result: ~0.454
KDist::dcor(x,y)     # Result: ~0.454

# Unbiased distance covariance
energy::dcovU(x,y)   # Result: ~0.263
KDist::dcov(x, y, u_center = TRUE)  # Result: ~0.263

# Multi-variable independence measures
set.seed(456)
u <- runif(500, -pi, pi)
x1 <- matrix(sin(u), ncol = 1)
x2 <- matrix(cos(u), ncol = 1)
x3 <- matrix(sin(u) * cos(u), ncol = 1)

# d-variable HSIC
dHSIC::dhsic(list(x1,x2,x3))      # Result: ~0.0538
KDist::dhsic(list(x1, x2, x3), type = "gaussian")  # Result: ~0.0538
```

### Speed Comparison

```r
library(microbenchmark)

# Compare speed for distance covariance
n <- 1000
x <- rnorm(n)
y <- x + rnorm(n)

microbenchmark(
  energy = energy::dcov(x, y)^2,
  KDist = KDist::dcov(x, y),
  times = 1000
)
# Result:
# Unit: milliseconds
#   expr       min        lq      mean    median        uq      max neval
# energy 26.389937 30.262490 37.462686 32.261280 34.303224 116.4194  1000
#  KDist  3.945307  5.009933  9.549976  5.135844  8.245859  77.6524  1000

# Compare speed for multi-variable independence measures
microbenchmark(
  dHSIC = dHSIC::dhsic(list(x1, x2, x3)),
  KDist = KDist::dhsic(list(x1, x2, x3), type = "gaussian"),
  times = 100
)
# Result:
# Unit: milliseconds
#  expr      min       lq      mean   median        uq      max neval
# dHSIC 7.540105 8.813565 10.421584 9.059462 11.489123 83.57329  1000
# KDist 6.176732 6.908684  7.819992 6.986277  7.158436 74.45932  1000
```

These benchmarks demonstrate that KDist implementations are typically 2-3× faster than comparable packages due to efficient C++ optimization.

KDist provides significant advantages:
- **Unified framework** for kernel and distance-based methods
- **Optimized C++ implementation** for computational efficiency
- **Automatic bandwidth selection** for kernel methods
- **High-dimensional extensions** specifically designed for p >> n settings
- **Comprehensive toolkit** with methods for various hypothesis testing scenarios

## 📊 Examples

### Two-Sample Testing

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

```r
# Generate dependent data
set.seed(123)
n <- 100
x <- matrix(rnorm(n*2), n, 2)
y <- x + matrix(rnorm(n*2, sd = 0.5), n, 2)

# HSIC test
test_result <- hsic_test(x, y, type = "gaussian", n_perm = 199)
print(test_result)
```

### Change-point Detection in High-dimension

```r
# Generate data with change point
set.seed(123)
p <- 200
cp <- 30   # True change point
data1 <- matrix(rnorm(cp*p), cp, p)
data2 <- matrix(rnorm(cp*p, mean=0.3), cp, p)
data <- rbind(data1, data2)

# Detect single change point
result <- kcpd_single(data, type="e-dist", method="permutation", B=99)
print(result)
```

*Additional examples available in the package documentation.*

## 📚 Available Functions

<details>
<summary>Click to expand function list</summary>

### Core Functions
- `KDist_matrix()`: Compute distance or kernel matrices
- `bw_optim()`: Automatic bandwidth selection for kernel methods

### Two-Sample Testing
- `mmd()`: Calculate Maximum Mean Discrepancy (MMD)
- `mmd_test()`: Permutation test based on MMD
- `ed()`: Calculate Energy Distance
- `ed_test()`: Energy distance-based test
- `hd_two_test()`: High-dimensional two-sample test

### Independence Testing
- `hsic()`: Calculate Hilbert-Schmidt Independence Criterion
- `hsic_test()`: HSIC-based independence test
- `hsic_cor()`: HSIC-based correlation coefficient
- `dcov()`: Calculate distance covariance
- `dcov_test()`: Distance covariance test
- `dcor()`: Calculate distance correlation
- `dhsic()`: Multi-variable independence measure
- `dhsic_test()`: Multi-variable independence test
- `jhsic()`: Joint independence measure
- `jhsic_test()`: Joint independence test
- `mhsic()`: Mutual independence test for high-dimensional data
- `hd_dep_test()`: High-dimensional dependence test

### Multivariate Analysis
- `dmanova()`: Distance-based MANOVA
- `dmanova2()`: Fully distance-based MANOVA

### Conditional Testing
- `kcid()`: Kernel conditional independence test
- `tcdt()`: Two-sample conditional distribution test
- `knn_conditional_sampling()`: KNN-based conditional sampling

### Change-point Detection
- `kcpd_single()`: Single change-point detection
- `kcpd_sbs()`: Multiple change-point detection with seeded binary segmentation
- `kcpd_wbs()`: Multiple change-point detection with wild binary segmentation
- `get_seeded_intervals()`: Generate intervals for binary segmentation

### Utilities
- `v_center()`: V-centering for matrices
- `u_center()`: U-centering for matrices
- Various plot and print methods

</details>

## 📝 Citation

If you use KDist in your research, please cite:

```
Zhang, X. (2025). KDist: A Collection of Kernel and Distance Methods for Statistical Inference. 
R package version 0.1.0.
```

## 📄 License

This package is available under the [MIT License](LICENSE).

## 📚 References

<details>
<summary>Click to expand references</summary>

1. Székely, G. J., & Rizzo, M. L. (2004). Testing for Equal Distributions in High Dimension. *InterStat*, Nov. (5).

2. Gretton, A., Fukumizu, K., Teo, C., Song, L., Schölkopf, B., & Smola, A. (2007). A kernel statistical test of independence. *Advances in Neural Information Processing Systems*, 20.

3. Székely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing dependence by correlation of distances. *The Annals of Statistics*, 35(6), 2769-2794.

4. Gretton, A., Borgwardt, K. M., Rasch, M. J., Schölkopf, B., & Smola, A. (2012). A kernel two-sample test. *Journal of Machine Learning Research*, 13, 723-773.

5. Székely, G. J., & Rizzo, M. L. (2013). Energy statistics: A class of statistics based on distances. *Journal of Statistical Planning and Inference*, 143(8), 1249-1272.

6. Pfister, N., Bühlmann, P., Schölkopf, B., & Peters, J. (2018). Kernel-based tests for joint independence. *Journal of the Royal Statistical Society: Series B*, 80(1), 5-31.

7. Yao, S., Zhang, X., & Shao, X. (2018). Testing Mutual Independence in High Dimension via Distance Covariance. *Journal of the Royal Statistical Society: Series B*, 80, 455-480.

8. Chakraborty, S., & Zhang, X. (2019). Distance Metrics for Measuring Joint Dependence with Application to Causal Inference. *Journal of the American Statistical Association*, 114, 1638-1650.

9. Zhu, C., Zhang, X., Yao, S., & Shao, X. (2020). Distance-based and RKHS-based Dependence Metrics in High-dimension. *Annals of Statistics*, 48, 3366-3394.

10. Chakraborty, S., & Zhang, X. (2021). A New Framework for Distance and Kernel-based Metrics in High-dimension. *Electronic Journal of Statistics*, 15, 5455-5522.

11. Chakraborty, S., & Zhang, X. (2021). High-dimensional Change-point Detection Using Generalized Homogeneity Metrics. arXiv:2105.08976.

12. Chen, J., & Zhang, X. (2022). D-MANOVA: Fast Distance-based Multivariate Analysis of Variance for Large-scale Microbiome Association Studies. *Bioinformatics*, 38, 286-288.

13. Zhang, K., Peters, J., Janzing, D., & Schölkopf, B. (2012). Kernel-based Conditional Independence Test and Application in Causal Discovery. arXiv preprint arXiv:1202.3775.

14. Yan, J., Li, Z., & Zhang, X. (2024). Distance and Kernel-Based Measures for Global and Local Two-Sample Conditional Distribution Testing. arXiv preprint arXiv:2210.08149.

</details>
