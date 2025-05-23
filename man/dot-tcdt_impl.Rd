% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conditional_tests.R
\name{.tcdt_impl}
\alias{.tcdt_impl}
\title{Implementation of Two-Sample Conditional Distribution Test}
\usage{
.tcdt_impl(
  X1,
  X2,
  Y1,
  Y2,
  x0 = NULL,
  B = 299,
  h = NULL,
  adj_bw = 0.1,
  const_factor = 1,
  h_boot = NULL,
  adj_bw_boot = 0,
  const_factor_boot = 1,
  scale_factor = 1,
  kern_smooth = "gaussian",
  nu = 2,
  stat = "cmmd",
  kern_mmd = "gaussian",
  sampling_method = "bootstrap",
  knn_k = 5,
  knn_weighted = FALSE,
  num_cores = 1
)
}
\arguments{
\item{X1}{First conditioning dataset (matrix, data frame, or vector)}

\item{X2}{Second conditioning dataset with the same number of columns as X1}

\item{Y1}{First response dataset (matrix, data frame, or vector)}

\item{Y2}{Second response dataset with the same number of columns as Y1}

\item{x0}{Local point for local testing (default: NULL for global testing)}

\item{B}{Number of bootstrap or resampling iterations (default: 299)}

\item{h}{Bandwidth for local smoothing of X (default: NULL for automatic selection)}

\item{adj_bw}{Adjustment factor for automatic bandwidth (default: 0.1)}

\item{const_factor}{Constant scaling factor for automatic bandwidth (default: 1)}

\item{h_boot}{Bandwidth for bootstrap resampling (default: NULL)}

\item{adj_bw_boot}{Adjustment factor for bootstrap bandwidth (default: 0)}

\item{const_factor_boot}{Constant scaling factor for bootstrap bandwidth (default: 1)}

\item{scale_factor}{Scaling factor for kernel bandwidth in MMD calculation (default: 1)}

\item{kern_smooth}{Kernel type for smoothing (default: "gaussian")}

\item{nu}{Order of smoothing kernel (default: 2)}

\item{stat}{Test statistic type (default: "cmmd")}

\item{kern_mmd}{Kernel type for MMD when stat="cmmd" (default: "gaussian")}

\item{sampling_method}{Resampling method (default: "bootstrap")}

\item{knn_k}{Number of nearest neighbors for KNN (default: 5)}

\item{knn_weighted}{Whether to use weighted KNN (default: FALSE)}

\item{num_cores}{Number of cores for parallel computation (default: 1)}
}
\value{
A list containing test results similar to the original tcdt function
}
\description{
Internal implementation function for the two-sample conditional distribution test.
This function should not be called directly by users unless they need
fine-grained control over all parameters.
}
\keyword{internal}
