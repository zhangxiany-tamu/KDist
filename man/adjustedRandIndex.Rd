% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/change_point.R
\name{adjustedRandIndex}
\alias{adjustedRandIndex}
\title{Adjusted Rand Index for Comparing Clusterings}
\usage{
adjustedRandIndex(x, y)
}
\arguments{
\item{x}{First clustering (vector of cluster assignments)}

\item{y}{Second clustering (vector of cluster assignments)}
}
\value{
A numeric value between -1 and 1:
\itemize{
\item 1: Perfect agreement between the two clusterings
\item ~0: Agreement equivalent to random chance
\item Negative values: Agreement less than random chance
}
}
\description{
Calculates the Adjusted Rand Index (ARI) between two clusterings or partitions of the same data.
The ARI is a measure of the similarity between two data clusterings, adjusted for chance.
}
\details{
The Adjusted Rand Index is a modified version of the Rand Index that is adjusted for the
chance groupings of elements. It has an expected value of 0 for random partitions and a
maximum value of 1 for identical clusterings.

The formula for the ARI is:

\deqn{ARI = \frac{\sum_{ij} \binom{n_{ij}}{2} - \big[\sum_i \binom{a_i}{2} \sum_j \binom{b_j}{2} \big] / \binom{n}{2}}
{\frac{1}{2}\big[\sum_i \binom{a_i}{2} + \sum_j \binom{b_j}{2}\big] - \big[\sum_i \binom{a_i}{2} \sum_j \binom{b_j}{2}\big] / \binom{n}{2}}}

where:
\itemize{
\item \eqn{n_{ij}} is the number of objects in both class \eqn{i} of clustering \eqn{x} and class \eqn{j} of clustering \eqn{y}
\item \eqn{a_i} is the number of objects in class \eqn{i} of clustering \eqn{x}
\item \eqn{b_j} is the number of objects in class \eqn{j} of clustering \eqn{y}
\item \eqn{n} is the total number of objects
}

The ARI is especially useful for evaluating change point detection algorithms by comparing the
detected segmentation with a ground truth segmentation.
}
\examples{
# Example 1: Identical clusterings
x <- c(1, 1, 1, 2, 2, 3, 3, 3)
y <- c(1, 1, 1, 2, 2, 3, 3, 3)
adjustedRandIndex(x, y)  # Should return 1

# Example 2: Different label names but same structure
x <- c(1, 1, 1, 2, 2, 3, 3, 3)
y <- c(5, 5, 5, 9, 9, 2, 2, 2)
adjustedRandIndex(x, y)  # Should return 1

# Example 3: Similar but not identical clusterings
x <- c(1, 1, 1, 2, 2, 3, 3, 3)
y <- c(1, 1, 2, 2, 2, 3, 3, 3)
adjustedRandIndex(x, y)  # Should return a value between 0 and 1

# Example 4: Very different clusterings
x <- c(1, 1, 1, 1, 2, 2, 2, 2)
y <- c(1, 2, 1, 2, 1, 2, 1, 2)
adjustedRandIndex(x, y)  # Should return a value close to 0

# Example 5: Comparing change point detection results with ground truth
true_segments <- c(rep(1, 50), rep(2, 30), rep(3, 40))
detected_segments <- c(rep(1, 45), rep(2, 35), rep(3, 40))
ari <- adjustedRandIndex(true_segments, detected_segments)
print(paste("ARI between true and detected segments:", ari))

}
\references{
Hubert, L., & Arabie, P. (1985). Comparing partitions. Journal of Classification, 2(1), 193-218.
}
\seealso{
\code{\link{kcpd_single}} and \code{\link{kcpd_sbs}} for change point detection methods
that can be evaluated using the Adjusted Rand Index
}
