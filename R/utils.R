#' Print a summary of a mmd_test object
#'
#' @param x A mmd_test object
#' @param ... Additional arguments to be passed to methods
#'
#' @return Invisibly returns the object
#' @export
print.mmd_test <- function(x, ...) {
  cat("Maximum Mean Discrepancy Test\n\n")
  cat("Statistic:", format(x$statistic, digits = 4), "\n")
  cat("P-value:", format(x$p.value, digits = 4), "\n")
  cat("\n")
  cat("Based on", length(x$permutation_values), "permutations\n")
  invisible(x)
}

#' Print a summary of a hsic_test object
#'
#' @param x A hsic_test object
#' @param ... Additional arguments to be passed to methods
#'
#' @return Invisibly returns the object
#' @export
print.hsic_test <- function(x, ...) {
  cat("Hilbert-Schmidt Independence Criterion Test\n\n")
  cat("Statistic:", format(x$statistic, digits = 4), "\n")
  cat("P-value:", format(x$p.value, digits = 4), "\n")
  cat("\n")
  cat("Based on", length(x$permutation_values), "permutations\n")
  invisible(x)
}

#' Print a summary of a dhsic_test object
#'
#' @param x A dhsic_test object
#' @param ... Additional arguments to be passed to methods
#'
#' @return Invisibly returns the object
#' @export
print.dhsic_test <- function(x, ...) {
  cat("Distance Hilbert-Schmidt Independence Criterion Test\n\n")
  cat("Statistic:", format(x$statistic, digits = 4), "\n")
  cat("P-value:", format(x$p.value, digits = 4), "\n")
  cat("\n")
  cat("Based on", length(x$permutation_values), "permutations\n")
  invisible(x)
}

#' Print a summary of a jhsic_test object
#'
#' @param x A jhsic_test object
#' @param ... Additional arguments to be passed to methods
#'
#' @return Invisibly returns the object
#' @export
print.jhsic_test <- function(x, ...) {
  cat("Joint Hilbert-Schmidt Independence Criterion Test\n\n")
  cat("Statistic type:", x$stat_description, "\n")
  cat("Statistic:", format(x$statistic, digits = 4), "\n")
  cat("P-value:", format(x$p.value, digits = 4), "\n")
  cat("\n")
  cat("Based on", length(x$permutation_values), "permutations\n")
  invisible(x)
}

#' Plot a histogram of permutation values for independence tests
#'
#' @param x A test object (mmd_test, hsic_test, dhsic_test, or jhsic_test)
#' @param ... Additional arguments to be passed to methods
#'
#' @return Invisibly returns the plot object
#' @export
plot.mmd_test <- function(x, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  hist(x$permutation_values, 
       main = "Histogram of Permutation Values",
       xlab = "MMD Statistic",
       freq = FALSE, 
       ...)
  
  # Draw a vertical line for the observed statistic
  abline(v = x$statistic, col = "red", lwd = 2)
  
  # Add a legend
  legend("topright", 
         legend = c(paste("Observed =", format(x$statistic, digits = 4)),
                    paste("P-value =", format(x$p.value, digits = 4))),
         bty = "n")
  
  invisible(x)
}

#' Plot method for hsic_test
#'
#' @param x A hsic_test object
#' @param ... Additional arguments to be passed to methods
#'
#' @return Invisibly returns the plot object
#' @export
plot.hsic_test <- function(x, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  hist(x$permutation_values, 
       main = "Histogram of Permutation Values",
       xlab = "HSIC Statistic",
       freq = FALSE, 
       ...)
  
  # Draw a vertical line for the observed statistic
  abline(v = x$statistic, col = "red", lwd = 2)
  
  # Add a legend
  legend("topright", 
         legend = c(paste("Observed =", format(x$statistic, digits = 4)),
                    paste("P-value =", format(x$p.value, digits = 4))),
         bty = "n")
  
  invisible(x)
}

#' Plot method for dhsic_test
#'
#' @param x A dhsic_test object
#' @param ... Additional arguments to be passed to methods
#'
#' @return Invisibly returns the plot object
#' @export
plot.dhsic_test <- function(x, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  hist(x$permutation_values, 
       main = "Histogram of Permutation Values",
       xlab = "dHSIC Statistic",
       freq = FALSE, 
       ...)
  
  # Draw a vertical line for the observed statistic
  abline(v = x$statistic, col = "red", lwd = 2)
  
  # Add a legend
  legend("topright", 
         legend = c(paste("Observed =", format(x$statistic, digits = 4)),
                    paste("P-value =", format(x$p.value, digits = 4))),
         bty = "n")
  
  invisible(x)
}

#' Plot method for jhsic_test
#'
#' @param x A jhsic_test object
#' @param ... Additional arguments to be passed to methods
#'
#' @return Invisibly returns the plot object
#' @export
plot.jhsic_test <- function(x, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  hist(x$permutation_values, 
       main = "Histogram of Permutation Values",
       xlab = "JHSIC Statistic",
       freq = FALSE, 
       ...)
  
  # Draw a vertical line for the observed statistic
  abline(v = x$statistic, col = "red", lwd = 2)
  
  # Add a legend
  legend("topright", 
         legend = c(paste("Test type:", x$stat_description),
                    paste("Observed =", format(x$statistic, digits = 4)),
                    paste("P-value =", format(x$p.value, digits = 4))),
         bty = "n")
  
  invisible(x)
}
