// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Calculate the bandwidth
// [[Rcpp::export]]
double bw_rcpp(NumericMatrix x, int sample_limit = 1000, double scale_factor = 0.5) {
  int n = x.nrow();
  int p = x.ncol();

  // Limit sample size if needed
  if(n > sample_limit) {
    n = sample_limit;
  }

  // Calculate total pairs and middle position
  int npairs = n * (n - 1) / 2;
  int middle = npairs / 2;

  // Allocate vector for squared distances
  NumericVector d_squared(npairs);

  // Calculate squared distances - optimized implementation
  int count = 0;
  for(int i = 0; i < n - 1; i++) {
    NumericMatrix::Row row_i = x.row(i); // Get reference to row i
    for(int j = i + 1; j < n; j++) {
      double sum_squared = 0.0;
      NumericMatrix::Row row_j = x.row(j); // Get reference to row j

      // Unrolled loop for better performance when p is small
      int k = 0;
      while(k + 3 < p) { // Process 4 elements at once
        double diff1 = row_i[k] - row_j[k];
        double diff2 = row_i[k+1] - row_j[k+1];
        double diff3 = row_i[k+2] - row_j[k+2];
        double diff4 = row_i[k+3] - row_j[k+3];

        sum_squared += diff1*diff1 + diff2*diff2 + diff3*diff3 + diff4*diff4;
        k += 4;
      }

      // Handle remaining elements
      while(k < p) {
        double diff = row_i[k] - row_j[k];
        sum_squared += diff * diff;
        k++;
      }

      d_squared[count++] = sum_squared;
    }
  }

  // Find median using nth_element (partial sort)
  NumericVector d_copy = clone(d_squared);
  std::nth_element(d_copy.begin(), d_copy.begin() + middle, d_copy.end());
  double median_val = d_copy[middle];

  // Apply transformation with customizable scale factor
  return sqrt(median_val * scale_factor);
}

// Euclidean distance with optional exponent
// [[Rcpp::export]]
NumericMatrix euclidean_dist_rcpp(NumericMatrix data, double expo = 1) {
  int n = data.nrow();
  int p = data.ncol();
  NumericMatrix result(n, n);

  for (int i = 0; i < n; i++) {
    result(i, i) = 0.0;  // Diagonal is always 0

    for (int j = i + 1; j < n; j++) {
      double sum_squared = 0.0;

      // Get references to rows for better performance
      NumericMatrix::Row row_i = data.row(i);
      NumericMatrix::Row row_j = data.row(j);

      // Unrolled loop processing 4 elements at a time
      int k = 0;
      while (k + 3 < p) {
        double diff1 = row_i[k] - row_j[k];
        double diff2 = row_i[k+1] - row_j[k+1];
        double diff3 = row_i[k+2] - row_j[k+2];
        double diff4 = row_i[k+3] - row_j[k+3];

        sum_squared += diff1 * diff1 + diff2 * diff2 + diff3 * diff3 + diff4 * diff4;
        k += 4;
      }

      // Handle remaining elements
      for (; k < p; k++) {
        double diff = row_i[k] - row_j[k];
        sum_squared += diff * diff;
      }

      // Calculate distance
      double dist;
      if (expo == 1.0) {
        dist = std::sqrt(sum_squared);
      } else if (expo == 2.0) {
        // Optimization for squared Euclidean distance
        dist = sum_squared;
      } else {
        dist = std::pow(std::sqrt(sum_squared), expo);
      }

      result(i, j) = dist;
      result(j, i) = dist;  // Matrix is symmetric
    }
  }

  return result;
}

// Gaussian kernel
// [[Rcpp::export]]
NumericMatrix gaussian_kernel_rcpp(NumericMatrix data, double bw) {
  int n = data.nrow();
  int p = data.ncol();
  NumericMatrix result(n, n);

  // Precompute 1/(2*bw^2) for efficiency
  double bw_factor = 1.0 / (2.0 * pow(bw, 2.0));

  for (int i = 0; i < n; i++) {
    result(i, i) = 1.0;  // Self-similarity is 1

    for (int j = i + 1; j < n; j++) {
      double sum_squared = 0.0;

      // Get references to rows for better performance
      NumericMatrix::Row row_i = data.row(i);
      NumericMatrix::Row row_j = data.row(j);

      // Unrolled loop processing 4 elements at a time
      int k = 0;
      while (k + 3 < p) {
        double diff1 = row_i[k] - row_j[k];
        double diff2 = row_i[k+1] - row_j[k+1];
        double diff3 = row_i[k+2] - row_j[k+2];
        double diff4 = row_i[k+3] - row_j[k+3];

        sum_squared += diff1 * diff1 + diff2 * diff2 + diff3 * diff3 + diff4 * diff4;
        k += 4;
      }

      // Handle remaining elements
      for (; k < p; k++) {
        double diff = row_i[k] - row_j[k];
        sum_squared += diff * diff;
      }

      // Apply Gaussian kernel transformation
      double val = exp(-sum_squared * bw_factor);

      result(i, j) = val;
      result(j, i) = val;  // Matrix is symmetric
    }
  }

  return result;
}

// laplacian kernel
// [[Rcpp::export]]
NumericMatrix laplacian_kernel_rcpp(NumericMatrix data, double bw) {
  int n = data.nrow();
  int p = data.ncol();
  NumericMatrix result(n, n);

  // Precompute 1/bw for efficiency
  double bw_inv = 1.0 / bw;

  for (int i = 0; i < n; i++) {
    result(i, i) = 1.0;  // Self-similarity is 1

    for (int j = i + 1; j < n; j++) {
      double sum_squared = 0.0;

      // Get references to rows for better performance
      NumericMatrix::Row row_i = data.row(i);
      NumericMatrix::Row row_j = data.row(j);

      // Unrolled loop processing 4 elements at a time
      int k = 0;
      while (k + 3 < p) {
        double diff1 = row_i[k] - row_j[k];
        double diff2 = row_i[k+1] - row_j[k+1];
        double diff3 = row_i[k+2] - row_j[k+2];
        double diff4 = row_i[k+3] - row_j[k+3];

        sum_squared += diff1 * diff1 + diff2 * diff2 + diff3 * diff3 + diff4 * diff4;
        k += 4;
      }

      // Handle remaining elements
      for (; k < p; k++) {
        double diff = row_i[k] - row_j[k];
        sum_squared += diff * diff;
      }

      // Calculate Euclidean distance
      double dist = std::sqrt(sum_squared);

      // Apply laplacian kernel transformation
      double val = std::exp(-dist * bw_inv);

      result(i, j) = val;
      result(j, i) = val;  // Matrix is symmetric
    }
  }

  return result;
}

// Polynomial kernel
// [[Rcpp::export]]
NumericMatrix polynomial_kernel_rcpp(NumericMatrix data, double degree) {
  int n = data.nrow();
  int p = data.ncol();
  NumericMatrix result(n, n);

  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      double dot_product = 0.0;

      // Get references to rows for better performance
      NumericMatrix::Row row_i = data.row(i);
      NumericMatrix::Row row_j = data.row(j);

      // Unrolled loop processing 4 elements at a time
      int k = 0;
      while (k + 3 < p) {
        dot_product += row_i[k] * row_j[k] +
          row_i[k+1] * row_j[k+1] +
          row_i[k+2] * row_j[k+2] +
          row_i[k+3] * row_j[k+3];
        k += 4;
      }

      // Handle remaining elements
      for (; k < p; k++) {
        dot_product += row_i[k] * row_j[k];
      }

      // Apply polynomial transformation: (1 + dot_product)^degree
      double val = std::pow(1.0 + dot_product, degree);

      result(i, j) = val;
      result(j, i) = val;  // Matrix is symmetric
    }
  }

  return result;
}

// Calculate group-wise Euclidean norm induced distance
// [[Rcpp::export]]
NumericMatrix e_dist_rcpp(NumericMatrix data, Nullable<IntegerVector> group_ = R_NilValue) {
  int n = data.nrow();
  int p = data.ncol();

  // Handle default group (each coordinate is its own group)
  IntegerVector group;
  if (group_.isNotNull()) {
    group = as<IntegerVector>(group_);
    // Check if group vector has correct length
    if (group.size() != p) {
      stop("Length of group vector must match number of columns in data");
    }
  } else {
    group = IntegerVector(p);
    for (int k = 0; k < p; k++) {
      group[k] = k + 1;  // Create groups 1 to p
    }
  }

  int num_groups = max(group);
  NumericMatrix result(n, n);

  // Pre-calculate group membership
  std::vector<std::vector<int> > group_members(num_groups + 1);
  for (int k = 0; k < p; k++) {
    group_members[group[k]].push_back(k);
  }

  for (int i = 0; i < n; i++) {
    result(i, i) = 0.0;  // Diagonal is always 0

    for (int j = i + 1; j < n; j++) {
      double total_dist = 0.0;

      // Calculate distance for each group
      for (int g = 1; g <= num_groups; g++) {
        double group_sum_squared = 0.0;

        // Only iterate through columns in this group
        for (size_t idx = 0; idx < group_members[g].size(); idx++) {
          int k = group_members[g][idx];
          double diff = data(i, k) - data(j, k);
          group_sum_squared += diff * diff;
        }

        // Add Euclidean norm for this group to total distance
        total_dist += std::sqrt(group_sum_squared);
      }

      // Apply final square root to match original function behavior
      result(i, j) = std::sqrt(total_dist);
      result(j, i) = result(i, j);  // Matrix is symmetric
    }
  }

  return result;
}

// Calculate group-wise laplacian and Gaussian kernel induced distance
// [[Rcpp::export]]
arma::mat kernel_induced_distance_matrix(
    const arma::mat& data,
    const arma::vec& g,
    const std::string& type,
    Nullable<IntegerVector> group_ = R_NilValue) {

  int n = data.n_rows;
  int p = data.n_cols;

  // Handle default group (each coordinate is its own group)
  IntegerVector group;
  if (group_.isNotNull()) {
    group = as<IntegerVector>(group_);
    // Check if group vector has correct length
    if (group.size() != p) {
      stop("Length of group vector must match number of columns in data");
    }
  } else {
    group = IntegerVector(p);
    for (int k = 0; k < p; k++) {
      group[k] = k + 1;  // Create groups 1 to p
    }
  }

  int num_groups = max(group);

  // Check if g has correct length (same as number of groups)
  if (g.n_elem != (unsigned)num_groups) {
    stop("Length of bandwidth vector must match number of groups");
  }

  // Check for non-positive bandwidth values for used groups
  for (unsigned int i = 0; i < g.n_elem; i++) {
    // Find if there are any columns in this group
    bool group_used = false;
    for (int k = 0; k < p; k++) {
      if (group[k] == (int)(i + 1)) {
        group_used = true;
        break;
      }
    }

    // Only check bandwidth for groups that are actually used
    if (group_used && g(i) <= 0) {
      stop("Bandwidth must be positive for used groups");
    }
  }

  arma::mat dist_matrix(n, n);

  // Pre-calculate group membership
  std::vector<std::vector<int> > group_members(num_groups + 1);
  for (int k = 0; k < p; k++) {
    group_members[group[k]].push_back(k);
  }

  for (int i = 0; i < n; i++) {
    dist_matrix(i, i) = 0.0;  // Diagonal is 0

    for (int j = i + 1; j < n; j++) {
      double sum_val = 0.0;

      // For each group
      for (int g_idx = 1; g_idx <= num_groups; g_idx++) {
        const std::vector<int>& cols = group_members[g_idx];

        // Skip empty groups
        if (cols.size() == 0) {
          continue;
        }

        // Use the g value for this group
        double bandwidth = g(g_idx - 1); // g is 0-indexed, groups are 1-indexed

        if (type == "l-dist") {
          // laplacian kernel induced distance
          double euclidean_dist = 0.0;

          // Calculate Euclidean distance within this group
          for (size_t idx = 0; idx < cols.size(); idx++) {
            int k = cols[idx];
            double diff = data(i, k) - data(j, k);
            euclidean_dist += std::pow(diff, 2);
          }
          euclidean_dist = std::sqrt(euclidean_dist);

          sum_val += 2.0 - 2.0 * std::exp(-euclidean_dist / bandwidth);

        } else {
          // Gaussian kernel induced distance
          double squared_euclidean_dist = 0.0;

          // Calculate squared Euclidean distance within this group
          for (size_t idx = 0; idx < cols.size(); idx++) {
            int k = cols[idx];
            double diff = data(i, k) - data(j, k);
            squared_euclidean_dist += std::pow(diff, 2);
          }

          sum_val += 2.0 - 2.0 * std::exp(-squared_euclidean_dist / (2.0 * std::pow(bandwidth, 2)));
        }
      }

      dist_matrix(i, j) = std::sqrt(sum_val);
      dist_matrix(j, i) = dist_matrix(i, j);  // Symmetry
    }
  }

  return dist_matrix;
}

// Calculate HSIC
// [[Rcpp::export]]
double hsic_cpp(NumericMatrix x, NumericMatrix y, std::string type = "euclidean",
                SEXP bw_sexp = R_NilValue, double expo = 1, double scale_factor = 0.5,
                Nullable<List> group_ = R_NilValue, bool u_center = false,
                bool is_distance = false) {
  // Get package environment
  Environment KDist = Environment::namespace_env("KDist");
  Function KDist_matrix = KDist["KDist_matrix"];

  // Check if x and y are the same matrix by comparing memory addresses
  bool same_matrices = (x.begin() == y.begin());

  // Extract group vectors from the list if provided
  Nullable<IntegerVector> group_x = R_NilValue;
  Nullable<IntegerVector> group_y = R_NilValue;

  if (group_.isNotNull()) {
    List group_list(group_);
    if (group_list.size() != 2) {
      stop("The 'group' parameter must be a list of size 2");
    }

    // Extract the two group vectors
    if (!Rf_isNull(group_list[0])) {
      IntegerVector g_x = as<IntegerVector>(group_list[0]);
      group_x = Nullable<IntegerVector>(g_x);
    }

    if (!Rf_isNull(group_list[1])) {
      IntegerVector g_y = as<IntegerVector>(group_list[1]);
      group_y = Nullable<IntegerVector>(g_y);
    }
  }

  // Compute distance matrices (only one if x=y) or use input directly if is_distance=TRUE
  NumericMatrix dist_x_rcpp;
  NumericMatrix dist_y_rcpp;

  if (is_distance) {
    // Use the input matrices directly as distance matrices
    dist_x_rcpp = x;
    if (same_matrices) {
      dist_y_rcpp = dist_x_rcpp;
    } else {
      dist_y_rcpp = y;
    }
  } else {
    // Calculate distance matrices using KDist_matrix with appropriate groups
    dist_x_rcpp = KDist_matrix(x, type, bw_sexp, expo, scale_factor, group_x);

    if (same_matrices) {
      dist_y_rcpp = dist_x_rcpp; // Reuse the same distance matrix
    } else {
      dist_y_rcpp = KDist_matrix(y, type, bw_sexp, expo, scale_factor, group_y);
    }
  }

  // Convert to Armadillo matrices for efficiency
  int n = x.nrow();
  arma::mat Dx(dist_x_rcpp.begin(), n, n, false);
  arma::mat Dy;

  if (same_matrices) {
    Dy = Dx; // Just use the same matrix reference
  } else {
    Dy = arma::mat(dist_y_rcpp.begin(), n, n, false);
  }

  // For Gaussian and laplacian kernel types with u_center=TRUE, ensure diagonal elements are zero
  if (u_center && (type == "gaussian" || type == "laplacian" || type == "polynomial")) {
    for (int i = 0; i < n; i++) {
      Dx(i, i) = 0.0;
      Dy(i, i) = 0.0;
    }
  }

  double sum = 0.0;

  if (u_center) {
    // U-statistic version

    // Calculate row sums
    arma::vec rowSums_x = arma::sum(Dx, 1);
    arma::vec rowSums_y;

    if (same_matrices) {
      rowSums_y = rowSums_x; // Reuse the same sums
    } else {
      rowSums_y = arma::sum(Dy, 1);
    }

    // Calculate total sums
    double totalSum_x = arma::accu(Dx);
    double totalSum_y = same_matrices ? totalSum_x : arma::accu(Dy);

    // Pre-calculate constants
    double rowDivisor = n - 2.0;
    double totalDivisor = (n - 1.0) * (n - 2.0);

    // Calculate U-centered product sum directly
    if (same_matrices) {
      // When x=y, U_x_ij = U_y_ij, so we can calculate U_x_ij^2
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          if (i != j) {  // Skip diagonal elements
            // Calculate U-centered value once
            double U_ij = Dx(i, j) +
              (totalSum_x / totalDivisor) -
              (rowSums_x(i) / rowDivisor) -
              (rowSums_x(j) / rowDivisor);

            // Square it instead of multiplying two identical values
            sum += U_ij * U_ij;
          }
        }
      }
    } else {
      // Regular case when x≠y
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          if (i != j) {  // Skip diagonal elements
            double Ux_ij = Dx(i, j) +
              (totalSum_x / totalDivisor) -
              (rowSums_x(i) / rowDivisor) -
              (rowSums_x(j) / rowDivisor);

            double Uy_ij = Dy(i, j) +
              (totalSum_y / totalDivisor) -
              (rowSums_y(i) / rowDivisor) -
              (rowSums_y(j) / rowDivisor);

            sum += Ux_ij * Uy_ij;
          }
        }
      }
    }

    // Return the final result for U-statistic
    return sum / n / (n - 3);

  } else {
    // Standard version (double-centering)

    // Calculate row means
    arma::rowvec means_x = arma::mean(Dx, 1).t();
    arma::rowvec means_y;

    if (same_matrices) {
      means_y = means_x; // Reuse the same means
    } else {
      means_y = arma::mean(Dy, 1).t();
    }

    // Calculate grand means
    double grand_mean_x = arma::mean(means_x);
    double grand_mean_y = same_matrices ? grand_mean_x : arma::mean(means_y);

    // Calculate double-centered product sum directly
    if (same_matrices) {
      // When x=y, we can calculate centered_x^2
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          double centered = Dx(i, j) - means_x(i) - means_x(j) + grand_mean_x;
          sum += centered * centered; // Square it
        }
      }
    } else {
      // Regular case
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          double centered_x = Dx(i, j) - means_x(i) - means_x(j) + grand_mean_x;
          double centered_y = Dy(i, j) - means_y(i) - means_y(j) + grand_mean_y;
          sum += centered_x * centered_y;
        }
      }
    }

    // Return the final result for standard version
    return sum / (n * n);
  }
}

// Calculate the cross sample distance covariance or HSIC
// [[Rcpp::export]]
double chsic_cpp(const arma::mat& D) {
  int n = D.n_rows;
  int m = D.n_cols;

  // Calculate row and column sums
  arma::vec rowSums = arma::sum(D, 1);  // Row sums
  arma::rowvec colSums = arma::sum(D, 0);  // Column sums
  double totalSum = arma::accu(D);  // Total sum

  // Calculate sum of squared centered values directly
  double sum_squared = 0.0;

  for (int i = 0; i < n; i++) {
    double row_term = rowSums(i) / m;

    for (int j = 0; j < m; j++) {
      double col_term = colSums(j) / n;
      double total_term = totalSum / (n * m);
      double centered_val = D(i, j) - row_term - col_term + total_term;

      sum_squared += centered_val * centered_val;
    }
  }

  return sum_squared / (n-1) / (m-1);
}

// Optimized Matrix V-centering function
// [[Rcpp::export]]
arma::mat matrix_v_center(arma::mat A) {
  int nrows = A.n_rows;
  int ncols = A.n_cols;

  // Calculate row and column sums
  arma::vec rowSums = arma::sum(A, 1);  // Row sums
  arma::rowvec colSums = arma::sum(A, 0);  // Column sums
  double totalSum = arma::accu(A);  // Total sum

  // Create V-centered matrix directly
  arma::mat VA(nrows, ncols);

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      // V-centering formula: (A - r - c + t)
      VA(i, j) = A(i, j) -
        (rowSums(i) / ncols) -
        (colSums(j) / nrows) +
        (totalSum / (nrows * ncols));
    }
  }

  return VA;
}

// J-HSIC V-statistic function
// [[Rcpp::export]]
double jhsic_v(Rcpp::List x, double cc = 1.0, std::string type = "euclidean",
               SEXP bw_sexp = R_NilValue, double expo = 1,
               double scale_factor = 0.5, Nullable<List> group_ = R_NilValue) {
  // Get environment
  Environment KDist = Environment::namespace_env("KDist");
  Function KDist_matrix = KDist["KDist_matrix"];
  Function matrix_v_center = KDist["matrix_v_center"];

  // Get dimensions
  NumericMatrix first_matrix = x[0];
  int n = first_matrix.nrow();
  int d = x.length();

  // Check that all matrices have the same number of rows
  for (int i = 1; i < d; i++) {
    NumericMatrix curr_matrix = x[i];
    if (curr_matrix.nrow() != n) {
      stop("Unequal Sample Sizes");
    }
  }

  // Extract group vectors from the list if provided
  List group_list;
  bool has_groups = false;

  if (group_.isNotNull()) {
    has_groups = true;
    group_list = List(group_);

    // Check that the group list has the right size
    if (group_list.size() != d) {
      stop("The 'group' parameter must be a list with one element per matrix in x");
    }
  }

  // Process first matrix
  NumericMatrix D_rcpp;
  if (has_groups) {
    D_rcpp = KDist_matrix(first_matrix, type, bw_sexp, expo, scale_factor, group_list[0]);
  } else {
    D_rcpp = KDist_matrix(first_matrix, type, bw_sexp, expo, scale_factor, R_NilValue);
  }

  arma::mat D(D_rcpp.begin(), n, n, false);

  // Transform kernel matrices to distance matrices if gaussian or laplacian
  if (type == "gaussian" || type == "laplacian" || type == "polynomial") {
    // For each element: 1 - original_value
    D = 1 - D;
  }

  arma::mat centered = as<arma::mat>(matrix_v_center(D));
  arma::mat A = cc - centered;

  // Process remaining matrices
  for (int i = 1; i < d; i++) {
    NumericMatrix curr_matrix = x[i];

    NumericMatrix D_curr_rcpp;
    if (has_groups) {
      D_curr_rcpp = KDist_matrix(curr_matrix, type, bw_sexp, expo, scale_factor, group_list[i]);
    } else {
      D_curr_rcpp = KDist_matrix(curr_matrix, type, bw_sexp, expo, scale_factor, R_NilValue);
    }

    arma::mat D_curr(D_curr_rcpp.begin(), n, n, false);

    // Transform kernel matrices to distance matrices if gaussian or laplacian
    if (type == "gaussian" || type == "laplacian" || type == "polynomial") {
      // For each element: 1 - original_value
      D_curr = 1 - D_curr;
    }

    arma::mat centered_curr = as<arma::mat>(matrix_v_center(D_curr));

    // Element-wise multiplication
    A = A % (cc - centered_curr);
  }

  // Calculate final result
  double result = accu(A) / (n * n) - pow(cc, d);
  return result;
}

// Matrix U-centering function for symmetric matrices
// [[Rcpp::export]]
arma::mat matrix_u_center(arma::mat A) {
  int n = A.n_rows;

  // Calculate row sums (no need for column sums as they are identical)
  arma::vec sums = arma::sum(A, 1);  // Row sums
  double total = arma::accu(A);  // Total sum

  // Define divisors for U-centering
  double rowDivisor = n - 2.0;
  double totalDivisor = (n - 1.0) * (n - 2.0);

  // Create U-centered matrix directly
  arma::mat UA(n, n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i != j) {  // Skip diagonal elements
        UA(i, j) = A(i, j) +
          (total / totalDivisor) -
          (sums(i) / rowDivisor) -
          (sums(j) / rowDivisor);
      } else {
        UA(i, j) = 0.0;  // Set diagonal to zero
      }
    }
  }

  return UA;
}

// J-HSIC U-statistic function with optional scaling
// [[Rcpp::export]]
double jhsic_u(Rcpp::List x, double cc = 1.0, std::string type = "euclidean",
               SEXP bw_sexp = R_NilValue, double expo = 1,
               bool scale = false, double scale_factor = 0.5,
               Nullable<List> group_ = R_NilValue) {
  // Get environment
  Environment KDist = Environment::namespace_env("KDist");
  Function KDist_matrix = KDist["KDist_matrix"];
  Function matrix_u_center = KDist["matrix_u_center"];
  Function hsic = KDist["hsic"];

  // Get dimensions
  NumericMatrix first_matrix = x[0];
  int n = first_matrix.nrow();
  int d = x.length();

  // Check that all matrices have the same number of rows
  for (int i = 1; i < d; i++) {
    NumericMatrix curr_matrix = x[i];
    if (curr_matrix.nrow() != n) {
      stop("Unequal Sample Sizes");
    }
  }

  // Extract group vectors from the list if provided
  List group_list;
  bool has_groups = false;

  if (group_.isNotNull()) {
    has_groups = true;
    group_list = List(group_);

    // Check that the group list has the right size
    if (group_list.size() != d) {
      stop("The 'group' parameter must be a list with one element per matrix in x");
    }
  }

  // Process first matrix
  NumericMatrix D_rcpp;
  if (has_groups) {
    D_rcpp = KDist_matrix(first_matrix, type, bw_sexp, expo, scale_factor, group_list[0]);
  } else {
    D_rcpp = KDist_matrix(first_matrix, type, bw_sexp, expo, scale_factor, R_NilValue);
  }

  arma::mat D(D_rcpp.begin(), n, n, false);

  // Transform kernel matrices to distance matrices if gaussian or laplacian
  if (type == "gaussian" || type == "laplacian" || type == "polynomial") {
    // For each element: 1 - original_value
    D = 1 - D;
  }

  arma::mat centered = as<arma::mat>(matrix_u_center(D));
  arma::mat A;

  if (scale) {
    // Calculate scale using hsic function
    double scale_factor_hsic;
    if (has_groups) {
      scale_factor_hsic = as<double>(hsic(first_matrix, first_matrix, type, bw_sexp, expo, scale_factor, group_list[0], true));
    } else {
      scale_factor_hsic = as<double>(hsic(first_matrix, first_matrix, type, bw_sexp, expo, scale_factor, R_NilValue, true));
    }
    // Apply scaling
    A = cc - centered/std::sqrt(scale_factor_hsic);
  } else {
    // No scaling (effectively scale = 1.0)
    A = cc - centered;
  }

  // Process remaining matrices
  for (int i = 1; i < d; i++) {
    NumericMatrix curr_matrix = x[i];

    NumericMatrix D_curr_rcpp;
    if (has_groups) {
      D_curr_rcpp = KDist_matrix(curr_matrix, type, bw_sexp, expo, scale_factor, group_list[i]);
    } else {
      D_curr_rcpp = KDist_matrix(curr_matrix, type, bw_sexp, expo, scale_factor, R_NilValue);
    }

    arma::mat D_curr(D_curr_rcpp.begin(), n, n, false);

    // Transform kernel matrices to distance matrices if gaussian or laplacian
    if (type == "gaussian" || type == "laplacian" || type == "polynomial") {
      // For each element: 1 - original_value
      D_curr = 1 - D_curr;
    }

    arma::mat centered_curr = as<arma::mat>(matrix_u_center(D_curr));

    if (scale) {
      // Calculate scale using hsic function for current matrix
      double scale_curr;
      if (has_groups) {
        scale_curr = as<double>(hsic(curr_matrix, curr_matrix, type, bw_sexp, expo, scale_factor, group_list[i], true));
      } else {
        scale_curr = as<double>(hsic(curr_matrix, curr_matrix, type, bw_sexp, expo, scale_factor, R_NilValue, true));
      }
      // Apply scaling
      A = A % (cc - centered_curr/std::sqrt(scale_curr));
    } else {
      // No scaling
      A = A % (cc - centered_curr);
    }
  }

  // Calculate final result
  double result = arma::accu(A) / n / (n-3) - std::pow(cc, d) * n / (n-3);
  return result;
}

// Function to rank a list of matrices
// [[Rcpp::export]]
List rank_list(List x) {
  // Create a copy of the input list
  List x_r = clone(x);
  int d = x.length();

  for (int j = 0; j < d; j++) {
    NumericMatrix mat = x[j];
    NumericMatrix result = x_r[j];
    int n_rows = mat.nrow();
    int n_cols = mat.ncol();

    for (int l = 0; l < n_cols; l++) {
      // Extract column
      NumericVector col = mat(_, l);

      // Create sorted copy of column for ECDF
      NumericVector sorted_col = clone(col);
      std::sort(sorted_col.begin(), sorted_col.end());

      // Calculate ECDF for each value
      for (int i = 0; i < n_rows; i++) {
        double value = col[i];

        // Find position in sorted array (binary search for efficiency)
        int count = 0;
        for (int k = 0; k < n_rows; k++) {
          if (sorted_col[k] <= value) {
            count++;
          }
        }

        // Compute empirical CDF
        result(i, l) = (double)count / n_rows;
      }
    }
  }

  return x_r;
}

// V-statistic estimate of the dHSIC
// [[Rcpp::export]]
double dhsic_fast_rcpp(List x_list,
                       String type = "gaussian",
                       SEXP bw = R_NilValue,
                       double expo = 1.0,
                       double scale_factor = 0.5,
                       Nullable<List> group_ = R_NilValue) {
  int d = x_list.size();
  // Get first matrix and dimensions
  NumericMatrix x1 = as<NumericMatrix>(x_list[0]);
  int n = x1.nrow();
  int p = x1.ncol();

  // Get environment and KDist function
  Environment KDist = Environment::namespace_env("KDist");
  Function KDist_matrix = KDist["KDist_matrix"];

  // Pre-allocate vectors for row sums and means
  NumericVector R(n, 0.0);
  double S = 0.0;

  // Extract group vectors from the list if provided
  List group_list;
  bool has_groups = false;

  if (group_.isNotNull()) {
    has_groups = true;
    group_list = List(group_);

    // Check that the group list has the right size
    if (group_list.size() != d) {
      stop("The 'group' parameter must be a list with one element per matrix in x_list");
    }
  }

  // Handle bandwidth parameter
  bool has_vector_bw = false;
  NumericVector bw_vector;
  if (bw != R_NilValue && Rf_isVector(bw) && Rf_length(bw) > 1) {
    has_vector_bw = true;
    bw_vector = as<NumericVector>(bw);
    if (bw_vector.size() != d) {
      stop("Length of bandwidth vector must match the number of matrices");
    }
  }

  // Check if we're using a distance-based type
  bool is_distance_based = (type == "euclidean" || type == "e-dist" ||
                            type == "g-dist" || type == "l-dist");

  // Calculate the first matrix with appropriate bandwidth
  SEXP bw_k = R_NilValue;
  if (has_vector_bw) {
    bw_k = wrap(bw_vector[0]); // Use first element of bw for first matrix
  } else {
    bw_k = bw; // Use the original bw (NULL or scalar)
  }

  // Get group for first matrix
  SEXP group_k = R_NilValue;
  if (has_groups) {
    group_k = group_list[0]; // Use first element of group list
  }

  NumericMatrix M;

  if (is_distance_based) {
    // Create y1 as a zero matrix with same dimensions as x1
    NumericMatrix y1(n, p);
    // Fill with zeros (already the default, but being explicit)
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < p; j++) {
        y1(i, j) = 0.0;
      }
    }

    // Create combined matrix [x1; y1]
    NumericMatrix combined(2 * n, p);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < p; j++) {
        combined(i, j) = x1(i, j);
        combined(i + n, j) = y1(i, j);
      }
    }

    // Calculate the full distance matrix
    NumericMatrix L = as<NumericMatrix>(KDist_matrix(combined, type, bw_k, expo, scale_factor, group_k));

    // Create the properly transformed matrix M
    M = NumericMatrix(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        // M = L[1:n,(n+1):(2n)] + L[(n+1):(2n),1:n] - L[1:n,1:n]
        M(i, j) = L(i, j + n) + L(i + n, j) - L(i, j);
      }
    }
  } else {
    // For kernel-based types, use the original approach
    M = as<NumericMatrix>(KDist_matrix(x1, type, bw_k, expo, scale_factor, group_k));
  }

  // Calculate row means of first matrix
  for (int i = 0; i < n; i++) {
    double row_sum = 0.0;
    for (int j = 0; j < n; j++) {
      row_sum += M(i, j);
    }
    R[i] = row_sum / n;
    S += row_sum;
  }
  S = S / (n * n);

  // Process remaining matrices
  for (int k = 1; k < d; k++) {
    // Select the appropriate bandwidth for this matrix
    if (has_vector_bw) {
      bw_k = wrap(bw_vector[k]); // Use k-th element of bw
    }

    // Select the appropriate group for this matrix
    if (has_groups) {
      group_k = group_list[k]; // Use k-th element of group list
    }

    NumericMatrix D;
    NumericMatrix xk = as<NumericMatrix>(x_list[k]);

    if (is_distance_based) {
      // Create yk as a zero matrix with same dimensions as xk
      int pk = xk.ncol();
      NumericMatrix yk(n, pk);
      // Fill with zeros
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < pk; j++) {
          yk(i, j) = 0.0;
        }
      }

      // Create combined matrix [xk; yk]
      NumericMatrix combined(2 * n, pk);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < pk; j++) {
          combined(i, j) = xk(i, j);
          combined(i + n, j) = yk(i, j);
        }
      }

      // Calculate the full distance matrix
      NumericMatrix L = as<NumericMatrix>(KDist_matrix(combined, type, bw_k, expo, scale_factor, group_k));

      // Create the properly transformed matrix D
      D = NumericMatrix(n, n);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          // D = L[1:n,(n+1):(2n)] + L[(n+1):(2n),1:n] - L[1:n,1:n]
          D(i, j) = L(i, j + n) + L(i + n, j) - L(i, j);
        }
      }
    } else {
      // For kernel-based types, use the original approach
      D = as<NumericMatrix>(KDist_matrix(xk, type, bw_k, expo, scale_factor, group_k));
    }

    // Update M and calculate row means efficiently
    double D_mean = 0.0;
    for (int i = 0; i < n; i++) {
      double row_sum = 0.0;
      for (int j = 0; j < n; j++) {
        M(i, j) *= D(i, j);
        row_sum += D(i, j);
      }
      row_sum /= n;
      D_mean += row_sum;
      R[i] *= row_sum;
    }

    // Update S with product of means
    S *= (D_mean / n);
  }

  // Calculate final terms
  double Term1 = 0.0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Term1 += M(i, j);
    }
  }
  Term1 /= (n * n);
  double Term2 = S;
  double Term3 = 0.0;
  for (int i = 0; i < n; i++) {
    Term3 += R[i];
  }
  Term3 /= n;

  return Term1 + Term2 - 2 * Term3;
}

// Mutual independence tests based on HSIC
// [[Rcpp::export]]
List mhsic_cpp(NumericMatrix x, std::string type = "euclidean",
          SEXP bw_sexp = R_NilValue, double expo = 1,
          double scale_factor = 0.5) {
  // Get environment
  Environment KDist = Environment::namespace_env("KDist");
  Function hsic_cpp_function = KDist["hsic_cpp"];

  int n = x.nrow();
  int p = x.ncol();

  // Pre-compute all self-HSIC values (diagonal elements)
  NumericVector hsic_diag(p);
  for (int i = 0; i < p; i++) {
    NumericMatrix col_i = NumericMatrix(n, 1, x(_, i).begin());
    hsic_diag[i] = as<double>(hsic_cpp_function(col_i, col_i, type, bw_sexp, expo,
                                                scale_factor, R_NilValue, true));
  }

  double num = 0.0;
  double den = 0.0;

  // Use lower triangular part only due to symmetry
  for (int i = 1; i < p; i++) {
    NumericMatrix col_i = NumericMatrix(n, 1, x(_, i).begin());
    for (int j = 0; j < i; j++) {
      NumericMatrix col_j = NumericMatrix(n, 1, x(_, j).begin());

      // Calculate HSIC between columns i and j
      double hsic_val = as<double>(hsic_cpp_function(col_i, col_j, type, bw_sexp, expo,
                                                     scale_factor, R_NilValue, true));

      // Update numerator
      num += hsic_val;

      // Update denominator using pre-computed diagonal elements
      den += hsic_diag[i] * hsic_diag[j];
    }
  }

  // Calculate the test statistic
  double stat = sqrt(n * (n - 1.0) / 2.0 / den) * num;

  // Calculate the p-value
  Environment stats = Environment::namespace_env("stats");
  Function pnorm = stats["pnorm"];
  double pval = 1.0 - as<double>(pnorm(stat));

  // Create and return a named list
  List result = List::create(
    Named("statistic") = stat,
    Named("p.value") = pval
  );

  return result;
}

// Calculating the global statistic for conditional two sample tests
// [[Rcpp::export]]
double teststatg(int n1, int n2,
                 const arma::mat& kY11, const arma::mat& kY22, const arma::mat& kY12,
                 const arma::mat& G1_11, const arma::mat& G1_12,
                 const arma::mat& G2_22, const arma::mat& G2_21,
                 const arma::vec& S1_1, const arma::vec& S1_2,
                 const arma::vec& S2_1, const arma::vec& S2_2) {

  double T1 = 0.0;
  double T2 = 0.0;
  double T3 = 0.0;
  double T = 0.0;

  double tmp = 0.0;

  for (int i = 0; i < n1; i++) {
    tmp = 0.0;
    for (int j = 0; j < n1; j++) {
      if(i != j) {
        tmp += kY11(i,j) * G1_11(i,j);
      }
    }
    T1 += tmp * S2_1(i);
  }

  for (int l = 0; l < n2; l++) {
    tmp = 0.0;
    for (int m = 0; m < n2; m++) {
      if(l != m) {
        tmp += kY22(l,m) * G2_22(l,m);
      }
    }
    T2 += tmp * S2_2(l);
  }

  for (int i = 0; i < n1; i++) {
    for (int l = 0; l < n2; l++) {
      T3 += kY12(i,l) * (G1_12(i,l) + G2_21(l,i)) * (S1_2(l) - G1_12(i,l)) * (S1_1(i) - G2_21(l,i));
    }
  }

  T = (T1 + T2 - T3) / (n1 * (n1 - 1)) / (n2 * (n2 - 1));

  return T;

}

// Calculating the local statistic for conditional two sample tests
// [[Rcpp::export]]
double teststatl(int n1, int n2,
                 const arma::mat& kY11, const arma::mat& kY22, const arma::mat& kY12,
                 const arma::vec& G1X, const arma::vec& G2X,
                 double S1_1, double S1_2, double S2_1, double S2_2) {

  double T1 = 0.0;
  double T2 = 0.0;
  double T3 = 0.0;
  double T = 0.0;

  for (int i = 0; i < (n1 - 1); i++) {
    for (int j = (i + 1); j < n1; j++) {
      T1 += kY11(i,j) * G1X(i) * G1X(j);
    }
  }
  T1 /= S2_1 / 2.0;

  for (int l = 0; l < (n2 - 1); l++) {
    for (int m = (l + 1); m < n2; m++) {
      T2 += kY22(l,m) * G2X(l) * G2X(m);
    }
  }
  T2 /= S2_2 / 2.0;

  for (int i = 0; i < n1; i++) {
    for (int l = 0; l < n2; l++) {
      T3 += kY12(i,l) * G1X(i) * G2X(l) * (S1_1 - G1X(i)) * (S1_2 - G2X(l));
    }
  }
  T3 /= S2_1 * S2_2;

  T = T1 + T2 - 2.0 * T3;

  return T;
}

// Rcpp implementation of calculate_partial_aves
// [[Rcpp::export]]
Rcpp::NumericVector calculate_partial_aves_cpp(const arma::mat& full_dist, int min_k = 4, int max_k = 0) {
  int n = full_dist.n_rows;

  // Set default max_k if not provided or if it is 0
  if (max_k <= 0) {
    max_k = n - 4;
  }

  // Validate inputs
  min_k = std::max(1, min_k);
  max_k = std::min(n, max_k);
  int num_k = max_k - min_k + 1;

  // Create a copy of the matrix with zeros on the diagonal
  arma::mat dist_copy = full_dist;
  for (int i = 0; i < n; i++) {
    dist_copy(i, i) = 0.0;
  }

  // Initialize result vector
  Rcpp::NumericVector partial_sums(num_k);

  // Calculate initial sum for min_k
  double initial_sum = 0.0;
  for (int i = 0; i < min_k; i++) {
    for (int j = 0; j < min_k; j++) {
      initial_sum += dist_copy(i, j);
    }
  }
  partial_sums[0] = initial_sum;

  // Efficiently calculate partial sums for each k
  for (int i = 1; i < num_k; i++) {
    int current_k = min_k + i;
    double additional_sum = 0.0;

    // Sum the new row and column elements
    for (int j = 0; j < current_k; j++) {
      additional_sum += dist_copy(current_k-1, j) + dist_copy(j, current_k-1);
    }

    // Subtract the diagonal element that was double-counted
    additional_sum -= dist_copy(current_k-1, current_k-1);

    // Add to previous sum
    partial_sums[i] = partial_sums[i-1] + additional_sum;

    // Normalize by k*(k-1)
    partial_sums[i-1] = partial_sums[i-1] / (min_k + i - 1) / (min_k + i - 2);
  }

  // Normalize the last element
  partial_sums[num_k-1] = partial_sums[num_k-1] / max_k / (max_k - 1);

  return partial_sums;
}

// Rcpp implementation of calculate_partial_aves_cross
// [[Rcpp::export]]
Rcpp::NumericVector calculate_partial_aves_cross_cpp(const arma::mat& full_dist, int min_k = 4, int max_k = 0) {
  int n = full_dist.n_rows;

  // Set default max_k if not provided or if it is 0
  if (max_k <= 0) {
    max_k = n - 4;
  }

  // Validate inputs
  min_k = std::max(1, min_k);
  max_k = std::min(n, max_k);
  int num_k = max_k - min_k + 1;

  // Create a copy of the matrix with zeros on the diagonal
  arma::mat dist_copy = full_dist;
  for (int i = 0; i < n; i++) {
    dist_copy(i, i) = 0.0;
  }

  // Initialize result vector
  Rcpp::NumericVector partial_sums(num_k);

  // Calculate initial sum for min_k
  double initial_sum = 0.0;
  for (int i = 0; i < min_k; i++) {
    for (int j = min_k; j < n; j++) {
      initial_sum += dist_copy(i, j);
    }
  }
  partial_sums[0] = initial_sum;

  // Efficiently calculate cross partial sums for each k
  for (int i = 1; i < num_k; i++) {
    int current_k = min_k + i;

    // Remove sum of elements in the row that is moving from right to left side
    double sum_to_remove = 0.0;
    for (int j = 0; j < current_k-1; j++) {
      sum_to_remove += dist_copy(j, current_k-1);
    }

    // Add sum of elements in the column that is moving from right to left side
    double sum_to_add = 0.0;
    for (int j = current_k; j < n; j++) {
      sum_to_add += dist_copy(current_k-1, j);
    }

    // Update partial sum
    partial_sums[i] = partial_sums[i-1] - sum_to_remove + sum_to_add;

    // Normalize previous element
    partial_sums[i-1] = partial_sums[i-1] / (min_k + i - 1) / (n - min_k - i + 1);
  }

  // Normalize the last element
  partial_sums[num_k-1] = partial_sums[num_k-1] / max_k / (n - max_k);

  return partial_sums;
}

// Rcpp implementation of hsic_recur
// [[Rcpp::export]]
Rcpp::NumericVector hsic_recur_cpp(const arma::mat& full_dist, int min_k = 4, int max_k = 0) {
  int n = full_dist.n_rows;

  // Set default max_k if not provided or if it is 0
  if (max_k <= 0) {
    max_k = n - 4;
  }

  // Validate inputs
  min_k = std::max(1, min_k);
  max_k = std::min(n, max_k);
  int num_k = max_k - min_k + 1;

  // Create a copy of the matrix with zeros on the diagonal
  arma::mat dist_copy = full_dist;
  for (int i = 0; i < n; i++) {
    dist_copy(i, i) = 0.0;
  }

  // Initialize result vectors
  Rcpp::NumericVector T1(num_k);
  Rcpp::NumericVector T2(num_k);
  Rcpp::NumericVector T3(num_k);
  Rcpp::NumericVector result(num_k);

  // Extract initial submatrix
  arma::mat A = dist_copy.submat(0, 0, min_k-1, min_k-1);
  arma::mat Asq = A * A;

  // Initial T values
  T1[0] = arma::trace(Asq);
  T2[0] = arma::accu(A);
  T3[0] = arma::accu(Asq);

  // Calculate rowsums of initial matrix
  arma::rowvec RS = arma::sum(A, 1).t();

  // Calculate values for each k
  for (int i = 1; i < num_k; i++) {
    int ii = min_k + i;

    // Get current column/row values
    arma::vec current_col = dist_copy.col(ii-1).subvec(0, ii-2);

    // Calculate intermediate values
    double s = arma::sum(current_col);
    double sq = arma::accu(arma::square(current_col));

    // Update T values
    T1[i] = T1[i-1] + 2.0 * sq;
    T2[i] = T2[i-1] + 2.0 * s;
    T3[i] = T3[i-1] + s*s + sq + 2.0 * arma::as_scalar(RS * current_col);

    // Update RS
    arma::rowvec RS_new(ii);
    for (int j = 0; j < ii-1; j++) {
      RS_new[j] = RS[j] + current_col[j];
    }
    RS_new[ii-1] = s;
    RS = RS_new;
  }

  // Calculate k_values
  arma::vec k_values = arma::linspace<arma::vec>(min_k, max_k, num_k);

  // Calculate final result
  for (int i = 0; i < num_k; i++) {
    double k = k_values[i];
    double hsic_unscale = T1[i] + T2[i]*T2[i]/((k-1)*(k-2)) - 2.0*T3[i]/(k-2);
    result[i] = hsic_unscale/(k*(k-3));
  }

  return result;
}

// Rcpp implementation of chsic_recur
// [[Rcpp::export]]
Rcpp::NumericVector chsic_recur_cpp(const arma::mat& full_dist, int min_k = 4, int max_k = 0) {
  int n = full_dist.n_rows;

  // Set default max_k if not provided
  if (max_k <= 0) {
    max_k = n - 4;
  }

  // Validate inputs
  min_k = std::max(1, min_k);
  max_k = std::min(n, max_k);
  int num_k = max_k - min_k + 1;

  // Create a copy of full_dist with zeros on diagonal
  arma::mat dist_copy = full_dist;
  for (int i = 0; i < n; i++) {
    dist_copy(i, i) = 0.0;
  }

  // Create vector equivalent to k_values <- min_k:max_k
  Rcpp::IntegerVector k_values(num_k);
  for (int i = 0; i < num_k; i++) {
    k_values[i] = min_k + i;
  }

  // Initialize result vectors
  Rcpp::NumericVector T1(num_k);
  Rcpp::NumericVector T2(num_k);
  Rcpp::NumericVector T3(num_k);
  Rcpp::NumericVector T4(num_k);

  // Extract the initial submatrix A (equivalent to R: A <- full_dist[1:min_k, -(1:min_k)])
  arma::mat A(min_k, n - min_k);
  for (int i = 0; i < min_k; i++) {
    for (int j = 0; j < n - min_k; j++) {
      A(i, j) = dist_copy(i, min_k + j);
    }
  }

  // Calculate initial squared matrices
  arma::mat Asq_1 = A * A.t();  // Equivalent to R: Asq_1 <- A %*% t(A)
  arma::mat Asq_2 = A.t() * A;  // Equivalent to R: Asq_2 <- t(A) %*% A

  // Set initial T values (R uses 1-based indexing, C++ uses 0-based)
  T1[0] = arma::trace(Asq_1);  // Equivalent to R: T1[1] <- sum(diag(Asq_1))
  T2[0] = arma::accu(A);       // Equivalent to R: T2[1] <- sum(A)
  T3[0] = arma::accu(Asq_1);   // Equivalent to R: T3[1] <- sum(Asq_1)
  T4[0] = arma::accu(Asq_2);   // Equivalent to R: T4[1] <- sum(Asq_2)

  // Initialize RS and CS vectors
  arma::rowvec RS;
  arma::rowvec CS;

  if (A.n_cols > 1) {
    // In R: RS <- rowSums(A[,-1])
    RS = arma::sum(A.cols(1, A.n_cols - 1), 1).t();
    // In R: CS <- colSums(A[,-1])
    CS = arma::sum(A.cols(1, A.n_cols - 1), 0);
  } else {
    RS = arma::zeros<arma::rowvec>(min_k);
    CS = arma::rowvec();
  }

  // Iterate through remaining k values
  for (int i = 1; i < num_k; i++) {
    int ii = k_values[i];

    // Calculate s_left: sum(full_dist[1:ii,ii])
    double s_left = 0.0;
    for (int j = 0; j < ii; j++) {
      s_left += dist_copy(j, ii - 1);
    }

    // Calculate sq_left: sum(full_dist[1:ii,ii]^2)
    double sq_left = 0.0;
    for (int j = 0; j < ii; j++) {
      sq_left += dist_copy(j, ii - 1) * dist_copy(j, ii - 1);
    }

    // Calculate s_low: sum(full_dist[ii,(ii+1):n])
    double s_low = 0.0;
    for (int j = ii; j < n; j++) {
      s_low += dist_copy(ii - 1, j);
    }

    // Calculate sq_low: sum(full_dist[ii,(ii+1):n]^2)
    double sq_low = 0.0;
    for (int j = ii; j < n; j++) {
      sq_low += dist_copy(ii - 1, j) * dist_copy(ii - 1, j);
    }

    // Update T values
    T1[i] = T1[i-1] - sq_left + sq_low;
    T2[i] = T2[i-1] - s_left + s_low;

    // Extract the vector full_dist[ii,(ii+1):n] for CS dot product
    arma::vec ii_to_n(n - ii);
    for (int j = 0; j < n - ii; j++) {
      ii_to_n(j) = dist_copy(ii - 1, ii + j);
    }

    // Calculate CS dot product (CS %*% full_dist[ii,(ii+1):n])
    double cs_dot_prod = 0.0;
    int min_len = std::min(CS.n_elem, ii_to_n.n_elem);
    for (int j = 0; j < min_len; j++) {
      cs_dot_prod += CS(j) * ii_to_n(j);
    }

    T3[i] = T3[i-1] - s_left * s_left + sq_low + 2 * cs_dot_prod;

    // Extract the vector full_dist[1:(ii-1),ii] for RS dot product
    arma::vec one_to_ii(ii - 1);
    for (int j = 0; j < ii - 1; j++) {
      one_to_ii(j) = dist_copy(j, ii - 1);
    }

    // Calculate RS dot product (RS %*% full_dist[1:(ii-1),ii])
    double rs_dot_prod = 0.0;
    min_len = std::min(RS.n_elem, one_to_ii.n_elem);
    for (int j = 0; j < min_len; j++) {
      rs_dot_prod += RS(j) * one_to_ii(j);
    }

    T4[i] = T4[i-1] - sq_left - 2 * rs_dot_prod + s_low * s_low;

    // Update RS and CS for next iteration
    if (i < num_k - 1) {
      // Update RS: RS <- c(RS - full_dist[1:(ii-1),(ii+1)], sum(full_dist[ii,(ii+2):n]))

      // Extract full_dist[1:(ii-1),(ii+1)]
      arma::vec ii_plus_1_col(ii - 1);
      for (int j = 0; j < ii - 1; j++) {
        ii_plus_1_col(j) = dist_copy(j, ii);
      }

      // Subtract from RS (only up to the minimum length)
      arma::rowvec rs_minus_col = RS;
      int min_rs_len = (rs_minus_col.n_elem < (size_t)(ii - 1)) ? rs_minus_col.n_elem : (ii - 1);
      for (int j = 0; j < min_rs_len; j++) {
        rs_minus_col(j) -= ii_plus_1_col(j);
      }

      // Calculate sum(full_dist[ii,(ii+2):n])
      double sum_ii_plus_2_to_n = 0.0;
      for (int j = ii + 1; j < n; j++) {
        sum_ii_plus_2_to_n += dist_copy(ii - 1, j);
      }

      // Combine into new RS
      arma::rowvec new_RS(ii);
      for (int j = 0; j < ii - 1; j++) {
        if (j < rs_minus_col.n_elem) {
          new_RS(j) = rs_minus_col(j);
        } else {
          new_RS(j) = 0.0;  // Pad with zeros if necessary
        }
      }
      new_RS(ii - 1) = sum_ii_plus_2_to_n;

      RS = new_RS;

      // Update CS: CS <- CS[-1] + full_dist[ii,(ii+2):n]
      // This reduces the length of CS by 1 in each iteration

      // Handle the case where CS is already empty
      if (CS.n_elem == 0) {
        // Keep CS as empty
        continue;
      }

      // Remove the first element of CS
      arma::rowvec cs_tail;
      if (CS.n_elem > 1) {
        cs_tail = CS.subvec(1, CS.n_elem - 1);
      } else {
        cs_tail = arma::rowvec();
      }

      // Extract the vector full_dist[ii,(ii+2):n]
      arma::rowvec vec_to_add;
      if (ii + 2 <= n) {
        // Get the elements from row ii, columns (ii+2) to n
        int n_elements = n - (ii + 2) + 1;
        vec_to_add.set_size(n_elements);
        for (int j = 0; j < n_elements; j++) {
          vec_to_add(j) = dist_copy(ii - 1, ii + 1 + j);
        }
      } else {
        // Handle the case where there are no elements
        vec_to_add = arma::rowvec();
      }

      // Perform element-wise addition with Rs recycling rules
      size_t result_len = std::max(cs_tail.n_elem, vec_to_add.n_elem);

      if (result_len == 0) {
        // If both vectors are empty, the result is empty
        CS = arma::rowvec();
      } else {
        // Allocate result vector
        CS.set_size(result_len);

        // Add with recycling
        for (size_t j = 0; j < result_len; j++) {
          double val_cs = (cs_tail.n_elem > 0) ? cs_tail(j % cs_tail.n_elem) : 0.0;
          double val_add = (vec_to_add.n_elem > 0) ? vec_to_add(j % vec_to_add.n_elem) : 0.0;
          CS(j) = val_cs + val_add;
        }
      }
    }
  }

  // Calculate final chsic_unscale and result
  Rcpp::NumericVector result(num_k);
  for (int i = 0; i < num_k; i++) {
    double k = k_values[i];
    double nk = n - k;
    double chsic_unscale = T1[i] + (T2[i] * T2[i]) / (k * nk) - T3[i] / k - T4[i] / nk;
    result[i] = chsic_unscale / ((k - 1) * (nk - 1));
  }

  return result;
}
