#include <Rcpp.h>
#include "Candidate.h"
#include "MultivariateInfo.h"
#include <vector>
#include <limits>

using namespace Rcpp;
using namespace changepoint;

//' Prune Candidates using Convex Hull in Multivariate Space
//'
//' This function applies the pruning algorithm from MultivariateInfo to a set of candidate
//' changepoints. The pruning uses Qhull to compute the convex hull of candidate points in
//' the space of (tau, st_1, st_2, ..., st_d), retaining only the extreme points.
//'
//' @param candidates_df A data.frame with columns:
//'   - `tau`: integer column with time indices of candidates
//'   - `dim_1`, `dim_2`, ..., `dim_d`: numeric columns for each dimension that needs to be considered for pruning. Note that columns without "dim_" will be excluded
//'
//' @return A data.frame with the same structure as `candidates_df` but containing only the
//'   pruned candidates (extreme points on the convex hull), sorted by tau.
//'
//' @details
//' The pruning algorithm works as follows:
//' 1. Constructs points in (tau, st_1, ..., st_d) space
//' 2. Identifies numerically independent dimensions using Modified Gram-Schmidt
//' 3. Computes the convex hull using Qhull with numerical stability guards
//' 4. Returns only the vertices of the convex hull, sorted by tau
//'
//' If fewer than `target_dim + 2` candidates are provided, all candidates are returned.
//'
//' For univariate analysis, the function still requires the multivariate format (e.g., `dim_1`).
//'
//' @examples
//' set.seed(123)
//'  candidates <- data.frame(
//'      tau = 1:500,
//'      dim_1 = cumsum(rnorm(500)),
//'      dim_2 = cumsum(rnorm(500))
//'  )
//'  pruned <- rcpp_prune_multivariate(candidates)
//'  nrow(pruned)
//' @export
// [[Rcpp::export]]
DataFrame rcpp_prune_multivariate(DataFrame candidates_df) {

  try {
    // Extract tau column
    IntegerVector tau_vec = candidates_df["tau"];
    int num_candidates = tau_vec.size();

    if (num_candidates == 0) {
      // Return empty dataframe with same structure
      return candidates_df;
    }

    // Determine dimensions (count of dim_* columns)
    CharacterVector col_names = candidates_df.names();
    std::vector<int> dim_col_indices;
    for (int i = 0; i < col_names.size(); ++i) {
      std::string name = as<std::string>(col_names[i]);
      if (name.substr(0, 4) == "dim_") {
        dim_col_indices.push_back(i);
      }
    }

    if (dim_col_indices.empty()) {
      Rcpp::stop("No 'dim_*' columns found in candidates_df. "
                 "Columns should be named 'dim_1', 'dim_2', etc.");
    }

    int num_dims = dim_col_indices.size();

    // Extract dimension columns
    std::vector<NumericVector> dim_vecs(num_dims);
    for (int d = 0; d < num_dims; ++d) {
      dim_vecs[d] = candidates_df[dim_col_indices[d]];
      if (dim_vecs[d].size() != num_candidates) {
        Rcpp::stop("All 'dim_*' columns must have the same length as tau");
      }
    }

    // Extract side column if present
    CharacterVector side_vec;
    bool has_side = false;
    for (int i = 0; i < col_names.size(); ++i) {
      std::string name = as<std::string>(col_names[i]);
      if (name == "side") {
        side_vec = candidates_df[i];
        has_side = true;
        break;
      }
    }

    // Convert R candidates to C++ Candidate objects
    std::vector<Candidate> cpp_candidates;
    cpp_candidates.reserve(num_candidates);

    for (int i = 0; i < num_candidates; ++i) {
      std::vector<double> dim_vals(num_dims);
      for (int d = 0; d < num_dims; ++d) {
        dim_vals[d] = dim_vecs[d][i];
      }

      std::string side_str = "right";
      if (has_side) {
        side_str = as<std::string>(side_vec[i]);
      }

      cpp_candidates.emplace_back(dim_vals, tau_vec[i], side_str);
    }

    // Create MultivariateInfo object with minimal parameters and prune immediately
    // sn is initialized with 0s (not used since anomaly_intensity is NaN by default)
    // n is set to 0 (not used)
    std::vector<double> sn_dummy(num_dims, 0.0);
    MultivariateInfo info(sn_dummy, 0);
    info.set_pruning_counter_to_zero();
    std::vector<Candidate> pruned_candidates = info.prune(cpp_candidates);

    // Convert pruned candidates back to R dataframe
    int pruned_size = pruned_candidates.size();

    IntegerVector pruned_tau(pruned_size);
    std::vector<NumericVector> pruned_dim(num_dims);
    for (int d = 0; d < num_dims; ++d) {
      pruned_dim[d] = NumericVector(pruned_size);
    }
    CharacterVector pruned_side(pruned_size);

    for (int i = 0; i < pruned_size; ++i) {
      pruned_tau[i] = pruned_candidates[i].tau;
      for (int d = 0; d < num_dims; ++d) {
        pruned_dim[d][i] = pruned_candidates[i].st[d];
      }
      pruned_side[i] = pruned_candidates[i].side;
    }

    // Build output dataframe
    DataFrame result = DataFrame::create(
      Named("tau") = pruned_tau
    );

    // Add dimension columns in order
    for (int d = 0; d < num_dims; ++d) {
      std::string col_name = "dim_" + std::to_string(d + 1);
      result.push_back(pruned_dim[d], col_name);
    }

    // Add side column if it was present in input
    if (has_side) {
      result.push_back(pruned_side, "side");
    }

    return result;

  } catch (const std::exception& e) {
    Rcpp::stop(std::string("Error in rcpp_prune_multivariate: ") + e.what());
  }
}
