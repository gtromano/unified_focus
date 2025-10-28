// focus_rcpp_module.cpp
#include <Rcpp.h>
#include "Info.h"
#include "MultivariateInfo.h"
#include "Detectors.h"
#include "Costs.h"

using namespace Rcpp;
using namespace changepoint;

// ------------------------
// Factories & wrappers
// ------------------------

// [[Rcpp::export]]
SEXP detector_create(std::string type,
                     Nullable<NumericVector> theta0 = R_NilValue,
                     Nullable<List> dim_indexes = R_NilValue,
                     int pruning_mult = 2,
                     int pruning_offset = 1,
                     std::string side = "right") {
  using namespace changepoint;

  std::shared_ptr<Info> cs;

  // ---- Prepare theta0 ----
  std::vector<double> theta0_vec;
  double theta0_scalar = std::numeric_limits<double>::quiet_NaN();

  if (!theta0.isNull()) {
    NumericVector th(theta0.get());
    if (th.size() == 1) {
      theta0_scalar = th[0];
      theta0_vec = {th[0]};
    } else if (th.size() > 1) {
      theta0_vec = as<std::vector<double>>(th);
      theta0_scalar = th[0];
    }
  }

  // ---- Build Info object ----
  if (type == "multivariate") {
    // Parse dim_indexes if provided
    std::vector<std::vector<int>> dim_idx; // NEW: arbitrary-length projections
    if (!dim_indexes.isNull()) {
      List l(dim_indexes);
      for (int i = 0; i < l.size(); ++i) {
        // Each element may be an integer vector of any length >= 1
        IntegerVector iv = as<IntegerVector>(l[i]);
        if (iv.size() == 0)
          stop("Each dim_indexes element must be an integer vector of length >= 1");

        // convert to std::vector<int>
        std::vector<int> proj;
        proj.reserve(iv.size());
        for (int j = 0; j < iv.size(); ++j) {
          proj.push_back(static_cast<int>(iv[j]));
        }
        // Optionally: validate indices are non-negative here or within bounds later
        dim_idx.push_back(std::move(proj));
      }
    }

    cs = std::make_shared<MultivariateInfo>(
      theta0_vec,                   // may be empty or contain NaN
      std::vector<double>{0.0},     // sn
      0,                            // n
      dim_idx,                      // now vector<vector<int>>
      pruning_mult,
      pruning_offset
    );

  } else if (type == "univariate") {
    // two-sided
    cs = std::make_shared<UnivariateInfo>(
      theta0_scalar,  // defaults to NaN if not provided
      0.0,            // sn
      0               // n
    );

  } else if (type == "univariate_one_sided") {
    // one-sided
    if (side != "right" && side != "left")
      stop("side must be 'right' or 'left'");
    cs = std::make_shared<OneSideUnivariateInfo>(
      theta0_scalar,
      0.0,
      0,
      side
    );

  } else {
    stop("type must be one of: 'multivariate', 'univariate', 'univariate_one_sided'");
  }

  // ---- Return Info pointer ----
  XPtr<std::shared_ptr<Info>> ptr(new std::shared_ptr<Info>(cs), true);
  return ptr;
}




// Update detector with new observation vector y
// [[Rcpp::export]]
void detector_update(SEXP info_ptr, NumericVector y) {
  XPtr<std::shared_ptr<Info>> ptr(info_ptr);
  if (!ptr || !(*ptr)) stop("Invalid info pointer");
  (*ptr)->update(as<std::vector<double>>(y));
}

// Get statistics (changepoint result) by computing costs
// [[Rcpp::export]]
List get_statistics(SEXP info_ptr,
                    std::string family,
                    Nullable<NumericVector> theta0 = R_NilValue) {
  XPtr<std::shared_ptr<Info>> ptr(info_ptr);
  if (!ptr || !(*ptr)) stop("Invalid info pointer");

  const Info& cs = **ptr;

  // ---- Prepare theta0 for cost function ----
  std::vector<double> theta0_vec;

  if (!theta0.isNull()) {
    NumericVector th(theta0.get());
    if (th.size() >= 1) {
      theta0_vec = as<std::vector<double>>(th);
    }
  }

  // If theta0 not provided, use the one from Info (if available)
  if (theta0_vec.empty() && cs.has_theta0()) {
    theta0_vec = cs.theta0();
  }

  // ---- Select and apply cost function ----
  ChangepointResult result;

  if (family == "gaussian") {
    // Check if we need two-sided version
    result = compute_costs_gaussian(cs, theta0_vec);
  } else if (family == "poisson") {
    result = compute_costs_poisson(cs, theta0_vec);
  } else {
    stop("Unknown family: must be 'gaussian' or 'poisson'");
  }

  // ---- Convert to R list ----
  RObject cp = R_NilValue;
  RObject st = R_NilValue;

  if (result.changepoint.has_value()) {
    cp = wrap(result.changepoint.value());
  }
  if (result.stat.has_value()) {
    st = wrap(result.stat.value());
  }

  return List::create(
    Named("stopping_time") = result.stopping_time,
    Named("changepoint")   = cp,
    Named("stat")          = st
  );
}

// Get number of candidates for inspection
// [[Rcpp::export]]
int detector_pieces_len(SEXP info_ptr) {
  XPtr<std::shared_ptr<Info>> ptr(info_ptr);
  if (!ptr || !(*ptr)) stop("Invalid info pointer");
  return static_cast<int>((*ptr)->candidates().size());
}

// Get info.n()
// [[Rcpp::export]]
int detector_info_n(SEXP info_ptr) {
  XPtr<std::shared_ptr<Info>> ptr(info_ptr);
  if (!ptr || !(*ptr)) stop("Invalid info pointer");
  return (*ptr)->n();
}

// Get info.sn()
// [[Rcpp::export]]
std::vector<double> detector_info_sn(SEXP info_ptr) {
  XPtr<std::shared_ptr<Info>> ptr(info_ptr);
  if (!ptr || !(*ptr)) stop("Invalid info pointer");
  return (*ptr)->sn();
}

// Get all candidates as a data frame
// [[Rcpp::export]]
List detector_candidates(SEXP info_ptr) {
  XPtr<std::shared_ptr<Info>> ptr(info_ptr);
  if (!ptr || !(*ptr)) stop("Invalid info pointer");

  const auto& candidates = (*ptr)->candidates();
  const size_t K = candidates.size();

  IntegerVector tau(K);
  CharacterVector side(K);
  List st_list(K);      // each element: NumericVector for candidate.st

  for (size_t i = 0; i < K; ++i) {
    const Candidate& c = candidates[i];

    // tau
    tau[i] = c.tau;

    // side
    side[i] = c.side;

    // st - always return as NumericVector (even scalar case)
    if (!c.st.empty()) {
      st_list[i] = NumericVector(c.st.begin(), c.st.end());
    } else {
      // defensive: empty st -> NULL
      st_list[i] = R_NilValue;
    }
  }

  // Construct a data-frame-like list (no theta0 column anymore)
  List out = List::create(
    Named("tau") = tau,
    Named("st")  = st_list,
    Named("side") = side
  );

  // Optionally give it class "data.frame" and set rownames
  out.attr("class") = CharacterVector::create("tbl_df", "tbl", "data.frame");
  IntegerVector rn = seq(1, static_cast<int>(K));
  out.attr("row.names") = rn;

  return out;
}

// [[Rcpp::export]]
std::vector<std::vector<int>> generate_projection_indexes(int D, int p) {
  std::vector<std::vector<int>> combs = generate_circular_combinations(D, p);
  return combs;
}

// Run complete offline detection in C++ for efficiency
// [[Rcpp::export]]
List focus_offline(SEXP Y,
                   double threshold,
                   std::string type,
                   std::string family,
                   Nullable<NumericVector> theta0 = R_NilValue,
                   Nullable<List> dim_indexes = R_NilValue,
                   int pruning_mult = 2,
                   int pruning_offset = 1,
                   std::string side = "right") {

  // ---- Parse input data Y ----
  NumericMatrix Y_mat;
  NumericVector Y_vec;
  bool is_matrix = false;
  int n_obs = 0;
  int p_dim = 1;

  if (Rf_isMatrix(Y)) {
    Y_mat = as<NumericMatrix>(Y);
    n_obs = Y_mat.nrow();
    p_dim = Y_mat.ncol();
    is_matrix = true;
  } else {
    Y_vec = as<NumericVector>(Y);
    n_obs = Y_vec.size();
    p_dim = 1;
    is_matrix = false;
  }

  // ---- Validate type against data dimensions (unchanged) ----
  if (type == "multivariate" && p_dim == 1) {
    warning("type='multivariate' specified but Y is univariate (vector or single column). Consider using type='univariate' or 'univariate_one_sided'.");
  }
  if ((type == "univariate" || type == "univariate_one_sided") && p_dim > 1) {
    warning("type='%s' specified but Y is multivariate (%d columns). Consider using type='multivariate'.",
            type.c_str(), p_dim);
  }

  // check to see if the projection dimension indexes are not out of range
  if (!dim_indexes.isNull() && type == "multivariate") {
    List l(dim_indexes);
    for (int i = 0; i < l.size(); ++i) {
      IntegerVector iv = as<IntegerVector>(l[i]);
      for (int j = 0; j < iv.size(); ++j) {
        int idx = static_cast<int>(iv[j]);
        if (idx < 0 || idx >= p_dim) {
          stop("dim_indexes contains index out of range [0, %d) for column count %d",
               p_dim, p_dim);
        }
      }
    }
  }

  // ---- Create detector (Info object) ----
  SEXP detector_ptr = detector_create(type, theta0, dim_indexes, pruning_mult, pruning_offset, side);
  XPtr<std::shared_ptr<Info>> ptr(detector_ptr);
  if (!ptr || !(*ptr)) stop("Failed to create detector");
  std::shared_ptr<Info>& info = *ptr;

  // ---- Prepare theta0 for cost function ----
  std::vector<double> theta0_vec;
  if (!theta0.isNull()) {
    NumericVector th(theta0.get());
    if (th.size() >= 1) theta0_vec = as<std::vector<double>>(th); // single copy
  }
  if (theta0_vec.empty() && info->has_theta0()) {
    theta0_vec = info->theta0(); // copy from Info if present
  }

  // ---- Select cost function ----
  CostFunction cost_fn;
  const auto* maybe_uni = dynamic_cast<const UnivariateInfo*>(info.get());
  bool is_univariate_info = (maybe_uni != nullptr);

  if (family == "gaussian") {
    cost_fn = compute_costs_gaussian;
  } else if (family == "poisson") {
    cost_fn = compute_costs_poisson;
  } else {
    stop("Unknown family: must be 'gaussian' or 'poisson'");
  }

  // ---- Prepare reusable containers ----
  std::vector<double> stats;
  std::vector<int> changepoints;
  stats.reserve(n_obs);
  changepoints.reserve(n_obs);

  // Reusable observation vector for updates (avoid per-iteration alloc)
  std::vector<double> y_t;
  y_t.resize(static_cast<size_t>(p_dim));

  // For matrix input: precompute column base pointers for faster indexing
  std::vector<const double*> col_ptrs;
  if (is_matrix) {
    col_ptrs.resize(static_cast<size_t>(p_dim));
    for (int j = 0; j < p_dim; ++j) {
      col_ptrs[static_cast<size_t>(j)] = &Y_mat(0, j); // pointer to column j
    }
  }

  int detection_time = NA_INTEGER;
  int detected_changepoint = NA_INTEGER;
  int actual_length = 0;

  // ---- Run online detection ----
  for (int t = 0; t < n_obs; ++t) {
    // Fill y_t in-place (no allocation)
    if (is_matrix) {
      // column-major: accessing by column pointer + t is cache-friendly for per-row extraction
      for (int j = 0; j < p_dim; ++j) {
        y_t[static_cast<size_t>(j)] = col_ptrs[static_cast<size_t>(j)][t];
      }
    } else {
      y_t[0] = Y_vec[t];
    }

    // Update detector
    info->update(y_t);

    // Compute statistics once
    ChangepointResult result = cost_fn(*info, theta0_vec);

    // extract stat (use 0.0 if no stat)
    double stat_val = result.stat.has_value() ? result.stat.value() : 0.0;
    stats.push_back(stat_val);

    if (result.changepoint.has_value()) {
      changepoints.push_back(result.changepoint.value());
    } else {
      changepoints.push_back(NA_INTEGER);
    }

    actual_length = t + 1;

    if (stat_val > threshold) {
      detection_time = t + 1; // 1-based for R
      if (result.changepoint.has_value()) detected_changepoint = result.changepoint.value();
      break;
    }
  }

  // ---- Convert to R vectors once ----
  NumericVector stat_vec = wrap(stats);
  IntegerVector changepoint_vec = wrap(changepoints);

  // ---- Return results ----
  return List::create(
    Named("stat") = stat_vec,
    Named("changepoint") = changepoint_vec,
    Named("detection_time") = detection_time == NA_INTEGER ? R_NilValue : wrap(detection_time),
    Named("detected_changepoint") = detected_changepoint == NA_INTEGER ? R_NilValue : wrap(detected_changepoint),
    Named("candidates") = detector_candidates(detector_ptr),
    Named("threshold") = threshold,
    Named("n") = actual_length,
    Named("type") = type,
    Named("family") = family
  );
}

