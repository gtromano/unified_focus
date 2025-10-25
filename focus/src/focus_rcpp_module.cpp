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
SEXP detector_create(std::string info,
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
  if (info == "multivariate") {
    // Parse dim_indexes if provided
    std::vector<std::pair<int,int>> dim_idx;
    if (!dim_indexes.isNull()) {
      List l(dim_indexes);
      for (int i = 0; i < l.size(); ++i) {
        IntegerVector pair = as<IntegerVector>(l[i]);
        if (pair.size() != 2)
          stop("Each dim_indexes element must be a vector of length 2");
        // assume 0-based; adjust if you want 1-based from R
        dim_idx.emplace_back(pair[0], pair[1]);
      }
    }

    cs = std::make_shared<MultivariateInfo>(
      theta0_vec,           // may be empty or contain NaN
      std::vector<double>{0.0}, // sn
      0,                    // n
      dim_idx,
      pruning_mult,
      pruning_offset
    );

  } else if (info == "univariate") {
    // two-sided
    cs = std::make_shared<UnivariateInfo>(
      theta0_scalar,  // defaults to NaN if not provided
      0.0,            // sn
      0               // n
    );

  } else if (info == "univariate_one_sided") {
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
    stop("info must be one of: 'multivariate', 'univariate', 'univariate_one_sided'");
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
    const auto* uni_cs = dynamic_cast<const UnivariateInfo*>(&cs);
    if (uni_cs) {
      // Two-sided univariate
      auto cost_fn = make_two_sided_gaussian();
      result = cost_fn(cs, theta0_vec);
    } else {
      // One-sided or multivariate
      result = compute_costs_gaussian(cs, theta0_vec);
    }
  } else if (family == "poisson") {
    // Check if we need two-sided version
    const auto* uni_cs = dynamic_cast<const UnivariateInfo*>(&cs);
    if (uni_cs) {
      // Two-sided univariate
      auto cost_fn = make_two_sided_poisson();
      result = cost_fn(cs, theta0_vec);
    } else {
      // One-sided or multivariate
      result = compute_costs_poisson(cs, theta0_vec);
    }
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
