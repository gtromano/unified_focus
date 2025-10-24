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
SEXP detector_create(std::string family,
                     std::string info,
                     Nullable<NumericVector> theta0 = R_NilValue,
                     Nullable<List> dim_indexes = R_NilValue,
                     int pruning_mult = 2,
                     int pruning_offset = 1,
                     std::string side = "right") {
  using namespace changepoint;

  std::shared_ptr<Info> cs;
  CostFunction cost_fn;

  // ---- Select cost function ----
  if (family == "gaussian") {
    cost_fn = compute_costs_gaussian;
  } else if (family == "poisson") {
    cost_fn = compute_costs_poisson;
  } else {
    stop("Unknown family: must be 'gaussian' or 'poisson'");
  }

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

    // make the cost function two sided with make_two_sided_cost_fn
    cost_fn = make_two_sided_cost_fn(cost_fn);

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

  // ---- Construct Detector ----
  auto detector = std::make_shared<Detector>(cs, cost_fn);
  XPtr<std::shared_ptr<Detector>> ptr(new std::shared_ptr<Detector>(detector), true);
  return ptr;
}


// Update detector with new observation vector y
// [[Rcpp::export]]
void detector_update(SEXP detector_ptr, NumericVector y) {
  XPtr<std::shared_ptr<Detector>> ptr(detector_ptr);
  if (!ptr || !(*ptr)) stop("Invalid detector pointer");
  (*ptr)->update(as<std::vector<double>>(y));
}

// Get statistic (double)
// [[Rcpp::export]]
double detector_statistic(SEXP detector_ptr) {
  XPtr<std::shared_ptr<Detector>> ptr(detector_ptr);
  if (!ptr || !(*ptr)) stop("Invalid detector pointer");
  return (*ptr)->statistic();
}

// Get changepoint result as an R list with NULL or scalar values.
// [[Rcpp::export]]
List detector_changepoint(SEXP detector_ptr) {
  XPtr<std::shared_ptr<Detector>> ptr(detector_ptr);
  if (!ptr || !(*ptr)) stop("Invalid detector pointer");

  auto result = (*ptr)->changepoint();

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

// Optional helper: extract number of pieces (size_t -> integer) for inspection.
// returns integer scalar
// [[Rcpp::export]]
int detector_pieces_len(SEXP detector_ptr) {
  XPtr<std::shared_ptr<Detector>> ptr(detector_ptr);
  if (!ptr || !(*ptr)) stop("Invalid detector pointer");
  return static_cast<int>((*ptr)->pieces().size());
}

// Optional helper: get info.n() from detector (delegates to Info::n())
// [[Rcpp::export]]
int detector_info_n(SEXP detector_ptr) {
  XPtr<std::shared_ptr<Detector>> ptr(detector_ptr);
  if (!ptr || !(*ptr)) stop("Invalid detector pointer");
  return (*ptr)->info().n();
}

// [[Rcpp::export]]
std::vector<double> detector_info_sn(SEXP detector_ptr) {
  XPtr<std::shared_ptr<Detector>> ptr(detector_ptr);
  if (!ptr || !(*ptr)) stop("Invalid detector pointer");
  return (*ptr)->info().sn();
}


// [[Rcpp::export]]
List detector_candidates(SEXP detector_ptr) {
  XPtr<std::shared_ptr<Detector>> ptr(detector_ptr);
  if (!ptr || !(*ptr)) stop("Invalid detector pointer");

  // Assume Detector::pieces() returns std::vector<Candidate>
  const auto& candidates = (*ptr)->pieces();
  const size_t K = candidates.size();

  IntegerVector tau(K);
  CharacterVector side(K);
  List st_list(K);      // each element: NumericVector for candidate.st
  List theta0_list(K);  // each element: NumericVector for candidate.theta0 or R_NilValue if empty

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

    // theta0 - return NULL if empty
    if (c.has_theta0()) {
      theta0_list[i] = NumericVector(c.theta0.begin(), c.theta0.end());
    } else {
      theta0_list[i] = R_NilValue;
    }
  }

  // Construct a data-frame-like list (list-columns are allowed)
  List out = List::create(
    Named("tau") = tau,
    Named("st")  = st_list,
    Named("theta0") = theta0_list,
    Named("side") = side
  );

  // Optionally give it class "data.frame" and set rownames if you want a proper data.frame:
  out.attr("class") = CharacterVector::create("tbl_df", "tbl", "data.frame"); // optional: tibble-like
  // set rownames
  IntegerVector rn = seq(1, static_cast<int>(K));
  out.attr("row.names") = rn;

  return out;
}
