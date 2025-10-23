// focus_rcpp_module.cpp
#include <Rcpp.h>
#include "Info.h"
#include "MultivariateInfo.h"
#include "Detectors.h"
#include "Costs.h"

using namespace Rcpp;
using namespace changepoint;

// // ------------------------
// // Rcpp Modules: Info types
// // ------------------------
// // Expose only R-friendly Info methods. Avoid exposing complex return types (like Candidate vectors)
// // here; use factory functions / wrappers for anything more complicated.
// RCPP_MODULE(InfoModule) {
//   class_<Info>("Info")
//   .method("update", &Info::update, "Update the cumulative state with new data")
//   .method("sn", &Info::sn, "Get the cumulative sum")
//   .method("n", &Info::n, "Get the number of observations")
//   ;
// }
//
// // Register MultivariateInfo as deriving from Info so R objects of class MultivariateInfo
// // will inherit Info methods. Do NOT expose a constructor with complex signature here.
// RCPP_MODULE(MultivariateInfoModule) {
//   class_<MultivariateInfo>("MultivariateInfo");
// }

// ------------------------
// Factories & wrappers
// ------------------------


// Detector factory that mirrors your earlier detector_create but returns an XPtr.
// family: "gaussian" or "poisson" (uses compute_costs_gaussian / compute_costs_poisson).
// Optionally, you could accept an XPtr to an Info-derived object instead of creating one.
// [[Rcpp::export]]
SEXP detector_create(std::string family, std::string info) {
  std::shared_ptr<Info> cs;
  CostFunction cost_fn;

  if (info == "multivariate") {
    cs = std::make_shared<MultivariateInfo>();
  } else if (info == "univariate") {
    cs = std::make_shared<UnivariateInfo>();
  } else {
    stop("Please either have univariate or multivariate");
  }

  if (family == "gaussian") {
    // default-constructed MultivariateInfo (uses default ctor in your header)
    cost_fn = compute_costs_gaussian;
  } else if (family == "poisson") {
    cost_fn = compute_costs_poisson;
  } else {
    stop("Unknown family: %s", family);
  }

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
