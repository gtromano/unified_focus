// focus_rcpp_module.cpp
//
// Provides R bindings for the high-performance C++ implementation of the FOCuS
// and md-FOCuS algorithms for online and offline changepoint detection.
//
// Main interfaces:
// - `detector_create()`: Create online/sequential detector
// - `detector_update()`: Add new observations step-by-step
// - `get_statistics()`: Compute current test statistic and detection result
// - `focus_offline()`: Full C++ offline detector (fast, batch mode)

#include <Rcpp.h>
#include "Info.h"
#include "MultivariateInfo.h"
#include "ChangepointResult.h"
#include "Costs.h"
#include "NonparametricInfo.h"
#include "CostsNonparametric.h"
#include "ARpInfo.h"
#include "CostsArp.h"

using namespace Rcpp;
using namespace changepoint;

// ------------------------
// Factories & wrappers
// ------------------------

//' Create a FOCuS changepoint detector
//'
//' Creates an online (sequential) changepoint detector object that provides
//' a step-by-step interface to the FOCuS algorithm. Each call to
//' \code{\link{detector_update}()} adds new data, and \code{\link{get_statistics}()} computes
//' the current test statistic and detection result.
//'
//' @param type Character string specifying detector type. One of:
//'   \itemize{
//'     \item \code{"univariate"}: Two-sided univariate detection
//'     \item \code{"univariate_one_sided"}: One-sided univariate detection
//'     \item \code{"multivariate"}: Multivariate detection with projections
//'     \item \code{"npfocus"}: Nonparametric detection (NP-FOCuS). See details.
//'     \item \code{"arp"}: AutoRegressive Process detection. Requires \code{rho} parameter.
//'   }
//' @param dim_indexes List of integer vectors specifying projection index sets
//'   for high-dimensional multivariate detectors. Not required for sequences of dimentions less than 5.
//'    Each element is a vector of 0-based column indices. Default is \code{NULL}.
//' @param quantiles Numeric vector of quantiles for nonparametric
//'   (\code{"npfocus"}) detectors. Required when \code{type = "npfocus"}.
//'   Default is \code{NULL}.
//' @param pruning_mult Integer. Candidate pruning multiplier parameter.
//'   Default is 2.
//' @param pruning_offset Integer. Candidate pruning offset parameter.
//'   Default is 1.
//' @param side Character string. For one-sided detectors, either \code{"right"}
//'   (detects increases) or \code{"left"} (detects decreases). Default is
//'   \code{"right"}.
//' @param anomaly_intensity Numeric scalar. Anomaly intensity threshold for
//'   pruning candidates. Only candidates with sufficient signal magnitude are
//'   retained. Default is \code{NULL} (disabled).
//' @param rho Numeric vector. AR coefficients for AutoRegressive Process (ARP)
//'   detectors. Required when \code{type = "arp"}. Default is \code{NULL}.
//' @param mu0_arp Numeric scalar. Pre-change mean for ARP detectors (optional).
//'   When provided, enables more efficient pruning by filtering candidates based on
//'   the known pre-change parameter. Only used when \code{type = "arp"}.
//'   Default is \code{NULL}.
//'
//' @return An external pointer (SEXP) to the detector object. This should be
//'   passed to other detector functions like \code{\link{detector_update}()} and
//'   \code{\link{get_statistics}()}.
//'
//' @details
//' The detector maintains sufficient statistics internally and uses pruning
//' to efficiently track candidate changepoints. The \code{pruning_mult} and
//' \code{pruning_offset} parameters control the pruning strategy.
//'
//' AutoRegressive Process (ARP):
//' When \code{type = "arp"}, the \code{rho} parameter must be provided as a numeric
//' vector of AR coefficients (lag-1, lag-2, ..., lag-p). The detector then computes
//' statistics optimal for detecting changepoints in AR(p) processes. Use
//' \code{get_statistics(family = "arp")} to retrieve the test statistics.
//' The optional \code{theta0} parameter specifies the pre-change mean and is tied
//' to the pruning logic: if provided, it enables more efficient pruning by allowing
//' the algorithm to filter candidates based on the known pre-change parameter.
//'
//' High-dimentional multivariate detectors:
//' For high-dimensional multivariate detection, computing the full hull would be too prohibitive.
//' Additionally, the complexity is expected to be log(n)^p, where n is the number of iterations and p is the
//' dimensions. So for short, high-dimentional sequences, it is possible to reconstruct the set of changepoint
//' locations by approximating the hull on projections in smaller dimentions.  \code{dim_indexes} specifies which dimensions
//' to use for each projection of the convex hull for pruning. Use \code{generate_projection_indexes()} to
//' generate systematic projection sets.
//'
//' NPFOCuS:
//'   For non-parametric detection one needs to set \code{detector(type = "npfocus")} and the cost can be computed as \code{get_statistics(family = "npfocus")}.
//'   Any other cost will not work with this detector type. The \code{quantiles} vector argument is required, see \code{\link{get_statistics}()} for details.
//'
//' @examples
//' # Univariate detector
//' det <- detector_create(type = "univariate")
//' detector_update(det, 0.5)
//' detector_update(det, 1.2)
//' r <- get_statistics(det, family = "gaussian")
//' print(r)
//'
//' ## Online (sequential) example
//' # Generate data with a changepoint
//' set.seed(123)
//' Y <- c(rnorm(500, mean = 0), rnorm(500, mean = 1))
//' det <- detector_create(type = "univariate")
//' stat_trace <- numeric(length(Y))
//' threshold <- 20
//' for (i in seq_along(Y)) {
//'   detector_update(det, Y[i])
//'   r <- get_statistics(det, family = "gaussian")
//'   stat_trace[i] <- r$stat
//'   if (!is.null(r$stat) && r$stat > threshold) {
//'     cat("Online detection at", i, "estimate tau =", r$changepoint, "\n")
//'     plot(stat_trace[1:i], type = "l", ylab = "Test Statistic", xlab = "Time")
//'     break
//'   }
//' }
//'
//' # Multivariate detector with projections
//' dim_indexes <- list(c(0,1), c(1,2), c(0,2))  # 0-based indices
//' det_mv <- detector_create(type = "multivariate", dim_indexes = dim_indexes)
//' detector_update(det_mv, c(0.5, 1.2, -0.3))
//'
//' # Nonparametric detector
//' quants <- qnorm(c(0.25, 0.5, 0.75))
//' det_np <- detector_create(type = "npfocus", quantiles = quants)
//'
//' # One-sided univariate detector
//' det_one_sided <- detector_create(type = "univariate_one_sided", side = "left")
//'
//' @export
// [[Rcpp::export]]
SEXP detector_create(std::string type,
                     Nullable<List> dim_indexes = R_NilValue,
                     Nullable<NumericVector> quantiles = R_NilValue,
                     int pruning_mult = 2,
                     int pruning_offset = 1,
                     std::string side = "right",
                     Nullable<NumericVector> anomaly_intensity = R_NilValue,
                     Nullable<NumericVector> rho = R_NilValue,
                     Nullable<NumericVector> mu0_arp = R_NilValue) {
  using namespace changepoint;

  std::shared_ptr<Info> cs;

  // If quantiles provided but type != "npfocus" -> error
  if (!quantiles.isNull() && type != "npfocus") {
    warning("`quantiles` parameter provided but will be ignored unless type == \"npfocus\".");
  }

  // ---- Parse anomaly_intensity (scalar) ----
  double anomaly_intensity_val = std::numeric_limits<double>::quiet_NaN();
  if (!anomaly_intensity.isNull()) {
    NumericVector ai(anomaly_intensity.get());
    if (ai.size() >= 1) anomaly_intensity_val = ai[0];
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
      std::vector<double>{0.0},     // sn
      0,                            // n
      dim_idx,                      // now vector<vector<int>>
      pruning_mult,
      pruning_offset,
      anomaly_intensity_val
    );

  } else if (type == "univariate") {
    // two-sided
    cs = std::make_shared<UnivariateInfo>(
      0.0,            // sn
      0,              // n
      anomaly_intensity_val
    );

  } else if (type == "univariate_one_sided") {
    // one-sided
    if (side != "right" && side != "left")
      stop("side must be 'right' or 'left'");
    cs = std::make_shared<OneSideUnivariateInfo>(
      0.0,
      0,
      side,
      anomaly_intensity_val
    );

  } else if (type == "npfocus") {
    // Nonparametric NP-FOCuS: require quantiles vector
    if (quantiles.isNull()) {
      stop("type == 'npfocus' requires a NumericVector `quantiles` argument.");
    }
    NumericVector qv(quantiles.get());
    if (qv.size() < 1) stop("`quantiles` must be a NumericVector of length >= 1");

    std::vector<double> quants = as<std::vector<double>>(qv);

    // NonparametricInfo constructor
    cs = std::make_shared<NonparametricInfo>(quants);

  } else if (type == "arp") {
    // ARP (AutoRegressive Process) detector: require rho vector
    if (rho.isNull()) {
      stop("type == 'arp' requires a NumericVector `rho` argument (AR coefficients).");
    }
    NumericVector rv(rho.get());
    if (rv.size() < 1) stop("`rho` must be a NumericVector of length >= 1 (AR order p)");

    std::vector<double> rho_vec = as<std::vector<double>>(rv);

    // Check mu0_arp: distinguish between NULL (unknown) and NA (user mistake)
    bool known_prechange = false;

    if (!mu0_arp.isNull()) {
      NumericVector mu0_vec(mu0_arp.get());

      if (mu0_vec.size() < 1 || NumericVector::is_na(mu0_vec[0])) {
        stop("`mu0_arp` is NA. If the pre-change mean is unknown, set `mu0_arp = NULL`.");
      }

      known_prechange = true;
      double mu0_val = mu0_vec[0];

      cs = std::make_shared<ARpInfo>(rho_vec, known_prechange, mu0_val);  // <-- pass mu0
    } else {
      //stop("The mu0_arp unknown case is still under development, please specify a parameter `mu0_arp` when using type='arp'.");
      cs = std::make_shared<ARpInfo>(rho_vec, known_prechange);  // <-- defaults to mu0 = 0.0
    }

  } else {
    stop("type must be one of: 'multivariate', 'univariate', 'univariate_one_sided', 'npfocus', 'arp'");
  }

  // ---- Return Info pointer ----
  XPtr<std::shared_ptr<Info>> ptr(new std::shared_ptr<Info>(cs), true);
  return ptr;
}


//' Update detector with new observation(s)
//'
//' Adds new observation(s) to the detector's internal state and updates
//' sufficient statistics.
//'
//' @param det_ptr External pointer to detector created by
//'   \code{\link{detector_create}()}.
//' @param y Numeric vector of new observation(s). For univariate detectors,
//'   this should be a scalar (length-1 vector). For multivariate detectors,
//'   this should be a vector matching the number of dimensions.
//' @param lambda Numeric scalar. Rate parameter for background process (default: 1.0).
//'   Allows for non-fixed background rate. For example, use \code{lambda_i} for observation-specific rates.
//'   Default is \code{1.0} (standard CUSUM, one observation per update).
//'
//' @return
//' An external pointer to the detector (the same object that was passed in).
//' The detector is updated *in place* — no copy is made — so this return value
//' is provided only for convenience, for example when using the native R pipe
//' operator (`|>`) e.g.,
//' \code{det |> detector_update(y) |> get_statistics()}.
//'
//' @examples
//' # Univariate example
//' det <- detector_create(type = "univariate")
//' detector_update(det, 0.5)
//' detector_update(det, 1.2)
//'
//' # Multivariate example
//' det_mv <- detector_create(type = "multivariate")
//' detector_update(det_mv, c(0.5, 1.2, -0.3))
//'
//' ## Online (sequential) example
//' # Generate data with a changepoint
//' set.seed(123)
//' Y <- c(rnorm(500, mean = 0), rnorm(500, mean = 1))
//' det <- detector_create(type = "univariate")
//' stat_trace <- numeric(length(Y))
//' threshold <- 20
//' for (i in seq_along(Y)) {
//'   detector_update(det, Y[i])
//'   r <- get_statistics(det, family = "gaussian")
//'   stat_trace[i] <- r$stat
//'   if (!is.null(r$stat) && r$stat > threshold) {
//'     cat("Online detection at", i, "estimate tau =", r$changepoint, "\n")
//'     plot(stat_trace[1:i], type = "l", ylab = "Test Statistic", xlab = "Time")
//'     break
//'   }
//' }
//'
//'
//' @export
// [[Rcpp::export]]
SEXP detector_update(SEXP det_ptr, NumericVector y, double lambda = 1.0) {
  XPtr<std::shared_ptr<Info>> ptr(det_ptr);
  if (!ptr || !(*ptr)) stop("Invalid info pointer");
  (*ptr)->update(as<std::vector<double>>(y), lambda);
  return det_ptr;
}

//' Compute current changepoint statistics
//'
//' Computes the current changepoint test statistic and detection result based
//' on all observations processed so far.
//'
//' @param det_ptr External pointer to detector created by
//'   \code{\link{detector_create}()}.
//' @param family Character string specifying the distribution family:
//'   \itemize{
//'     \item \code{"gaussian"}: Gaussian (normal) distribution
//'     \item \code{"poisson"}: Poisson distribution
//'     \item \code{"bernoulli"}: Bernoulli (binary) distribution
//'     \item \code{"gamma"}: Gamma distribution (requires \code{shape} parameter)
//'     \item \code{"npfocus"}: Nonparametric detection (see details)
//'     \item \code{"arp"}: AutoRegressive Process detection (requires detector created with \code{type = "arp"})
//'   }
//' @param theta0 Numeric vector specifying the null hypothesis parameter.
//'   For univariate detectors: scalar (length-1 vector).
//'   For multivariate detectors: vector matching the number of dimensions.
//'   Default is \code{NULL}.
//' @param shape Numeric scalar. Shape parameter for gamma distribution.
//'   Required and must be positive when \code{family = "gamma"}.
//'   Default is \code{NULL}.
//'
//' @return A list with components:
//'   \item{stopping_time}{Integer. Current time index (number of observations
//'     processed).}
//'   \item{changepoint}{Integer or NULL. Detected changepoint location
//'     (1-based index), or NULL if no changepoint detected.}
//'   \item{stat}{Numeric scalar or vector. Test statistic(s). For univariate
//'     detectors, a scalar. For multivariate detectors, a vector of statistics
//'     for each projection.}
//'
//' @details
//' The function computes a log-likelihood ratio test statistic comparing the
//' null hypothesis (no change) against the alternative (a change at the optimal
//' location). The statistic is typically compared against a threshold to
//' determine if a changepoint should be declared.
//'
//' Gamma family:
//' When \code{family = "gamma"} a positive \code{shape} parameter must
//' be provided; otherwise an error is raised. Passing \code{shape} for a
//' non-gamma family will be ignored (and in some interfaces will trigger
//'                                     a warning).
//'
//' NPFOCuS:
//'   For non-parametric detection one needs to set \code{detector(type = "npfocus")} and the cost can be computed as \code{get_statistics(family = "npfocus")}.
//'   The \code{quantiles} vector argument is required. NPFOCuS returns two
//' statistics (sum and max over quantiles) as a vector; in the offline interface
//' \code{stat} will be a matrix with two columns.
//'
//' AutoRegressive Process (ARP):
//' For ARP detection, use \code{family = "arp"} with a detector created via
//' \code{detector_create(type = "arp", rho = ...)}. The AR coefficients (rho)
//' are already built into the detector at creation time. The optional \code{theta0}
//' parameter specifies the pre-change mean (if known) and is used for pruning logic;
//' if not provided (\code{NULL}), pruning operates without this information.
//' Returns a scalar test statistic optimized for detecting changepoints in AR processes.
//'
//' @examples
//'
//' ## Online (sequential) example
//' # Generate data with a changepoint
//' set.seed(123)
//' Y <- c(rnorm(500, mean = 0), rnorm(500, mean = 1))
//' det <- detector_create(type = "univariate")
//' stat_trace <- numeric(length(Y))
//' threshold <- 20
//' for (i in seq_along(Y)) {
//'   detector_update(det, Y[i])
//'   r <- get_statistics(det, family = "gaussian")
//'   stat_trace[i] <- r$stat
//'   if (!is.null(r$stat) && r$stat > threshold) {
//'     cat("Online detection at", i, "estimate tau =", r$changepoint, "\n")
//'     plot(stat_trace[1:i], type = "l", ylab = "Test Statistic", xlab = "Time")
//'     break
//'   }
//' }
//'
//' ## Note that multiple models can be tested simultaneously on the same detector
//' # as the statistics is independent of the detector state.
//' # For example, testing both Gaussian and Poisson costs
//' set.seed(2024)
//' # Generate Poisson count data with a rate change
//' Y_counts <- c(rpois(500, lambda = 10), rpois(500, lambda = 15))
//' # Compute full trajectories for comparison
//' det2 <- detector_create(type = "univariate")
//' stat_gaussian <- numeric(length(Y_counts))
//' stat_poisson <- numeric(length(Y_counts))
//'
//' for (i in seq_along(Y_counts)) {
//'   detector_update(det2, Y_counts[i])
//'   stat_gaussian[i] <- get_statistics(det2, family = "gaussian")$stat
//'   stat_poisson[i] <- get_statistics(det2, family = "poisson", theta0 = 10)$stat
//' }
//'
//' # Plot comparison
//' oldpar <- par(mfrow = c(1, 2))
//' plot(stat_gaussian, type = "l", main = "Gaussian Statistic on Poisson Data",
//'      xlab = "Time", ylab = "Statistic", lwd = 2, col = "blue")
//' abline(v = 500, col = "green", lty = 3, lwd = 2)
//'
//' plot(stat_poisson, type = "l", main = "Poisson Statistic on Poisson Data",
//'      xlab = "Time", ylab = "Statistic", lwd = 2, col = "red")
//' abline(v = 500, col = "green", lty = 3, lwd = 2)
//' par(oldpar)
//'
//'
//' @export
// [[Rcpp::export]]
List get_statistics(SEXP det_ptr,
                    std::string family,
                    Nullable<NumericVector> theta0 = R_NilValue,
                    Nullable<NumericVector> shape = R_NilValue) {   // <-- added shape

  XPtr<std::shared_ptr<Info>> ptr(det_ptr);
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


  // ---- Prepare shape (scalar) ----
  double shape_scalar = std::numeric_limits<double>::quiet_NaN();
  if (!shape.isNull()) {
    NumericVector sh(shape.get());
    if (sh.size() >= 1) shape_scalar = sh[0];
  }

  // If shape provided but family is not gamma -> warn that it will be ignored
  if (!std::isnan(shape_scalar) && family != "gamma") {
    warning("`shape` parameter provided but will be ignored unless family == \"gamma\".");
  }

  // If family is gamma, shape must be provided and > 0
  if (family == "gamma") {
    if (std::isnan(shape_scalar) || !(shape_scalar > 0.0)) {
      stop("family == 'gamma' requires a positive scalar `shape` parameter.");
    }
  }

  // ---- Select and apply cost function ----
  ChangepointResult result;

  if (family == "gaussian") {
    result = compute_costs_gaussian(cs, theta0_vec);
  } else if (family == "poisson") {
    result = compute_costs_poisson(cs, theta0_vec);
  } else if (family == "bernoulli") {
    result = compute_costs_bernoulli(cs, theta0_vec);
  } else if (family == "gamma") {
    // bind shape into a lambda that matches CostFunction signature
    auto gamma_fn = [shape_scalar](const Info& cs_inner, const std::vector<double>& th0) -> ChangepointResult {
      return compute_costs_gamma(cs_inner, th0, shape_scalar);
    };
    result = gamma_fn(cs, theta0_vec);
  } else if (family == "npfocus") {
    // typed npfocus wrapper will check the Info type and throw if mismatch
    result = compute_costs_npfocus(cs, theta0_vec);
  } else if (family == "arp") {
    // Special warning for ARP: theta0 should be set at detector creation (mu0_arp)
    if (!theta0.isNull()) {
      warning("For ARP detection, the pre-change mean should be specified at detector creation time (mu0_arp parameter in detector_create), not as theta0 in get_statistics. The theta0 parameter is ignored for ARP family.");
    }
    // typed arp wrapper will check the Info type and throw if mismatch
    result = compute_costs_arp(cs, theta0_vec);
  } else {
    stop("Unknown family: must be 'gaussian', 'poisson', 'bernoulli', 'gamma', 'npfocus' or 'arp'");
  }

  // ---- Convert to R list ----
  RObject cp = R_NilValue;
  RObject st = R_NilValue;

  if (result.changepoint.has_value()) {
    cp = wrap(result.changepoint.value());
  }
  if (result.stat.has_value()) {
    st = std::visit([](auto&& x) -> RObject { return wrap(x); }, *result.stat);
  }

  return List::create(
    Named("stopping_time") = result.stopping_time,
    Named("changepoint")   = cp,
    Named("stat")          = st
  );
}


//' Get number of candidate segments
//'
//' Returns the number of candidate changepoint segments currently tracked
//' by the detector.
//'
//' @param det_ptr External pointer to detector created by
//'   \code{\link{detector_create}()}.
//'
//' @return Integer. Number of candidate segments.
//'
//' @details
//' The FOCuS algorithm maintains a set of candidate segments that could
//' potentially contain changepoints. This number grows with time but is
//' controlled by the pruning parameters.
//'
//' @export
// [[Rcpp::export]]
int detector_cands_len(SEXP det_ptr) {
  XPtr<std::shared_ptr<Info>> ptr(det_ptr);
  if (!ptr || !(*ptr)) stop("Invalid info pointer");
  return static_cast<int>((*ptr)->candidates().size());
}

//' Get number of observations processed
//'
//' Returns the total number of observations processed by the detector.
//'
//' @param det_ptr External pointer to detector created by
//'   \code{\link{detector_create}()}.
//'
//' @return Integer. Number of observations processed (current time index).
//'
//' @export
// [[Rcpp::export]]
int detector_info_n(SEXP det_ptr) {
  XPtr<std::shared_ptr<Info>> ptr(det_ptr);
  if (!ptr || !(*ptr)) stop("Invalid info pointer");
  return (*ptr)->n();
}

//' Get cumulative sum statistic
//'
//' Returns the current cumulative sum statistic maintained by the detector.
//'
//' @param det_ptr External pointer to detector created by
//'   \code{\link{detector_create}()}.
//'
//' @return Numeric vector. Cumulative sum statistic. For univariate detectors,
//'   a scalar (length-1 vector). For multivariate detectors, a vector of
//'   length equal to the number of dimensions.
//'
//' @export
// [[Rcpp::export]]
std::vector<double> detector_info_sn(SEXP det_ptr) {
  XPtr<std::shared_ptr<Info>> ptr(det_ptr);
  if (!ptr || !(*ptr)) stop("Invalid info pointer");
  return (*ptr)->sn();
}

//' Get candidate segments
//'
//' Returns detailed information about all candidate changepoint segments
//' currently tracked by the detector.
//'
//' @param det_ptr External pointer to detector created by
//'   \code{\link{detector_create}()}.
//'
//' @return A data frame (tibble) with columns:
//'   \item{tau}{Integer vector. Candidate changepoint locations (0-based indices).}
//'   \item{st}{List of numeric vectors. Sufficient statistics for each
//'     candidate segment (e.g., cumulative sums of the data).}
//'   \item{side}{Character vector. Side indicator for each candidate
//'     (relevant for one-sided detectors).}
//'
//' @details
//' Each row represents a candidate segment from time \code{tau} to the current
//' time. The sufficient statistics in \code{st} are used to efficiently compute
//' test statistics without reprocessing the data.
//'
//' @export
// [[Rcpp::export]]
List detector_candidates(SEXP det_ptr) {
  XPtr<std::shared_ptr<Info>> ptr(det_ptr);
  if (!ptr || !(*ptr)) stop("Invalid info pointer");

  const auto& candidates = (*ptr)->candidates();
  const size_t K = candidates.size();

  // if size is 0, return empty list (as this is arp case)
  if (K == 0) {
    return List::create(
      Named("tau") = IntegerVector(0),
      Named("st") = List(0),
      Named("side") = CharacterVector(0)
    );
  }

  NumericVector tau(K);
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

//' Generate projection index sets
//'
//' Generates projection index sets for high-dimensional multivariate detectors
//' using circular combinations.
//'
//' @param d Integer. Total number of dimensions.
//' @param p Integer. Projection subset size (number of dimensions per projection).
//'
//' @return A list of integer vectors. Each element is a vector of 0-based
//'   column indices representing one projection.
//'
//' @details
//' This function generates systematic projection sets for use with multivariate
//' detectors. The circular combination approach ensures good coverage of the
//' dimensional space while keeping the number of projections manageable.
//'
//' @examples
//' \donttest{
//' # Generate 2-dimensional projections from 5 dimensions
//' proj <- generate_projection_indexes(d = 5, p = 2)
//' print(proj)
//'
//' # Use with multivariate detector
//' det <- detector_create(type = "multivariate", dim_indexes = proj)
//' set.seed(42)
//' d <- 5
//'
//' # Create data: changepoint at t=1000
//' Y_multi <- rbind(
//'     matrix(rnorm(1000 * d, mean = -1, 1), ncol = d),
//'     matrix(rnorm(500 * d, mean = 1.2), ncol = d)
//' )
//'
//' # Full multivariate detection
//' system.time(
//' res_multi <- focus_offline(Y_multi, threshold = Inf,
//'                            type = "multivariate", family = "gaussian")
//' )
//'
//' # Low-dimensional projection approximation
//' dim_indexes <- generate_projection_indexes(5, 2)
//' system.time(
//' res_multi_approx <- focus_offline(Y_multi, threshold = Inf,
//'                                   type = "multivariate", family = "gaussian",
//'                                   dim_indexes = dim_indexes)
//' )
//'
//' # Verify similarity
//' all.equal(res_multi$stat, res_multi_approx$stat)
//' }
//'
//' @export
// [[Rcpp::export]]
std::vector<std::vector<int>> generate_projection_indexes(int d, int p) {
  std::vector<std::vector<int>> combs = generate_circular_combinations(d, p);
  return combs;
}


//' @name focus_offline
//' @title Run FOCuS detector in offline batch mode
//'
//' @description
//' Processes all data at once and returns detection results and trajectories.
//' This is the most efficient way to run changepoint detection when all data
//' is available upfront.
//'
//' @param Y Numeric vector or matrix. Data array. For univariate detection,
//'   a numeric vector. For multivariate detection, a matrix with observations
//'   in rows and dimensions in columns.
//' @param threshold Numeric scalar or vector. Detection threshold(s). Can be:
//'   \itemize{
//'     \item A scalar: applied to all test statistics (default behavior for most cost functions)
//'     \item A vector: length must match the number of test statistics (for example, np-focus
//'       returns both the sum and max statistics, so threshold should be length 2).
//'     \item \code{Inf}: no thresholding (returns all statistics)
//'   }
//' @param type Character string specifying detector type. See
//'   \code{\link{detector_create}()} for options. Defaults to \code{"univariate"}.
//' @param family Character string specifying distribution family. See
//'   \code{\link{get_statistics}()} for options. Defaults to \code{"gaussian"}.
//' @param theta0 Numeric vector. Null hypothesis parameter. Default is \code{NULL}.
//' @param dim_indexes List of integer vectors. Projection index sets for
//'   multivariate detectors. Default is \code{NULL}.
//' @param quantiles Numeric vector. Quantiles for nonparametric detectors.
//'   Default is \code{NULL}.
//' @param pruning_mult Integer. Pruning multiplier parameter. Default is 2.
//' @param pruning_offset Integer. Pruning offset parameter. Default is 1.
//' @param side Character string. For one-sided detectors: \code{"right"} or
//'   \code{"left"}. Default is \code{"right"}.
//' @param shape Numeric scalar. Shape parameter for gamma distribution.
//'   Default is \code{NULL}.
//' @param anomaly_intensity Numeric scalar. Anomaly intensity threshold for
//'   pruning candidates. Only candidates with sufficient signal magnitude are
//'   retained. Default is \code{NULL} (disabled).
//' @param rho Numeric vector. AR coefficients for AutoRegressive Process (ARP)
//'   detectors. Required when \code{type = "arp"}. Default is \code{NULL}.
//' @param mu0_arp Numeric scalar. Pre-change mean for ARP detectors (optional).
//'   Only used when \code{type = "arp"}.
//'   Default is \code{NULL}.
//'
//' @return A list with components:
//'   \item{stat}{Numeric matrix. Test statistics over time (n_obs × n_stats).
//'     Each row corresponds to one time point, each column to one statistic.}
//'   \item{changepoint}{Integer vector. Detected changepoints at each time
//'     point (1-based indices), or NA if no changepoint detected at that time.}
//'   \item{detection_time}{Integer or NULL. Time of first detection (1-based),
//'     or NULL if no detection occurred.}
//'   \item{detected_changepoint}{Integer or NULL. Changepoint location at
//'     detection time (1-based), or NULL if no detection occurred.}
//'   \item{candidates}{Data frame. Final candidate segments (see
//'     \code{detector_candidates()}).}
//'   \item{threshold}{Numeric vector. Threshold(s) used for detection.}
//'   \item{n}{Integer. Number of observations processed.}
//'   \item{type}{Character. Detector type used.}
//'   \item{family}{Character. Distribution family used.}
//'   \item{shape}{Numeric or NULL. Shape parameter (for gamma family).}
//'
//' @details
//' This function runs the complete detection algorithm in C++ for maximum
//' efficiency. It processes observations sequentially and stops at the first
//' detection (when any statistic exceeds its threshold).
//'
//' For multivariate data, the algorithm computes multiple statistics (one per
//' projection). Detection occurs when ANY statistic exceeds its threshold.
//'
//' @examples
//' # Univariate Gaussian detection
//' set.seed(123)
//' Y <- c(rnorm(100, mean = 0), rnorm(100, mean = 2))
//' result <- focus_offline(Y, threshold = 10, type = "univariate",
//'                        family = "gaussian")
//' cat("Detection at time:", result$detection_time, "\n")
//' cat("Changepoint at:", result$detected_changepoint, "\n")
//'
//' # Plot statistics
//' plot(result$stat, type = "l", ylab = "Test Statistic", xlab = "Time")
//' abline(h = result$threshold, col = "red", lty = 2)
//' if (!is.null(result$detection_time)) {
//'   abline(v = result$detection_time, col = "blue", lty = 2)
//' }
//'
//' # Poisson detection
//' Y_poisson <- c(rpois(100, lambda = 2), rpois(100, lambda = 5))
//' result_poisson <- focus_offline(Y_poisson, threshold = 10,
//'                                 type = "univariate",
//'                                 family = "poisson")
//'
//' @export
// [[Rcpp::export]]
List focus_offline(SEXP Y,
                   SEXP threshold,
                   std::string type = "univariate",
                   std::string family = "gaussian",
                   Nullable<NumericVector> theta0 = R_NilValue,
                   Nullable<List> dim_indexes = R_NilValue,
                   Nullable<NumericVector> quantiles = R_NilValue,
                   int pruning_mult = 2,
                   int pruning_offset = 1,
                   std::string side = "right",
                   Nullable<NumericVector> shape = R_NilValue,
                   Nullable<NumericVector> anomaly_intensity = R_NilValue,
                   Nullable<NumericVector> rho = R_NilValue,
                   Nullable<NumericVector> mu0_arp = R_NilValue) {

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

  // ---- Parse threshold (scalar or vector) ----
  std::vector<double> threshold_vec;
  if (Rf_isNumeric(threshold)) {
    NumericVector th = as<NumericVector>(threshold);
    if (th.size() == 0) {
      stop("threshold must be a non-empty numeric vector or scalar");
    }
    threshold_vec = as<std::vector<double>>(th);
  } else {
    stop("threshold must be numeric");
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

  // If quantiles provided but type != "npfocus" -> error (require npfocus)
  if (!quantiles.isNull() && type != "npfocus") {
    stop("`quantiles` parameter provided but will be ignored unless type == \"npfocus\".");
  }

  // Auto-set family to "arp" when type == "arp"
  if (type == "arp" && family == "gaussian") {
    family = "arp";
  }

  // Warning for ARP: if both theta0 and mu0_arp are specified, or one is missing
  if (type == "arp" && (!theta0.isNull() || !mu0_arp.isNull())) {
    if (!theta0.isNull() && !mu0_arp.isNull()) {
      warning("For ARP detector: both theta0 and mu0_arp specified. Using mu0_arp (mu0_arp takes precedence for ARP initialization). theta0 will be ignored.");
    } else if (!theta0.isNull() && mu0_arp.isNull()) {
      warning("For ARP detector: theta0 provided but mu0_arp is NULL. Using theta0 as pre-change mean for ARP initialization.");
      // In this case, we need to pass theta0 to mu0_arp for detector creation
      mu0_arp = theta0;
    }
  }

  // ---- Create detector (Info object) ----
  SEXP detector_ptr = detector_create(type, dim_indexes, quantiles, pruning_mult, pruning_offset, side, anomaly_intensity, rho, mu0_arp);
  XPtr<std::shared_ptr<Info>> ptr(detector_ptr);
  if (!ptr || !(*ptr)) stop("Failed to create detector");
  std::shared_ptr<Info>& info = *ptr;

  // ---- Prepare theta0 for cost function ----
  std::vector<double> theta0_vec;
  if (!theta0.isNull()) {
    NumericVector th(theta0.get());
    if (th.size() >= 1) theta0_vec = as<std::vector<double>>(th);
  }

  // ---- Prepare shape scalar ----
  double shape_scalar = std::numeric_limits<double>::quiet_NaN();
  if (!shape.isNull()) {
    NumericVector sh(shape.get());
    if (sh.size() >= 1) shape_scalar = sh[0];
  }

  // If shape provided but family is not gamma -> warn that it will be ignored
  if (!std::isnan(shape_scalar) && family != "gamma") {
    warning("`shape` parameter provided but will be ignored unless family == \"gamma\".");
  }

  // If family is gamma, shape must be provided and > 0
  if (family == "gamma") {
    if (std::isnan(shape_scalar) || !(shape_scalar > 0.0)) {
      stop("family == 'gamma' requires a positive scalar `shape` parameter.");
    }
  }

  // ---- Select cost function ----
  CostFunction cost_fn;

  if (family == "gaussian") {
    cost_fn = compute_costs_gaussian;
  } else if (family == "poisson") {
    cost_fn = compute_costs_poisson;
  } else if (family == "bernoulli") {
    cost_fn = compute_costs_bernoulli;
  } else if (family == "gamma") {
    // bind shape into a lambda matching CostFunction
    cost_fn = [shape_scalar](const Info& cs_inner, const std::vector<double>& th0) -> ChangepointResult {
      return compute_costs_gamma(cs_inner, th0, shape_scalar);
    };
  } else if (family == "npfocus") {
    cost_fn = compute_costs_npfocus;
  } else if (family == "arp") {
    cost_fn = compute_costs_arp;
  } else {
    stop("Unknown family: must be 'gaussian', 'poisson', 'bernoulli', 'gamma', 'npfocus' or 'arp'");
  }

  // ---- Prepare reusable containers ----
  std::vector<std::vector<double>> stats_per_time;  // Each element is a vector of stats for that time
  std::vector<int> changepoints;
  stats_per_time.reserve(n_obs);
  changepoints.reserve(n_obs);

  // Track the number of statistics returned (determined from first result)
  int n_stats = -1;
  bool threshold_warning_issued = false;

  // Reusable observation vector for updates (avoid per-iteration alloc)
  std::vector<double> y_t;
  y_t.resize(static_cast<size_t>(p_dim));

  // For matrix input: precompute column base pointers for faster indexing
  std::vector<const double*> col_ptrs;
  if (is_matrix) {
    col_ptrs.resize(static_cast<size_t>(p_dim));
    for (int j = 0; j < p_dim; ++j) {
      col_ptrs[static_cast<size_t>(j)] = &Y_mat(0, j);
    }
  }

  int detection_time = NA_INTEGER;
  int detected_changepoint = NA_INTEGER;
  int actual_length = 0;

  // ---- Run online detection ----
  for (int t = 0; t < n_obs; ++t) {
    // Fill y_t in-place (no allocation)
    if (is_matrix) {
      for (int j = 0; j < p_dim; ++j) {
        y_t[static_cast<size_t>(j)] = col_ptrs[static_cast<size_t>(j)][t];
      }
    } else {
      y_t[0] = Y_vec[t];
    }

    // Update detector
    info->update(y_t);

    // Compute statistics
    ChangepointResult result = cost_fn(*info, theta0_vec);

    // Extract stats as vector
    std::vector<double> stat_vec;
    if (result.stat.has_value()) {
      std::visit([&stat_vec](auto&& arg) {
        using T = std::decay_t<decltype(arg)>;
        if constexpr (std::is_same_v<T, double>) {
          stat_vec.push_back(arg);
        } else if constexpr (std::is_same_v<T, std::vector<double>>) {
          stat_vec = arg;
        }
      }, *result.stat);
    }

    // On first iteration, determine dimensionality and validate threshold
    if (n_stats == -1) {
      n_stats = stat_vec.empty() ? 1 : static_cast<int>(stat_vec.size());

      // Validate threshold dimensions
      if (threshold_vec.size() != 1 && static_cast<int>(threshold_vec.size()) != n_stats) {
        stop("threshold must be a scalar or a vector of length %d (matching number of statistics)", n_stats);
      }

      // Issue warning if single threshold used for multiple statistics
      if (threshold_vec.size() == 1 && n_stats > 1 && !threshold_warning_issued) {
        warning("Single threshold provided for %d statistics. Using threshold = %g for all statistics.",
                n_stats, threshold_vec[0]);
        threshold_warning_issued = true;
      }
    }

    // Ensure consistent dimensionality
    if (static_cast<int>(stat_vec.size()) != n_stats) {
      if (stat_vec.empty()) {
        stat_vec.resize(n_stats, 0.0);
      } else {
        stop("Inconsistent number of statistics returned across time steps");
      }
    }

    stats_per_time.push_back(stat_vec);

    if (result.changepoint.has_value()) {
      changepoints.push_back(result.changepoint.value());
    } else {
      changepoints.push_back(NA_INTEGER);
    }

    actual_length = t + 1;

    // Check for detection: any statistic exceeding its threshold
    bool detected = false;
    for (int s = 0; s < n_stats; ++s) {
      double thresh = (threshold_vec.size() == 1) ? threshold_vec[0] : threshold_vec[s];
      if (stat_vec[s] > thresh) {
        detected = true;
        break;
      }
    }

    if (detected) {
      detection_time = t + 1; // 1-based for R
      if (result.changepoint.has_value()) detected_changepoint = result.changepoint.value();
      break;
    }
  }

  // For ARP, remove the first observation (initialization only, no stat computed)
  if (type == "arp" && !stats_per_time.empty()) {
    stats_per_time.erase(stats_per_time.begin());
    changepoints.erase(changepoints.begin());
    actual_length--;
    if (detection_time != NA_INTEGER) detection_time--;
  }

  // ---- Convert stats to matrix ----
  NumericMatrix stat_mat(actual_length, n_stats);
  for (int t = 0; t < actual_length; ++t) {
    for (int s = 0; s < n_stats; ++s) {
      stat_mat(t, s) = stats_per_time[t][s];
    }
  }

  // ---- Convert changepoints to R vector ----
  IntegerVector changepoint_vec = wrap(changepoints);

  // ---- Get candidates ----
  // For ARP, skip candidates extraction since it doesn't use the candidate structure
  List candidates_list;
  candidates_list = detector_candidates(detector_ptr);

  // ---- Return results ----
  return List::create(
    Named("stat") = stat_mat,
    Named("changepoint") = changepoint_vec,
    Named("detection_time") = detection_time == NA_INTEGER ? R_NilValue : wrap(detection_time),
    Named("detected_changepoint") = detected_changepoint == NA_INTEGER ? R_NilValue : wrap(detected_changepoint),
    Named("candidates") = candidates_list,
    Named("threshold") = wrap(threshold_vec),
    Named("n") = actual_length,
    Named("type") = type,
    Named("family") = family,
    Named("shape") = (family == "gamma" ? wrap(shape_scalar) : R_NilValue)
  );
}
