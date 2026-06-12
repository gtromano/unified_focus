#pragma once

#include "Info.h"
#include "Candidate.h"
#include "ChangepointResult.h"

#include <vector>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <optional>

namespace changepoint {

// Type alias for cost functions - returns ChangepointResult
using CostFunction = std::function<ChangepointResult(const Info&, const std::vector<double>&)>;

////////////////////////////////////////////////////////////////////////////////
// Generic helper to construct a ChangepointResult from a helper that computes
// per-candidate costs. 'HelperFn' must have signature:
//   std::vector<double> helper(const std::vector<Candidate>&, const std::vector<double>& S_n, double n, const std::vector<double>& theta0)
template <typename HelperFn>
inline ChangepointResult compute_from_helper(const Info& cs,
                                             const std::vector<double>& theta0,
                                             HelperFn helper) {
  ChangepointResult result;
  result.stopping_time = cs.n();

  const auto& candidates = cs.candidates();
  if (candidates.empty()) {
    result.changepoint = std::nullopt;
    result.stat = std::nullopt;
    return result;
  }

  const std::vector<double> vals = helper(candidates, cs.sn(), cs.n(), theta0);

  if (vals.empty()) {
    result.changepoint = std::nullopt;
    result.stat = std::nullopt;
    return result;
  }

  auto max_it = std::max_element(vals.begin(), vals.end());
  int max_idx = static_cast<int>(std::distance(vals.begin(), max_it));

  result.changepoint = candidates[max_idx].tau;
  result.stat = *max_it;
  return result;
}


////////////////////////////////////////////////////////////////////////////////
// Gaussian per-candidate cost helper (keeps math localized)
inline std::vector<double> compute_costs_gaussian_helper(
    const std::vector<Candidate>& candidates,
    const std::vector<double>& S_n,
    double n,
    const std::vector<double>& theta0) {

  int K = static_cast<int>(candidates.size());
  std::vector<double> costs(K, -1e300);

  bool use_theta0 = !theta0.empty() && !std::isnan(theta0[0]);

  // Precompute term3 when no null hypothesis theta0 is provided
  double term3 = 0.0;
  if (!use_theta0 && n > 0) {
    for (double v : S_n) term3 += (v * v) / static_cast<double>(n);
  }

  for (int i = 0; i < K; ++i) {
    const Candidate& c = candidates[i];
    double tau = c.tau;
    const auto& S_i = c.st;
    auto right_len = n - tau;

    if (tau <= 0 || right_len <= 0 || n <= 0) {
      // print tau and right_len for debugging
      // Rcpp::Rcout << "Skipping candidate with tau=" << tau << " and right_len=" << right_len << std::endl;
      costs[i] = 0.0;
      continue;
    }

    double cost = 0.0;
    if (!use_theta0) {
      double term1 = 0.0;
      for (double v : S_i) term1 += v * v;
      term1 /= static_cast<double>(tau);

      double term2 = 0.0;
      for (size_t j = 0; j < S_n.size(); ++j) {
        double d = S_n[j] - S_i[j];
        term2 += d * d;
      }
      term2 /= static_cast<double>(right_len);

      cost = term1 + term2 - term3;
    } else {
      double sum = 0.0;
      for (size_t j = 0; j < S_n.size(); ++j) {
        double d = S_n[j] - S_i[j] - right_len * theta0[j];
        sum += d * d;
      }
      cost = sum / static_cast<double>(right_len);
    }

    costs[i] = std::isnan(cost) ? 0.0 : cost;
  }

  return costs;
}

////////////////////////////////////////////////////////////////////////////////
// Top-level wrapper for Gaussian
inline ChangepointResult compute_costs_gaussian(const Info& cs, const std::vector<double>& theta0) {
  return compute_from_helper(cs, theta0, compute_costs_gaussian_helper);
}

////////////////////////////////////////////////////////////////////////////////
// Poisson per-candidate cost helper
inline std::vector<double> compute_costs_poisson_helper(
    const std::vector<Candidate>& candidates,
    const std::vector<double>& S_n,
    double n,
    const std::vector<double>& theta0) {

  auto max_l = [](const std::vector<double>& st, double tau) -> double {
    if (tau <= 0) return 0.0;
    double acc = 0.0;
    bool any_nan = false;
    for (double v : st) {
      double rate = v / static_cast<double>(tau);
      if (rate <= 0.0) { // log undefined or meaningless -> mark NaN
        any_nan = true;
        continue;
      }
      double term = -v + v * std::log(rate);
      if (std::isnan(term)) { any_nan = true; }
      else acc += term;
    }
    return any_nan ? std::numeric_limits<double>::quiet_NaN() : acc;
  };

  int K = static_cast<int>(candidates.size());
  std::vector<double> costs(K, -1e300);

  bool use_theta0 = !theta0.empty() && !std::isnan(theta0[0]);

  double term3 = 0.0;
  if (!use_theta0 && n > 0) {
    term3 = max_l(S_n, n);
  }

  for (int i = 0; i < K; ++i) {
    const Candidate& c = candidates[i];
    double tau = c.tau;
    const auto& S_i = c.st;
    auto right_len = n - tau;

    if (right_len <= 0) {
      costs[i] = 0.0;
      continue;
    }

    double cost = 0.0;
    if (!use_theta0) {
      double term1 = max_l(S_i, tau);

      std::vector<double> diff(S_n.size());
      for (size_t j = 0; j < S_n.size(); ++j) diff[j] = S_n[j] - S_i[j];

      double term2 = max_l(diff, right_len);

      if (std::isnan(term1) || std::isnan(term2) || std::isnan(term3)) {
        cost = std::numeric_limits<double>::quiet_NaN();
      } else {
        cost = term1 + term2 - term3;
      }
    } else {
      std::vector<double> diff(S_n.size());
      for (size_t j = 0; j < S_n.size(); ++j) diff[j] = S_n[j] - S_i[j];

      double term2 = max_l(diff, right_len);

      bool any_nan = false;
      double null_val = 0.0;
      for (size_t j = 0; j < S_n.size(); ++j) {
        double th = (j < theta0.size()) ? theta0[j] : std::numeric_limits<double>::quiet_NaN();
        if (!(th > 0.0)) { any_nan = true; continue; } // log undefined
        double contrib = -static_cast<double>(right_len) * th + diff[j] * std::log(th);
        if (std::isnan(contrib)) any_nan = true;
        else null_val += contrib;
      }

      if (std::isnan(term2) || any_nan) {
        cost = std::numeric_limits<double>::quiet_NaN();
      } else {
        cost = term2 - null_val;
      }
    }

    costs[i] = std::isnan(cost) ? 0.0 : cost;
  }

  return costs;
}

////////////////////////////////////////////////////////////////////////////////
// Top-level wrapper for Poisson
inline ChangepointResult compute_costs_poisson(const Info& cs, const std::vector<double>& theta0) {
  return compute_from_helper(cs, theta0, compute_costs_poisson_helper);
}


////////////////////////////////////////////////////////////////////////////////
// Bernoulli per-candidate cost helper
inline std::vector<double> compute_costs_bernoulli_helper(
    const std::vector<Candidate>& candidates,
    const std::vector<double>& S_n,
    double n,
    const std::vector<double>& theta0) {

  // helper: compute max log-likelihood for a vector of successes `st`
  // where each element is number of successes in a block of length `tau`.
  auto max_l = [](const std::vector<double>& st, double tau) -> double {
    if (tau <= 0) return 0.0;
    double acc = 0.0;
    bool any_nan = false;
    for (double v : st) {
      // invalid counts
      if (v < 0.0 || v > static_cast<double>(tau)) { any_nan = true; break; }

      // handle edge cases where p_hat == 0 or p_hat == 1 to avoid log(0)
      if (v == 0.0) {
        // contribution = 0*log(0) + tau*log(1) -> 0 in the limit
        continue;
      } else if (v == static_cast<double>(tau)) {
        // contribution = tau*log(1) + 0*log(0) -> 0 in the limit
        continue;
      } else {
        double p = v / static_cast<double>(tau);
        double term = v * std::log(p) + (static_cast<double>(tau) - v) * std::log(1.0 - p);
        if (std::isnan(term)) { any_nan = true; break; }
        acc += term;
      }
    }
    return any_nan ? std::numeric_limits<double>::quiet_NaN() : acc;
  };

  int K = static_cast<int>(candidates.size());
  std::vector<double> costs(K, -1e300);

  bool use_theta0 = !theta0.empty() && !std::isnan(theta0[0]);

  // Precompute term3 (MLE log-lik for whole sample) when no null provided
  double term3 = 0.0;
  if (!use_theta0 && n > 0) {
    term3 = max_l(S_n, n);
  }

  for (int i = 0; i < K; ++i) {
    const Candidate& c = candidates[i];
    double tau = c.tau;
    const auto& S_i = c.st;
    auto right_len = n - tau;

    if (right_len <= 0) {
      costs[i] = 0.0;
      continue;
    }

    double cost = 0.0;
    if (!use_theta0) {
      double term1 = max_l(S_i, tau);

      std::vector<double> diff(S_n.size());
      for (size_t j = 0; j < S_n.size(); ++j) diff[j] = S_n[j] - S_i[j];

      double term2 = max_l(diff, right_len);

      if (std::isnan(term1) || std::isnan(term2) || std::isnan(term3)) {
        cost = std::numeric_limits<double>::quiet_NaN();
      } else {
        cost = term1 + term2 - term3;
      }
    } else {
      // compute alternative MLE log-lik (term2) and null log-lik (null_val)
      std::vector<double> diff(S_n.size());
      for (size_t j = 0; j < S_n.size(); ++j) diff[j] = S_n[j] - S_i[j];

      double term2 = max_l(diff, right_len);

      bool any_nan = false;
      double null_val = 0.0;
      for (size_t j = 0; j < S_n.size(); ++j) {
        double th = (j < theta0.size()) ? theta0[j] : std::numeric_limits<double>::quiet_NaN();
        // theta must be strictly between 0 and 1 for log to be defined in general
        if (!(th > 0.0 && th < 1.0)) { any_nan = true; continue; }

        double k = diff[j];
        double contrib = 0.0;
        // handle k == 0 or k == right_len to avoid log(0) but those contributions
        // reduce to the finite pieces:
        if (k == 0.0) {
          // contrib = 0*log(th) + right_len*log(1-th)
          contrib = static_cast<double>(right_len) * std::log(1.0 - th);
        } else if (k == static_cast<double>(right_len)) {
          // contrib = right_len*log(th) + 0*log(1-th)
          contrib = static_cast<double>(right_len) * std::log(th);
        } else {
          contrib = k * std::log(th) + (static_cast<double>(right_len) - k) * std::log(1.0 - th);
        }

        if (std::isnan(contrib)) any_nan = true;
        else null_val += contrib;
      }

      if (std::isnan(term2) || any_nan) {
        cost = std::numeric_limits<double>::quiet_NaN();
      } else {
        cost = term2 - null_val;
      }
    }

    costs[i] = std::isnan(cost) ? 0.0 : cost;
  }

  return costs;
}

////////////////////////////////////////////////////////////////////////////////
// Top-level wrapper for Bernoulli
inline ChangepointResult compute_costs_bernoulli(const Info& cs, const std::vector<double>& theta0) {
  return compute_from_helper(cs, theta0, compute_costs_bernoulli_helper);
}

////////////////////////////////////////////////////////////////////////////////
// Gamma per-candidate cost helper
// shape: user-specified known shape parameter (k > 0)
// Log-likelihood contribution used:  -right_len * shape * log(theta) - diff * (1/theta)
static inline std::vector<double> compute_costs_gamma_impl(
    const std::vector<Candidate>& candidates,
    const std::vector<double>& S_n,
    double n,
    const std::vector<double>& theta0,
    double shape) {

  if (!(shape > 0.0)) {
    throw std::invalid_argument("Gamma shape parameter must be > 0");
  }

  // helper: compute max log-likelihood (using the MLE for theta) for a vector of sums `st`
  // where each element is the sum of observations in a block of length `tau`.
  auto max_l = [shape](const std::vector<double>& st, double tau) -> double {
    if (tau <= 0) return 0.0;
    double acc = 0.0;
    bool any_nan = false;
    for (double v : st) {
      // For the parameterization assumed by the user, v must be > 0 to form theta_hat
      if (!(v > 0.0)) { any_nan = true; break; }

      // MLE for theta under the model used: theta_hat = v / (tau * shape)
      double theta_hat = v / (static_cast<double>(tau) * shape);

      if (!(theta_hat > 0.0)) { any_nan = true; break; }

      // contribution: -tau * shape * log(theta_hat) - v * (1/theta_hat)
      double term = -static_cast<double>(tau) * shape * std::log(theta_hat) - v * (1.0 / theta_hat);

      if (std::isnan(term)) { any_nan = true; break; }
      acc += term;
    }
    return any_nan ? std::numeric_limits<double>::quiet_NaN() : acc;
  };

  int K = static_cast<int>(candidates.size());
  std::vector<double> costs(K, -1e300);

  bool use_theta0 = !theta0.empty() && !std::isnan(theta0[0]);

  double term3 = 0.0;
  if (!use_theta0 && n > 0) {
    term3 = max_l(S_n, n);
  }

  for (int i = 0; i < K; ++i) {
    const Candidate& c = candidates[i];
    double tau = c.tau;
    const auto& S_i = c.st;
    auto right_len = n - tau;

    if (right_len <= 0) {
      costs[i] = 0.0;
      continue;
    }

    double cost = 0.0;
    if (!use_theta0) {
      double term1 = max_l(S_i, tau);

      std::vector<double> diff(S_n.size());
      for (size_t j = 0; j < S_n.size(); ++j) diff[j] = S_n[j] - S_i[j];

      double term2 = max_l(diff, right_len);

      if (std::isnan(term1) || std::isnan(term2) || std::isnan(term3)) {
        cost = std::numeric_limits<double>::quiet_NaN();
      } else {
        cost = term1 + term2 - term3;
      }
    } else {
      std::vector<double> diff(S_n.size());
      for (size_t j = 0; j < S_n.size(); ++j) diff[j] = S_n[j] - S_i[j];

      double term2 = max_l(diff, right_len);

      bool any_nan = false;
      double null_val = 0.0;
      for (size_t j = 0; j < S_n.size(); ++j) {
        double th = (j < theta0.size()) ? theta0[j] : std::numeric_limits<double>::quiet_NaN();

        // theta must be > 0 for the log and reciprocal to be defined
        if (!(th > 0.0)) { any_nan = true; continue; }

        double ksum = diff[j];
        // null contribution: -right_len * shape * log(th) - diff * (1/th)
        double contrib = -static_cast<double>(right_len) * shape * std::log(th) - ksum * (1.0 / th);

        if (std::isnan(contrib)) any_nan = true;
        else null_val += contrib;
      }

      if (std::isnan(term2) || any_nan) {
        cost = std::numeric_limits<double>::quiet_NaN();
      } else {
        cost = term2 - null_val;
      }
    }

    costs[i] = std::isnan(cost) ? 0.0 : cost;
  }

  return costs;
}

////////////////////////////////////////////////////////////////////////////////
// Top-level wrapper for Gamma: binds the shape parameter and calls compute_from_helper
inline ChangepointResult compute_costs_gamma(const Info& cs,
                                             const std::vector<double>& theta0,
                                             double shape) {
  // lambda captures the shape and matches the HelperFn signature expected by compute_from_helper
  auto helper = [shape](const std::vector<Candidate>& candidates,
                        const std::vector<double>& S_n,
                        double n,
                        const std::vector<double>& theta0_inner) {
    return compute_costs_gamma_impl(candidates, S_n, n, theta0_inner, shape);
  };
  return compute_from_helper(cs, theta0, helper);
}


} // namespace changepoint
