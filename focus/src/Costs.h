#pragma once

#include "Info.h"
#include "Candidate.h"
#include "Detectors.h"

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
//   std::vector<double> helper(const std::vector<Candidate>&, const std::vector<double>& S_n, int n, const std::vector<double>& theta0)
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
    int n,
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
    int tau = c.tau;
    const auto& S_i = c.st;
    int right_len = n - tau;

    if (tau <= 0 || right_len <= 0 || n <= 0) {
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
// Poisson per-candidate cost helper
inline std::vector<double> compute_costs_poisson_helper(
    const std::vector<Candidate>& candidates,
    const std::vector<double>& S_n,
    int n,
    const std::vector<double>& theta0) {

  auto max_l = [](const std::vector<double>& st, int tau) -> double {
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
    int tau = c.tau;
    const auto& S_i = c.st;
    int right_len = n - tau;

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
// Top-level cost functions using compute_from_helper
inline ChangepointResult compute_costs_gaussian(const Info& cs, const std::vector<double>& theta0) {
  return compute_from_helper(cs, theta0, compute_costs_gaussian_helper);
}

inline ChangepointResult compute_costs_poisson(const Info& cs, const std::vector<double>& theta0) {
  return compute_from_helper(cs, theta0, compute_costs_poisson_helper);
}

////////////////////////////////////////////////////////////////////////////////
// Two-sided cost wrapper (reuses the same helpers)
inline CostFunction make_two_sided_cost_fn(
    const std::function<std::vector<double>(const std::vector<Candidate>&,
                                            const std::vector<double>&,
                                            int,
                                            const std::vector<double>&)>& base_cost_helper) {

  return [base_cost_helper](const Info& cs, const std::vector<double>& theta0) -> ChangepointResult {
    ChangepointResult result;
    result.stopping_time = cs.n();

    const auto& candidates = cs.candidates();
    if (candidates.empty()) {
      result.changepoint = std::nullopt;
      result.stat = std::nullopt;
      return result;
    }

    const auto* uni_cs = dynamic_cast<const UnivariateInfo*>(&cs);
    if (!uni_cs) {
      throw std::runtime_error("make_two_sided_cost_fn requires UnivariateInfo");
    }

    std::vector<double> costs(candidates.size(), -1e300);

    for (size_t i = 0; i < candidates.size(); ++i) {
      const auto& c = candidates[i];
      std::vector<Candidate> single = { c };

      if (c.side == "left") {
        auto v = base_cost_helper(single, uni_cs->left()->sn(), uni_cs->left()->n(), theta0);
        costs[i] = (v.empty() ? 0.0 : v[0]);
      } else {
        auto v = base_cost_helper(single, uni_cs->right()->sn(), uni_cs->right()->n(), theta0);
        costs[i] = (v.empty() ? 0.0 : v[0]);
      }
    }

    auto max_it = std::max_element(costs.begin(), costs.end());
    int max_idx = static_cast<int>(std::distance(costs.begin(), max_it));

    result.changepoint = candidates[max_idx].tau;
    result.stat = *max_it;
    return result;
  };
}


inline CostFunction make_two_sided_gaussian() {
  return make_two_sided_cost_fn(compute_costs_gaussian_helper);
}

inline CostFunction make_two_sided_poisson() {
  return make_two_sided_cost_fn(compute_costs_poisson_helper);
}

} // namespace changepoint
