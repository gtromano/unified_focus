#pragma once

#include <vector>
#include <stdexcept>
#include <cmath>
#include <optional>

#include "Info.h"
#include "NonparametricInfo.h"
#include "Costs.h" // for compute_costs_bernoulli, ChangepointResult

namespace changepoint {

/*
  Strongly-typed NP-FOCuS cost function.

  - compute_costs_npfocus_typed(const NonparametricInfo&, const std::vector<double>&)
      The real implementation: accepts only NonparametricInfo (typed), computes
      per-sub-detector bernoulli costs and returns the SUM of per-detector max-statistics.

  - compute_costs_npfocus(const Info&, const std::vector<double>&)
      Thin wrapper matching the CostFunction alias so you can assign it into
      your existing dispatcher. The wrapper WILL THROW std::invalid_argument if
      the provided Info is not a NonparametricInfo. This prevents accidental
      invocation on the parent concatenated candidates (which would mix quantiles).
*/

inline ChangepointResult compute_costs_npfocus_typed(const NonparametricInfo& npinfo,
                                                     const std::vector<double>& theta0_vec) {
  ChangepointResult out;
  out.stopping_time = npinfo.n();

  const size_t K = npinfo.n_detectors();
  if (K == 0) {
    out.stat = std::nullopt;
    out.changepoint = std::nullopt;
    return out;
  }

  double sum_stats = 0.0;
  double max_stat = -1e300;

  bool use_provided_theta0 = !theta0_vec.empty();
  if (use_provided_theta0 && theta0_vec.size() != K) {
    throw std::invalid_argument("compute_costs_npfocus: theta0 vector must be empty or have length equal to n_detectors");
  }

  for (size_t i = 0; i < K; ++i) {
    const UnivariateInfo* sub = npinfo.sub_detector(i);
    if (!sub) continue;

    std::vector<double> theta_sub;
    if (use_provided_theta0 && !std::isnan(theta0_vec[i]))
      theta_sub.push_back(theta0_vec[i]);

    // Call bernoulli cost *on the sub-detector* (uses sub->candidates(), sub->sn(), sub->n()).
    ChangepointResult r = compute_costs_bernoulli(*sub, theta_sub);
    double s = r.stat.has_value() && std::holds_alternative<double>(*r.stat)
      ? std::get<double>(*r.stat)
        : 0.0;
    sum_stats += s;
    max_stat = std::max(max_stat, s);
  }

  out.stat = std::vector<double>{sum_stats, max_stat};
  out.changepoint = std::nullopt;
  return out;
}


// Wrapper matching the generic CostFunction signature that enforces the typed input.
// This is the function to assign into your CostFunction variable when you want NP-FOCuS,
// e.g. cost_fn = compute_costs_npfocus;
inline ChangepointResult compute_costs_npfocus(const Info& cs,
                                               const std::vector<double>& theta0) {
  const NonparametricInfo* np = dynamic_cast<const NonparametricInfo*>(&cs);
  if (!np) {
    throw std::invalid_argument("compute_costs_npfocus: Info must be a NonparametricInfo (use family='npfocus').");
  }
  return compute_costs_npfocus_typed(*np, theta0);
}

} // namespace changepoint
