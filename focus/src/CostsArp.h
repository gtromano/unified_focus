#pragma once

#include <vector>
#include <stdexcept>
#include <cmath>
#include <optional>

#include "Info.h"
#include "ARpInfo.h"
#include "ChangepointResult.h"

namespace changepoint {

/*
  CostsArp:
  
  Strongly-typed ARP cost function that computes the maximum test statistic
  across all four detection directions (right_pos, left_pos, right_neg, left_neg).

  - compute_costs_arp_typed(const ARpInfo&)
      The real implementation: accepts only ARpInfo (typed), extracts the
      maximum test statistic computed during the last update step.

  - compute_costs_arp(const Info&, const std::vector<double>&)
      Thin wrapper matching the CostFunction alias that enforces typed input.
      Will throw std::invalid_argument if the provided Info is not an ARpInfo.
*/

inline ChangepointResult compute_costs_arp_typed(const ARpInfo& arp_info) {
  ChangepointResult out;
  out.stopping_time = arp_info.n();
  
  // Extract the max statistic and changepoint computed during update
  double stat = arp_info.max_stat();
  int cpt = arp_info.cpt();
  
  out.stat = stat;
  out.changepoint = (cpt < 0) ? std::nullopt : std::optional<int>(cpt);
  
  return out;
}

// Wrapper matching the generic CostFunction signature that enforces the typed input.
// NOTE: theta0 parameter passed here is actually mu0_arp (pre-change mean for ARP).
// It is used by the pruning logic: if provided at detector creation time via mu0_arp,
// it enables more efficient candidate filtering. Unlike gaussian/poisson/bernoulli
// families which use theta0 for null hypothesis parameters, ARP's pre-change parameter
// is tied to pruning and specified at detector creation time as mu0_arp.
inline ChangepointResult compute_costs_arp(const Info& cs,
                                           const std::vector<double>& /* theta0/mu0_arp unused */) {
  const ARpInfo* arp = dynamic_cast<const ARpInfo*>(&cs);
  if (!arp) {
    throw std::invalid_argument("compute_costs_arp: Info must be an ARpInfo (use family='arp').");
  }
  return compute_costs_arp_typed(*arp);
}

} // namespace changepoint
