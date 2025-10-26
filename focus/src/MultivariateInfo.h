#pragma once

#include "Info.h"
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullVertexSet.h>
#include <set>
#include <utility>
#include <vector>
#include <stdexcept>
#include <algorithm>

namespace changepoint {

class MultivariateInfo : public Info {
private:
    std::vector<std::pair<int, int>> dim_indexes_;  // Pairs of dimension indices for 2D projections
    int pruning_params_[2];  // (multiplier, offset)
    mutable int pruning_in_;          // Counter for pruning frequency

public:
    MultivariateInfo(const std::vector<double>& theta0 = {},
                     const std::vector<double>& sn = {0.0},
                     int n = 0,
                     const std::vector<std::pair<int, int>>& dim_indexes = {},
                     int pruning_mult = 2,
                     int pruning_offset = 1)
        : Info(sn, n, theta0),
          dim_indexes_(dim_indexes),
          pruning_in_(5) {
        pruning_params_[0] = pruning_mult;
        pruning_params_[1] = pruning_offset;
        
        // Initialize with first candidate
        auto initial = new_candidate();
        candidates_.insert(candidates_.end(), initial.begin(), initial.end());
    }

  std::vector<Candidate> prune(const std::vector<Candidate>& candidates) const override {
    int K = static_cast<int>(candidates.size());
    if (K <= 1 || pruning_in_ > 0) {
      pruning_in_--;
      return candidates;
    }

    int target_dim = static_cast<int>(sn_.size());
    if (target_dim == 0) throw std::runtime_error("sn_ must have at least one dimension.");
    if (K < target_dim + 2) return candidates; // not enough to form a full hull

    std::set<int> hull_indices;

    // Build flat_points: each point is (tau, st_0, st_1, ..., st_{d-1})
    const int point_dim = 1 + target_dim;
    std::vector<double> flat_points;
    flat_points.reserve(static_cast<size_t>(K * point_dim));
    for (int i = 0; i < K; ++i) {
      const auto& c = candidates[i];
      flat_points.push_back(static_cast<double>(c.tau));
      // write st values (padding with zeros if necessary)
      const auto& st_vec = c.st;
      for (int j = 0; j < target_dim; ++j) {
        double v = (j < static_cast<int>(st_vec.size())) ? st_vec[j] : 0.0;
        flat_points.push_back(v);
      }
    }

    if (dim_indexes_.empty()) {
      // Full-dimensional hull via Qhull
      try {
        orgQhull::Qhull qhull;
        qhull.runQhull("", point_dim, K, flat_points.data(), "");
        for (const auto& vertex : qhull.vertexList()) {
          hull_indices.insert(vertex.point().id());
        }
      } catch (...) {
        for (int i = 0; i < K; ++i) hull_indices.insert(i);
      }
    } else {
      // For each 2D projection, build a small flat array (tau + two dims) and run Qhull
      for (const auto& pr : dim_indexes_) {
        int d1 = pr.first;
        int d2 = pr.second;
        std::vector<double> flat_proj;
        flat_proj.reserve(static_cast<size_t>(K * 3));
        for (int i = 0; i < K; ++i) {
          // idx into flat_points for this row = i * point_dim
          size_t base = static_cast<size_t>(i) * point_dim;
          double tau_v = flat_points[base];                   // tau
          double v1 = flat_points[base + 1 + d1];             // st[d1]
          double v2 = flat_points[base + 1 + d2];             // st[d2]
          flat_proj.push_back(tau_v);
          flat_proj.push_back(v1);
          flat_proj.push_back(v2);
        }
        try {
          orgQhull::Qhull qhull;
          qhull.runQhull("", 3, K, flat_proj.data(), "");
          for (const auto& vertex : qhull.vertexList()) {
            hull_indices.insert(vertex.point().id());
          }
        } catch (...) {
          for (int i = 0; i < K; ++i) hull_indices.insert(i);
        }
      }
    }

    // Build pruned list in tau order
    std::vector<Candidate> pruned;
    pruned.reserve(hull_indices.size());
    for (int idx : hull_indices) pruned.push_back(candidates[idx]);
    std::sort(pruned.begin(), pruned.end(), [](const Candidate& a, const Candidate& b){ return a.tau < b.tau; });

    int pruned_size = static_cast<int>(pruned.size());
    pruning_in_ = pruned_size * pruning_params_[0] + pruning_params_[1];
    return pruned;
  }
};

} // namespace changepoint