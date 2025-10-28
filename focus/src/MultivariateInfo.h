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

class MultivariateInfo : public CandidateListInfo {
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
        : CandidateListInfo(sn, n, theta0),
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
      const auto& st_vec = c.st;
      for (int j = 0; j < target_dim; ++j) {
        double v = (j < static_cast<int>(st_vec.size())) ? st_vec[j] : 0.0;
        flat_points.push_back(v);
      }
    }

    // cheap helper: count how many coordinates actually vary above `eps`
    auto count_varying_dims = [&](const double* data, int dims, int rows, double eps) -> int {
      for (int d = 0; d < dims; ++d) {
        double minv = std::numeric_limits<double>::infinity();
        double maxv = -std::numeric_limits<double>::infinity();
        for (int r = 0; r < rows; ++r) {
          double v = data[r * dims + d];
          if (v < minv) minv = v;
          if (v > maxv) maxv = v;
        }
        // store min/max per dim externally if you want; here we just test
      }
      // We'll compute count outside because we need per-dim min/max; do the explicit version below.
      return 0; // not used
    };

    // Use a small eps relative to the data scale to detect "constant" coordinates.
    // For tau we expect it to vary; choose eps = 1e-12 * (1 + span)
    auto compute_span_per_dim = [&](const double* data, int dims, int rows, std::vector<double>& spans) {
      spans.assign(dims, 0.0);
      for (int d = 0; d < dims; ++d) {
        double minv = std::numeric_limits<double>::infinity();
        double maxv = -std::numeric_limits<double>::infinity();
        for (int r = 0; r < rows; ++r) {
          double v = data[r * dims + d];
          if (v < minv) minv = v;
          if (v > maxv) maxv = v;
        }
        spans[d] = maxv - minv;
      }
    };

    // --- Full-dimensional hull via Qhull, guarded for degeneracy ---
    if (dim_indexes_.empty()) {
      // compute spans for all coords (tau + st_i)
      std::vector<double> spans;
      compute_span_per_dim(flat_points.data(), point_dim, K, spans);

      // set eps per-dim relative to span (tiny absolute floor too)
      std::vector<double> epss(point_dim);
      for (int d = 0; d < point_dim; ++d) {
        epss[d] = std::max(1e-12, 1e-12 * (1.0 + spans[d]));
      }
      int varying = 0;
      for (int d = 0; d < point_dim; ++d) if (spans[d] > epss[d]) ++varying;

      if (varying <= 1) {
        // effectively 1D: only one coordinate varies (likely only tau). Qhull would fail.
        // Cheap and safe fallback: consider all candidates as "on the hull"
        for (int i = 0; i < K; ++i) hull_indices.insert(i);
      } else if (varying == 2) {
        // effectively 2D: compress to 2D and run 2D Qhull (safer)
        // find the two dims that vary
        int d0 = -1, d1 = -1;
        for (int d = 0; d < point_dim; ++d) {
          if (spans[d] > epss[d]) {
            if (d0 == -1) d0 = d; else d1 = d;
          }
        }
        std::vector<double> flat2;
        flat2.reserve(static_cast<size_t>(K * 2));
        for (int i = 0; i < K; ++i) {
          flat2.push_back(flat_points[i * point_dim + d0]);
          flat2.push_back(flat_points[i * point_dim + d1]);
        }
        try {
          orgQhull::Qhull qhull;
          qhull.runQhull("", 2, K, flat2.data(), "");
          for (const auto& vertex : qhull.vertexList()) hull_indices.insert(vertex.point().id());
        } catch (...) {
          for (int i = 0; i < K; ++i) hull_indices.insert(i);
        }
      } else {
        // full or >=3 effective dims -> safe to call Qhull in the full dimension
        try {
          orgQhull::Qhull qhull;
          qhull.runQhull("", point_dim, K, flat_points.data(), "");
          for (const auto& vertex : qhull.vertexList()) hull_indices.insert(vertex.point().id());
        } catch (...) {
          for (int i = 0; i < K; ++i) hull_indices.insert(i);
        }
      }
    } else {
      // For each 2D projection, build a small flat array (tau + two dims) and run Qhull, but guard degeneracy
      for (const auto& pr : dim_indexes_) {
        int d1 = pr.first;
        int d2 = pr.second;
        // projection dims are indices into st: resulting projection dims (in flat_proj) are:
        // col0 = tau (flat_points[base + 0])
        // col1 = st[d1] (flat_points[base + 1 + d1])
        // col2 = st[d2] (flat_points[base + 1 + d2])
        std::vector<double> flat_proj;
        flat_proj.reserve(static_cast<size_t>(K * 3));
        for (int i = 0; i < K; ++i) {
          size_t base = static_cast<size_t>(i) * point_dim;
          double tau_v = flat_points[base + 0];
          double v1 = flat_points[base + 1 + d1];
          double v2 = flat_points[base + 1 + d2];
          flat_proj.push_back(tau_v);
          flat_proj.push_back(v1);
          flat_proj.push_back(v2);
        }

        // compute spans for the 3 coords
        std::vector<double> spans3(3, 0.0);
        for (int d = 0; d < 3; ++d) {
          double minv = std::numeric_limits<double>::infinity();
          double maxv = -std::numeric_limits<double>::infinity();
          for (int r = 0; r < K; ++r) {
            double v = flat_proj[r * 3 + d];
            if (v < minv) minv = v;
            if (v > maxv) maxv = v;
          }
          spans3[d] = maxv - minv;
        }
        // eps per dim
        std::vector<double> eps3(3);
        for (int d = 0; d < 3; ++d) eps3[d] = std::max(1e-12, 1e-12 * (1.0 + spans3[d]));
        int varying3 = 0;
        for (int d = 0; d < 3; ++d) if (spans3[d] > eps3[d]) ++varying3;

        if (varying3 <= 1) {
          // effectively 1D -> skip qhull, include all
          for (int i = 0; i < K; ++i) hull_indices.insert(i);
          continue;
        } else if (varying3 == 2) {
          // find which two dims vary and build a 2D array
          int a = -1, b = -1;
          for (int d = 0; d < 3; ++d) if (spans3[d] > eps3[d]) { if (a == -1) a = d; else b = d; }
          std::vector<double> flat2;
          flat2.reserve(static_cast<size_t>(K * 2));
          for (int i = 0; i < K; ++i) {
            flat2.push_back(flat_proj[i * 3 + a]);
            flat2.push_back(flat_proj[i * 3 + b]);
          }
          try {
            orgQhull::Qhull qhull;
            qhull.runQhull("", 2, K, flat2.data(), "");
            for (const auto& vertex : qhull.vertexList()) hull_indices.insert(vertex.point().id());
          } catch (...) {
            for (int i = 0; i < K; ++i) hull_indices.insert(i);
          }
        } else {
          // full 3D projection
          try {
            orgQhull::Qhull qhull;
            qhull.runQhull("", 3, K, flat_proj.data(), "");
            for (const auto& vertex : qhull.vertexList()) hull_indices.insert(vertex.point().id());
          } catch (...) {
            for (int i = 0; i < K; ++i) hull_indices.insert(i);
          }
        }
      }
    }

    // Build pruned list in tau order
    std::vector<Candidate> pruned;
    pruned.reserve(hull_indices.size());
    for (int idx : hull_indices) pruned.push_back(candidates[idx]);
    std::sort(pruned.begin(), pruned.end(),
              [](const Candidate& a, const Candidate& b){ return a.tau < b.tau; });

    int pruned_size = static_cast<int>(pruned.size());
    pruning_in_ = pruned_size * pruning_params_[0] + pruning_params_[1];
    return pruned;
  }
};

} // namespace changepoint
