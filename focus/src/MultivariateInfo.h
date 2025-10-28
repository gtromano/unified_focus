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
  std::vector<std::vector<int>> dim_indexes_;  // Projections of arbitrary p size
  int pruning_params_[2];  // (multiplier, offset)
    mutable int pruning_in_;          // Counter for pruning frequency

public:
    MultivariateInfo(const std::vector<double>& theta0 = {},
                     const std::vector<double>& sn = {0.0},
                     int n = 0,
                     const std::vector<std::vector<int>>& dim_indexes = {},
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

    // helper: compute spans (max - min) per column
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
        // effectively 1D: only one coordinate varies (likely only tau). Qhull would fail. Retain all.
        for (int i = 0; i < K; ++i) hull_indices.insert(i);
      } else if (varying == 2) {
        // compress to 2D and run 2D Qhull (find the two varying dims)
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
        // full or >=3 effective dims -> call Qhull in the full dimension
        try {
          orgQhull::Qhull qhull;
          qhull.runQhull("", point_dim, K, flat_points.data(), "");
          for (const auto& vertex : qhull.vertexList()) hull_indices.insert(vertex.point().id());
        } catch (...) {
          for (int i = 0; i < K; ++i) hull_indices.insert(i);
        }
      }
    } else {
      // Generalised projection loop: each 'pr' is a vector<int> of indices into st (0-based)
      for (const auto& pr : dim_indexes_) {
        const int p = static_cast<int>(pr.size());
        if (p <= 0) continue; // skip empty projections

        const int proj_cols = 1 + p; // tau + p selected st coords
        std::vector<double> flat_proj;
        flat_proj.reserve(static_cast<size_t>(K * proj_cols));

        // Build projected points (preserve order => Qhull vertex ids map to candidate indices)
        for (int i = 0; i < K; ++i) {
          const auto& st_vec = candidates[i].st;
          flat_proj.push_back(flat_points[i * point_dim + 0]); // tau
          for (int pi = 0; pi < p; ++pi) {
            int st_index = pr[pi];
            double v = (st_index >= 0 && st_index < static_cast<int>(st_vec.size())) ? st_vec[st_index] : 0.0;
            flat_proj.push_back(v);
          }
        }

        // compute spans and eps per column for this projection
        std::vector<double> spans_proj;
        compute_span_per_dim(flat_proj.data(), proj_cols, K, spans_proj);

        std::vector<double> eps_proj(proj_cols);
        for (int d = 0; d < proj_cols; ++d) eps_proj[d] = std::max(1e-12, 1e-12 * (1.0 + spans_proj[d]));

        // detect varying dimensions (indices into projection columns)
        std::vector<int> varying_dims;
        for (int d = 0; d < proj_cols; ++d) if (spans_proj[d] > eps_proj[d]) varying_dims.push_back(d);

        if (static_cast<int>(varying_dims.size()) <= 1) {
          // effectively 1D -> Qhull not meaningful here; include all
          for (int i = 0; i < K; ++i) hull_indices.insert(i);
          continue;
        }

        try {
          if (static_cast<int>(varying_dims.size()) == 2) {
            // run 2D Qhull on the two varying dims
            std::vector<double> flat2;
            flat2.reserve(static_cast<size_t>(K * 2));
            int a = varying_dims[0], b = varying_dims[1];
            for (int i = 0; i < K; ++i) {
              flat2.push_back(flat_proj[i * proj_cols + a]);
              flat2.push_back(flat_proj[i * proj_cols + b]);
            }
            orgQhull::Qhull qhull;
            qhull.runQhull("", 2, K, flat2.data(), "");
            for (const auto& vertex : qhull.vertexList()) hull_indices.insert(vertex.point().id());
          } else {
            // >= 3 varying dims: build compact array of only varying columns and run Qhull in eff_dim
            const int eff_dim = static_cast<int>(varying_dims.size());
            std::vector<double> flat_eff;
            flat_eff.reserve(static_cast<size_t>(K * eff_dim));
            for (int i = 0; i < K; ++i) {
              for (int vd : varying_dims) {
                flat_eff.push_back(flat_proj[i * proj_cols + vd]);
              }
            }
            orgQhull::Qhull qhull;
            qhull.runQhull("", eff_dim, K, flat_eff.data(), "");
            for (const auto& vertex : qhull.vertexList()) hull_indices.insert(vertex.point().id());
          }
        } catch (...) {
          // conservative fallback: include all points for this projection
          for (int i = 0; i < K; ++i) hull_indices.insert(i);
        }
      } // end projection loop
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
