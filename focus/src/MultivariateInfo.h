#pragma once

#include "Info.h"
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullVertexSet.h"
#include <set>
#include <utility>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <limits>


namespace changepoint {

// data: K x m (row-major) matrix (data[r*m + c])
static std::vector<int> find_independent_columns_mgs(const double* data, int K, int m, double tol = 1e-12) {
  if (K <= 0 || m <= 0) return {};
  // precompute column 2-norms and scale
  std::vector<double> col_norm(m, 0.0);
  double max_col_norm = 0.0;
  for (int c = 0; c < m; ++c) {
    double s = 0.0;
    const double* cp = data + c;
    for (int r = 0; r < K; ++r) { double v = *(cp + static_cast<size_t>(r) * m); s += v*v; }
    col_norm[c] = std::sqrt(s);
    if (col_norm[c] > max_col_norm) max_col_norm = col_norm[c];
  }
  if (max_col_norm == 0.0) return {}; // all-zero matrix

  double scaled_tol = tol * max_col_norm;

  std::vector<std::vector<double>> Q; // orthonormal basis columns
  std::vector<int> keep;
  Q.reserve(m);

  for (int c = 0; c < m; ++c) {
    // copy column c
    std::vector<double> v(K);
    const double* cp = data + c;
    for (int r = 0; r < K; ++r) v[r] = *(cp + static_cast<size_t>(r) * m);

    // subtract projection on existing Q
    for (const auto& qcol : Q) {
      double dot = 0.0;
      for (int r = 0; r < K; ++r) dot += qcol[r] * v[r];
      for (int r = 0; r < K; ++r) v[r] -= dot * qcol[r];
    }
    // residual norm
    double n2 = 0.0;
    for (int r = 0; r < K; ++r) n2 += v[r]*v[r];
    double n = std::sqrt(n2);
    if (n > scaled_tol) {
      // accept basis vector
      for (int r = 0; r < K; ++r) v[r] /= n;
      Q.push_back(std::move(v));
      keep.push_back(c);
    }
  }
  return keep;
}

inline std::vector<std::vector<int>> generate_circular_combinations(int D, int p) {
    std::vector<std::vector<int>> out;
    if (p <= 0 || p > D) return out;

    for (int start = 0; start < D; ++start) {
        std::vector<int> comb(p);
        for (int i = 0; i < p; ++i) {
            comb[i] = (start + i) % D;  // wrap around
        }
        out.push_back(comb);
    }

    return out;
}



class MultivariateInfo : public CandidateListInfo {
private:
  std::vector<std::vector<int>> dim_indexes_;  // Projections of arbitrary p size
  int pruning_params_[2];  // (multiplier, offset)
  mutable int pruning_in_;          // Counter for pruning frequency
  double anomaly_intensity_;

public:
  MultivariateInfo(const std::vector<double>& sn = {0.0},
                   int n = 0,
                   const std::vector<std::vector<int>>& dim_indexes = {},
                   int pruning_mult = 2,
                   int pruning_offset = 1,
                   double anomaly_intensity = std::numeric_limits<double>::quiet_NaN())
    : CandidateListInfo(sn, n),
      dim_indexes_(dim_indexes),
      pruning_in_(5),
      anomaly_intensity_(anomaly_intensity) {
    pruning_params_[0] = pruning_mult;
    pruning_params_[1] = pruning_offset;

    // Initialize with first candidate
    auto initial = new_candidate();
    candidates_.insert(candidates_.end(), initial.begin(), initial.end());
  }

  void set_pruning_counter_to_zero() {
    pruning_in_ = 0;
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
          double v = data[static_cast<size_t>(r) * dims + d];
          if (v < minv) minv = v;
          if (v > maxv) maxv = v;
        }
        spans[d] = maxv - minv;
      }
    };

    // tolerance used by MGS independence test (tune if needed)
    const double mgs_tol = 1e-12;

    // --- Full-dimensional hull via Qhull, guarded for degeneracy & numeric rank ---
    if (dim_indexes_.empty()) {
      // compute spans for all coords (tau + st_i)
      std::vector<double> spans;
      compute_span_per_dim(flat_points.data(), point_dim, K, spans);

      // set eps per-dim relative to span (tiny absolute floor too)
      std::vector<double> epss(point_dim);
      for (int d = 0; d < point_dim; ++d) {
        epss[d] = std::max(1e-12, 1e-12 * (1.0 + spans[d]));
      }

      // select columns that appear to vary
      std::vector<int> varying_cols;
      for (int d = 0; d < point_dim; ++d) if (spans[d] > epss[d]) varying_cols.push_back(d);

      if (static_cast<int>(varying_cols.size()) <= 1) {
        // effectively 1D: keep all candidates (conservative)
        for (int i = 0; i < K; ++i) hull_indices.insert(i);
      } else {
        // build K x m matrix with only the varying columns, row-major
        const int m = static_cast<int>(varying_cols.size());
        std::vector<double> flat_eff;
        flat_eff.reserve(static_cast<size_t>(K * m));
        for (int r = 0; r < K; ++r) {
          size_t base = static_cast<size_t>(r) * point_dim;
          for (int cidx = 0; cidx < m; ++cidx) {
            int col = varying_cols[cidx];
            flat_eff.push_back(flat_points[base + col]);
          }
        }

        // numeric independence test (select independent columns from flat_eff)
        std::vector<int> kept = find_independent_columns_mgs(flat_eff.data(), K, m, mgs_tol);

        if (kept.empty() || static_cast<int>(kept.size()) <= 1) {
          // numerically rank <= 1 -> degenerate: keep all
          for (int i = 0; i < K; ++i) hull_indices.insert(i);
        } else if (static_cast<int>(kept.size()) == 2) {
          // run Qhull in 2D on the two independent columns
          std::vector<double> flat2;
          flat2.reserve(static_cast<size_t>(K * 2));
          for (int r = 0; r < K; ++r) {
            flat2.push_back(flat_eff[static_cast<size_t>(r) * m + kept[0]]);
            flat2.push_back(flat_eff[static_cast<size_t>(r) * m + kept[1]]);
          }
          try {
            orgQhull::Qhull qhull;
            qhull.runQhull("", 2, K, flat2.data(), "");
            for (const auto& vertex : qhull.vertexList()) hull_indices.insert(vertex.point().id());
          } catch (...) {
            for (int i = 0; i < K; ++i) hull_indices.insert(i);
          }
        } else {
          // rank >= 3: build compact data with only the independent columns and call Qhull
          const int eff_rank = static_cast<int>(kept.size());
          std::vector<double> flat_rank;
          flat_rank.reserve(static_cast<size_t>(K * eff_rank));
          for (int r = 0; r < K; ++r) {
            for (int j = 0; j < eff_rank; ++j) {
              int col = kept[j];
              flat_rank.push_back(flat_eff[static_cast<size_t>(r) * m + col]);
            }
          }
          try {
            orgQhull::Qhull qhull;
            qhull.runQhull("", eff_rank, K, flat_rank.data(), "");
            for (const auto& vertex : qhull.vertexList()) hull_indices.insert(vertex.point().id());
          } catch (...) {
            for (int i = 0; i < K; ++i) hull_indices.insert(i);
          }
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
          flat_proj.push_back(flat_points[static_cast<size_t>(i) * point_dim + 0]); // tau
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
          // effectively 1D -> include all
          for (int i = 0; i < K; ++i) hull_indices.insert(i);
          continue;
        }

        // Build K x m_effective array with only the varying columns (row-major)
        const int m_eff = static_cast<int>(varying_dims.size());
        std::vector<double> flat_eff;
        flat_eff.reserve(static_cast<size_t>(K * m_eff));
        for (int r = 0; r < K; ++r) {
          for (int vd = 0; vd < m_eff; ++vd) {
            int col = varying_dims[vd];
            flat_eff.push_back(flat_proj[static_cast<size_t>(r) * proj_cols + col]);
          }
        }

        // numeric independence selection
        std::vector<int> kept = find_independent_columns_mgs(flat_eff.data(), K, m_eff, mgs_tol);

        if (kept.empty() || static_cast<int>(kept.size()) <= 1) {
          // degenerate -> include all
          for (int i = 0; i < K; ++i) hull_indices.insert(i);
          continue;
        }

        if (static_cast<int>(kept.size()) == 2) {
          // run Qhull in 2D on the two independent columns
          std::vector<double> flat2;
          flat2.reserve(static_cast<size_t>(K * 2));
          for (int r = 0; r < K; ++r) {
            flat2.push_back(flat_eff[static_cast<size_t>(r) * m_eff + kept[0]]);
            flat2.push_back(flat_eff[static_cast<size_t>(r) * m_eff + kept[1]]);
          }
          try {
            orgQhull::Qhull qhull;
            qhull.runQhull("", 2, K, flat2.data(), "");
            for (const auto& vertex : qhull.vertexList()) hull_indices.insert(vertex.point().id());
          } catch (...) {
            for (int i = 0; i < K; ++i) hull_indices.insert(i);
          }
        } else {
          // higher rank -> Qhull with reduced basis
          const int eff_rank = static_cast<int>(kept.size());
          std::vector<double> flat_rank;
          flat_rank.reserve(static_cast<size_t>(K * eff_rank));
          for (int r = 0; r < K; ++r) {
            for (int j = 0; j < eff_rank; ++j) {
              int col = kept[j];
              flat_rank.push_back(flat_eff[static_cast<size_t>(r) * m_eff + col]);
            }
          }
          try {
            orgQhull::Qhull qhull;
            qhull.runQhull("", eff_rank, K, flat_rank.data(), "");
            for (const auto& vertex : qhull.vertexList()) hull_indices.insert(vertex.point().id());
          } catch (...) {
            for (int i = 0; i < K; ++i) hull_indices.insert(i);
          }
        }
      } // end projection loop
    }

    // Apply anomaly_intensity pruning if enabled
    if (!std::isnan(anomaly_intensity_) && anomaly_intensity_ > 0.0) {
      std::set<int> intensity_filtered;

      for (int idx : hull_indices) {
        const auto& c = candidates[idx];
        int denom = n_ - c.tau;

        // Calculate infinity norm: max absolute ratio across all dimensions
        double max_abs_ratio = 0.0;
        bool has_valid_ratio = false;

        for (int dim = 0; dim < target_dim; ++dim) {
          if (denom > 0) {
            double st_val = (dim < static_cast<int>(c.st.size())) ? c.st[dim] : 0.0;
            double num = sn_[dim] - st_val;
            double ratio = num / denom;
            double abs_ratio = std::abs(ratio);

            if (abs_ratio > max_abs_ratio) {
              max_abs_ratio = abs_ratio;
            }
            has_valid_ratio = true;
          }
        }

        // Keep candidate if infinity norm >= anomaly_intensity
        // (i.e., at least one dimension has strong signal)
        if (!has_valid_ratio || max_abs_ratio >= anomaly_intensity_) {
          intensity_filtered.insert(idx);
        }
      }

      hull_indices = intensity_filtered;
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
