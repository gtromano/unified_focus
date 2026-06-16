#pragma once

#include "Info.h"
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullVertexSet.h"
#include <set>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <limits>

namespace changepoint {

// Center K x m matrix in-place (row-major)
static void center_points(std::vector<double>& data, int K, int m) {
  std::vector<double> mean(m, 0.0);
  for (int r = 0; r < K; ++r)
    for (int c = 0; c < m; ++c)
      mean[c] += data[r * m + c];

  for (int c = 0; c < m; ++c)
    mean[c] /= static_cast<double>(K);

  for (int r = 0; r < K; ++r)
    for (int c = 0; c < m; ++c)
      data[r * m + c] -= mean[c];
}

// Modified Gram-Schmidt rank (on centered data)
static int compute_affine_rank_mgs(const double* data, int K, int m, double tol = 1e-12) {
  std::vector<std::vector<double>> Q;

  for (int c = 0; c < m; ++c) {
    std::vector<double> v(K);
    for (int r = 0; r < K; ++r)
      v[r] = data[r * m + c];

    for (const auto& q : Q) {
      double dot = 0.0;
      for (int r = 0; r < K; ++r)
        dot += q[r] * v[r];
      for (int r = 0; r < K; ++r)
        v[r] -= dot * q[r];
    }

    double norm = 0.0;
    for (double x : v) norm += x * x;
    norm = std::sqrt(norm);

    if (norm > tol) {
      for (double& x : v) x /= norm;
      Q.push_back(std::move(v));
    }
  }
  return static_cast<int>(Q.size());
}

// Remove near-duplicate rows
static std::vector<double> deduplicate_points(const std::vector<double>& data,
                                              int K, int m, int& K_out,
                                              double tol = 1e-12) {
  std::vector<double> out;
  out.reserve(data.size());

  for (int i = 0; i < K; ++i) {
    const double* pi = &data[i * m];
    bool duplicate = false;

    for (int j = 0; j < K_out; ++j) {
      const double* pj = &out[j * m];
      double dist2 = 0.0;
      for (int d = 0; d < m; ++d) {
        double diff = pi[d] - pj[d];
        dist2 += diff * diff;
      }
      if (dist2 < tol * tol) {
        duplicate = true;
        break;
      }
    }

    if (!duplicate) {
      out.insert(out.end(), pi, pi + m);
      K_out++;
    }
  }
  return out;
}

static std::vector<std::vector<int>> generate_circular_combinations(int d, int p) {
    std::vector<std::vector<int>> out;
    if (p <= 0 || p > d) return out;

    for (int start = 0; start < d; ++start) {
        std::vector<int> comb(p);
        for (int i = 0; i < p; ++i) {
            comb[i] = (start + i) % d;
        }
        out.push_back(comb);
    }

    return out;
}
class MultivariateInfo : public CandidateListInfo {
private:
  std::vector<std::vector<int>> dim_indexes_;
  int pruning_params_[2];
  mutable int pruning_in_;
  double anomaly_intensity_;

public:
  MultivariateInfo(const std::vector<double>& sn = {0.0},
                   double n = 0.0,
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
    if (K < target_dim + 2) return candidates;

    std::set<int> hull_indices;

    const int point_dim = 1 + target_dim;
    std::vector<double> flat_points;
    flat_points.reserve(K * point_dim);

    for (int i = 0; i < K; ++i) {
      flat_points.push_back(static_cast<double>(candidates[i].tau));
      for (int j = 0; j < target_dim; ++j) {
        double v = (j < (int)candidates[i].st.size()) ? candidates[i].st[j] : 0.0;
        flat_points.push_back(v);
      }
    }

    auto run_qhull_safe = [&](std::vector<double> data, int dims) {
      int K_eff = 0;
      data = deduplicate_points(data, K, dims, K_eff);

      if (K_eff < 3) return;

      center_points(data, K_eff, dims);

      int rank = compute_affine_rank_mgs(data.data(), K_eff, dims);
      if (rank < 2 || K_eff < rank + 1) return;

      try {
        // first try to run qhull normally
        orgQhull::Qhull qh;
        qh.runQhull("", dims, K_eff, data.data(), "");
        for (const auto& v : qh.vertexList())
          hull_indices.insert(v.point().id());
      } catch (...) {
        try {
          // if this fails try to run qhull with "QJ" option (joggle input)
          orgQhull::Qhull qh2;
          qh2.runQhull("", dims, K_eff, data.data(), "QJ");
          for (const auto& v : qh2.vertexList())
            hull_indices.insert(v.point().id());
        } catch (...) {
          // if this fails again, just return all indices (no pruning)
          for (int i = 0; i < K; ++i)
            hull_indices.insert(i);
        }
      }
    };

    if (dim_indexes_.empty()) {
      run_qhull_safe(flat_points, point_dim);
    } else {
      for (const auto& pr : dim_indexes_) {
        int p = pr.size();
        if (p <= 0) continue;

        int dims = 1 + p;
        std::vector<double> proj;
        proj.reserve(K * dims);

        for (int i = 0; i < K; ++i) {
          proj.push_back(flat_points[i * point_dim]);
          for (int j = 0; j < p; ++j) {
            int idx = pr[j];
            double v = (idx < (int)candidates[i].st.size()) ? candidates[i].st[idx] : 0.0;
            proj.push_back(v);
          }
        }

        run_qhull_safe(proj, dims);
      }
    }

    if (hull_indices.empty()) {
      for (int i = 0; i < K; ++i)
        hull_indices.insert(i);
    }

    // Apply anomaly_intensity pruning if enabled
    if (!std::isnan(anomaly_intensity_) && anomaly_intensity_ > 0.0) {
      std::set<int> intensity_filtered;

      for (int idx : hull_indices) {
        const auto& c = candidates[idx];

        const double denom = static_cast<double>(n_) - c.tau;

        // Skip unstable ratios when tau is numerically equal to n_
        if (denom <= 1e-12) {
          intensity_filtered.insert(idx);
          continue;
        }

        double max_abs_ratio = 0.0;

        for (int dim = 0; dim < target_dim; ++dim) {
          const double st_val =
            (dim < static_cast<int>(c.st.size())) ? c.st[dim] : 0.0;

          const double num = sn_[dim] - st_val;
          const double abs_ratio = std::abs(num / denom);

          max_abs_ratio = std::max(max_abs_ratio, abs_ratio);
        }

        // Keep candidate if infinity norm exceeds threshold
        if (max_abs_ratio >= anomaly_intensity_) {
          intensity_filtered.insert(idx);
        }
      }

      hull_indices = std::move(intensity_filtered);
    }

    std::vector<Candidate> pruned;
    for (int idx : hull_indices)
      pruned.push_back(candidates[idx]);

    std::sort(pruned.begin(), pruned.end(),
              [](const Candidate& a, const Candidate& b){ return a.tau < b.tau; });

    int pruned_size = pruned.size();
    pruning_in_ = pruned_size * pruning_params_[0] + pruning_params_[1];

    return pruned;
  }
};

} // namespace changepoint
