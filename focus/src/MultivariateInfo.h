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
    }

    std::vector<Candidate> prune(const std::vector<Candidate>& candidates) const override {
        int K = candidates.size();
        if (K <= 1 || pruning_in_ > 0) {
            pruning_in_--;
            return candidates;
        }

        // Ensure sn_ is at least 1D
        int target_dim = sn_.size();
        if (target_dim == 0) {
            throw std::runtime_error("sn_ must have at least one dimension.");
        }

        // Build matrix: [tau, st_0, st_1, ..., st_{d-1}]
        std::vector<std::vector<double>> points(K, std::vector<double>(1 + target_dim));
        for (int i = 0; i < K; ++i) {
            const auto& c = candidates[i];
            points[i][0] = static_cast<double>(c.tau);

            const std::vector<double>& st_vec = c.st;
            if (st_vec.size() != target_dim) {
                throw std::runtime_error("Candidate st dimension mismatch.");
            }

            for (size_t j = 0; j < st_vec.size(); ++j) {
                points[i][1 + j] = st_vec[j];
            }
        }

        std::set<int> hull_indices;

        if (dim_indexes_.empty()) {
            // Full-dimensional hull
            try {
                orgQhull::Qhull qhull;
                std::vector<double> flat_points;
                for (const auto& row : points) {
                    flat_points.insert(flat_points.end(), row.begin(), row.end());
                }
                qhull.runQhull("", 1 + target_dim, K, flat_points.data(), "");

                for (const auto& vertex : qhull.vertexList()) {
                    hull_indices.insert(vertex.point().id());
                }
            } catch (...) {
                // If hull fails, include all indices
                for (int i = 0; i < K; ++i) {
                    hull_indices.insert(i);
                }
            }
        } else {
            // Project to 2D subspaces (tau + pair of dims)
            for (const auto& pair : dim_indexes_) {
                std::vector<std::vector<double>> proj(K, std::vector<double>(3));  // tau + 2 dims
                for (int i = 0; i < K; ++i) {
                    proj[i][0] = points[i][0];  // tau
                    proj[i][1] = points[i][1 + pair.first];
                    proj[i][2] = points[i][1 + pair.second];
                }

                try {
                    orgQhull::Qhull qhull;
                    std::vector<double> flat_proj;
                    for (const auto& row : proj) {
                        flat_proj.insert(flat_proj.end(), row.begin(), row.end());
                    }
                    qhull.runQhull("", 3, K, flat_proj.data(), "");

                    for (const auto& vertex : qhull.vertexList()) {
                        hull_indices.insert(vertex.point().id());
                    }
                } catch (...) {
                    // If hull fails, include all indices
                    for (int i = 0; i < K; ++i) {
                        hull_indices.insert(i);
                    }
                }
            }
        }

        // Build pruned list
        std::vector<Candidate> pruned;
        for (int idx : hull_indices) {
            pruned.push_back(candidates[idx]);
        }

        std::sort(pruned.begin(), pruned.end(),
                 [](const Candidate& a, const Candidate& b) {
                     return a.tau < b.tau;
                 });

        // Update pruning counter
        int pruned_size = pruned.size();
        pruning_in_ = pruned_size * pruning_params_[0] + pruning_params_[1];

        return pruned;
    }
};

} // namespace changepoint