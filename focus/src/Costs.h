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

namespace changepoint {

// Type alias for cost functions - now returns ChangepointResult
using CostFunction = std::function<ChangepointResult(const Info&, const std::vector<double>&)>;

// Helper function to compute costs for all candidates
inline std::vector<double> compute_costs_gaussian_helper(
    const std::vector<Candidate>& candidates, 
    const std::vector<double>& S_n,
    int n,
    const std::vector<double>& theta0) {
    
    int K = candidates.size();
    std::vector<double> costs(K, -1e300);

    bool use_theta0 = !theta0.empty() && !std::isnan(theta0[0]);

    double term3 = 0.0;
    if (!use_theta0) {
        for (double val : S_n) {
            term3 += (val * val) / static_cast<double>(n);
        }
    }

    for (int i = 0; i < K; ++i) {
        const auto& c = candidates[i];
        int tau = c.tau;
        const std::vector<double>& S_i = c.st;
        int right_len = n - tau;

        if (tau <= 0 || right_len <= 0 || n <= 0) {
            costs[i] = 0.0;
            continue;
        }

        double cost;
        if (!use_theta0) {
            double term1 = 0.0;
            for (double val : S_i) {
                term1 += val * val;
            }
            term1 /= static_cast<double>(tau);

            double term2 = 0.0;
            for (size_t j = 0; j < S_n.size(); ++j) {
                double diff = S_n[j] - S_i[j];
                term2 += diff * diff;
            }
            term2 /= static_cast<double>(right_len);

            cost = term1 + term2 - term3;
        } else {
            double shifted_sum = 0.0;
            for (size_t j = 0; j < S_n.size(); ++j) {
                double diff = S_n[j] - S_i[j] - right_len * theta0[j];
                shifted_sum += diff * diff;
            }
            cost = shifted_sum / static_cast<double>(right_len);
        }

        costs[i] = std::isnan(cost) ? 0.0 : cost;
    }

    return costs;
}

// Gaussian cost function - returns ChangepointResult
inline ChangepointResult compute_costs_gaussian(const Info& cs, const std::vector<double>& theta0) {
    ChangepointResult result;
    result.stopping_time = cs.n();
    
    const auto& candidates = cs.candidates();
    
    if (candidates.size() <= 1) {
        result.changepoint = std::nullopt;
        result.stat = std::nullopt;
        return result;
    }
    
    // Exclude last candidate(s) (dummy for current time)
    auto last_cands = cs.new_candidate();
    int exclude_count = last_cands.size();
    
    if (candidates.size() <= static_cast<size_t>(exclude_count)) {
        result.changepoint = std::nullopt;
        result.stat = std::nullopt;
        return result;
    }
    
    std::vector<Candidate> considered(candidates.begin(), 
                                     candidates.end() - exclude_count);
    
    if (considered.empty()) {
        result.changepoint = std::nullopt;
        result.stat = std::nullopt;
        return result;
    }
    
    // Compute costs
    std::vector<double> vals = compute_costs_gaussian_helper(
        considered, cs.sn(), cs.n(), theta0);
    
    // Find maximum
    auto max_it = std::max_element(vals.begin(), vals.end());
    int max_idx = std::distance(vals.begin(), max_it);
    
    result.changepoint = considered[max_idx].tau;
    result.stat = *max_it;
    
    return result;
}

// Helper function to compute Poisson costs for all candidates
inline std::vector<double> compute_costs_poisson_helper(
    const std::vector<Candidate>& candidates,
    const std::vector<double>& S_n,
    int n,
    const std::vector<double>& theta0) {
    
    int K = static_cast<int>(candidates.size());
    std::vector<double> costs(K, -1e300);

    bool use_theta0 = !theta0.empty() && !std::isnan(theta0[0]);

    auto max_l = [](const std::vector<double>& st, int tau) -> double {
        if (tau <= 0) return 0.0;
        double result = 0.0;
        bool any_nan = false;
        for (double val : st) {
            double ratio = val / static_cast<double>(tau);
            double log_ratio = std::log(ratio);
            double term2 = val * log_ratio;
            double contrib = -val + term2;

            if (std::isnan(contrib)) {
                any_nan = true;
            } else {
                result += contrib;
            }
        }
        return any_nan ? std::numeric_limits<double>::quiet_NaN() : result;
    };

    double term3 = 0.0;
    if (!use_theta0) {
        term3 = max_l(S_n, n);
    }

    for (int i = 0; i < K; ++i) {
        const auto& c = candidates[i];
        int tau = c.tau;
        const std::vector<double>& S_i = c.st;
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
                double log_th = std::log(th);
                double contrib = - static_cast<double>(right_len) * th + diff[j] * log_th;
                if (std::isnan(contrib)) {
                    any_nan = true;
                } else {
                    null_val += contrib;
                }
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

// Poisson cost function - returns ChangepointResult
inline ChangepointResult compute_costs_poisson(const Info& cs, const std::vector<double>& theta0) {
    ChangepointResult result;
    result.stopping_time = cs.n();
    
    const auto& candidates = cs.candidates();
    
    if (candidates.size() <= 1) {
        result.changepoint = std::nullopt;
        result.stat = std::nullopt;
        return result;
    }
    
    // Exclude last candidate(s) (dummy for current time)
    auto last_cands = cs.new_candidate();
    int exclude_count = last_cands.size();
    
    if (candidates.size() <= static_cast<size_t>(exclude_count)) {
        result.changepoint = std::nullopt;
        result.stat = std::nullopt;
        return result;
    }
    
    std::vector<Candidate> considered(candidates.begin(), 
                                     candidates.end() - exclude_count);
    
    if (considered.empty()) {
        result.changepoint = std::nullopt;
        result.stat = std::nullopt;
        return result;
    }
    
    // Compute costs
    std::vector<double> vals = compute_costs_poisson_helper(
        considered, cs.sn(), cs.n(), theta0);
    
    // Find maximum
    auto max_it = std::max_element(vals.begin(), vals.end());
    int max_idx = std::distance(vals.begin(), max_it);
    
    result.changepoint = considered[max_idx].tau;
    result.stat = *max_it;
    
    return result;
}

// Two-sided cost wrapper
inline CostFunction make_two_sided_cost_fn(
    const std::function<std::vector<double>(const std::vector<Candidate>&, 
                                            const std::vector<double>&,
                                            int,
                                            const std::vector<double>&)>& base_cost_helper) {
    
    return [base_cost_helper](const Info& cs, const std::vector<double>& theta0) -> ChangepointResult {
        ChangepointResult result;
        result.stopping_time = cs.n();
        
        const auto& candidates = cs.candidates();
        
        if (candidates.size() <= 1) {
            result.changepoint = std::nullopt;
            result.stat = std::nullopt;
            return result;
        }
        
        // Downcast to UnivariateInfo to access left/right
        const auto* uni_cs = dynamic_cast<const UnivariateInfo*>(&cs);
        if (!uni_cs) {
            throw std::runtime_error("make_two_sided_cost_fn requires UnivariateInfo");
        }
        
        // Exclude last candidate(s) (dummy for current time)
        auto last_cands = cs.new_candidate();
        int exclude_count = last_cands.size();
        
        if (candidates.size() <= static_cast<size_t>(exclude_count)) {
            result.changepoint = std::nullopt;
            result.stat = std::nullopt;
            return result;
        }
        
        std::vector<Candidate> considered(candidates.begin(), 
                                         candidates.end() - exclude_count);
        
        if (considered.empty()) {
            result.changepoint = std::nullopt;
            result.stat = std::nullopt;
            return result;
        }
        
        // Compute costs for each side separately
        int K = considered.size();
        std::vector<double> costs(K, -1e300);
        
        for (int i = 0; i < K; ++i) {
            const auto& c = considered[i];
            std::vector<Candidate> single = {c};
            
            std::vector<double> single_cost;
            if (c.side == "left") {
                single_cost = base_cost_helper(single, uni_cs->left()->sn(), 
                                              uni_cs->left()->n(), theta0);
            } else {
                single_cost = base_cost_helper(single, uni_cs->right()->sn(), 
                                              uni_cs->right()->n(), theta0);
            }
            costs[i] = single_cost[0];
        }
        
        // Find maximum
        auto max_it = std::max_element(costs.begin(), costs.end());
        int max_idx = std::distance(costs.begin(), max_it);
        
        result.changepoint = considered[max_idx].tau;
        result.stat = *max_it;
        
        return result;
    };
}

// Convenience function to create two-sided Gaussian cost function
inline CostFunction make_two_sided_gaussian() {
    return make_two_sided_cost_fn(compute_costs_gaussian_helper);
}

// Convenience function to create two-sided Poisson cost function
inline CostFunction make_two_sided_poisson() {
    return make_two_sided_cost_fn(compute_costs_poisson_helper);
}

} // namespace changepoint