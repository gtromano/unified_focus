#pragma once

#include "Info.h"
#include "Candidate.h"
#include <vector>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <limits>

namespace changepoint {

// Type alias for cost functions
using CostFunction = std::function<std::vector<double>(const std::vector<Candidate>&, const Info&)>;

// Gaussian cost function
inline std::vector<double> compute_costs_gaussian(const std::vector<Candidate>& candidates, const Info& cs) {

    // print number of candidates and n for debugging
    // std::cout << "Number of candidates: " << candidates.size() << ", n: " << cs.n() << std::endl;

    // print all candidates
        // for (size_t idx = 0; idx < candidates.size(); ++idx) {
        //     const auto& c = candidates[idx];
        //     std::cout << "Candidate " << idx << ": tau=" << c.tau << ", st=[";
        //     for (size_t j = 0; j < c.st.size(); ++j) {
        //         std::cout << c.st[j];
        //         if (j < c.st.size() - 1) std::cout << ", ";
        //     }
        //     std::cout << "], theta0=[";
        //     for (size_t j = 0; j < c.theta0.size(); ++j) {
        //         std::cout << c.theta0[j];
        //         if (j < c.theta0.size() - 1) std::cout << ", ";
        //     }
        //     std::cout << "], side=" << c.side << std::endl;
        // }

    int K = candidates.size();
    std::vector<double> costs(K, -1e300);

    const std::vector<double>& S_n = cs.sn();
    int n = cs.n();


    double term3 = 0.0;
    for (double val : S_n) {
        term3 += (val * val) / static_cast<double>(n);
    }

    for (int i = 0; i < K; ++i) {

        // print candidate index and candidates for debugging
        // std::cout << "Processing candidate index: " << i << std::endl;
        // std::cout << "Candidate tau: " << candidates[i].tau << std::endl;


        const auto& c = candidates[i];
        int tau = c.tau;
        const std::vector<double>& S_i = c.st;
        int right_len = n - tau;

        if (tau <= 0 || right_len <= 0 || n <= 0) {
            costs[i] = 0.0;
            continue;
        }
        // print checkpoint
        // std::cout << "Checkpoint 1" << std::endl;
        double cost;
        if (!c.has_theta0()) {
            // std::cout << "Checkpoint 2" << std::endl;
            double term1 = 0.0;
            for (double val : S_i) {
                // print S_i values for debugging
                // std::cout << "S_i value: " << val << std::endl;
                term1 += val * val;
            }
            // print checkpoint
            // std::cout << "Checkpoint 3" << std::endl;
            term1 /= static_cast<double>(tau);

            double term2 = 0.0;
            for (size_t j = 0; j < S_n.size(); ++j) {
                // print S_n and S_i values for debugging
                // std::cout << "S_n[" << j << "]: " << S_n[j] << ", S_i[" << j << "]: " << S_i[j] << std::endl;
                double diff = S_n[j] - S_i[j];
                term2 += diff * diff;
            }
            term2 /= static_cast<double>(right_len);

            cost = term1 + term2 - term3;
            // print cost for debugging
            // std::cout << "Cost (no theta0): " << cost << std::endl;
        } else {
            const std::vector<double>& theta0_vec = c.theta0;
            double shifted_sum = 0.0;
            for (size_t j = 0; j < S_n.size(); ++j) {
                double shifted = S_n[j] - S_i[j] - right_len * theta0_vec[j];
                shifted_sum += shifted * shifted;
            }
            cost = shifted_sum / static_cast<double>(right_len);
        }

        costs[i] = std::isnan(cost) ? 0.0 : cost;
    }

    return costs;
}

// Poisson cost function
inline std::vector<double> compute_costs_poisson(const std::vector<Candidate>& candidates, const Info& cs) {
    int K = candidates.size();
    std::vector<double> costs(K, -1e300);

    const std::vector<double>& S_n = cs.sn();
    int n = cs.n();

    auto max_l = [](const std::vector<double>& st, int tau) -> double {
        double result = 0.0;
        for (double val : st) {
            if (val > 0 && tau > 0) {
                result += -val + val * std::log(val / tau);
            }
        }
        return result;
    };

    double term3 = max_l(S_n, n);

    for (int i = 0; i < K; ++i) {
        const auto& c = candidates[i];
        int tau = c.tau;
        const std::vector<double>& S_i = c.st;
        int right_len = n - tau;

        if (right_len <= 0) {
            costs[i] = 0.0;
            continue;
        }

        double cost;
        if (c.has_theta0()) {
            double term1 = max_l(S_i, tau);
            std::vector<double> diff(S_n.size());
            for (size_t j = 0; j < S_n.size(); ++j) {
                diff[j] = S_n[j] - S_i[j];
            }
            double term2 = max_l(diff, right_len);
            cost = term1 + term2 - term3;
        } else {
            std::vector<double> diff(S_n.size());
            for (size_t j = 0; j < S_n.size(); ++j) {
                diff[j] = S_n[j] - S_i[j];
            }
            double term2 = max_l(diff, right_len);
            const std::vector<double>& theta0_vec = c.theta0;
            double null_val = 0.0;
            for (size_t j = 0; j < S_n.size(); ++j) {
                if (theta0_vec[j] > 0) {
                    null_val += -right_len * theta0_vec[j] + diff[j] * std::log(theta0_vec[j]);
                }
            }
            cost = term2 - null_val;
        }

        costs[i] = std::isnan(cost) ? 0.0 : cost;
    }

    return costs;
}

// Two-sided cost wrapper
inline CostFunction make_two_sided_cost_fn(const CostFunction& base_cost_fn) {
    return [base_cost_fn](const std::vector<Candidate>& candidates, const Info& cs) -> std::vector<double> {
        int K = candidates.size();
        std::vector<double> costs(K, -1e300);

        // Downcast to UnivariateInfo to access left/right
        const auto* uni_cs = dynamic_cast<const UnivariateInfo*>(&cs);
        if (!uni_cs) {
            throw std::runtime_error("make_two_sided_cost_fn requires UnivariateInfo");
        }

        for (int i = 0; i < K; ++i) {
            const auto& c = candidates[i];
            std::vector<Candidate> single = {c};

            if (c.side == "left") {
                std::vector<double> single_cost = base_cost_fn(single, *uni_cs->left());
                costs[i] = !single_cost.empty() ? single_cost[0] : -1e300;
            } else {
                std::vector<double> single_cost = base_cost_fn(single, *uni_cs->right());
                costs[i] = !single_cost.empty() ? single_cost[0] : -1e300;
            }
        }

        return costs;
    };
}

} // namespace changepoint
