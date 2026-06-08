#pragma once

#include <vector>
#include <string>

namespace changepoint {

struct Candidate {
    std::vector<double> st;                // Cumulative sum at tau (can be scalar or vector)
    double tau;                            // Time index of the candidate (effective observations accounting for lambda)
    std::string side;                      // "right" or "left" for univariate two-sided

    Candidate() : tau(0.0), side("right") {
        st = {0.0};
    }

    Candidate(const std::vector<double>& st_, double tau_,
              const std::string& side_ = "right")
        : st(st_), tau(tau_), side(side_) {}

    // Constructor for scalar case
    Candidate(double st_scalar, double tau_,
              const std::string& side_ = "right")
        : tau(tau_), side(side_) {
        st = {st_scalar};
    }

    double scalar_st() const { return st[0]; }
};

} // namespace changepoint