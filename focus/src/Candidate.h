#pragma once

#include <vector>
#include <string>
#include <limits>
#include <cmath> // For std::isnan

namespace changepoint {

struct Candidate {
    std::vector<double> st;                // Cumulative sum at tau (can be scalar or vector)
    int tau;                               // Time index of the candidate
    std::vector<double> theta0;            // Null hypothesis parameter (empty or NaN if missing)
    std::string side;                      // "right" or "left" for univariate two-sided

    Candidate() : tau(0), side("right") {
        st = {0.0};
    }

    Candidate(const std::vector<double>& st_, int tau_,
              const std::vector<double>& theta0_ = {},
              const std::string& side_ = "right")
        : st(st_), tau(tau_), theta0(theta0_), side(side_) {}

    // Constructor for scalar case
    Candidate(double st_scalar, int tau_,
              double theta0_scalar = std::numeric_limits<double>::quiet_NaN(),
              const std::string& side_ = "right")
        : tau(tau_), side(side_) {
        st = {st_scalar};
        if (!std::isnan(theta0_scalar)) {
            theta0 = {theta0_scalar};
        }
    }

    double scalar_st() const { return st[0]; }

    double scalar_theta0() const {
        return theta0.empty() ? std::numeric_limits<double>::quiet_NaN() : theta0[0];
    }

    bool has_theta0() const {
        // print "Checking if candidate has theta0" for debugging
        std::cout << "Checking if candidate has theta0" << std::endl;
        // print theta0 size for debugging
        std::cout << "theta0 size: " << theta0.size() << std::endl;
        return theta0.size() > 0;// && !std::isnan(theta0[0]);
    }
};

} // namespace changepoint
