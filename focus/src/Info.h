#pragma once

#include "Candidate.h"
#include <vector>
#include <memory>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <cmath>

namespace changepoint {

// Abstract base class for Info (cumulative state with behavior)
class Info {
protected:
    std::vector<double> sn_;                    // Cumulative sum
    int n_;                                     // Number of observations
    std::vector<double> theta0_;                // Null hypothesis parameter (empty or NaN if missing)
    std::vector<Candidate> candidates_;         // Candidate changepoints

public:
    Info(const std::vector<double>& sn = {0.0},
         int n = 0,
         const std::vector<double>& theta0 = {})
        : sn_(sn), n_(n), theta0_(theta0) {}

    virtual ~Info() = default;

    // Create new candidate(s) representing current state
    virtual std::vector<Candidate> new_candidate() const {
        return {Candidate(sn_, n_, "right")};
    }

    // Add new observation point (renamed from old update)
    virtual void add_new_point(const std::vector<double>& y) {
        n_++;
        if (sn_.size() != y.size()) {
            sn_.resize(y.size(), 0.0);
        }
        for (size_t i = 0; i < sn_.size(); ++i) {
            sn_[i] += y[i];
        }
    }

    // New update method that replaces Detector.update()
    virtual void update(const std::vector<double>& y) {
        // Update cumulative state
        add_new_point(y);
        
        // Prune candidates
        candidates_ = prune(candidates_);
        
        // Append new candidate(s) for current time
        auto new_cands = new_candidate();
        candidates_.insert(candidates_.end(), new_cands.begin(), new_cands.end());
    }

    virtual std::vector<Candidate> prune(const std::vector<Candidate>& candidates) const {
        return candidates;
    }

    const std::vector<double>& sn() const { return sn_; }
    int n() const { return n_; }
    const std::vector<double>& theta0() const { return theta0_; }
    bool has_theta0() const { return !theta0_.empty() && !std::isnan(theta0_[0]); }
    const std::vector<Candidate>& candidates() const { return candidates_; }
};

// One-sided univariate Info with monotone pruning
class OneSideUnivariateInfo : public Info {
private:
    std::string side_;

public:
    OneSideUnivariateInfo(double theta0 = std::numeric_limits<double>::quiet_NaN(),
                          double sn = 0.0, int n = 0,
                          const std::string& side = "right")
        : Info({sn}, n, std::isnan(theta0) ? std::vector<double>() : std::vector<double>({theta0})),
          side_(side) {
        if (side != "right" && side != "left") {
            throw std::invalid_argument("side must be 'right' or 'left'");
        }
        
        // Initialize with first candidate
        auto initial = new_candidate();
        candidates_.insert(candidates_.end(), initial.begin(), initial.end());
    }

    std::vector<Candidate> new_candidate() const override {
        return {Candidate(sn_, n_, side_)};
    }

    std::vector<Candidate> prune(const std::vector<Candidate>& candidates) const override {
        // Filter candidates for this side
        std::vector<Candidate> side_candidates;
        for (const auto& c : candidates) {
            if (c.side == side_) {
                side_candidates.push_back(c);
            }
        }

        int K = side_candidates.size();
        if (K <= 1) {
            return candidates;
        }

        // Monotone pruning
        int i = K;
        while (i > 1) {
            const auto& c1 = side_candidates[i - 1];
            const auto& c0 = side_candidates[i - 2];

            int tau1 = c1.tau;
            int tau0 = c0.tau;
            int denom1 = n_ - tau1;
            int denom0 = n_ - tau0;

            double num1 = sn_[0] - c1.scalar_st();
            double num0 = sn_[0] - c0.scalar_st();

            double ratio1 = (denom1 > 0) ? (num1 / denom1) : std::numeric_limits<double>::infinity();
            double ratio0 = (denom0 > 0) ? (num0 / denom0) : std::numeric_limits<double>::infinity();

            bool cond = (c1.side == "right") ? (ratio1 <= ratio0) : (ratio1 >= ratio0);

            if (cond) {
                i--;
                if (i == 1) break;
            } else {
                break;
            }
        }

        // Return pruned list
        std::vector<Candidate> pruned(side_candidates.begin(), side_candidates.begin() + i);

        if (side_candidates.size() == candidates.size()) {
            return pruned;
        } else {
            // Rebuild with non-side candidates
            std::vector<Candidate> result;
            for (const auto& c : candidates) {
                if (c.side != side_) {
                    result.push_back(c);
                }
            }
            result.insert(result.end(), pruned.begin(), pruned.end());
            std::sort(result.begin(), result.end(),
                     [](const Candidate& a, const Candidate& b) {
                         return a.tau < b.tau || (a.tau == b.tau && a.side < b.side);
                     });
            return result;
        }
    }

    const std::string& side() const { return side_; }
};

// Two-sided univariate Info
class UnivariateInfo : public Info {
private:
    std::unique_ptr<OneSideUnivariateInfo> right_;
    std::unique_ptr<OneSideUnivariateInfo> left_;

public:
    UnivariateInfo(double theta0 = std::numeric_limits<double>::quiet_NaN(),
                   double sn = 0.0, int n = 0)
        : Info({sn}, n, std::isnan(theta0) ? std::vector<double>() : std::vector<double>({theta0})) {

        right_ = std::make_unique<OneSideUnivariateInfo>(theta0, sn, n, "right");
        left_ = std::make_unique<OneSideUnivariateInfo>(theta0, sn, n, "left");
        
        // Initialize candidates from both sides
        auto right_cand = right_->new_candidate();
        auto left_cand = left_->new_candidate();
        candidates_.insert(candidates_.end(), right_cand.begin(), right_cand.end());
        candidates_.insert(candidates_.end(), left_cand.begin(), left_cand.end());
    }

    std::vector<Candidate> new_candidate() const override {
        auto right_cand = right_->new_candidate();
        auto left_cand = left_->new_candidate();
        std::vector<Candidate> result;
        result.insert(result.end(), right_cand.begin(), right_cand.end());
        result.insert(result.end(), left_cand.begin(), left_cand.end());
        return result;
    }

    void add_new_point(const std::vector<double>& y) override {
        // Update right with y
        right_->add_new_point(y);
        // Update left with y
        left_->add_new_point(y);

        // Keep top-level consistent
        n_ = right_->n();
        sn_ = right_->sn();
    }

    void update(const std::vector<double>& y) override {
        // Update both sides
        right_->update(y);
        left_->update(y);

        // Keep top-level consistent
        n_ = right_->n();
        sn_ = right_->sn();
        
        // Combine candidates from both sides
        candidates_.clear();
        const auto& right_cands = right_->candidates();
        const auto& left_cands = left_->candidates();
        candidates_.insert(candidates_.end(), right_cands.begin(), right_cands.end());
        candidates_.insert(candidates_.end(), left_cands.begin(), left_cands.end());
        
        std::sort(candidates_.begin(), candidates_.end(),
                 [](const Candidate& a, const Candidate& b) {
                     return a.tau < b.tau || (a.tau == b.tau && a.side < b.side);
                 });
    }

    std::vector<Candidate> prune(const std::vector<Candidate>& candidates) const override {
        // Split by side
        std::vector<Candidate> right_cands, left_cands;
        for (const auto& c : candidates) {
            if (c.side == "right") {
                right_cands.push_back(c);
            } else {
                left_cands.push_back(c);
            }
        }

        auto pr_right = right_->prune(right_cands);
        auto pr_left = left_->prune(left_cands);

        std::vector<Candidate> combined;
        combined.insert(combined.end(), pr_right.begin(), pr_right.end());
        combined.insert(combined.end(), pr_left.begin(), pr_left.end());

        std::sort(combined.begin(), combined.end(),
                 [](const Candidate& a, const Candidate& b) {
                     return a.tau < b.tau || (a.tau == b.tau && a.side < b.side);
                 });
        return combined;
    }

    const OneSideUnivariateInfo* right() const { return right_.get(); }
    const OneSideUnivariateInfo* left() const { return left_.get(); }
};

} // namespace changepoint