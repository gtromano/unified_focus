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
  std::vector<Candidate> candidates_;         // Candidate changepoints (pre-allocated)
  size_t k_;                                  // Active candidate count (index-based)

public:
  Info(const std::vector<double>& sn = {0.0},
       int n = 0,
       const std::vector<double>& theta0 = {})
    : sn_(sn), n_(n), theta0_(theta0), k_(0) {
    // Pre-allocate candidate storage
    candidates_.reserve(30);
  }

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

  // Add one observation
  virtual void update(const std::vector<double>& y) {
    // Update cumulative state
    add_new_point(y);

    // Prune candidates (in-place)
    candidates_ = prune(candidates_);

    // Append new candidate(s) for current time
    append_new_candidate();
    // std::cout << k_ << std::endl;
  }

  virtual void prune_inplace() {
    std::cout << "Running prune in place" << std::endl;
    // Default: no pruning
  }

  virtual void append_new_candidate() {
    candidates_.push_back(Candidate(sn_, n_, "right"));
  }

  virtual std::vector<Candidate> prune(const std::vector<Candidate>& candidates) const {
    return candidates;
  }

  const std::vector<double>& sn() const { return sn_; }
  int n() const { return n_; }
  const std::vector<double>& theta0() const { return theta0_; }
  bool has_theta0() const { return !theta0_.empty() && !std::isnan(theta0_[0]); }

  // Default: return the internal vector (may include inactive pre-allocated slots).
  // Derived classes (e.g. OneSideUnivariateInfo) can override to return only active entries.
  virtual const std::vector<Candidate>& candidates() const {
    return candidates_;
  }

  size_t active_candidate_count() const { return k_; }
};

// One-sided univariate Info with monotone pruning
class OneSideUnivariateInfo : public Info {
private:
  std::string side_;

  // Cache to return only active candidates from candidates_ (first k_ elements).
  // Mutable so it can be updated in const candidates() method.
  mutable std::vector<Candidate> active_cache_;
  mutable bool cache_valid_ = false;

  // Helper to invalidate cache whenever active set changes
  void invalidate_cache() {
    cache_valid_ = false;
  }

  // Helper to rebuild active_cache_ from candidates_ and k_
  void rebuild_active_cache() const {
    active_cache_.clear();
    active_cache_.reserve(k_);
    for (size_t i = 0; i < k_; ++i) {
      active_cache_.push_back(candidates_[i]);
    }
    cache_valid_ = true;
  }

public:
  OneSideUnivariateInfo(double theta0 = std::numeric_limits<double>::quiet_NaN(),
                        double sn = 0.0, int n = 0,
                        const std::string& side = "right")
    : Info({sn}, n, std::isnan(theta0) ? std::vector<double>() : std::vector<double>({theta0})),
      side_(side) {
    if (side != "right" && side != "left") {
      throw std::invalid_argument("side must be 'right' or 'left'");
    }

    // Pre-allocate candidate storage
    candidates_.reserve(30);
    for (int i = 0; i < 30; i++) {
      candidates_.push_back(Candidate({0.0}, 0, side_));
    }

    // Initialize with first candidate
    candidates_[0] = Candidate(sn_, n_, side_);
    k_ = 1;

    // Initialize cache
    invalidate_cache();
  }

  std::vector<Candidate> new_candidate() const override {
    return {Candidate(sn_, n_, side_)};
  }

  void update(const std::vector<double>& y) override {
    // Update cumulative state
    add_new_point(y);

    // Prune candidates (in-place)
    prune_inplace();

    // Append new candidate(s) for current time
    append_new_candidate();
  }

  void prune_inplace() override {
    if (k_ <= 1) {
      return;
    }

    // Monotone pruning - work backwards from end
    while (k_ > 1) {
      const auto& c1 = candidates_[k_ - 1];
      const auto& c0 = candidates_[k_ - 2];

      int tau1 = c1.tau;
      int tau0 = c0.tau;
      int denom1 = n_ - tau1;
      int denom0 = n_ - tau0;

      double num1 = sn_[0] - c1.scalar_st();
      double num0 = sn_[0] - c0.scalar_st();

      double ratio1 = (denom1 > 0) ? (num1 / denom1) : std::numeric_limits<double>::infinity();
      double ratio0 = (denom0 > 0) ? (num0 / denom0) : std::numeric_limits<double>::infinity();

      bool cond = (side_ == "right") ? (ratio1 <= ratio0) : (ratio1 >= ratio0);

      if (cond) {
        k_--;
        invalidate_cache();
        if (k_ == 1) break;
      } else {
        break;
      }
    }
  }

  void append_new_candidate() override {
    // Reuse existing slots or expand if needed
    if (k_ < candidates_.size()) {
      candidates_[k_].st = sn_;
      candidates_[k_].tau = n_;
      candidates_[k_].side = side_;
    } else {
      candidates_.push_back(Candidate(sn_, n_, side_));
    }
    k_++;
    invalidate_cache();
  }

  // Override to return only active candidates (first k_ elements)
  const std::vector<Candidate>& candidates() const override {
    if (!cache_valid_) {
      rebuild_active_cache();
    }
    return active_cache_;
  }

  const std::string& side() const { return side_; }
};

// Two-sided univariate Info
class UnivariateInfo : public Info {
private:
  std::unique_ptr<OneSideUnivariateInfo> right_;
  std::unique_ptr<OneSideUnivariateInfo> left_;

  // -----------------------------
  // Combined cache for active candidates
  // -----------------------------
  mutable std::vector<Candidate> combined_cache_;  // holds concatenated right+left active candidates
  mutable bool cache_valid_ = false;              // true if combined_cache_ is up-to-date

public:
  UnivariateInfo(double theta0 = std::numeric_limits<double>::quiet_NaN(),
                 double sn = 0.0, int n = 0)
    : Info({sn}, n, std::isnan(theta0) ? std::vector<double>() : std::vector<double>({theta0})) {

    right_ = std::make_unique<OneSideUnivariateInfo>(theta0, sn, n, "right");
    left_ = std::make_unique<OneSideUnivariateInfo>(theta0, sn, n, "left");

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

    // Invalidate cache since active candidates changed
    cache_valid_ = false;
  }

  const std::vector<Candidate>& candidates() const override {
    if (!cache_valid_) {
      combined_cache_.clear();
      const auto& right_cands = right_->candidates();
      const auto& left_cands  = left_->candidates();
      combined_cache_.reserve(right_cands.size() + left_cands.size());
      combined_cache_.insert(combined_cache_.end(), right_cands.begin(), right_cands.end());
      combined_cache_.insert(combined_cache_.end(), left_cands.begin(),  left_cands.end());
      cache_valid_ = true;
    }
    return combined_cache_;
  }

  const OneSideUnivariateInfo* right() const { return right_.get(); }
  const OneSideUnivariateInfo* left() const { return left_.get(); }
};

} // namespace changepoint
