#pragma once

#include "Candidate.h"
#include <vector>
#include <memory>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <cmath>

namespace changepoint {

// Minimal base class for Info types
class Info {
protected:
  std::vector<double> sn_;      // Cumulative sum
  double n_;                     // Effective number of observations (accounting for lambda)

public:
  Info(const std::vector<double>& sn = {0.0}, double n = 0.0)
    : sn_(sn), n_(n) {}

  virtual ~Info() = default;

  // Core interface methods
  virtual std::vector<Candidate> new_candidate() const {
    return {Candidate(sn_, n_, "right")};
  }

  virtual void add_new_point(const std::vector<double>& y, double lambda = 1.0) {
    n_ += lambda;
    if (sn_.size() != y.size()) {
      sn_.resize(y.size(), 0.0);
    }
    for (size_t i = 0; i < sn_.size(); ++i) {
      sn_[i] += y[i];
    }
  }

  virtual void update(const std::vector<double>& y, double lambda = 1.0) = 0;
  virtual const std::vector<Candidate>& candidates() const = 0;

  // Accessors for common state
  const std::vector<double>& sn() const { return sn_; }
  double n() const { return n_; }
};

// ---------------------------------------------------------------------------
// Base class for Info types that maintain a candidate list
// ---------------------------------------------------------------------------
class CandidateListInfo : public Info {
protected:
  std::vector<Candidate> candidates_;

public:
  CandidateListInfo(const std::vector<double>& sn = {0.0}, double n = 0.0)
    : Info(sn, n) {
    candidates_.reserve(30);
  }

  virtual ~CandidateListInfo() = default;

  const std::vector<Candidate>& candidates() const override {
    return candidates_;
  }

  virtual std::vector<Candidate> prune(const std::vector<Candidate>& candidates) const {
    return candidates;
  }

  void update(const std::vector<double>& y, double lambda = 1.0) override {
    add_new_point(y, lambda);
    candidates_ = prune(candidates_);
    append_new_candidate();
  }

protected:
  virtual void append_new_candidate() {
    candidates_.push_back(Candidate(sn_, n_, "right"));
  }
};

// ---------------------------------------------------------------------------
// One-sided univariate Info with in-place pruning and index-based management
// ---------------------------------------------------------------------------
class OneSideUnivariateInfo : public Info {
private:
  std::string side_;
  std::vector<Candidate> candidates_;  // Pre-allocated storage
  size_t left_;                       // Start of active candidates
  size_t k_;                          // End of active candidates (exclusive)
  double anomaly_intensity_;

  mutable std::vector<Candidate> active_cache_;
  mutable bool cache_valid_;

  void invalidate_cache() {
    cache_valid_ = false;
  }

  void rebuild_active_cache() const {
    active_cache_.clear();
    active_cache_.reserve(k_ - left_);
    for (size_t i = left_; i < k_; ++i) {
      active_cache_.push_back(candidates_[i]);
    }
    cache_valid_ = true;
  }

  void compact_if_needed() {
    // If we're about to exceed capacity, compact the array
    if (k_ >= candidates_.size() && left_ > 0) {
      size_t active_count = k_ - left_;
      for (size_t i = 0; i < active_count; ++i) {
        candidates_[i] = candidates_[left_ + i];
      }
      left_ = 0;
      k_ = active_count;
    }
  }

public:
  OneSideUnivariateInfo(double sn = 0.0, double n = 0.0,
                        const std::string& side = "right",
                        double anomaly_intensity = std::numeric_limits<double>::quiet_NaN())
    : Info({sn}, n), side_(side), left_(0), k_(0), 
      anomaly_intensity_(anomaly_intensity), cache_valid_(false) {
    
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
  }

  std::vector<Candidate> new_candidate() const override {
    return {Candidate(sn_, n_, side_)};
  }

  void update(const std::vector<double>& y, double lambda = 1.0) override {
    add_new_point(y, lambda);
    prune_inplace();
    append_new_candidate();
  }

  const std::vector<Candidate>& candidates() const override {
    if (!cache_valid_) {
      rebuild_active_cache();
    }
    return active_cache_;
  }

  const std::string& side() const { return side_; }
  size_t active_candidate_count() const { return k_ - left_; }

private:
  void prune_inplace() {
    // First, prune from the left based on anomaly_intensity
    prune_left_by_intensity();
    
    // Then, prune from the right based on monotonicity
    prune_right_by_monotonicity();
  }

  void prune_left_by_intensity() {
    if (std::isnan(anomaly_intensity_) || anomaly_intensity_ <= 0.0) {
      return;  // No intensity-based pruning
    }

    while (left_ < k_) {
      const auto& c = candidates_[left_];
      double tau = c.tau;
      double denom = n_ - tau;
      double num = sn_[0] - c.scalar_st();
      
      double ratio = (denom > 0) ? (num / denom) : std::numeric_limits<double>::infinity();
      double abs_ratio = std::abs(ratio);

      if (abs_ratio < anomaly_intensity_) {
        left_++;
        invalidate_cache();
      } else {
        break;  // Stop at first candidate that meets threshold
      }
    }
  }

  void prune_right_by_monotonicity() {
    if (k_ - left_ <= 1) {
      return;
    }

    while (k_ - left_ > 1) {
      const auto& c1 = candidates_[k_ - 1];
      const auto& c0 = candidates_[k_ - 2];

      double tau1 = c1.tau;
      double tau0 = c0.tau;
      double denom1 = n_ - tau1;
      double denom0 = n_ - tau0;

      double num1 = sn_[0] - c1.scalar_st();
      double num0 = sn_[0] - c0.scalar_st();

      double ratio1 = (denom1 > 0) ? (num1 / denom1) : std::numeric_limits<double>::infinity();
      double ratio0 = (denom0 > 0) ? (num0 / denom0) : std::numeric_limits<double>::infinity();

      bool cond = (side_ == "right") ? (ratio1 <= ratio0) : (ratio1 >= ratio0);

      if (cond) {
        k_--;
        invalidate_cache();
        if (k_ - left_ == 1) break;
      } else {
        break;
      }
    }
  }

  void append_new_candidate() {
    compact_if_needed();
    
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
};

// ---------------------------------------------------------------------------
// Two-sided univariate Info
// ---------------------------------------------------------------------------
class UnivariateInfo : public Info {
private:
  std::unique_ptr<OneSideUnivariateInfo> right_;
  std::unique_ptr<OneSideUnivariateInfo> left_;

  mutable std::vector<Candidate> combined_cache_;
  mutable bool cache_valid_;

public:
  UnivariateInfo(double sn = 0.0, double n = 0.0, 
                 double anomaly_intensity = std::numeric_limits<double>::quiet_NaN())
    : Info({sn}, n), cache_valid_(false) {

    right_ = std::make_unique<OneSideUnivariateInfo>(sn, n, "right", anomaly_intensity);
    left_ = std::make_unique<OneSideUnivariateInfo>(sn, n, "left", anomaly_intensity);
  }

  std::vector<Candidate> new_candidate() const override {
    auto right_cand = right_->new_candidate();
    auto left_cand = left_->new_candidate();
    std::vector<Candidate> result;
    result.insert(result.end(), right_cand.begin(), right_cand.end());
    result.insert(result.end(), left_cand.begin(), left_cand.end());
    return result;
  }

  void add_new_point(const std::vector<double>& y, double lambda = 1.0) override {
    right_->add_new_point(y, lambda);
    left_->add_new_point(y, lambda);
    n_ = right_->n();
    sn_ = right_->sn();
  }

  void update(const std::vector<double>& y, double lambda = 1.0) override {
    right_->update(y, lambda);
    left_->update(y, lambda);
    n_ = right_->n();
    sn_ = right_->sn();
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