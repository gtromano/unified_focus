#pragma once

#include "Info.h"
#include <vector>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace changepoint {

/*
  NonparametricInfo:
  - Constructs one UnivariateInfo (two-sided Bernoulli detector) per quantile.
  - update(y): requires y.size() == 1. For each quantile q_i compute indicator
    v_i = (y[0] <= q_i) ? 1.0 : 0.0 and call sub_i->update({v_i}).
  - Throws invalid_argument if y is multivariate (size != 1).
  - Maintains sn_ as vector of cumulative successes (one entry per quantile).
  - candidates(): concatenates candidates from all sub-detectors (pruning is local).
*/
class NonparametricInfo : public Info {
public:
  // Construct from quantiles and optional per-quantile theta0 vector.
  // If theta0_vec is empty or element is NaN -> sub-detector theta0 will be NaN.
  explicit NonparametricInfo(const std::vector<double>& quants,
                             const std::vector<double>& theta0_vec = {})
    : Info(), quants_(quants) {

    if (quants_.empty()) {
      throw std::invalid_argument("NonparametricInfo requires at least one quantile");
    }

    // create one UnivariateInfo per quantile
    sub_infos_.reserve(quants_.size());
    for (size_t i = 0; i < quants_.size(); ++i) {
      double th = std::numeric_limits<double>::quiet_NaN();
      if (i < theta0_vec.size()) th = theta0_vec[i];
      sub_infos_.push_back(std::make_unique<UnivariateInfo>(th, 0.0, 0));
    }

    // initialize parent sn_ to zeros (one entry per quantile)
    sn_.assign(quants_.size(), 0.0);
    n_ = 0;
    // store theta0_ vector aggregated from sub-infos (NaNs if missing)
    theta0_.clear();
    theta0_.reserve(quants_.size());
    for (const auto& s : sub_infos_) {
      if (s->has_theta0()) theta0_.push_back(s->theta0()[0]);
      else theta0_.push_back(std::numeric_limits<double>::quiet_NaN());
    }

    combined_cache_valid_ = false;
  }

  virtual ~NonparametricInfo() = default;

  // forbid copy to avoid copying unique_ptrs
  NonparametricInfo(const NonparametricInfo&) = delete;
  NonparametricInfo& operator=(const NonparametricInfo&) = delete;

  // number of sub-detectors (equals number of quantiles)
  size_t n_detectors() const noexcept { return sub_infos_.size(); }

  // Access a sub-detector (const)
  const UnivariateInfo* sub_detector(size_t idx) const {
    if (idx >= sub_infos_.size()) throw std::out_of_range("sub_detector index out of range");
    return sub_infos_[idx].get();
  }

  // Access a sub-detector (non-const)
  UnivariateInfo* sub_detector(size_t idx) {
    if (idx >= sub_infos_.size()) throw std::out_of_range("sub_detector index out of range");
    return sub_infos_[idx].get();
  }

  // update: preprocess y into indicators y <= q_i and forward to each UnivariateInfo
  void update(const std::vector<double>& y) override {
    if (y.size() != 1) {
      throw std::invalid_argument("NonparametricInfo::update requires univariate observation (length 1).");
    }
    double raw = y[0];

    for (size_t i = 0; i < sub_infos_.size(); ++i) {
      double indicator = (raw <= quants_[i]) ? 1.0 : 0.0;
      std::vector<double> single = { indicator };
      sub_infos_[i]->update(single);
    }

    // sync parent sn_ and n_ with sub-detectors (they should all have same n)
    n_ = sub_infos_[0]->n();
    sn_.clear();
    sn_.reserve(sub_infos_.size());
    for (const auto& s : sub_infos_) {
      const auto& s_sn = s->sn();
      // UnivariateInfo maintains sn as size-1 vector
      if (!s_sn.empty()) sn_.push_back(s_sn[0]);
      else sn_.push_back(0.0);
    }

    // invalidate combined candidates cache
    combined_cache_valid_ = false;
  }

  // return concatenated candidates from all sub-detectors
  const std::vector<Candidate>& candidates() const override {
    if (!combined_cache_valid_) {
      combined_cache_.clear();
      // reserve approximate size
      size_t total = 0;
      for (const auto& s : sub_infos_) total += s->candidates().size();
      combined_cache_.reserve(total);

      for (const auto& s : sub_infos_) {
        const auto& c = s->candidates();
        combined_cache_.insert(combined_cache_.end(), c.begin(), c.end());
      }
      combined_cache_valid_ = true;
    }
    return combined_cache_;
  }

  std::vector<Candidate> new_candidate() const override {
    // Not typically used: return a concatenation of each sub new_candidate
    std::vector<Candidate> out;
    for (const auto& s : sub_infos_) {
      auto nc = s->new_candidate();
      out.insert(out.end(), nc.begin(), nc.end());
    }
    return out;
  }

  // // aggregate helpers
  // bool has_theta0() const override {
  //   for (const auto& s : sub_infos_) if (s->has_theta0()) return true;
  //   return false;
  // }
  //
  // std::vector<double> theta0() const override {
  //   std::vector<double> out;
  //   out.reserve(sub_infos_.size());
  //   for (const auto& s : sub_infos_) {
  //     if (s->has_theta0()) out.push_back(s->theta0()[0]);
  //     else out.push_back(std::numeric_limits<double>::quiet_NaN());
  //   }
  //   return out;
  // }

  // statistic: sum of sub-detector statistics (convenience)
  double statistic() const {
    double acc = 0.0;
    for (const auto& s : sub_infos_) acc += s->right()->n() ? s->right()->n() : 0.0; // not meaningful generally
    // prefer costs to compute the real NP-FOCuS stat (sum of bernoulli max stats)
    return acc;
  }

private:
  std::vector<double> quants_;
  std::vector<std::unique_ptr<UnivariateInfo>> sub_infos_;

  // cache for concatenated candidates
  mutable std::vector<Candidate> combined_cache_;
  mutable bool combined_cache_valid_;
};

} // namespace changepoint
