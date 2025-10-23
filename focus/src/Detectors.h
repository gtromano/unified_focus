#pragma once

#include "Info.h"
#include "Candidate.h"
#include "Costs.h"
#include <memory>
#include <vector>
#include <optional>
#include <algorithm>
#include <stdexcept>

namespace changepoint {

struct ChangepointResult {
    int stopping_time;
    std::optional<int> changepoint;
    std::optional<double> stat;
};

class Detector {
private:
    std::shared_ptr<Info> cs_;
    CostFunction compute_costs_fn_;
    std::vector<Candidate> pieces_;

public:
    Detector(std::shared_ptr<Info> cs, CostFunction compute_costs_fn)
        : cs_(std::move(cs)), compute_costs_fn_(std::move(compute_costs_fn)) {
        
        if (!cs_) {
            throw std::invalid_argument("cs must be a valid Info instance");
        }
        
        // Initialize with new candidate(s) from cs
        auto initial = cs_->new_candidate();
        pieces_.insert(pieces_.end(), initial.begin(), initial.end());
    }

    void update(const std::vector<double>& y) {
        // Update cumulative state
        cs_->update(y);
        
        // Prune candidates
        pieces_ = cs_->prune(pieces_);
        
        // Append new candidate(s) for current time
        auto new_cands = cs_->new_candidate();
        pieces_.insert(pieces_.end(), new_cands.begin(), new_cands.end());
    }

    double statistic() const {
        if (pieces_.empty()) {
            return 0.0;
        }
        
        std::vector<double> vals = compute_costs_fn_(pieces_, *cs_);
        return *std::max_element(vals.begin(), vals.end());
    }

    ChangepointResult changepoint() const {
        ChangepointResult result;
        result.stopping_time = cs_->n();
        
        if (pieces_.size() <= 1) {
            result.changepoint = std::nullopt;
            result.stat = std::nullopt;
            return result;
        }
        
        // Exclude last candidate(s) (dummy for current time)
        auto last_cands = cs_->new_candidate();
        int exclude_count = last_cands.size();
        
        if (pieces_.size() <= static_cast<size_t>(exclude_count)) {
            result.changepoint = std::nullopt;
            result.stat = std::nullopt;
            return result;
        }
        
        std::vector<Candidate> considered(pieces_.begin(), 
                                         pieces_.end() - exclude_count);
        
        if (considered.empty()) {
            result.changepoint = std::nullopt;
            result.stat = std::nullopt;
            return result;
        }
        
        std::vector<double> vals = compute_costs_fn_(considered, *cs_);
        auto max_it = std::max_element(vals.begin(), vals.end());
        int max_idx = std::distance(vals.begin(), max_it);
        
        result.changepoint = considered[max_idx].tau;
        result.stat = *max_it;
        
        return result;
    }

    // Getters for testing/inspection
    const std::vector<Candidate>& pieces() const { return pieces_; }
    const Info& info() const { return *cs_; }
};

} // namespace changepoint