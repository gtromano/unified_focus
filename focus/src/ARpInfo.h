#pragma once

#include "Info.h"
#include <vector>
#include <memory>
#include <cmath>
#include <stdexcept>

namespace changepoint {

/*
  ARpInfo:
  Online (sequential) interface for ARP (AutoRegressive Process) changepoint detection.

  - Stores series data, ARP parameters (rho, p), and current statistics.
  - update(y): Takes a scalar observation and triggers an update via external function.
  - Maintains the maximum test statistic and corresponding changepoint across all four
    detection directions (right_pos, left_pos, right_neg, left_neg).
  - candidates() returns a dummy/empty candidate list (the real computation happens in cost function).
*/

// Forward declaration from focus_ARp.cpp
// This function will be implemented in focus_ARp.cpp and handles the actual state updates
void arp_detector_update_impl(double obs,
                               const std::vector<double>& rho,
                               int p,
                               int buf_max,
                               bool known_prechange,
                               double n,
                               void*& opaque_states,  // Opaque pointer to hold the four State objects
                               double& out_max_stat,
                               int& out_cpt);

// Cleanup helper for opaque states
void cleanup_arp_states(void* opaque_states);

class ARpInfo : public Info {
public:
  // Constructor: Initialize with ARP parameters
  ARpInfo(const std::vector<double>& rho,
          bool known_prechange = false,
          double mu0 = 0.0)              // <-- added
    : Info(std::vector<double>{0.0}, 0),
      rho_(normalise_rho(rho)),
      known_prechange_(known_prechange),
      mu0_(mu0),                         // <-- stored
      p_((int)rho_.size()),
      buf_max_(std::max(2 * p_, p_ + 1)),
      max_stat_(-1.0),
      cpt_(-1),
      cumsum_(0.0),
      opaque_states_(nullptr) {
  }

  virtual ~ARpInfo() {
    // Clean up opaque state if it exists
    if (opaque_states_) {
      changepoint::cleanup_arp_states(opaque_states_);
      opaque_states_ = nullptr;
    }
  }

  // Forbid copy/move to avoid copying State objects
  ARpInfo(const ARpInfo&) = delete;
  ARpInfo& operator=(const ARpInfo&) = delete;

  // update: Process a new observation (univariate scalar)
  void update(const std::vector<double>& y, double lambda = 1.0) override {
    if (y.size() != 1) {
      throw std::invalid_argument("ARpInfo::update requires univariate observation (length 1).");
    }

    double obs = y[0] - mu0_;   // <-- shift by known pre-change mean

    n_ += lambda;
    cumsum_ += obs;

    // Call implementation for all observations starting from n==1
    changepoint::arp_detector_update_impl(obs, rho_, p_, buf_max_,
                                          known_prechange_, n_,
                                          opaque_states_,
                                          max_stat_, cpt_);

    // Sync parent sn_ with cumulative sum
    sn_.assign(1, cumsum_);
  }

  // candidates: Return dummy candidates (actual computation in cost function)
  const std::vector<Candidate>& candidates() const override {
    return dummy_candidates_;
  }

  // Accessors for external cost function
  double max_stat()        const { return max_stat_; }
  int cpt()                const { return cpt_; }
  const std::vector<double>& rho() const { return rho_; }
  int p()                  const { return p_; }
  bool known_prechange()   const { return known_prechange_; }
  double mu0()             const { return mu0_; }   // <-- accessor if needed

private:

  static std::vector<double> normalise_rho(const std::vector<double>& rho) {
    if (rho.size() == 1) {
      return std::vector<double>{rho[0], 0.0};
    }
    return rho;
  }

  // ARP parameters
  std::vector<double> rho_;       // AR coefficients
  bool known_prechange_;          // Whether pre-change mean is known
  double mu0_;                    // Known pre-change mean (if applicable)
  int p_;                         // AR order
  int buf_max_;                   // Buffer size

  // Current statistics
  double max_stat_;               // Maximum test statistic across all four states
  int cpt_;                       // Corresponding changepoint

  // Data
  double cumsum_;                 // Cumulative sum of all observations

  // Dummy candidates (required by interface but not used for ARP)
  mutable std::vector<Candidate> dummy_candidates_;

  // Opaque pointer to hold the four State objects (implementation detail in focus_ARp.cpp)
  void* opaque_states_;
};

} // namespace changepoint
