#include "ARpInfo.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace changepoint {
namespace {

/*
  AR(p)-focus, as in Algorithms 1-2 and Appendix C of the AR(p)-focus paper.

  The first p observations are only used to whiten the others, y_t = x_t - sum_j rho_j x_(t-j). Write v_m = 1 -
  sum_(j<m) rho_j and g_m = sum_(j>=m) rho_j (the paper's u_m) for m = 1..p, and Phi = v_(p+1) = 1 - sum_j rho_j. A
  change after y_tau from mu0 to mu1 gives the log-likelihood, up to a constant,
    Cost_(tau,n)(mu0, mu1) = A mu0^2 + B mu1^2 + C mu0 mu1 + D mu0 + E mu1 + F,
  with the coefficients of Appendix C. The changes whose p transitional observations have been seen are kept in two
  stacks, for mu1 > mu0 and for mu1 < mu0. An element holds its change, its coefficients and l, the Inter of its
  curve with that of the element below it: a later change at tau' beats the one at tau for mu1 above (below)
    Inter = B_(tau,tau') / (h Phi^2),  h = tau' - tau,
    B_(tau,tau') = 2 (sum_(t=tau+1)^(tau'+p) v_(t-tau) y_t - sum_(t=tau'+1)^(tau'+p) v_(t-tau') y_t),
  with v_m = Phi for m > p. A new change removes the elements from the top while its Inter with them is at most (at
  least) their l. With the pre-change mean known, l is floored at 0, since a change must also beat no change; with it
  unknown, the first change has l = -inf (+inf) and is never removed.

  At time n (Algorithm 1), the change at n - p - 1 joins both stacks, every coefficient is updated with y_n (D.3), and
  the statistic is the largest of the stacked quadratics and of the last p changes, from (D.2). Changes are counted
  on the whitened observations, from 1; CostsArp adds p.

  Each of the last p changes keeps its sums over its transitional observations seen so far, which y_n extends by one
  term, so that an observation costs O(p) besides the stacks (whitening it is O(p) too). The change at n - p - 1 joins
  the stacks with its completed sums.
*/

struct Quad {
  int tau;
  double W;                        // sum_m v_m y_(tau+m)
  double M;                        // M_(tau+p), the sum of y up to tau + p
  double S;                        // M_tau, the sum of y up to tau
  double A, B, C, D, E, F;         // the coefficients of Cost_(tau,n)
  double l;                        // its Inter with the element below
};

// One of the last p changes, at tau = n - r: its sums over m = 1..r
struct Pending {
  double E;                        // sum_m v_m y_(tau+m)
  double gy;                       // sum_m g_m y_(tau+m)
  double window;                   // sum_m y_(tau+m)
};

struct ArpState {
  ArpState(const std::vector<double>& rho_, bool known_)
      : p(static_cast<int>(rho_.size())), known(known_), rho(rho_), v(p), g(p),
        cv2(p + 1, 0.0), cg2(p + 1, 0.0), cgv(p + 1, 0.0) {
    double sum_rho = 0.0;
    for (int j = 0; j < p; ++j) sum_rho += rho[j];
    Phi = 1.0 - sum_rho;
    double before = 0.0;           // sum_(j<m) rho_j
    for (int m = 1; m <= p; ++m) {
      v[m - 1] = 1.0 - before;
      g[m - 1] = sum_rho - before;
      cv2[m] = cv2[m - 1] + v[m - 1] * v[m - 1];
      cg2[m] = cg2[m - 1] + g[m - 1] * g[m - 1];
      cgv[m] = cgv[m - 1] + g[m - 1] * v[m - 1];
      before += rho[m - 1];
    }
    x_last.reserve(p);
    pending.assign(p, Pending{0.0, 0.0, 0.0});
  }

  int p;
  bool known;                      // pre-change mean known (and the data centred on it)
  std::vector<double> rho;
  double Phi = 1.0;
  std::vector<double> v, g;        // v_m and g_m, m = 1..p
  std::vector<double> cv2, cg2, cgv;   // cumulative sums of v_m^2, g_m^2 and g_m v_m, c[r] = sum_(m<=r)

  std::vector<double> x_last;      // last p observations, a ring: x_head is the slot of the oldest
  int x_head = 0;
  std::vector<Pending> pending;    // the last p changes, the one at tau in slot tau mod p
  int n = 0;                       // whitened observations so far
  double M = 0.0;                  // their sum, M_n
  double sumsq = 0.0;              // and the sum of their squares

  std::vector<Quad> pos, neg;      // the stacks for mu1 > mu0 and mu1 < mu0
};

// Inter: the non-zero mu1 where the curves of the change at Q.tau and of the later one at N.tau cross, for mu0 = 0
inline double inter(const ArpState& st, const Quad& Q, const Quad& N) {
  const double denom = static_cast<double>(N.tau - Q.tau) * st.Phi * st.Phi;
  if (denom == 0.0) return 0.0;
  return 2.0 * (Q.W + st.Phi * (N.M - Q.M) - N.W) / denom;
}

// Algorithm 2: prune the stack with the new change q, then push it
inline void insert(const ArpState& st, std::vector<Quad>& stack, Quad q, bool increase) {
  double c = 0.0;
  while (!stack.empty()) {
    c = inter(st, stack.back(), q);
    if (increase ? (c <= stack.back().l) : (c >= stack.back().l)) {
      stack.pop_back();
    } else {
      break;
    }
  }
  const double inf = std::numeric_limits<double>::infinity();
  if (stack.empty()) {
    q.l = st.known ? 0.0 : (increase ? -inf : inf);
  } else {
    q.l = st.known ? (increase ? std::max(0.0, c) : std::min(0.0, c)) : c;
  }
  stack.push_back(q);
}

// twice the gain of the quadratic over no change: maximised over (mu0, mu1), or over mu1 with mu0 = 0 when known;
// null_max is the maximum of the no-change cost
inline double lr(const ArpState& st, double null_max, double A, double B, double C, double D, double E, double F) {
  if (st.known) return 2.0 * (F - E * E / (4.0 * B) - null_max);
  const double inv = 1.0 / (4.0 * A * B - C * C);
  const double mu0 = (C * E - 2.0 * B * D) * inv;
  const double mu1 = (C * D - 2.0 * A * E) * inv;
  return 2.0 * (F + 0.5 * (D * mu0 + E * mu1) - null_max);
}

inline void consider(double value, int tau, double& best, int& best_tau) {
  if (value > best) {
    best = value;
    best_tau = tau;
  }
}

}  // namespace

void arp_detector_update_impl(double obs,
                              const std::vector<double>& rho,
                              bool known_prechange,
                              void*& opaque_states,
                              double& out_max_stat,
                              int& out_cpt) {
  out_max_stat = -1.0;
  out_cpt = -1;
  if (rho.empty()) return;

  ArpState* st = reinterpret_cast<ArpState*>(opaque_states);
  if (st == nullptr) {
    st = new ArpState(rho, known_prechange);
    opaque_states = st;
  }
  const int p = st->p;

  // the first p observations only start the whitening
  if (static_cast<int>(st->x_last.size()) < p) {
    st->x_last.push_back(obs);
    return;
  }
  double y = obs;
  int k = st->x_head;                                // x_(t-j) is in slot x_head - j, mod p
  for (int j = 1; j <= p; ++j) {
    k = (k == 0 ? p : k) - 1;
    y -= st->rho[j - 1] * st->x_last[k];
  }
  st->x_last[st->x_head] = obs;
  st->x_head = (st->x_head + 1 == p) ? 0 : st->x_head + 1;

  const int n = ++st->n;
  st->M += y;
  st->sumsq += y * y;

  // Algorithm 2, lines 1-10: the change at n - p - 1, whose p transitional observations end at n - 1, joins both stacks.
  // It is in the slot that the change at n - 1 takes next.
  Pending& slot = st->pending[(n - 1) % p];
  if (n - p - 1 >= 1) {
    Quad q;
    q.tau = n - p - 1;
    q.W = slot.E;
    q.M = st->M - y;                                 // M_(tau+p) = M_(n-1)
    const double M_tau = q.M - slot.window;
    q.S = M_tau;
    // (D.2) at time tau + p
    q.A = -0.5 * (q.tau * st->Phi * st->Phi + st->cg2[p]);
    q.B = -0.5 * st->cv2[p];
    q.C = st->cgv[p];
    q.D = st->Phi * M_tau - slot.gy;
    q.E = q.W;
    q.F = -0.5 * (st->sumsq - y * y);
    insert(*st, st->pos, q, true);
    insert(*st, st->neg, q, false);
  }
  if (n < 2) return;
  slot = Pending{0.0, 0.0, 0.0};                    // the change at n - 1

  double best = -1.0;
  int best_tau = -1;
  const double null_max = (st->known ? 0.0 : st->M * st->M / (2.0 * n)) - 0.5 * st->sumsq;

  // Algorithm 2, line 11: update the coefficients with y_n (D.3); Algorithm 1, line 3: maximise the stacked quadratics
  const double half_phi2 = 0.5 * st->Phi * st->Phi;
  for (std::vector<Quad>* stack : {&st->pos, &st->neg}) {
    for (Quad& q : *stack) {
      q.B -= half_phi2;
      q.E += st->Phi * y;
      q.F -= 0.5 * y * y;
      consider(lr(*st, null_max, q.A, q.B, q.C, q.D, q.E, q.F), q.tau, best, best_tau);
    }
  }

  // Algorithm 1, line 4: the last p changes, at n - r for r = 1..p, from (D.2), their sums extended with y_n = y_(tau+r)
  int s = n % p;
  for (int r = 1; r <= std::min(p, n - 1); ++r) {
    s = (s == 0 ? p : s) - 1;                        // the slot of n - r
    Pending& c = st->pending[s];
    c.E += st->v[r - 1] * y;
    c.gy += st->g[r - 1] * y;
    c.window += y;
    const int tau = n - r;
    const double A = -0.5 * (tau * st->Phi * st->Phi + st->cg2[r]);
    const double D = st->Phi * (st->M - c.window) - c.gy;
    consider(lr(*st, null_max, A, -0.5 * st->cv2[r], st->cgv[r], D, c.E, -0.5 * st->sumsq), tau, best, best_tau);
  }

  out_max_stat = best;
  out_cpt = best_tau;
}

void cleanup_arp_states(void* opaque_states) {
  delete reinterpret_cast<ArpState*>(opaque_states);
}

// The changes in the stack for increases (side "right"), then in the one for decreases ("left"). tau is counted on
// the original data, as the changepoint in CostsArp, and st is M_tau.
void arp_detector_candidates(const void* opaque_states, std::vector<Candidate>& out) {
  out.clear();
  const ArpState* st = reinterpret_cast<const ArpState*>(opaque_states);
  if (st == nullptr) return;
  out.reserve(st->pos.size() + st->neg.size());
  for (const Quad& q : st->pos) out.emplace_back(q.S, static_cast<double>(q.tau + st->p), "right");
  for (const Quad& q : st->neg) out.emplace_back(q.S, static_cast<double>(q.tau + st->p), "left");
}

}  // namespace changepoint
