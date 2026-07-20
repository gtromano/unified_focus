#include "ARpInfo.h"

#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <limits>

namespace changepoint {
namespace {

struct Triple {
  int tau = 0;
  std::vector<double> M_tau;
  double l = std::numeric_limits<double>::quiet_NaN();
  double A = std::numeric_limits<double>::quiet_NaN();
  double B = std::numeric_limits<double>::quiet_NaN();
  double C = std::numeric_limits<double>::quiet_NaN();
  double D = std::numeric_limits<double>::quiet_NaN();
  double E = std::numeric_limits<double>::quiet_NaN();
  double f = std::numeric_limits<double>::quiet_NaN();
};

struct State {
  std::vector<Triple> triples;
  std::vector<Triple> triples_out;
  double max_val = -1.0;
  int cpt = -1;
};

struct TwoSidedState {
  State state_right;
  State state_left;
};

struct ArpState {
  explicit ArpState(int p_, bool known_prechange_)
      : p(p_), known_prechange(known_prechange_) {}

  int p;
  bool known_prechange;
  bool initialized = false;
  int current_i = 2;
  std::vector<double> data_old;
  std::vector<double> x_buf;
  std::vector<double> y_buf;
  double M_n = 0.0;
  double sum_square = 0.0;
  State state_right_pos;
  State state_left_pos;
  State state_right_neg;
  State state_left_neg;
};

inline double dot_product(const std::vector<double>& a, const std::vector<double>& b) {
  double out = 0.0;
  const std::size_t n = std::min(a.size(), b.size());
  for (std::size_t i = 0; i < n; ++i) {
    out += a[i] * b[i];
  }
  return out;
}

inline double sum_vector(const std::vector<double>& values) {
  return std::accumulate(values.begin(), values.end(), 0.0);
}

inline std::vector<double> reverse_copy(const std::vector<double>& values) {
  return std::vector<double>(values.rbegin(), values.rend());
}

inline std::vector<double> make_u(const std::vector<double>& rho) {
  std::vector<double> u(1, 1.0);
  for (int i = 1; i < static_cast<int>(rho.size()); ++i) {
    double s = 0.0;
    for (int j = 0; j < i; ++j) {
      s += rho[j];
    }
    u.push_back(1.0 - s);
  }
  return u;
}

inline std::vector<double> make_v(const std::vector<double>& rho) {
  std::vector<double> v;
  v.reserve(rho.size());
  for (std::size_t i = 0; i < rho.size(); ++i) {
    double s = 0.0;
    for (std::size_t j = i; j < rho.size(); ++j) {
      s += rho[j];
    }
    v.push_back(s);
  }
  return v;
}

inline Triple create_triple(int tau, const std::vector<double>& M_tau, double l = std::numeric_limits<double>::quiet_NaN(),
                            double A = std::numeric_limits<double>::quiet_NaN(),
                            double B = std::numeric_limits<double>::quiet_NaN(),
                            double C = std::numeric_limits<double>::quiet_NaN(),
                            double D = std::numeric_limits<double>::quiet_NaN(),
                            double E = std::numeric_limits<double>::quiet_NaN(),
                            double f = std::numeric_limits<double>::quiet_NaN()) {
  Triple out;
  out.tau = tau;
  out.M_tau = M_tau;
  out.l = l;
  out.A = A;
  out.B = B;
  out.C = C;
  out.D = D;
  out.E = E;
  out.f = f;
  return out;
}

inline std::vector<Triple> exact_form_out_start(const std::vector<double>& rho, int current_n, const std::vector<double>& y_buf) {
  const int p = static_cast<int>(rho.size());
  std::vector<Triple> out;
  out.reserve(std::max(0, current_n - 1));
  const std::vector<double> u = make_u(rho);
  const std::vector<double> v = make_v(rho);
  const double sum_sq = dot_product(y_buf, y_buf);

  for (int tau_out = 1; tau_out <= current_n - 1; ++tau_out) {
    const int width = current_n - tau_out;
    std::vector<double> u_new(u.begin(), u.begin() + width);
    std::vector<double> v_new(v.begin(), v.begin() + width);
    std::vector<double> y_tau(y_buf.begin() + tau_out, y_buf.begin() + current_n);
    double M_tau = 0.0;
    for (int j = 0; j < tau_out; ++j) {
      M_tau += y_buf[j];
    }

    const double A = -0.5 * (tau_out * std::pow(1.0 - sum_vector(rho), 2.0) + dot_product(v_new, v_new));
    const double B = -0.5 * dot_product(u_new, u_new);
    const double C = dot_product(u_new, v_new);
    const double D = -0.5 * (-2.0 * (1.0 - sum_vector(rho)) * M_tau + 2.0 * dot_product(y_tau, v_new));
    const double E = -0.5 * (-2.0 * dot_product(y_tau, u_new));
    const double f = -0.5 * sum_sq;
    out.push_back(create_triple(tau_out, {}, std::numeric_limits<double>::quiet_NaN(), A, B, C, D, E, f));
  }
  return out;
}

inline std::vector<Triple> exact_form_out_start_pre0(const std::vector<double>& rho, int current_n, const std::vector<double>& y_buf) {
  std::vector<Triple> out;
  out.reserve(std::max(0, current_n - 1));
  const std::vector<double> u = make_u(rho);

  for (int tau_out = 1; tau_out <= current_n - 1; ++tau_out) {
    const int width = current_n - tau_out;
    std::vector<double> v_new(u.begin(), u.begin() + width);
    std::vector<double> y_tau(y_buf.begin() + tau_out, y_buf.begin() + current_n);
    const double A = -0.5 * dot_product(v_new, v_new);
    const double B = -0.5 * (-2.0 * dot_product(v_new, y_tau));
    out.push_back(create_triple(tau_out, {}, std::numeric_limits<double>::quiet_NaN(), A, B));
  }
  return out;
}

inline Triple coef_introduce(const Triple& q_new, const std::vector<double>& rho, const std::vector<double>& y_buf, double sum_square) {
  const int p = static_cast<int>(rho.size());
  const std::vector<double> u = make_u(rho);
  const std::vector<double> v = make_v(rho);
  Triple out = q_new;
  const double M_tau = q_new.M_tau[p];
  const std::vector<double> y_tau(y_buf.end() - p, y_buf.end());
  const double A = -0.5 * (q_new.tau * std::pow(1.0 - sum_vector(rho), 2.0) + dot_product(v, v));
  const double B = -0.5 * dot_product(u, u);
  const double C = dot_product(u, v);
  const double D = -0.5 * (-2.0 * (1.0 - sum_vector(rho)) * M_tau + 2.0 * dot_product(y_tau, v));
  const double E = -0.5 * (-2.0 * dot_product(y_tau, u));
  const double f = -0.5 * sum_square;
  out.A = A; out.B = B; out.C = C; out.D = D; out.E = E; out.f = f;
  return out;
}

inline Triple coef_introduce_pre0(const Triple& q_new, const std::vector<double>& rho, const std::vector<double>& y_buf) {
  const int p = static_cast<int>(rho.size());
  const std::vector<double> u = make_u(rho);
  const std::vector<double> y_tau(y_buf.end() - p, y_buf.end());
  Triple out = q_new;
  out.A = -0.5 * dot_product(u, u);
  out.B = dot_product(y_tau, u);
  return out;
}

inline std::vector<Triple> coef_introduce_outre(const Triple& q_new, const std::vector<double>& rho, double sum_square, double M_n, const std::vector<double>& y_buf) {
  const int p = static_cast<int>(rho.size());
  const std::vector<double> u = make_u(rho);
  const std::vector<double> v = make_v(rho);
  std::vector<Triple> out;
  const std::vector<double> y_tau(y_buf.end() - p, y_buf.end());
  for (int tau_out = q_new.tau + 1; tau_out <= q_new.tau + p - 1; ++tau_out) {
    const int index_out = tau_out - q_new.tau;
    const int width = p - index_out;
    std::vector<double> u_new(u.begin(), u.begin() + width);
    std::vector<double> v_new(v.begin(), v.begin() + width);
    const std::vector<double> y_tau_new(y_tau.begin() + index_out, y_tau.end());
    const double M_tau = M_n - sum_vector(y_tau_new);
    const double A = -0.5 * (tau_out * std::pow(1.0 - sum_vector(rho), 2.0) + dot_product(v_new, v_new));
    const double B = -0.5 * dot_product(u_new, u_new);
    const double C = dot_product(u_new, v_new);
    const double D = -0.5 * (-2.0 * (1.0 - sum_vector(rho)) * M_tau + 2.0 * dot_product(y_tau_new, v_new));
    const double E = -0.5 * (-2.0 * dot_product(y_tau_new, u_new));
    const double f = -0.5 * sum_square;
    out.push_back(create_triple(tau_out, {}, std::numeric_limits<double>::quiet_NaN(), A, B, C, D, E, f));
  }
  return out;
}

inline std::vector<Triple> coef_introduce_outre_pre0(const Triple& q_new, const std::vector<double>& rho, const std::vector<double>& y_buf) {
  const int p = static_cast<int>(rho.size());
  const std::vector<double> u = make_u(rho);
  std::vector<Triple> out;
  const std::vector<double> y_tau(y_buf.end() - p, y_buf.end());
  for (int tau_out = q_new.tau + 1; tau_out <= q_new.tau + p - 1; ++tau_out) {
    const int index_out = tau_out - q_new.tau;
    const int width = p - index_out;
    std::vector<double> u_new(u.begin(), u.begin() + width);
    const std::vector<double> y_tau_new(y_tau.begin() + index_out, y_tau.end());
    const double A = -0.5 * dot_product(u_new, u_new);
    const double B = dot_product(y_tau_new, u_new);
    out.push_back(create_triple(tau_out, {}, std::numeric_limits<double>::quiet_NaN(), A, B));
  }
  return out;
}

inline std::vector<double> make_initial_M_tau(const std::vector<double>& y_buf, int p) {
  std::vector<double> M_tau(2 * p + 1, 0.0);
  for (int k = 0; k <= p; ++k) {
    double sum_val = 0.0;
    for (int j = 0; j <= k; ++j) {
      sum_val += y_buf[j];
    }
    M_tau[p + k] = sum_val;
  }
  return M_tau;
}

inline double intersec_point_newcurve(const std::vector<double>& M_new, const std::vector<double>& rho) {
  const int p = static_cast<int>(rho.size());
  std::vector<double> y_new;
  y_new.reserve(2 * p);
  for (int i = 1; i < 2 * p + 1; ++i) {
    y_new.push_back(M_new[i] - M_new[i - 1]);
  }
  std::vector<double> rho_star(1, 1.0);
  for (int i = 1; i < p; ++i) {
    double s = 0.0;
    for (int j = 0; j < i; ++j) {
      s += rho[j];
    }
    rho_star.push_back(1.0 - s);
  }
  std::vector<double> y_new_vec(y_new.begin() + p, y_new.begin() + 2 * p);
  const double denom = dot_product(rho_star, rho_star);
  if (denom == 0.0) {
    return 0.0;
  }
  return 2.0 * dot_product(rho_star, y_new_vec) / denom;
}

inline double intersec_point(const Triple& input_triple, int current_n, const std::vector<double>& rho, double M_n, const std::vector<double>& y_buf) {
  const int p = static_cast<int>(rho.size());
  const int tau_new = current_n - p;
  std::vector<double> rho_star(1, 1.0);
  for (int i = 1; i < p; ++i) {
    double s = 0.0;
    for (int j = 0; j < i; ++j) {
      s += rho[j];
    }
    rho_star.push_back(1.0 - s);
  }
  rho_star.push_back(1.0 - sum_vector(rho));

  std::vector<double> y_tau;
  y_tau.reserve(2 * p);
  for (int i = 1; i < 2 * p + 1; ++i) {
    y_tau.push_back(input_triple.M_tau[i] - input_triple.M_tau[i - 1]);
  }
  std::vector<double> m;
  m.reserve(p + 1);
  const int buf_len = static_cast<int>(y_buf.size());
  for (int i = 0; i < p; ++i) {
    m.push_back(y_tau[p + i] - y_buf[buf_len - p + i]);
  }
  m.push_back(M_n - input_triple.M_tau[2 * p]);

  const double numer = 2.0 * dot_product(rho_star, m);
  const double denom = std::pow(1.0 - sum_vector(rho), 2.0) * (tau_new - input_triple.tau);
  if (denom == 0.0) {
    return 0.0;
  }
  return std::max(0.0, numer / denom);
}

inline std::vector<Triple> coef_update_arp(std::vector<Triple> new_triples, const std::vector<double>& rho, double y_new, bool known_prechange) {
  if (new_triples.size() <= 1) {
    return new_triples;
  }
  const double rho_sum = sum_vector(rho);
  for (std::size_t i = 0; i + 1 < new_triples.size(); ++i) {
    if (!known_prechange) {
      new_triples[i].B -= 0.5 * std::pow(1.0 - rho_sum, 2.0);
      new_triples[i].E += (1.0 - rho_sum) * y_new;
      new_triples[i].f -= 0.5 * y_new * y_new;
    } else {
      new_triples[i].A -= 0.5 * std::pow(1.0 - rho_sum, 2.0);
      new_triples[i].B += (1.0 - rho_sum) * y_new;
    }
  }
  return new_triples;
}

// Only the triples visited while
// popping from the top of the stack get their `l` (intersection point)
// recomputed. Triples never reached by the walk keep whatever `l` they were
// last assigned.
inline std::vector<Triple> prune_triples(std::vector<Triple> triples, int current_n, const std::vector<double>& rho,
                                         double M_n, const std::vector<double>& y_buf, bool side, bool prune) {
  if (triples.empty()) {
    return triples;
  }

  int i = static_cast<int>(triples.size());
  triples[static_cast<std::size_t>(i - 1)].l = intersec_point(triples[static_cast<std::size_t>(i - 1)], current_n, rho, M_n, y_buf);
  const double l_0 = side ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity();

  bool stop = false;
  while (!stop && i >= 1) {
    const double l_i = triples[static_cast<std::size_t>(i - 1)].l;
    const double l_i_1 = (i - 1 == 0) ? l_0 : triples[static_cast<std::size_t>(i - 2)].l;
    const bool cond = side ? (l_i <= l_i_1) : (l_i >= l_i_1);
    if (cond) {
      --i;
      if (i > 0) {
        triples[static_cast<std::size_t>(i - 1)].l = intersec_point(triples[static_cast<std::size_t>(i - 1)], current_n, rho, M_n, y_buf);
      }
    } else {
      stop = true;
    }
  }

  if (prune && i < static_cast<int>(triples.size())) {
    triples.resize(static_cast<std::size_t>(i));
  }
  return triples;
}

inline double optimize_quadratic(const Triple& tr) {
  const double denom = 4.0 * tr.A * tr.B - tr.C * tr.C;
  if (std::abs(denom) < 1e-12) {
    return std::numeric_limits<double>::lowest();
  }
  const double inv = 1.0 / denom;
  const double mu0 = -inv * (2.0 * tr.B * tr.D - tr.C * tr.E);
  const double mu1 = -inv * (-tr.C * tr.D + 2.0 * tr.A * tr.E);
  return tr.A * mu0 * mu0 + tr.B * mu1 * mu1 + tr.C * mu0 * mu1 + tr.D * mu0 + tr.E * mu1 + tr.f;
}

inline double optimize_pre0(const Triple& tr) {
  const double a = tr.A;
  const double b = tr.B;
  const double opt_mu = -b / (2.0 * a);
  return a * opt_mu * opt_mu + b * opt_mu;
}

inline State max_val_compute_arp_out_start(const std::vector<Triple>& triples, int current_n, const std::vector<double>& rho, double sum_square, double M_n) {
  State out;
  const double sum_rho = sum_vector(rho);
  const double a_no = -0.5 * current_n * std::pow(1.0 - sum_rho, 2.0);
  const double b_no = (1.0 - sum_rho) * M_n;
  const double c_no = -0.5 * sum_square;
  const double opt_no = -(b_no / (2.0 * a_no));
  const double max_val_no = a_no * opt_no * opt_no + b_no * opt_no + c_no;
  for (const Triple& tr : triples) {
    const double max_val_change = optimize_quadratic(tr);
    const double value = 2.0 * (max_val_change - max_val_no);
    if (value > out.max_val) {
      out.max_val = value;
      out.cpt = tr.tau;
    }
  }
  return out;
}

inline State max_val_compute_arp_out(const std::vector<Triple>& triples, int current_n, const std::vector<double>& rho, double sum_square, double M_n) {
  return max_val_compute_arp_out_start(triples, current_n, rho, sum_square, M_n);
}

inline State max_val_compute_arp(const std::vector<Triple>& triples, int current_n, const std::vector<double>& rho, double sum_square, double M_n) {
  return max_val_compute_arp_out_start(triples, current_n, rho, sum_square, M_n);
}

inline State max_val_compute_pre0(const std::vector<Triple>& triples) {
  State out;
  for (const Triple& tr : triples) {
    const double value = 2.0 * optimize_pre0(tr);
    if (value > out.max_val) {
      out.max_val = value;
      out.cpt = tr.tau;
    }
  }
  return out;
}

inline State focus_arp_one_iter(const std::vector<double>& y_buf, int i, const std::vector<Triple>& triples_in, double M_n, double sum_square,
                                const std::vector<double>& rho, bool known_prechange, bool side, bool prune, int n) {
  const int p = static_cast<int>(rho.size());
  State out;
  out.triples = triples_in;
  out.triples_out.clear();

  if ((i >= 2) && (i <= (p + 1))) {
    std::vector<Triple> triples_out = known_prechange ? exact_form_out_start_pre0(rho, i, y_buf) : exact_form_out_start(rho, i, y_buf);
    State max_info_out = known_prechange ? max_val_compute_pre0(triples_out) : max_val_compute_arp_out_start(triples_out, i, rho, sum_square, M_n);
    out.max_val = max_info_out.max_val;
    out.cpt = max_info_out.cpt;
    out.triples_out = triples_out;
  }

  bool enter_recursion = (i >= (p + 2));
  if ((i == (p + 2)) && (i <= n)) {
    std::vector<double> M_tau = make_initial_M_tau(y_buf, p);
    Triple triple_1 = create_triple(1, M_tau);
    if (!known_prechange) {
      std::vector<double> y_init(y_buf.begin(), y_buf.begin() + p + 1);
      triple_1 = coef_introduce(triple_1, rho, y_init, sum_square - y_buf[p + 1] * y_buf[p + 1]);
    } else {
      std::vector<double> y_init(y_buf.begin(), y_buf.begin() + p + 1);
      triple_1 = coef_introduce_pre0(triple_1, rho, y_init);
    }
    out.triples = {triple_1};
    enter_recursion = true;
  }

  if (enter_recursion && i <= n) {
    std::vector<Triple> triples = out.triples;
    if (!triples.empty()) {
      triples = prune_triples(std::move(triples), i, rho, M_n, y_buf, side, prune);
    }

    std::vector<double> M_tau_new(1, M_n);
    const int buf_len = static_cast<int>(y_buf.size());
    for (int idx = buf_len; idx >= 1; --idx) {
      double prefix = 0.0;
      for (int j = idx - 1; j < buf_len; ++j) {
        prefix += y_buf[j];
      }
      M_tau_new.insert(M_tau_new.begin(), M_n - prefix);
    }
    while (static_cast<int>(M_tau_new.size()) < 2 * p + 1) {
      M_tau_new.insert(M_tau_new.begin(), 0.0);
    }
    while (static_cast<int>(M_tau_new.size()) > 2 * p + 1) {
      M_tau_new.erase(M_tau_new.begin());
    }

    const double new_l = intersec_point_newcurve(M_tau_new, rho);
    Triple q_new = create_triple(i - p, M_tau_new, new_l);
    if (!known_prechange) {
      out.triples_out = coef_introduce_outre(q_new, rho, sum_square, M_n, y_buf);
      q_new = coef_introduce(q_new, rho, y_buf, sum_square);
    } else {
      out.triples_out = coef_introduce_outre_pre0(q_new, rho, y_buf);
      q_new = coef_introduce_pre0(q_new, rho, y_buf);
    }
    triples.push_back(q_new);
    out.triples = coef_update_arp(triples, rho, y_buf.back(), known_prechange);
  }

  if (!out.triples.empty()) {
    State max_info = known_prechange ? max_val_compute_pre0(out.triples) : max_val_compute_arp(out.triples, i, rho, sum_square, M_n);
    State max_info_out = known_prechange ? max_val_compute_pre0(out.triples_out) : max_val_compute_arp_out(out.triples_out, i, rho, sum_square, M_n);
    if (max_info.max_val > max_info_out.max_val) {
      out.max_val = max_info.max_val;
      out.cpt = max_info.cpt;
    } else {
      out.max_val = max_info_out.max_val;
      out.cpt = max_info_out.cpt;
    }
  } else {
    State max_info_out = known_prechange ? max_val_compute_pre0(out.triples_out) : max_val_compute_arp_out(out.triples_out, i, rho, sum_square, M_n);
    out.max_val = max_info_out.max_val;
    out.cpt = max_info_out.cpt;
  }

  return out;
}

inline TwoSidedState focus_arp_two_sides_update(const std::vector<double>& y_buf, int i, const State& state_right_in, const State& state_left_in,
                                                const std::vector<double>& rho, double M_n, double sum_square, bool known_prechange, bool prune, int n) {
  TwoSidedState out;
  out.state_right = focus_arp_one_iter(y_buf, i, state_right_in.triples, M_n, sum_square, rho, known_prechange, true, prune, n);
  out.state_left = focus_arp_one_iter(y_buf, i, state_left_in.triples, M_n, sum_square, rho, known_prechange, false, prune, n);
  return out;
}

inline std::pair<double, int> compute_focus_arp(const std::vector<double>& data_old, const std::vector<double>& rho, bool known_prechange, bool prune) {
  const int p = static_cast<int>(rho.size());
  if (data_old.size() <= static_cast<std::size_t>(p)) {
    return { -1.0, -1 };
  }

  return { -1.0, -1 };
}

}  // namespace

void arp_detector_update_impl(double obs,
                              const std::vector<double>& rho,
                              int p,
                              int /* buf_max */,
                              bool known_prechange,
                              double /* n */,
                              void*& opaque_states,
                              double& out_max_stat,
                              int& out_cpt) {
  if (p <= 0) {
    out_max_stat = -1.0;
    out_cpt = -1;
    return;
  }

  ArpState* state = reinterpret_cast<ArpState*>(opaque_states);
  if (state == nullptr) {
    state = new ArpState(p, known_prechange);
    opaque_states = state;
  }

  state->data_old.push_back(obs);

  if (!state->initialized) {
    if (static_cast<int>(state->data_old.size()) < p + 1) {
      out_max_stat = -1.0;
      out_cpt = -1;
      return;
    }

    state->x_buf.assign(state->data_old.begin(), state->data_old.begin() + p);
    state->y_buf.push_back(state->data_old[static_cast<std::size_t>(p)] - dot_product(rho, reverse_copy(state->x_buf)));
    state->M_n = sum_vector(state->y_buf);
    state->sum_square = dot_product(state->y_buf, state->y_buf);
    state->initialized = true;
    state->current_i = 2;
    out_max_stat = -1.0;
    out_cpt = -1;
    return;
  }

  if (static_cast<int>(state->data_old.size()) >= p + 2) {
    const int current_i = state->current_i;
    state->x_buf.erase(state->x_buf.begin());
    state->x_buf.push_back(state->data_old[static_cast<std::size_t>(p + current_i - 2)]);

    const int buf_max_y = 2 * p;
    const double new_y = state->data_old[static_cast<std::size_t>(p + current_i - 1)] - dot_product(rho, reverse_copy(state->x_buf));
    if (static_cast<int>(state->y_buf.size()) < buf_max_y) {
      state->y_buf.push_back(new_y);
    } else {
      state->y_buf.erase(state->y_buf.begin());
      state->y_buf.push_back(new_y);
    }

    state->M_n += new_y;
    state->sum_square += new_y * new_y;

    TwoSidedState res_pos = focus_arp_two_sides_update(state->y_buf, state->current_i, state->state_right_pos, state->state_left_pos,
                                                       rho, state->M_n, state->sum_square, known_prechange, true, state->data_old.size() - p);
    state->state_right_pos = res_pos.state_right;
    state->state_left_pos = res_pos.state_left;

    std::vector<double> y_neg;
    y_neg.reserve(state->y_buf.size());
    for (double v : state->y_buf) {
      y_neg.push_back(-v);
    }
    TwoSidedState res_neg = focus_arp_two_sides_update(y_neg, state->current_i, state->state_right_neg, state->state_left_neg,
                                                       rho, -state->M_n, state->sum_square, known_prechange, true, state->data_old.size() - p);
    state->state_right_neg = res_neg.state_right;
    state->state_left_neg = res_neg.state_left;

    std::vector<double> max_vals = {
        state->state_right_pos.max_val,
        state->state_left_pos.max_val,
        state->state_right_neg.max_val,
        state->state_left_neg.max_val
    };
    double best_val = -std::numeric_limits<double>::infinity();
    int best_cpt = -1;
    for (std::size_t idx = 0; idx < max_vals.size(); ++idx) {
      if (max_vals[idx] > best_val) {
        best_val = max_vals[idx];
        if (idx == 0) {
          best_cpt = state->state_right_pos.cpt;
        } else if (idx == 1) {
          best_cpt = state->state_left_pos.cpt;
        } else if (idx == 2) {
          best_cpt = state->state_right_neg.cpt;
        } else {
          best_cpt = state->state_left_neg.cpt;
        }
      }
    }
    out_max_stat = best_val;
    out_cpt = best_cpt;
    ++state->current_i;
  } else {
    out_max_stat = -1.0;
    out_cpt = -1;
  }
}

void cleanup_arp_states(void* opaque_states) {
  delete reinterpret_cast<ArpState*>(opaque_states);
}

}  // namespace changepoint
