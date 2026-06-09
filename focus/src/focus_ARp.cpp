// ARP (AutoRegressive Process) implementation for FOCuS C++ backend
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <limits>

// ---------- helpers ----------
inline double dot_vec(const std::vector<double>& a, const std::vector<double>& b){
  double s = 0.0;
  size_t n = std::min(a.size(), b.size());
  for (size_t i = 0; i < n; ++i) s += a[i] * b[i];
  return s;
}

inline double sum_vec(const std::vector<double>& v){
  return std::accumulate(v.begin(), v.end(), 0.0);
}

// a function to build y_tau vector
inline std::vector<double> build_y_tau(const std::vector<double>& x_tau, const std::vector<double>& rho, int p) {
  std::vector<double> y_tau(p, 0.0);
  for (int j = 1; j <= p; ++j){
    int start_idx = p + j - 2;
    int end_idx = j - 1;
    std::vector<double> tmp;
    tmp.reserve(p);
    for (int k = start_idx; k >= end_idx; --k) {
        tmp.push_back(x_tau[k]);
    }
    y_tau[j - 1] = x_tau[p + j - 1] - dot_vec(rho, tmp);
  }
  return y_tau;
}

// ---------- Result struct to replace Rcpp::List ----------
struct MaxValResult {
  int cpt;
  double opt_max_val;
};

// ---------- Triple struct ----------
struct Triple {
  int tau;
  std::vector<double> S_tau;
  double l;
  double A, B, C, D, E, f;
  Triple(): tau(0), l(std::nan("")), A(std::nan("")), B(std::nan("")), C(std::nan("")), D(std::nan("")), E(std::nan("")), f(std::nan("")) {}
  Triple(int tau_, const std::vector<double>& S_tau_, double l_ = std::nan(""))
    : tau(tau_), S_tau(S_tau_), l(l_), A(std::nan("")), B(std::nan("")), C(std::nan("")), D(std::nan("")), E(std::nan("")), f(std::nan("")) {}
};

// ---------- return type for Q_n_mu_arp_unified ----------
struct QnResult {
  double              S_n;
  std::vector<Triple> triples_out;  // the p-1 "out" triples produced this step
  Triple              q_new;        // the newly introduced triple (needed by max_val_compute_arp_out)
};

// ---------- forward declarations ----------
double intersec_point_newcurve(const std::vector<double>& S_new_padded, const std::vector<double>& rho, int current_n);
double intersec_point(const Triple& input_triple, const std::vector<double>& S_new, const std::vector<double>& rho, int current_n);
Triple coef_introduce(const Triple& q_new, const std::vector<double>& y_start, const std::vector<double>& rho, double sum_square);
Triple coef_introduce_pre0(const Triple& q_new, const std::vector<double>& rho);
std::vector<Triple> coef_introduce_outre(const Triple& q_new, const std::vector<double>& y_start, const std::vector<double>& rho, double sum_square);
std::vector<Triple> coef_introduce_outre_pre0(const Triple& q_new, const std::vector<double>& rho);
QnResult Q_n_mu_arp_unified(std::vector<Triple>& triples,
                               const std::vector<double>& buf,
                               int buf_start,
                               double buf_sum_offset,
                               double xn,
                               int current_n,
                               double S_n_1,
                               const std::vector<double>& rho,
                               bool known_prechange,
                               double sum_square,
                               const std::vector<double>& y_start,
                               bool right_side);
std::vector<Triple> coef_update_arp(std::vector<Triple> new_triples, const std::vector<double>& rho, double y_new, bool known_prechange);
MaxValResult max_val_compute_arp(const std::vector<Triple>& triples, int current_n, const std::vector<double>& rho, double sum_square, const std::vector<double>& y_start);
MaxValResult max_val_compute_arp_out(const std::vector<Triple>& triples_out, int current_n, const std::vector<double>& rho, double sum_square, const std::vector<double>& y_start, const Triple& q_new_pre);
MaxValResult max_val_compute_pre0(const std::vector<Triple>& triples);
std::vector<Triple> exact_form_out_start(const std::vector<double>& rho, const std::vector<double>& y_current);
std::vector<Triple> exact_form_out_start_pre0(const std::vector<double>& rho, const std::vector<double>& y_current);
MaxValResult max_val_compute_arp_out_start(const std::vector<Triple>& triples, int current_n, const std::vector<double>& rho, const std::vector<double>& y_current);

struct State; // forward-declare struct if defined later

State init_state(double first_value);

void focus_arp_one_iter_cpp(double x_new,
                            int i,
                            State& st,
                            const std::vector<double>& rho,
                            int p,
                            int n,
                            int buf_max,
                            bool known_prechange,
                            bool right_side);

// ---------- implementations ----------
double intersec_point_newcurve(const std::vector<double>& S_new,
                               const std::vector<double>& rho,
                               int current_n)
{
  int p = (int)rho.size();

  // build x_new (length 2*p)
  std::vector<double> x_new(2*p);
  for (int k = 0; k < 2*p; ++k) x_new[k] = S_new[k+1] - S_new[k];

  // build rho_star (length p)
  std::vector<double> rho_star;
  rho_star.reserve(p);
  rho_star.push_back(1.0);
  for (int i = 1; i <= p-1; ++i){
    double s = 0.0;
    for (int j = 0; j <= i-1; ++j) s += rho[j];
    rho_star.push_back(1.0 - s);
  }

  // x_new_vec = x_new[(p)..(2*p-1)]
  std::vector<double> x_new_vec;
  x_new_vec.reserve(p);
  for (int idx = p; idx <= 2*p - 1; ++idx) x_new_vec.push_back(x_new[idx]);

  // compute m %*% rho  (fix: use -1 for 0-based indices)
  std::vector<double> m_times_rho(p, 0.0);
  for (int i = 0; i < p; ++i){
    double accum = 0.0;
    for (int j = 0; j < p; ++j){
      int idx = p + i - j - 1;            // <-- FIX: -1 for 0-based
      if (idx >= 0 && idx < 2*p) accum += x_new[idx] * rho[j];
    }
    m_times_rho[i] = accum;
  }

  // tmp = x_new_vec - m_times_rho
  std::vector<double> tmp(p, 0.0);
  for (int i = 0; i < p; ++i) tmp[i] = x_new_vec[i] - m_times_rho[i];

  // numerator = rho_star' * tmp
  double numer = 0.0;
  for (int i = 0; i < p; ++i) numer += rho_star[i] * tmp[i];

  // denom = rho_star' * rho_star
  double denom = 0.0;
  for (int i = 0; i < p; ++i) denom += rho_star[i] * rho_star[i];

  double mu_inter = 0.0;
  if (denom != 0.0) mu_inter = 2.0 * numer / denom;

  return mu_inter;
}



double intersec_point(const Triple& input_triple,
                      const std::vector<double>& S_new,
                      const std::vector<double>& rho,
                      int current_n) {
    int p = (int)rho.size();

    std::vector<double> x_tau(2*p);
    for (int k = 0; k < 2*p; ++k)
        x_tau[k] = input_triple.S_tau[k+1] - input_triple.S_tau[k];

    double S_tau_p = input_triple.S_tau[2*p];

    std::vector<double> S_tau_vec;
    S_tau_vec.reserve(p);
    for (int idx = 2*p - 1; idx >= p; --idx)
        S_tau_vec.push_back(input_triple.S_tau[idx]); // R order: descending

    std::vector<double> x_new(2*p);
    for (int k = 0; k < 2*p; ++k)
        x_new[k] = S_new[k+1] - S_new[k];

    double S_n = S_new[2*p];

    std::vector<double> S_new_vec;
    S_new_vec.reserve(p);
    for (int idx = 2*p - 1; idx >= p; --idx)
        S_new_vec.push_back(S_new[idx]); // match R's descending slice

    std::vector<double> rho_star;
    rho_star.reserve(p);
    rho_star.push_back(1.0);
    for (int i = 1; i <= p-1; ++i) {
        double s = 0.0;
        for (int j = 0; j <= i-1; ++j) s += rho[j];
        rho_star.push_back(1.0 - s);
    }

    std::vector<double> x_tau_vec;
    x_tau_vec.reserve(p);
    for (int idx = p; idx <= 2*p - 1; ++idx)
        x_tau_vec.push_back(x_tau[idx]);

    std::vector<double> x_new_vec;
    x_new_vec.reserve(p);
    for (int idx = p; idx <= 2*p - 1; ++idx)
        x_new_vec.push_back(x_new[idx]);

    std::vector<double> m_times_rho(p, 0.0);
    for (int i = 0; i < p; ++i) {
        double accum = 0.0;
        for (int j = 0; j < p; ++j) {
            int idx = p + i - j - 1; // FIXED: minus 1 to match R indexing
            double val = 0.0;
            if (idx >= 0 && idx < 2*p) val += x_tau[idx];
            if (idx >= 0 && idx < 2*p) val -= x_new[idx];
            accum += val * rho[j];
        }
        m_times_rho[i] = accum;
    }

    std::vector<double> tmp(p, 0.0);
    for (int i = 0; i < p; ++i)
        tmp[i] = x_tau_vec[i] - x_new_vec[i] - m_times_rho[i];

    double part1 = 0.0;
    for (int i = 0; i < p; ++i)
        part1 += rho_star[i] * tmp[i];

    double sum_rho = sum_vec(rho);
    double inner = 0.0;
    for (int i = 0; i < p; ++i)
        inner += rho[i] * (S_new_vec[i] - S_tau_vec[i]);

    double numer = 2.0 * (part1 + (1.0 - sum_rho) * (S_n - S_tau_p - inner));
    double denom = (1.0 - sum_rho) * (1.0 - sum_rho) *
                   ( (double)(current_n - p - input_triple.tau) );

    double mu_inter = 0.0;
    if (denom != 0.0) mu_inter = numer / denom;
    if (mu_inter < 0.0) mu_inter = 0.0;

    return mu_inter;
}


Triple coef_introduce(const Triple& q_new, const std::vector<double>& y_start, const std::vector<double>& rho, double sum_square){
  int p = (int)rho.size();
  Triple q = q_new;

  std::vector<double> x_tau(2*p);
  for (int k = 0; k < 2*p; ++k) x_tau[k] = q_new.S_tau[k+1] - q_new.S_tau[k];

  // compute y_tau vector using the helper function
  std::vector<double> y_tau = build_y_tau(x_tau, rho, p);

  // // print y_tau for debugging
  // Rcout << "y_tau: ";
  // for (auto v : y_tau) Rcout << v << " ";
  // Rcout << "\n";

  if (q_new.tau >= p){
    std::vector<double> u; u.reserve(p);
    u.push_back(1.0);
    if (p >= 2){
      for (int i = 1; i <= p-1; ++i){
        double s = 0.0;
        for (int j = 0; j <= i-1; ++j) s += rho[j];
        u.push_back(1.0 - s);
      }
    }
    std::vector<double> v; v.reserve(p);
    for (int i = 1; i <= p; ++i){
      double s = 0.0;
      for (int j = i-1; j <= p-1; ++j) s += rho[j];
      v.push_back(s);
    }
    std::vector<double> S_tau_mid(p+1);
    for (int i = 0; i <= p; ++i) S_tau_mid[i] = q_new.S_tau[i];
    double last_S = S_tau_mid[p];
    std::vector<double> Ssub; Ssub.reserve(p);
    for (int idx = p-1; idx >= 0; --idx) Ssub.push_back(S_tau_mid[idx]);
    while ((int)Ssub.size() < p) Ssub.push_back(0.0);
    double M_tau = last_S - dot_vec(rho, Ssub);
    double M_p = sum_vec(y_start);

    double u_dot_u = dot_vec(u,u);
    double v_dot_v = dot_vec(v,v);
    double u_dot_v = dot_vec(u,v);
    double sum_rho = sum_vec(rho);

    double A = -0.5 * ( u_dot_u + ( (double)(q_new.tau - p) ) * (1.0 - sum_rho) * (1.0 - sum_rho) + v_dot_v );
    double B = -0.5 * u_dot_u;
    double C = u_dot_v;
    double D = -0.5 * ( -2.0 * dot_vec(y_start, u) - 2.0 * (1.0 - sum_rho) * (M_tau - M_p) + 2.0 * dot_vec(v, y_tau) );
    double E = -0.5 * ( -2.0 * dot_vec(y_tau, u) );
    double f = -0.5 * sum_square;

    q.A = A; q.B = B; q.C = C; q.D = D; q.E = E; q.f = f;
    return q;
  } else {
    int tau = q_new.tau;
    std::vector<double> y_start_trunc(y_start.begin(), y_start.begin() + tau);

    std::vector<double> u1;
    u1.push_back(1.0);
    if ((tau - 1) != 0){
      for (int i = 1; i <= tau - 1; ++i){
        double s = 0.0;
        for (int j = 0; j <= i-1; ++j) s += rho[j];
        u1.push_back(1.0 - s);
      }
    }

    std::vector<double> u2_1;
    for (int i = 1; i <= (p - tau + 1); ++i){
      double s = 0.0;
      for (int j = i-1; j <= i-1 + (tau - 1); ++j) s += rho[j];
      u2_1.push_back(s);
    }
    std::vector<double> u2_2;
    if ((p - tau + 1) != p){
      for (int i = p - tau + 2; i <= p; ++i){
        double s = 0.0;
        for (int j = i-1; j <= p - 1; ++j) s += rho[j];
        u2_2.push_back(s);
      }
    }
    std::vector<double> u2 = u2_1;
    u2.insert(u2.end(), u2_2.begin(), u2_2.end());

    // // print u1 and u2 for debugging
    // if (q_new.tau <= 7) {
    //   Rcout << "u1: ";
    //   for (auto v : u1) Rcout << v << " ";
    //   Rcout << "\nu2: ";
    //   for (auto v : u2) Rcout << v << " ";
    //   Rcout << "\n";
    // }

    std::vector<double> v;
    v.push_back(1.0);
    if (p >= 2) {
      for (int i = 1; i <= p-1; ++i){
        double s = 0.0;
        for (int j = 0; j <= i-1; ++j) s += rho[j];
        v.push_back(1.0 - s);
      }
    }

    double A = -0.5 * ( dot_vec(u1, u1) + dot_vec(u2, u2) );
    double B = -0.5 * dot_vec(v, v);
    double C = dot_vec(u2, v);
    double D = -0.5 * ( -2.0 * dot_vec(y_start_trunc, u1) + 2.0 * dot_vec(y_tau, u2) );
    double E = -0.5 * ( -2.0 * dot_vec(y_tau, v) );
    double f = -0.5 * sum_square;

    q.A = A; q.B = B; q.C = C; q.D = D; q.E = E; q.f = f;
    return q;
  }
}

Triple coef_introduce_pre0(const Triple& q_new, const std::vector<double>& rho){
  int p = (int)rho.size();
  std::vector<double> u; u.reserve(p);
  u.push_back(1.0);
  if (p >= 2){
    for (int i = 1; i <= p-1; ++i){
      double s = 0.0;
      for (int j = 0; j <= i-1; ++j) s += rho[j];
      u.push_back(1.0 - s);
    }
  }

  std::vector<double> x_tau(2*p);
  for (int k = 0; k < 2*p; ++k) x_tau[k] = q_new.S_tau[k+1] - q_new.S_tau[k];

  std::vector<double> y_tau = build_y_tau(x_tau, rho, p);

  Triple q = q_new;
  q.A = -0.5 * dot_vec(u,u);
  q.B = dot_vec(y_tau, u);
  return q;
}

// ---------- coef_introduce_outre (unknown pre-change) ----------
// For q_new (tau >= 1), produces p-1 "out" triples for tau_out in (tau+1)..(tau+p-1).
// Mirrors R: coef_introduce_outre
std::vector<Triple> coef_introduce_outre(const Triple& q_new,
                                          const std::vector<double>& y_start,
                                          const std::vector<double>& rho,
                                          double sum_square)
{
  const int p   = (int)rho.size();
  const int tau = q_new.tau;

  // x_tau = differences of S_tau (length 2*p)
  std::vector<double> x_tau(2 * p);
  for (int k = 0; k < 2 * p; ++k) x_tau[k] = q_new.S_tau[k + 1] - q_new.S_tau[k];

  // y_tau[j] = x_{tau+j} - rho' * x_{tau+j-1..tau+j-p}   for j = 1..p
  std::vector<double> y_tau = build_y_tau(x_tau, rho, p);

  // Build u (length p): u[0]=1, u[i] = 1 - sum(rho[0..i-1])
  std::vector<double> u; u.reserve(p);
  u.push_back(1.0);
  for (int i = 1; i <= p - 1; ++i) {
    double s = 0.0;
    for (int j = 0; j < i; ++j) s += rho[j];
    u.push_back(1.0 - s);
  }

  // Build v (length p): v[i] = sum(rho[i..p-1])
  std::vector<double> v_full; v_full.reserve(p);
  for (int i = 0; i < p; ++i) {
    double s = 0.0;
    for (int j = i; j < p; ++j) s += rho[j];
    v_full.push_back(s);
  }

  double sum_rho = sum_vec(rho);

  std::vector<Triple> new_triples_out;
  new_triples_out.reserve(p - 1);

  for (int tau_out = tau + 1; tau_out <= tau + p - 1; ++tau_out) {
    const int index_out = tau_out - tau;  // 1 .. p-1

    // truncated y_tau: y_tau[index_out .. p-1]  (0-based: indices index_out .. p-1)
    std::vector<double> y_tau_new(y_tau.begin() + index_out, y_tau.end());
    const int len_new = (int)y_tau_new.size();  // = p - index_out

    double A, B, C, D, E, f_val;

    if (tau_out >= p) {
      // ---- tau >= p branch ----
      // u_new = u[0..p-index_out-1], v_new = v_full[0..p-index_out-1]
      std::vector<double> u_new(u.begin(), u.begin() + len_new);
      std::vector<double> v_new(v_full.begin(), v_full.begin() + len_new);

      // M_tau from S_tau[0..p]: S_tau_mid[p] - rho' * (S_tau_mid[p-1..0])
      std::vector<double> S_tau_mid(q_new.S_tau.begin(), q_new.S_tau.begin() + p + 1);
      std::vector<double> Ssub; Ssub.reserve(p);
      for (int idx = p - 1; idx >= 0; --idx) Ssub.push_back(S_tau_mid[idx]);
      double M_tau = S_tau_mid[p] - dot_vec(rho, Ssub);
      double M_p   = sum_vec(y_start);

      A = -0.5 * (dot_vec(u, u) + (double)(tau - p) * (1.0 - sum_rho) * (1.0 - sum_rho) + dot_vec(v_new, v_new));
      B = -0.5 * dot_vec(u_new, u_new);
      C =  dot_vec(u_new, v_new);
      D = -0.5 * (-2.0 * dot_vec(y_start, u)
                  - 2.0 * (1.0 - sum_rho) * (M_tau - M_p)
                  + 2.0 * dot_vec(v_new, y_tau_new));
      E = -0.5 * (-2.0 * dot_vec(y_tau_new, u_new));
      f_val = -0.5 * sum_square;
    } else {
      // ---- tau < p branch ----
      // y_start2 = y_start[0..tau-1]
      std::vector<double> y_start2(y_start.begin(), y_start.begin() + tau);

      // u1: length tau
      std::vector<double> u1; u1.reserve(tau);
      u1.push_back(1.0);
      for (int i = 1; i <= tau - 1; ++i) {
        double s = 0.0;
        for (int j = 0; j < i; ++j) s += rho[j];
        u1.push_back(1.0 - s);
      }

      // u2: length p (u2_1 then u2_2)
      std::vector<double> u2;
      for (int i = 1; i <= p - tau + 1; ++i) {
        double s = 0.0;
        for (int j = i - 1; j <= i - 1 + (tau - 1); ++j) s += rho[j];
        u2.push_back(s);
      }
      if ((p - tau + 1) != p) {
        for (int i = p - tau + 2; i <= p; ++i) {
          double s = 0.0;
          for (int j = i - 1; j < p; ++j) s += rho[j];
          u2.push_back(s);
        }
      }

      // v (same as u above, length p): already built as u; reuse
      // NOTE: in the tau<p branch of coef_introduce, v = u (same construction)
      // here we call it v to match R naming
      std::vector<double>& v_ref = u;  // alias: same vector

      // truncated u2 and v: first (p - index_out) elements
      std::vector<double> u2_new(u2.begin(), u2.begin() + len_new);
      std::vector<double> v_new(v_ref.begin(), v_ref.begin() + len_new);

      A = -0.5 * (dot_vec(u1, u1) + dot_vec(u2_new, u2_new));
      B = -0.5 * dot_vec(v_new, v_new);
      C =  dot_vec(u2_new, v_new);
      D = -0.5 * (-2.0 * dot_vec(y_start2, u1) + 2.0 * dot_vec(y_tau_new, u2_new));
      E = -0.5 * (-2.0 * dot_vec(y_tau_new, v_new));
      f_val = -0.5 * sum_square;
    }

    Triple t;
    t.tau = tau_out;
    // S_tau is not needed for out-triples (set empty)
    t.A = A; t.B = B; t.C = C; t.D = D; t.E = E; t.f = f_val;
    new_triples_out.push_back(t);
  }

  return new_triples_out;
}

// ---------- coef_introduce_outre_pre0 (known pre-change) ----------
// Mirrors R: coef_introduce_outre_pre0
std::vector<Triple> coef_introduce_outre_pre0(const Triple& q_new,
                                               const std::vector<double>& rho)
{
  const int p   = (int)rho.size();
  const int tau = q_new.tau;

  // Build u (length p)
  std::vector<double> u; u.reserve(p);
  u.push_back(1.0);
  for (int i = 1; i <= p - 1; ++i) {
    double s = 0.0;
    for (int j = 0; j < i; ++j) s += rho[j];
    u.push_back(1.0 - s);
  }

  // x_tau and y_tau
  std::vector<double> x_tau(2 * p);
  for (int k = 0; k < 2 * p; ++k) x_tau[k] = q_new.S_tau[k + 1] - q_new.S_tau[k];
  std::vector<double> y_tau = build_y_tau(x_tau, rho, p);

  std::vector<Triple> new_triples_out;
  new_triples_out.reserve(p - 1);

  for (int tau_out = tau + 1; tau_out <= tau + p - 1; ++tau_out) {
    const int index_out = tau_out - tau;      // 1 .. p-1
    const int len_new   = p - index_out;      // p - index_out

    std::vector<double> u_new(u.begin(), u.begin() + len_new);
    std::vector<double> y_tau_new(y_tau.begin() + index_out, y_tau.end());

    Triple t;
    t.tau = tau_out;
    t.A   = -0.5 * dot_vec(u_new, u_new);
    t.B   =  dot_vec(y_tau_new, u_new);
    new_triples_out.push_back(t);
  }

  return new_triples_out;
}

// ---------- exact_form_out_start (unknown pre-change, n <= p+1) ----------
// Mirrors R: exact_form_out_start
std::vector<Triple> exact_form_out_start(const std::vector<double>& rho,
                                          const std::vector<double>& y_current)
{
  const int p          = (int)rho.size();
  const int current_n  = (int)y_current.size();

  // vec_all[i] = 1 - sum(rho[0..i-1]),  vec_all[0] = 1   (length p)
  std::vector<double> vec_all; vec_all.reserve(p);
  vec_all.push_back(1.0);
  for (int i = 1; i <= p - 1; ++i) {
    double s = 0.0;
    for (int j = 0; j < i; ++j) s += rho[j];
    vec_all.push_back(1.0 - s);
  }

  double sum_square = dot_vec(y_current, y_current);

  std::vector<Triple> new_triples_out;
  new_triples_out.reserve(current_n - 1);

  for (int tau_out = 1; tau_out <= current_n - 1; ++tau_out) {
    const int len_u1   = tau_out;
    const int len_v    = current_n - tau_out;
    const int len_u2   = current_n - tau_out;   // same as len_v

    // u1 = vec_all[0..tau_out-1]
    std::vector<double> u1(vec_all.begin(), vec_all.begin() + len_u1);

    // v = vec_all[0..len_v-1]
    std::vector<double> v_vec(vec_all.begin(), vec_all.begin() + len_v);

    // u2[i] = sum(rho[i..i+tau_out-1])  for i = 0..len_u2-1   (0-based)
    std::vector<double> u2; u2.reserve(len_u2);
    for (int i = 0; i < len_u2; ++i) {
      double s = 0.0;
      for (int j = i; j <= i + tau_out - 1; ++j) s += rho[j];
      u2.push_back(s);
    }

    // y_start = y_current[0..tau_out-1], y_tau = y_current[tau_out..current_n-1]
    std::vector<double> y_start_loc(y_current.begin(), y_current.begin() + tau_out);
    std::vector<double> y_tau_loc(y_current.begin() + tau_out, y_current.end());

    double A     = -0.5 * (dot_vec(u1, u1) + dot_vec(u2, u2));
    double B     = -0.5 * dot_vec(v_vec, v_vec);
    double C     =  dot_vec(u2, v_vec);
    double D     = -0.5 * (-2.0 * dot_vec(y_start_loc, u1) + 2.0 * dot_vec(y_tau_loc, u2));
    double E     = -0.5 * (-2.0 * dot_vec(y_tau_loc, v_vec));
    double f_val = -0.5 * sum_square;

    Triple t;
    t.tau = tau_out;
    t.A = A; t.B = B; t.C = C; t.D = D; t.E = E; t.f = f_val;
    new_triples_out.push_back(t);
  }

  return new_triples_out;
}

// ---------- exact_form_out_start_pre0 (known pre-change, n <= p+1) ----------
// Mirrors R: exact_form_out_start_pre0
std::vector<Triple> exact_form_out_start_pre0(const std::vector<double>& rho,
                                               const std::vector<double>& y_current)
{
  const int p         = (int)rho.size();
  const int current_n = (int)y_current.size();

  // vec_all[i] = 1 - sum(rho[0..i-1]),  vec_all[0] = 1   (length p)
  std::vector<double> vec_all; vec_all.reserve(p);
  vec_all.push_back(1.0);
  for (int i = 1; i <= p - 1; ++i) {
    double s = 0.0;
    for (int j = 0; j < i; ++j) s += rho[j];
    vec_all.push_back(1.0 - s);
  }

  std::vector<Triple> new_triples_out;
  new_triples_out.reserve(current_n - 1);

  for (int tau_out = 1; tau_out <= current_n - 1; ++tau_out) {
    const int len_v = current_n - tau_out;

    // v = vec_all[0..len_v-1]
    std::vector<double> v_vec(vec_all.begin(), vec_all.begin() + len_v);

    // y_tau = y_current[tau_out..current_n-1]
    std::vector<double> y_tau_loc(y_current.begin() + tau_out, y_current.end());

    Triple t;
    t.tau = tau_out;
    t.A   = -0.5 * dot_vec(v_vec, v_vec);
    t.B   =  dot_vec(v_vec, y_tau_loc);   // R: B <- -0.5 * (-2 * (v %*% y_tau))
    new_triples_out.push_back(t);
  }

  return new_triples_out;
}

// ---------- max_val_compute_arp_out_start (unknown pre-change, n <= p+1) ----------
// Mirrors R: max_val_compute_arp_out_start
MaxValResult max_val_compute_arp_out_start(const std::vector<Triple>& triples,
                                            int current_n,
                                            const std::vector<double>& rho,
                                            const std::vector<double>& y_current)
{
  const int p = (int)rho.size();

  // u_whole[i] = 1 - sum(rho[0..i-1]),  u_whole[0]=1  (length p+1 in R, we use current_n)
  // R: u_whole <- c(1); for i in 1:p u_whole <- c(u_whole, 1-sum(rho[1:i]))
  // then u <- u_whole[1:current_n]
  std::vector<double> u_whole; u_whole.reserve(p + 1);
  u_whole.push_back(1.0);
  for (int i = 1; i <= p; ++i) {
    double s = 0.0;
    for (int j = 0; j < i; ++j) s += rho[j];
    u_whole.push_back(1.0 - s);
  }
  std::vector<double> u(u_whole.begin(), u_whole.begin() + current_n);

  double a_no    = -0.5 * dot_vec(u, u);
  double b_no    =  dot_vec(y_current, u);
  double c_no    = -0.5 * dot_vec(y_current, y_current);
  double opt_no  = (a_no != 0.0) ? -(b_no / (2.0 * a_no)) : 0.0;
  double max_val_no = a_no * opt_no * opt_no + b_no * opt_no + c_no;

  const int K = (int)triples.size();
  std::vector<double> max_vals(K, -std::numeric_limits<double>::infinity());
  std::vector<int>    cpts(K, -1);

  for (int k = 0; k < K; ++k) {
    const Triple& tr = triples[k];
    double A = tr.A, B = tr.B, C = tr.C, D = tr.D, E = tr.E, f = tr.f;
    double denom = 4.0 * A * B - C * C;
    double max_val_change = -std::numeric_limits<double>::infinity();
    if (denom != 0.0) {
      double inv  = 1.0 / denom;
      double mu0  = -(2.0 * B * inv * D + (-C * inv) * E);
      double mu1  = -( (-C * inv) * D + 2.0 * A * inv * E);
      max_val_change = A*mu0*mu0 + B*mu1*mu1 + C*mu0*mu1 + D*mu0 + E*mu1 + f;
    }
    max_vals[k] = 2.0 * (max_val_change - max_val_no);
    cpts[k]     = tr.tau;
  }

  int    opt_index = 0;
  double best      = max_vals[0];
  for (int k = 1; k < K; ++k) {
    if (max_vals[k] > best) { best = max_vals[k]; opt_index = k; }
  }
  return MaxValResult{cpts[opt_index], best};
}

// ---------- max_val_compute_arp_out (unknown pre-change, recursion phase) ----------
// Mirrors R: max_val_compute_arp_out
// q_new_pre is the last triple introduced (holds S_tau needed for M_n).
MaxValResult max_val_compute_arp_out(const std::vector<Triple>& triples_out,
                                      int current_n,
                                      const std::vector<double>& rho,
                                      double sum_square,
                                      const std::vector<double>& y_start,
                                      const Triple& q_new_pre)
{
  const int p = (int)rho.size();

  // M_n from q_new_pre.S_tau
  const std::vector<double>& S_new = q_new_pre.S_tau;
  std::vector<double> Ssub; Ssub.reserve(p);
  for (int idx = 2 * p - 1; idx >= p; --idx) Ssub.push_back(S_new[idx]);
  double M_n = S_new[2 * p] - dot_vec(rho, Ssub);

  // u (length p)
  std::vector<double> u; u.reserve(p);
  u.push_back(1.0);
  for (int i = 1; i <= p - 1; ++i) {
    double s = 0.0;
    for (int j = 0; j < i; ++j) s += rho[j];
    u.push_back(1.0 - s);
  }
  double M_p     = sum_vec(y_start);
  double sum_rho = sum_vec(rho);

  double a_no   = -0.5 * ((double)(current_n - p) * (1.0 - sum_rho) * (1.0 - sum_rho) + dot_vec(u, u));
  double b_no   = (1.0 - sum_rho) * (M_n - M_p) + dot_vec(y_start, u);
  double c_no   = -0.5 * sum_square;
  double opt_no = (a_no != 0.0) ? -(b_no / (2.0 * a_no)) : 0.0;
  double max_val_no = a_no * opt_no * opt_no + b_no * opt_no + c_no;

  const int K = (int)triples_out.size();
  std::vector<double> max_vals(K, -std::numeric_limits<double>::infinity());
  std::vector<int>    cpts(K, -1);

  for (int k = 0; k < K; ++k) {
    const Triple& tr = triples_out[k];
    double A = tr.A, B = tr.B, C = tr.C, D = tr.D, E = tr.E, f = tr.f;
    double denom = 4.0 * A * B - C * C;
    double max_val_change = -std::numeric_limits<double>::infinity();
    if (denom != 0.0) {
      double inv  = 1.0 / denom;
      double mu0  = -(2.0 * B * inv * D + (-C * inv) * E);
      double mu1  = -( (-C * inv) * D + 2.0 * A * inv * E);
      max_val_change = A*mu0*mu0 + B*mu1*mu1 + C*mu0*mu1 + D*mu0 + E*mu1 + f;
    }
    max_vals[k] = 2.0 * (max_val_change - max_val_no);
    cpts[k]     = tr.tau;
  }

  int    opt_index = 0;
  double best      = max_vals[0];
  for (int k = 1; k < K; ++k) {
    if (max_vals[k] > best) { best = max_vals[k]; opt_index = k; }
  }
  return MaxValResult{cpts[opt_index], best};
}

// NEW overload with side-control (right=true, left=false). Keeps existing API intact.
QnResult Q_n_mu_arp_unified(std::vector<Triple>& triples,
                               const std::vector<double>& buf,
                               int buf_start,
                               double buf_sum_offset,
                               double xn,
                               int current_n,
                               double S_n_1,
                               const std::vector<double>& rho,
                               bool known_prechange,
                               double sum_square,
                               const std::vector<double>& y_start,
                               bool right_side) {
  const int p = (int)rho.size();
  const int len_buf = (int)buf.size();

  // --- Construct backshift (last up to 2*p-1 previous observations) ---
  std::vector<double> backshift;
  if (len_buf <= (2 * p - 1)) {
    backshift.insert(backshift.end(), 2 * p - len_buf, 0.0);
    // first len_buf-1 elements
    for (int j = 0; j < len_buf - 1; ++j) backshift.push_back(buf[j]);
  } else {
    // buf[(len_buf+1-2*p):(len_buf-1)] in R -> [len_buf-2*p, len_buf-2] in C++
    for (int j = len_buf - 2 * p; j <= len_buf - 2; ++j) backshift.push_back(buf[j]);
  }

  // --- S_n and S_new construction ---
  const double S_n = S_n_1 + xn;

  // prepare S_n_minus: (length 2*p)
  std::vector<double> S_n_minus;
  S_n_minus.reserve(2 * p);
  for (int i = 1; i <= (2 * p - 1); ++i) {
    double sum_part = 0.0;
    for (int j = i - 1; j <= 2 * p - 2; ++j) {
      if (j >= 0 && j < (int)backshift.size()) sum_part += backshift[j];
    }
    S_n_minus.push_back(xn + sum_part);
  }
  S_n_minus.push_back(xn);

  // S_n_back = S_n - S_n_minus
  std::vector<double> S_n_back;
  S_n_back.reserve((int)S_n_minus.size());
  for (double v : S_n_minus) S_n_back.push_back(S_n - v);

  // S_new = c(S_n_back, S_n)
  std::vector<double> S_new = S_n_back;
  S_new.push_back(S_n);

  // --- Update l for newest triple and prune with side-specific rule ---
  int i = (int)triples.size();
  if (!triples.empty()) {
    triples[i - 1].l = intersec_point(triples[i - 1], S_new, rho, current_n);

    // side-dependent l_0 and monotonicity condition
    const double l_0 = right_side ? -std::numeric_limits<double>::infinity()
                                  :  std::numeric_limits<double>::infinity();

    int ind = 0;
    while (ind == 0 && i >= 1) {
      const double l_i   = triples[i - 1].l;
      const double l_i_1 = (i - 2 < 0) ? l_0 : triples[i - 2].l;

      const bool cond = right_side ? (l_i <= l_i_1) : (l_i >= l_i_1);

      if (cond) {
        i = i - 1;
        if (i > 0) triples[i - 1].l = intersec_point(triples[i - 1], S_new, rho, current_n);
      } else {
        ind = 1;
      }
    }
  } else {
    i = 0;
  }

  // keep only first i triples
  std::vector<Triple> new_triples;
  if (i < (int)triples.size()) {
    if (i > 0) new_triples.assign(triples.begin(), triples.begin() + i);
  } else {
    new_triples = triples;
  }

  // --- Compute S_tau_new for new intro curve (length 2*p+1) ---
  std::vector<double> S_tau_new;
  S_tau_new.reserve(2 * p + 1);
  for (int j = 2 * p; j >= 0; --j) {
    const int k = current_n - j;
    double s_k = 0.0;
    if (k > 0) {
      if (k <= (buf_start - 1)) {
        s_k = buf_sum_offset;
      } else {
        const int idx_in_buf = k - buf_start + 1; // 1-based count of items inside buf included
        s_k = buf_sum_offset;
        // this is the r equivalent of s_k <- buf_sum_offset + sum(buf[1:idx_in_buf]). Could use accumulate?
        for (int t = 0; t < idx_in_buf && t < (int)buf.size(); ++t) s_k += buf[t];
      }
    } else {
      s_k = 0.0;
    }
    S_tau_new.push_back(s_k);
  }

  // --- Intersection for new curve: pad S_new with p zeros to length 2p+1 ---
  std::vector<double> S_new_padded;
  S_new_padded.reserve(2 * p + 1);
  for (int z = 0; z < p; ++z) S_new_padded.push_back(0.0);
  for (double v : S_new) S_new_padded.push_back(v);

  const double new_l = intersec_point_newcurve(S_new_padded, rho, current_n);

  Triple q_new(current_n - p, S_tau_new, new_l);

  std::vector<Triple> new_triples_out;
  if (!known_prechange) {
    new_triples_out = coef_introduce_outre(q_new, y_start, rho, sum_square);
    q_new = coef_introduce(q_new, y_start, rho, sum_square);
  } else {
    new_triples_out = coef_introduce_outre_pre0(q_new, rho);
    q_new = coef_introduce_pre0(q_new, rho);
  }

  new_triples.push_back(q_new);
  triples.swap(new_triples);

  return QnResult{S_n, new_triples_out, q_new};
}


std::vector<Triple> coef_update_arp(std::vector<Triple> new_triples, const std::vector<double>& rho, double y_new, bool known_prechange){
  if ((int)new_triples.size() <= 1) return new_triples;
  int upto = (int)new_triples.size() - 1;
  double one_minus_sumrho = 1.0 - sum_vec(rho);

  for (int i = 0; i < upto; ++i){
    if (!known_prechange){
      new_triples[i].B = new_triples[i].B - 0.5 * (one_minus_sumrho * one_minus_sumrho);
      new_triples[i].E = new_triples[i].E + (one_minus_sumrho) * y_new;
      new_triples[i].f = new_triples[i].f - 0.5 * (y_new * y_new);
    } else {
      new_triples[i].A = new_triples[i].A - 0.5 * (one_minus_sumrho * one_minus_sumrho);
      new_triples[i].B = new_triples[i].B + (one_minus_sumrho) * y_new;
    }
  }
  return new_triples;
}

MaxValResult max_val_compute_arp(const std::vector<Triple>& triples, int current_n, const std::vector<double>& rho, double sum_square, const std::vector<double>& y_start){
  int p = (int)rho.size();
  if (triples.empty()) return MaxValResult{-1, -std::numeric_limits<double>::infinity()};

  const std::vector<double>& S_new = triples.back().S_tau;
  double S_new_last = S_new[2*p];
  std::vector<double> Ssub;
  Ssub.reserve(p);
  for (int idx = 2*p - 1; idx >= p; --idx) Ssub.push_back(S_new[idx]);
  double M_n = S_new_last - dot_vec(rho, Ssub);

  std::vector<double> u; u.reserve(p);
  u.push_back(1.0);
  if (p >= 2){
    for (int i = 1; i <= p-1; ++i){
      double s = 0.0;
      for (int j = 0; j <= i-1; ++j) s += rho[j];
      u.push_back(1.0 - s);
    }
  }
  double M_p = sum_vec(y_start);
  double a_no = -0.5 * ( (double)(current_n - p) * (1.0 - sum_vec(rho)) * (1.0 - sum_vec(rho)) + dot_vec(u,u) );
  double b_no = (1.0 - sum_vec(rho)) * (M_n - M_p) + dot_vec(y_start, u);
  double c_no = -0.5 * sum_square;

  double opt_no = 0.0;
  if (a_no != 0.0) opt_no = - (b_no / (2.0 * a_no));
  double max_val_no = a_no * opt_no * opt_no + b_no * opt_no + c_no;

  int K = (int)triples.size();
  std::vector<double> max_vals(K, -std::numeric_limits<double>::infinity());
  std::vector<int> cpts(K, -1);
  for (int k = 0; k < K; ++k){
    const Triple& tr = triples[k];
    double A = tr.A, B = tr.B, C = tr.C, D = tr.D, E = tr.E, f = tr.f;
    double denom = (4.0 * A * B - C * C);
    double max_val_change = -std::numeric_limits<double>::infinity();
    if (denom != 0.0){
      double inv = 1.0 / denom;
      double m00 = 2.0 * B * inv;
      double m01 = -C * inv;
      double m10 = -C * inv;
      double m11 = 2.0 * A * inv;
      double mu0 = - (m00 * D + m01 * E);
      double mu1 = - (m10 * D + m11 * E);
      max_val_change = A * mu0 * mu0 + B * mu1 * mu1 + C * mu0 * mu1 + D * mu0 + E * mu1 + f;
    }
    max_vals[k] = 2.0 * (max_val_change - max_val_no);
    cpts[k] = tr.tau;
  }

  int opt_index = 0;
  double best = max_vals[0];
  for (int k = 1; k < K; ++k){
    if (max_vals[k] > best){
      best = max_vals[k];
      opt_index = k;
    }
  }

  return MaxValResult{cpts[opt_index], best};
}

MaxValResult max_val_compute_pre0(const std::vector<Triple>& triples){
  int K = (int)triples.size();
  if (K == 0) return MaxValResult{-1, -std::numeric_limits<double>::infinity()};

  std::vector<double> max_vals(K, -std::numeric_limits<double>::infinity());
  std::vector<int> cpts(K, -1);
  for (int k = 0; k < K; ++k){
    const Triple& tr = triples[k];
    double a = tr.A;
    double b = tr.B;
    double opt_mu = 0.0;
    if (a != 0.0) opt_mu = - b / (2.0 * a);
    max_vals[k] = a * opt_mu * opt_mu + b * opt_mu;
    cpts[k] = tr.tau;
  }
  int opt_index = 0;
  double best = max_vals[0];
  for (int k = 1; k < K; ++k){
    if (max_vals[k] > best){
      best = max_vals[k];
      opt_index = k;
    }
  }
  return MaxValResult{cpts[opt_index], best};
}

struct State {
  std::vector<Triple> triples;
  std::vector<Triple> triples_out;  // "out" triples from last recursion step
  Triple              q_new_pre;    // last introduced triple (for max_val_compute_arp_out)
  std::vector<double> buf;
  int    buf_start = 1;
  double buf_sum_offset = 0.0;
  double S_n_1 = 0.0;

  // Only used when !known_prechange
  double sum_square = 0.0;
  std::vector<double> y_start;

  // outputs each iter
  double max_val = -1.0;
  int    cpt     = -1;
};

void focus_arp_one_iter_cpp(double x_new,
                            int i,
                            State& st,
                            const std::vector<double>& rho,
                            int p,
                            int n,
                            int buf_max,
                            bool known_prechange,
                            bool right_side) {
  // Append new observation to buffer
  st.buf.push_back(x_new);
  if ((int)st.buf.size() > buf_max) {
    st.buf_sum_offset += st.buf.front();
    st.buf.erase(st.buf.begin());
    st.buf_start += 1;
  }

  bool enter_recursion = (i >= (p + 2));

  // ---- early phase: n <= p+1 ----
  // Build y_current from the buffer (same logic as the R prototype's buf-based construction)
  if ((i >= 2) && (i <= (p + 1))) {
    std::vector<double> y_current;
    y_current.push_back(st.buf[0]);
    if (i >= 2) {
      for (int j = 2; j <= i; ++j) {
        std::vector<double> lag_vec;
        for (int k = j - 2; k >= 0; --k) lag_vec.push_back(st.buf[k]);
        while ((int)lag_vec.size() < p) lag_vec.push_back(0.0);
        y_current.push_back(st.buf[j - 1] - dot_vec(rho, lag_vec));
      }
    }

    if (!known_prechange) {
      st.triples_out = exact_form_out_start(rho, y_current);
    } else {
      st.triples_out = exact_form_out_start_pre0(rho, y_current);
    }

    MaxValResult max_info_out;
    if (!known_prechange) {
      max_info_out = max_val_compute_arp_out_start(st.triples_out, i, rho, y_current);
    } else {
      max_info_out = max_val_compute_pre0(st.triples_out);
    }
    st.max_val = max_info_out.opt_max_val;
    st.cpt     = max_info_out.cpt;
  }

  // Initial entry when i == p+2
  if ((i == (p + 2)) && (i <= n)) {
    if (!known_prechange) {
      st.y_start.clear();
      st.y_start.push_back(st.buf[0]);

      if (p >= 2) {
        for (int j = 2; j <= p; ++j) {
          std::vector<double> lag_vec;
          for (int k = j - 2; k >= 0; --k) lag_vec.push_back(st.buf[k]);
          while ((int)lag_vec.size() < p) lag_vec.push_back(0.0);
          st.y_start.push_back(st.buf[j - 1] - dot_vec(rho, lag_vec));
        }
      }

      // y_tau_p1 = (y_start, buf[p+1] - rho %*% buf[p:1])
      std::vector<double> rev_buf;
      for (int k = p - 1; k >= 0; --k) rev_buf.push_back(st.buf[k]);
      std::vector<double> take_for_rho;
      for (int t = 0; t < p; ++t)
        take_for_rho.push_back(t < (int)rev_buf.size() ? rev_buf[t] : 0.0);

      std::vector<double> y_tau_p1 = st.y_start;
      double lastterm = st.buf[p] - dot_vec(rho, take_for_rho);
      y_tau_p1.push_back(lastterm);

      st.sum_square = 0.0;
      for (double v : y_tau_p1) st.sum_square += v * v;
    }

    // Build S_tau for tau = 1
    // At i == p+2, buf contains observations 0..p+1 (p+2 elements)
    std::vector<double> S_tau; S_tau.reserve(2 * p + 1);
    for (int z = 0; z < p; ++z) S_tau.push_back(0.0);
    
    // Build partial sums from buf[0] to buf[p]
    double cumsum = 0.0;
    for (int k = 0; k <= p && k < (int)st.buf.size(); ++k) {
      cumsum += st.buf[k];
      S_tau.push_back(cumsum);
    }
    // If buffer is shorter than expected, pad with last value
    while ((int)S_tau.size() < 2 * p + 1) {
      S_tau.push_back(cumsum);
    }

    Triple triple_1(1, S_tau);
    if (!known_prechange)
      triple_1 = coef_introduce(triple_1, st.y_start, rho, st.sum_square);
    else
      triple_1 = coef_introduce_pre0(triple_1, rho);

    st.triples.clear();
    st.triples.push_back(triple_1);

    // Compute and store S_n_1 (cumulative sum from 0 to p)
    st.S_n_1 = 0.0;
    for (int t = 0; t <= p && t < (int)st.buf.size(); ++t) {
      st.S_n_1 += st.buf[t];
    }

    enter_recursion = true;
  }

  // Recursion
  if (enter_recursion && i <= n) {
    int len_buf = (int)st.buf.size();
    int len_prev = std::max(0, len_buf - 1);

    std::vector<double> lag_vec;
    if (len_prev >= p) {
      for (int idx = len_prev - p; idx <= len_prev - 1; ++idx) lag_vec.push_back(st.buf[idx]);
      std::reverse(lag_vec.begin(), lag_vec.end());
    } else if (len_prev > 0) {
      for (int idx = len_prev - 1; idx >= 0; --idx) lag_vec.push_back(st.buf[idx]);
      while ((int)lag_vec.size() < p) lag_vec.push_back(0.0);
    } else {
      lag_vec.assign(p, 0.0);
    }

    double y_next_n = x_new - dot_vec(rho, lag_vec);

    if (!known_prechange) {
      st.sum_square += y_next_n * y_next_n;
      QnResult res = Q_n_mu_arp_unified(st.triples, st.buf, st.buf_start, st.buf_sum_offset,
                                         x_new, i, st.S_n_1, rho,
                                         known_prechange, st.sum_square, st.y_start,
                                         right_side);
      st.S_n_1      = res.S_n;
      st.triples_out = res.triples_out;
      st.q_new_pre   = res.q_new;
      st.triples = coef_update_arp(st.triples, rho, y_next_n, known_prechange);
    } else {
      QnResult res = Q_n_mu_arp_unified(st.triples, st.buf, st.buf_start, st.buf_sum_offset,
                                         x_new, i, st.S_n_1, rho,
                                         known_prechange, 0.0, std::vector<double>(),
                                         right_side);
      st.S_n_1      = res.S_n;
      st.triples_out = res.triples_out;
      st.q_new_pre   = res.q_new;
      st.triples = coef_update_arp(st.triples, rho, y_next_n, known_prechange);
    }
  }

  // Evaluate statistic (max over triples and triples_out)
  if (!st.triples.empty()) {
    if (!known_prechange) {
      MaxValResult max_info     = max_val_compute_arp(st.triples, i, rho, st.sum_square, st.y_start);
      MaxValResult max_info_out = max_val_compute_arp_out(st.triples_out, i, rho, st.sum_square, st.y_start, st.q_new_pre);
      if (max_info.opt_max_val >= max_info_out.opt_max_val) {
        st.max_val = max_info.opt_max_val;
        st.cpt     = max_info.cpt;
      } else {
        st.max_val = max_info_out.opt_max_val;
        st.cpt     = max_info_out.cpt;
      }
    } else {
      MaxValResult max_info     = max_val_compute_pre0(st.triples);
      MaxValResult max_info_out = max_val_compute_pre0(st.triples_out);
      if (max_info.opt_max_val >= max_info_out.opt_max_val) {
        st.max_val = max_info.opt_max_val;
        st.cpt     = max_info.cpt;
      } else {
        st.max_val = max_info_out.opt_max_val;
        st.cpt     = max_info_out.cpt;
      }
    }
  } else {
    // triples is empty but triples_out may have been set by the early-phase block above
    // (the early-phase block already wrote max_val/cpt directly, so nothing to do here)
  }
}


// Helper: initialise state with first observation (positive or negative version)
State init_state(double first_value) {
  State s;
  s.triples.clear();
  s.buf.clear();
  s.buf.push_back(first_value);
  s.buf_start = 1;
  s.buf_sum_offset = 0.0;
  s.S_n_1 = 0.0;
  s.sum_square = 0.0;
  s.y_start.clear();
  s.max_val = -1.0;
  s.cpt     = -1;
  return s;
}

// ---------------------------------------------------------------------------
// Online/Sequential Interface for ARpInfo
// ---------------------------------------------------------------------------

// Structure to hold all four states (stored as opaque pointer in ARpInfo)
struct ARpStates {
  State right_pos;
  State left_pos;
  State right_neg;
  State left_neg;
};

namespace changepoint {

// Implementation function called from ARpInfo::update
void arp_detector_update_impl(double obs,
                               const std::vector<double>& rho,
                               int p,
                               int buf_max,
                               bool known_prechange,
                               double n,
                               void*& opaque_states,
                               double& out_max_stat,
                               int& out_cpt) {
  double i = n;  // Current iteration index
  
  // Initialize opaque_states on first call (when n == 1)
  if (!opaque_states) {
    ARpStates* arp_states = new ARpStates();
    arp_states->right_pos = init_state(obs);  // Initialize with first observation
    arp_states->left_pos  = init_state(obs);
    arp_states->right_neg = init_state(-obs);
    arp_states->left_neg  = init_state(-obs);
    
    opaque_states = static_cast<void*>(arp_states);
    // Return the initial stat (-1.0) to match offline behavior
    out_max_stat = -1.0;
    out_cpt = -1;
    return;
  }
  
  // Update and compute stats starting from n >= 2
  ARpStates* arp_states = static_cast<ARpStates*>(opaque_states);
  
  // Update all four states using the current observation
  focus_arp_one_iter_cpp(obs, i, arp_states->right_pos, rho, p, n, buf_max, known_prechange, true);
  focus_arp_one_iter_cpp(obs, i, arp_states->left_pos,  rho, p, n, buf_max, known_prechange, false);
  focus_arp_one_iter_cpp(-obs, i, arp_states->right_neg, rho, p, n, buf_max, known_prechange, true);
  focus_arp_one_iter_cpp(-obs, i, arp_states->left_neg,  rho, p, n, buf_max, known_prechange, false);
  
  // Take max among all four states
  const double max_vals[4] = {
    arp_states->right_pos.max_val, arp_states->left_pos.max_val,
    arp_states->right_neg.max_val, arp_states->left_neg.max_val
  };
  const int cpt_vals[4] = {
    arp_states->right_pos.cpt, arp_states->left_pos.cpt,
    arp_states->right_neg.cpt, arp_states->left_neg.cpt
  };
  
  int argmax = 0;
  double best = max_vals[0];
  for (int k = 1; k < 4; ++k) {
    if (max_vals[k] > best) { best = max_vals[k]; argmax = k; }
  }
  
  out_max_stat = best;
  out_cpt = cpt_vals[argmax];
}

// Cleanup function for ARpInfo destructor
void cleanup_arp_states(void* opaque_states) {
  if (opaque_states) {
    ARpStates* arp_states = static_cast<ARpStates*>(opaque_states);
    delete arp_states;
  }
}

} // namespace changepoint
