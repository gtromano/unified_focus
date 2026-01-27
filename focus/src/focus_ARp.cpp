// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

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

// ---------- Triple struct ----------
struct Triple {
  int tau;
  std::vector<double> S_tau;
  double l;
  double A, B, C, D, E, f;
  Triple(): tau(0), l(NA_REAL), A(NA_REAL), B(NA_REAL), C(NA_REAL), D(NA_REAL), E(NA_REAL), f(NA_REAL) {}
  Triple(int tau_, const std::vector<double>& S_tau_, double l_ = NA_REAL)
    : tau(tau_), S_tau(S_tau_), l(l_), A(NA_REAL), B(NA_REAL), C(NA_REAL), D(NA_REAL), E(NA_REAL), f(NA_REAL) {}
};

// ---------- forward declarations ----------
double intersec_point_newcurve(const std::vector<double>& S_new_padded, const std::vector<double>& rho, int current_n);
double intersec_point(const Triple& input_triple, const std::vector<double>& S_new, const std::vector<double>& rho, int current_n);
Triple coef_introduce(const Triple& q_new, const std::vector<double>& y_start, const std::vector<double>& rho, double sum_square);
Triple coef_introduce_pre0(const Triple& q_new, const std::vector<double>& rho);
double Q_n_mu_arp_unified(std::vector<Triple>& triples,
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
List max_val_compute_arp(const std::vector<Triple>& triples, int current_n, const std::vector<double>& rho, double sum_square, const std::vector<double>& y_start);
List max_val_compute_pre0(const std::vector<Triple>& triples);

struct State; // forward-declare struct if defined later

State init_state(double first_value);

void focus_arp_one_iter_cpp(const std::vector<double>& series,
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

// NEW overload with side-control (right=true, left=false). Keeps existing API intact.
double Q_n_mu_arp_unified(std::vector<Triple>& triples,
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
  if (!known_prechange) q_new = coef_introduce(q_new, y_start, rho, sum_square);
  else                  q_new = coef_introduce_pre0(q_new, rho);

  new_triples.push_back(q_new);
  triples.swap(new_triples);

  return S_n;
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

List max_val_compute_arp(const std::vector<Triple>& triples, int current_n, const std::vector<double>& rho, double sum_square, const std::vector<double>& y_start){
  int p = (int)rho.size();
  if (triples.empty()) return List::create(Named("cpt") = -1, Named("opt_max_val") = -INFINITY);

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
  std::vector<double> max_vals(K, -INFINITY);
  std::vector<int> cpts(K, -1);
  for (int k = 0; k < K; ++k){
    const Triple& tr = triples[k];
    double A = tr.A, B = tr.B, C = tr.C, D = tr.D, E = tr.E, f = tr.f;
    double denom = (4.0 * A * B - C * C);
    double max_val_change = -INFINITY;
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

  return List::create(Named("cpt") = cpts[opt_index], Named("opt_max_val") = best);
}

List max_val_compute_pre0(const std::vector<Triple>& triples){
  int K = (int)triples.size();
  if (K == 0) return List::create(Named("cpt") = -1, Named("opt_max_val") = -INFINITY);

  std::vector<double> max_vals(K, -INFINITY);
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
  return List::create(Named("cpt") = cpts[opt_index], Named("opt_max_val") = best);
}

struct State {
  std::vector<Triple> triples;
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

void focus_arp_one_iter_cpp(const std::vector<double>& series,
                            int i,
                            State& st,
                            const std::vector<double>& rho,
                            int p,
                            int n,
                            int buf_max,
                            bool known_prechange,
                            bool right_side) {
  // Append new observation to buffer
  const double x_new = series[i - 1];
  st.buf.push_back(x_new);
  if ((int)st.buf.size() > buf_max) {
    st.buf_sum_offset += st.buf.front();
    st.buf.erase(st.buf.begin());
    st.buf_start += 1;
  }

  bool enter_recursion = (i >= (p + 2));

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
    std::vector<double> S_tau; S_tau.reserve(2 * p + 1);
    for (int z = 0; z < p; ++z) S_tau.push_back(0.0);
    for (int k = 0; k <= p; ++k) {
      double s = 0.0;
      for (int t = 0; t <= k; ++t) s += series[t];
      S_tau.push_back(s);
    }

    Triple triple_1(1, S_tau);
    if (!known_prechange)
      triple_1 = coef_introduce(triple_1, st.y_start, rho, st.sum_square);
    else
      triple_1 = coef_introduce_pre0(triple_1, rho);

    st.triples.clear();
    st.triples.push_back(triple_1);

    st.S_n_1 = 0.0;
    for (int t = 0; t <= p; ++t) st.S_n_1 += series[t];

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

    double y_next_n = series[i - 1] - dot_vec(rho, lag_vec);

    if (!known_prechange) {
      st.sum_square += y_next_n * y_next_n;
      st.S_n_1 = Q_n_mu_arp_unified(st.triples, st.buf, st.buf_start, st.buf_sum_offset,
                                         series[i - 1], i, st.S_n_1, rho,
                                         known_prechange, st.sum_square, st.y_start,
                                         right_side);
      st.triples = coef_update_arp(st.triples, rho, y_next_n, known_prechange);
    } else {
      st.S_n_1 = Q_n_mu_arp_unified(st.triples, st.buf, st.buf_start, st.buf_sum_offset,
                                         series[i - 1], i, st.S_n_1, rho,
                                         known_prechange, 0.0, std::vector<double>(),
                                         right_side);
      st.triples = coef_update_arp(st.triples, rho, y_next_n, known_prechange);
    }
  }

  // Evaluate statistic
  if (!st.triples.empty()) {
    if (!known_prechange) {
      Rcpp::List max_info = max_val_compute_arp(st.triples, i, rho, st.sum_square, st.y_start);
      st.max_val = Rcpp::as<double>(max_info["opt_max_val"]);
      st.cpt     = Rcpp::as<int>(max_info["cpt"]);
    } else {
      Rcpp::List max_info = max_val_compute_pre0(st.triples);
      st.max_val = Rcpp::as<double>(max_info["opt_max_val"]);
      st.cpt     = Rcpp::as<int>(max_info["cpt"]);
    }
  } else {
    st.max_val = -1.0;
    st.cpt     = -1;
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
  std::vector<double> data_neg;  // Negative version of the data
};

namespace changepoint {

// Implementation function called from ARpInfo::update
void arp_detector_update_impl(const std::vector<double>& series,
                               const std::vector<double>& rho,
                               int p,
                               int buf_max,
                               bool known_prechange,
                               int n,
                               void*& opaque_states,
                               double& out_max_stat,
                               int& out_cpt) {
  int i = n;  // Current iteration index
  
  // Initialize opaque_states on first call (when n == 2)
  if (!opaque_states) {
    ARpStates* arp_states = new ARpStates();
    arp_states->right_pos = init_state(series[0]);
    arp_states->left_pos  = init_state(series[0]);
    
    // Prepare negative data from series
    arp_states->data_neg.assign(series.begin(), series.end());
    for (double &v : arp_states->data_neg) v = -v;
    
    arp_states->right_neg = init_state(arp_states->data_neg[0]);
    arp_states->left_neg  = init_state(arp_states->data_neg[0]);
    
    opaque_states = static_cast<void*>(arp_states);
  }
  
  ARpStates* arp_states = static_cast<ARpStates*>(opaque_states);
  
  // Ensure data_neg is kept in sync with negative of series
  if (arp_states->data_neg.size() != series.size()) {
    arp_states->data_neg.assign(series.begin(), series.end());
    for (double &v : arp_states->data_neg) v = -v;
  }
  
  // Update all four states using the full series (allows index access needed by focus_arp_one_iter_cpp)
  focus_arp_one_iter_cpp(series, i, arp_states->right_pos, rho, p, n, buf_max, known_prechange, true);
  focus_arp_one_iter_cpp(series, i, arp_states->left_pos,  rho, p, n, buf_max, known_prechange, false);
  focus_arp_one_iter_cpp(arp_states->data_neg, i, arp_states->right_neg, rho, p, n, buf_max, known_prechange, true);
  focus_arp_one_iter_cpp(arp_states->data_neg, i, arp_states->left_neg,  rho, p, n, buf_max, known_prechange, false);
  
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

// [[Rcpp::export]]
Rcpp::List Focus_arp_rcpp(Rcpp::NumericVector data_point_rcpp,
                          Rcpp::NumericVector rho_rcpp,
                          double lambda,
                          Rcpp::Nullable<Rcpp::NumericVector> pre_change_mean = R_NilValue) {

  // --- Inputs & pre-change mean handling ---
  std::vector<double> data_point = Rcpp::as<std::vector<double>>(data_point_rcpp);
  std::vector<double> rho        = Rcpp::as<std::vector<double>>(rho_rcpp);

  bool   known_prechange = false;
  double pre_mean = 0.0;
  if (pre_change_mean.isNotNull()) {
    Rcpp::NumericVector tmp(pre_change_mean.get());
    if (tmp.size() > 0) {
      pre_mean = tmp[0];
      known_prechange = true;
    }
  }
  if (known_prechange) {
    for (double &v : data_point) v -= pre_mean;
  }

  const int n = (int)data_point.size();
  const int p = (int)rho.size();
  const int buf_max = std::max(2 * p, p + 1);

  // Negative version of the data
  std::vector<double> data_neg(data_point.size());
  for (size_t ii = 0; ii < data_point.size(); ++ii) data_neg[ii] = -data_point[ii];

  // Initialise states (right/left × pos/neg)
  State state_right_pos = init_state(data_point.empty() ? 0.0 : data_point[0]);
  State state_left_pos  = init_state(data_point.empty() ? 0.0 : data_point[0]);
  State state_right_neg = init_state(data_neg.empty()   ? 0.0 : data_neg[0]);
  State state_left_neg  = init_state(data_neg.empty()   ? 0.0 : data_neg[0]);

  // Iteration loop
  int i = 2;
  bool no_detect = true;
  int stop_point = (n >= 2 ? 2 : n);
  int cpt = -1;
  std::vector<double> stat_history;

  while (no_detect && i <= n) {
    // Positive updates
    focus_arp_one_iter_cpp(data_point, i, state_right_pos, rho, p, n, buf_max, known_prechange, true);
    focus_arp_one_iter_cpp(data_point, i, state_left_pos,  rho, p, n, buf_max, known_prechange, false);

    // Negative updates
    focus_arp_one_iter_cpp(data_neg,   i, state_right_neg, rho, p, n, buf_max, known_prechange, true);
    focus_arp_one_iter_cpp(data_neg,   i, state_left_neg,  rho, p, n, buf_max, known_prechange, false);

    // // print the max_vals for debugging to output
    // Rcout << "i=" << i << ": "
    //       << "right_pos=" << state_right_pos.max_val << ", "
    //       << "left_pos="  << state_left_pos.max_val  << ", "
    //       << "right_neg=" << state_right_neg.max_val << ", "
    //       << "left_neg="  << state_left_neg.max_val  << "\n";

    // // print the number of triples for debugging to output
    // Rcout << "  #triples: "
    //       << "right_pos=" << state_right_pos.triples.size() << ", "
    //       << "left_pos="  << state_left_pos.triples.size()  << ", "
    //       << "right_neg=" << state_right_neg.triples.size() << ", "
    //       << "left_neg="  << state_left_neg.triples.size()  << "\n";

    // Take max among all four states
    const double max_vals[4] = {
      state_right_pos.max_val, state_left_pos.max_val,
      state_right_neg.max_val, state_left_neg.max_val
    };
    const int cpt_vals[4] = {
      state_right_pos.cpt, state_left_pos.cpt,
      state_right_neg.cpt, state_left_neg.cpt
    };

    int argmax = 0;
    double best = max_vals[0];
    for (int k = 1; k < 4; ++k) {
      if (max_vals[k] > best) { best = max_vals[k]; argmax = k; }
    }
    double max_val = best;
    int    cpt_local = cpt_vals[argmax];

    stat_history.push_back(max_val);

    if (max_val >= lambda) {
      no_detect = false;
      stop_point = i;
      cpt = cpt_local;
    } else {
      ++i;
      stop_point = i;
    }
  }

  if (no_detect) cpt = -1;

  return Rcpp::List::create(
    Rcpp::Named("cpt") = cpt,
    Rcpp::Named("stop_point") = stop_point,
    Rcpp::Named("stat_history") = Rcpp::wrap(stat_history)
  );
}

// [[Rcpp::export]]
Rcpp::List Focus_arp_rcpp_up_only(Rcpp::NumericVector data_point_rcpp,
                          Rcpp::NumericVector rho_rcpp,
                          double lambda,
                          Rcpp::Nullable<Rcpp::NumericVector> pre_change_mean = R_NilValue) {

  // --- Inputs & pre-change mean handling ---
  std::vector<double> data_point = Rcpp::as<std::vector<double>>(data_point_rcpp);
  std::vector<double> rho        = Rcpp::as<std::vector<double>>(rho_rcpp);

  bool   known_prechange = false;
  double pre_mean = 0.0;
  if (pre_change_mean.isNotNull()) {
    Rcpp::NumericVector tmp(pre_change_mean.get());
    if (tmp.size() > 0) {
      pre_mean = tmp[0];
      known_prechange = true;
    }
  }
  if (known_prechange) {
    for (double &v : data_point) v -= pre_mean;
  }

  const int n = (int)data_point.size();
  const int p = (int)rho.size();
  const int buf_max = std::max(2 * p, p + 1);

  // Negative version of the data
  std::vector<double> data_neg(data_point.size());
  for (size_t ii = 0; ii < data_point.size(); ++ii) data_neg[ii] = -data_point[ii];

  // Initialise states (right/left × pos/neg)
  State state_right_pos = init_state(data_point.empty() ? 0.0 : data_point[0]);
  State state_left_pos  = init_state(data_point.empty() ? 0.0 : data_point[0]);
  State state_right_neg = init_state(data_neg.empty()   ? 0.0 : data_neg[0]);
  State state_left_neg  = init_state(data_neg.empty()   ? 0.0 : data_neg[0]);

  // Iteration loop
  int i = 2;
  bool no_detect = true;
  int stop_point = (n >= 2 ? 2 : n);
  int cpt = -1;
  std::vector<double> stat_history;

  while (no_detect && i <= n) {
    // Positive updates
    focus_arp_one_iter_cpp(data_point, i, state_right_pos, rho, p, n, buf_max, known_prechange, true);
    //focus_arp_one_iter_cpp(data_point, i, state_left_pos,  rho, p, n, buf_max, known_prechange, false);

    // Negative updates
    //focus_arp_one_iter_cpp(data_neg,   i, state_right_neg, rho, p, n, buf_max, known_prechange, true);
    //focus_arp_one_iter_cpp(data_neg,   i, state_left_neg,  rho, p, n, buf_max, known_prechange, false);

    // // print the max_vals for debugging to output
    // Rcout << "i=" << i << ": "
    //       << "right_pos=" << state_right_pos.max_val << ", "
    //       << "left_pos="  << state_left_pos.max_val  << ", "
    //       << "right_neg=" << state_right_neg.max_val << ", "
    //       << "left_neg="  << state_left_neg.max_val  << "\n";

    // // print the number of triples for debugging to output
    // Rcout << "  #triples: "
    //       << "right_pos=" << state_right_pos.triples.size() << ", "
    //       << "left_pos="  << state_left_pos.triples.size()  << ", "
    //       << "right_neg=" << state_right_neg.triples.size() << ", "
    //       << "left_neg="  << state_left_neg.triples.size()  << "\n";

    // Take max among all four states
    const double max_vals[4] = {
      state_right_pos.max_val, state_left_pos.max_val,
      state_right_neg.max_val, state_left_neg.max_val
    };
    const int cpt_vals[4] = {
      state_right_pos.cpt, state_left_pos.cpt,
      state_right_neg.cpt, state_left_neg.cpt
    };

    int argmax = 0;
    double best = max_vals[0];
    for (int k = 1; k < 4; ++k) {
      if (max_vals[k] > best) { best = max_vals[k]; argmax = k; }
    }
    double max_val = best;
    int    cpt_local = cpt_vals[argmax];

    stat_history.push_back(max_val);

    if (max_val >= lambda) {
      no_detect = false;
      stop_point = i;
      cpt = cpt_local;
    } else {
      ++i;
      stop_point = i;
    }
  }

  if (no_detect) cpt = -1;

  return Rcpp::List::create(
    Rcpp::Named("cpt") = cpt,
    Rcpp::Named("stop_point") = stop_point,
    Rcpp::Named("stat_history") = Rcpp::wrap(stat_history)
  );
}
