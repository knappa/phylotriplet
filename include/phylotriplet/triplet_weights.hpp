#pragma once

#include <LBFGSB.h>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#include "sequence.hpp"
#include "warn.hpp"

namespace phylotriplet {

inline Eigen::Matrix<double, 4, 4> symmetrized_gtr4_q_matrix(double pi_a,
                                                             double pi_c,
                                                             double pi_g,
                                                             double pi_t,
                                                             double s_ac,
                                                             double s_ag,
                                                             double s_at,
                                                             double s_cg,
                                                             double s_ct,
                                                             double s_gt) {
  Eigen::Matrix<double, 4, 4> q;
  q << -(pi_c * s_ac + pi_g * s_ag + pi_t * s_at), sqrt(pi_a * pi_c) * s_ac,
      sqrt(pi_a * pi_g) * s_ag, sqrt(pi_a * pi_t) * s_at,
      sqrt(pi_a * pi_c) * s_ac, -(pi_a * s_ac + pi_g * s_cg + pi_t * s_ct),
      sqrt(pi_c * pi_g) * s_cg, sqrt(pi_c * pi_t) * s_ct,
      sqrt(pi_a * pi_g) * s_ag, sqrt(pi_c * pi_g) * s_cg,
      -(pi_a * s_ag + pi_c * s_cg + pi_t * s_gt), sqrt(pi_g * pi_t) * s_gt,
      sqrt(pi_a * pi_t) * s_at, sqrt(pi_c * pi_t) * s_ct,
      sqrt(pi_g * pi_t) * s_gt, -(pi_a * s_at + pi_c * s_ct + pi_g * s_gt);

  return q;
}

inline Eigen::Matrix<double, 4, 4> symmetrized_gtr4_q_matrix(
    const std::array<double, 4>& pi,
    const std::array<double, 6>& s) {
  return symmetrized_gtr4_q_matrix(
      pi[0], pi[1], pi[2], pi[3], s[0], s[1], s[2], s[3], s[4], s[5]);
}

inline double gtr4_rate(double pi_a,
                        double pi_c,
                        double pi_g,
                        double pi_t,
                        double s_ac,
                        double s_ag,
                        double s_at,
                        double s_cg,
                        double s_ct,
                        double s_gt) {
  return 2.0 * (pi_a * pi_c * s_ac + pi_a * pi_g * s_ag + pi_a * pi_t * s_at +
                pi_c * pi_g * s_cg + pi_c * pi_t * s_ct + pi_g * pi_t * s_gt);
}

inline double gtr4_rate(const std::array<double, 4>& pi,
                        const std::array<double, 6>& s) {
  return gtr4_rate(
      pi[0], pi[1], pi[2], pi[3], s[0], s[1], s[2], s[3], s[4], s[5]);
}

inline std::pair<double, std::array<double, 6>> normalize_gtr_params(
    const std::array<double, 4>& pi,
    const std::array<double, 6>& s) {
  const double rate = gtr4_rate(pi, s);
  if (rate == 0.0)
    throw std::invalid_argument("Zero rate GTR parameters");

  return std::make_pair(rate,
                        std::array<double, 6>{s[0] / rate,
                                              s[1] / rate,
                                              s[2] / rate,
                                              s[3] / rate,
                                              s[4] / rate,
                                              s[5] / rate});
}

inline std::tuple<Eigen::Matrix4d, Eigen::Matrix4d, Eigen::Vector4d>
compute_decomposition(const std::array<double, 4>& pis,
                      const std::array<double, 6>& ss) {
  Eigen::Matrix4d sym_q = symmetrized_gtr4_q_matrix(pis, ss);

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> eigen_solver(sym_q);

  auto evals = eigen_solver.eigenvalues();
  auto evecs = eigen_solver.eigenvectors();

  auto pi_sqrt =
      Eigen::Vector4d(sqrt(pis[0]), sqrt(pis[1]), sqrt(pis[2]), sqrt(pis[3]))
          .asDiagonal();
  auto pi_sqrt_inv = Eigen::Vector4d(1.0 / sqrt(pis[0]),
                                     1.0 / sqrt(pis[1]),
                                     1.0 / sqrt(pis[2]),
                                     1.0 / sqrt(pis[3]))
                         .asDiagonal();

  auto left = pi_sqrt_inv * evecs;
  auto right = evecs.transpose() * pi_sqrt;

  return std::make_tuple(left, right, evals);
}

inline std::array<Eigen::Matrix4d, 4> compute_transition_matrices(
    Eigen::Matrix4d left,
    Eigen::Matrix4d right,
    Eigen::Vector4d evals,
    const std::array<double, 4>& branch_lengths) {
  std::array<Eigen::Matrix4d, 4> transition_matrices;
  for (size_t i = 0; i < branch_lengths.size(); ++i) {
    transition_matrices[i] =
        left * (evals.array() * branch_lengths[i]).exp().matrix().asDiagonal() *
        right;
  }
  return transition_matrices;
}

inline std::array<Eigen::Matrix4d, 4> compute_transition_matrices(
    const std::array<double, 4>& pis,
    const std::array<double, 6>& ss,
    const std::array<double, 4>& branch_lengths) {
  Eigen::Matrix4d left, right;
  Eigen::Vector4d evals;
  std::tie(left, right, evals) = compute_decomposition(pis, ss);
  return compute_transition_matrices(left, right, evals, branch_lengths);
}

template <typename T, size_t N>
inline T log_sum_exp(std::array<T, N> values) {
  std::sort(values.begin(), values.end(), std::less<T>());
  // skip -infinity values to avoid NaN from (-inf) - (-inf)
  size_t start = 0;
  while (start < N && values[start] == -std::numeric_limits<T>::infinity()) {
    ++start;
  }
  if (start == N) {
    return -std::numeric_limits<T>::infinity();
  }
  T accum = values[start];
  for (size_t i = start + 1; i < N; ++i) {
    accum = values[i] + log1p(exp(accum - values[i]));
  }
  return accum;
}

inline double log_diff_exp(double log_a, double log_b) {
  if (log_a < log_b) {
    throw std::invalid_argument(
        "log_diff_exp: log_a must be >= log_b, got log_a=" +
        std::to_string(log_a) + " log_b=" + std::to_string(log_b));
  }
  return log_a + log1p(-exp(log_b - log_a));
}

inline double compute_log_likelihood(const std::array<double, 4>& pis,
                                     Eigen::Matrix4d p_t1,
                                     Eigen::Matrix4d p_t2,
                                     Eigen::Matrix4d p_t3,
                                     Eigen::Matrix4d p_t4,
                                     const TuplePatterns& site_patterns,
                                     const std::array<int, 3>& taxon_order) {
  // branch_lengths are (a,i), (b,i), (c,r), (i,r)
  // taxon_order identifies which taxon is a,b,c in site_patterns (which is
  // sorted by name) should be 0,1,2 in some order. will fail in interesting
  // ways if not.
  auto log_p_t1 = p_t1.array().max(1e-10).log().matrix();
  auto log_p_t2 = p_t2.array().max(1e-10).log().matrix();
  auto log_p_t3 = p_t3.array().max(1e-10).log().matrix();
  auto log_p_t4 = p_t4.array().max(1e-10).log().matrix();

  auto log_pis = Eigen::Vector4d(pis[0], pis[1], pis[2], pis[3])
                     .array()
                     .max(1e-10)
                     .log()
                     .matrix();

  double log_likelihood = 0.0;

  const size_t a_coeff = pow(4, 2 - taxon_order[0]);
  const size_t b_coeff = pow(4, 2 - taxon_order[1]);
  const size_t c_coeff = pow(4, 2 - taxon_order[2]);

  // Iterate over all 4^3=64 site patterns.
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = 0; b < 4; ++b) {
      for (size_t c = 0; c < 4; ++c) {
        const auto site_pattern_count =
            site_patterns.site_pattern[a * a_coeff + b * b_coeff + c * c_coeff];
        if (site_pattern_count == 0.0)
          continue;

        // for each site pattern, there are 4^2 = 16 possible internal states.
        // these are combined with a log-sum-exp. (sorted version for accuracy)
        std::array<double, 16> lse_summands = {};
        size_t i = 0;
        for (size_t x = 0; x < 4; ++x) {
          for (size_t y = 0; y < 4; ++y) {
            lse_summands[i++] = log_p_t1(x, a) + log_p_t2(x, b) +
                                log_p_t3(y, c) + log_p_t4(y, x) + log_pis(y);
          }
        }

        log_likelihood += site_pattern_count * log_sum_exp(lse_summands);
      }
    }
  }

  return log_likelihood;
}

inline double compute_log_likelihood(
    const std::array<double, 4>& pis,
    Eigen::Matrix4d left,
    Eigen::Matrix4d right,
    Eigen::Vector4d evals,
    const std::array<double, 4>& branch_lengths,
    const TuplePatterns& site_patterns,
    const std::array<int, 3>& taxon_order) {
  const auto mats =
      compute_transition_matrices(left, right, evals, branch_lengths);
  return compute_log_likelihood(
      pis, mats[0], mats[1], mats[2], mats[3], site_patterns, taxon_order);
}

inline double compute_log_likelihood(
    const std::array<double, 4>& pis,
    const std::array<double, 6>& ss,
    const std::array<double, 4>& branch_lengths,
    const TuplePatterns& site_patterns,
    const std::array<int, 3>& taxon_order) {
  const auto mats = compute_transition_matrices(pis, ss, branch_lengths);
  return compute_log_likelihood(
      pis, mats[0], mats[1], mats[2], mats[3], site_patterns, taxon_order);
}

inline std::array<double, 4> compute_jc69_site_pattern_log_probs(
    const double tau1,
    const double tau2,
    const double theta,
    const double mu = 4.0 / 3.0) {
  // this function does not check the validity of the parameters but
  // expects that tau2 >= tau1 >= 0, theta > 0, mu > 0

  const double tau1t = tau1 / mu;
  const double tau2t = tau2 / mu;
  const double beta = 1 / (mu * theta);

  const double log_beta = log(beta);
  const double log1p_beta = log1p(beta);
  const double log1p_2beta = log1p(2.0 * beta);
  const double log1p_3beta = log1p(3.0 * beta);

  const double log_s_xxx = -2 * tau1 - 2 * tau2 - 6 * log(2) -
                           log(2 * beta + 1) - log(beta + 1) +
                           log_sum_exp(std::array<double, 8>{
                               2 * tau1 + 2 * tau2 + log(2) + 2 * log(beta),
                               2 * tau1 + 2 * tau2 + log(3) + log(beta),
                               2 * tau1 + log(12) + 2 * log(beta),
                               2 * tau2 + log(6) + 2 * log(beta),
                               2 * tau1 + 2 * tau2,
                               tau1 + log(12) + 2 * log(beta),
                               2 * tau1 + log(6) + log(beta),
                               2 * tau2 + log(3) + log(beta)});

  const double log_s_xxy =
      -2 * (beta + 1) * tau2 - 2 * tau1 - 6 * log(2) - log(3 * beta + 1) -
      log(2 * beta + 1) - log(beta + 1) +
      log_diff_exp(
          log_sum_exp(std::array<double, 9>{
              2 * (beta + 1) * tau2 + 2 * tau1 + log(6) + 3 * log(beta),
              2 * (beta + 1) * tau2 + 2 * tau1 + log(11) + 2 * log(beta),
              2 * (beta + 1) * tau2 + log(18) + 3 * log(beta),
              2 * (beta + 1) * tau2 + 2 * tau1 + log(6) + log(beta),
              2 * (beta + 1) * tau1 + 4 * log(2) + 2 * log(beta),
              2 * (beta + 1) * tau2 + log(15) + 2 * log(beta),
              2 * (beta + 1) * tau2 + 2 * tau1,
              2 * (beta + 1) * tau1 + 3 * log(2) + log(beta),
              2 * (beta + 1) * tau2 + log(3) + log(beta)}),
          log_sum_exp(std::array<double, 5>{
              2 * beta * tau2 + 2 * tau1 + log(12) + 3 * log(beta),
              2 * beta * tau2 + tau1 + log(12) + 3 * log(beta),
              2 * beta * tau2 + 2 * tau1 + log(10) + 2 * log(beta),
              2 * beta * tau2 + tau1 + 2 * log(2) + 2 * log(beta),
              2 * beta * tau2 + 2 * tau1 + log(2) + log(beta)}));

  const double log_s_xyx =
      -2 * (beta + 1) * tau2 - 2 * tau1 - 6 * log(2) - log(3 * beta + 1) -
      log(2 * beta + 1) - log(beta + 1) +
      log_diff_exp(
          log_sum_exp(std::array<double, 7>{
              2 * beta * tau2 + 2 * tau1 + log(12) + 3 * log(beta),
              2 * (beta + 1) * tau2 + 2 * tau1 + log(6) + 3 * log(beta),
              2 * beta * tau2 + 2 * tau1 + log(10) + 2 * log(beta),
              2 * (beta + 1) * tau2 + 2 * tau1 + log(11) + 2 * log(beta),
              2 * beta * tau2 + 2 * tau1 + log(2) + log(beta),
              2 * (beta + 1) * tau2 + 2 * tau1 + log(6) + log(beta),
              2 * (beta + 1) * tau2 + 2 * tau1}),
          log_sum_exp(std::array<double, 7>{
              2 * beta * tau2 + tau1 + log(12) + 3 * log(beta),
              2 * (beta + 1) * tau2 + log(6) + 3 * log(beta),
              2 * beta * tau2 + tau1 + 2 * log(2) + 2 * log(beta),
              2 * (beta + 1) * tau1 + 3 * log(2) + 2 * log(beta),
              2 * (beta + 1) * tau2 + log(5) + 2 * log(beta),
              2 * (beta + 1) * tau1 + 2 * log(2) + log(beta),
              2 * (beta + 1) * tau2 + log(beta)}));

  const double log_s_xyz =
      -2 * (beta + 1) * tau2 - 2 * tau1 - 6 * log(2) - log(3 * beta + 1) -
      log(2 * beta + 1) - log(beta + 1) +
      log_diff_exp(log_sum_exp(std::array<double, 4>{
                       2 * tau1 + 2 * tau2 + log(2) + 2 * log(beta),
                       2 * tau1 + 2 * tau2 + log(3) + log(beta),
                       2 * tau1 + 2 * tau2,
                       tau1 + 2 * log(2) + 2 * log(beta)}),
                   log_sum_exp(std::array<double, 4>{
                       2 * tau1 + 2 * log(2) + 2 * log(beta),
                       2 * tau2 + log(2) + 2 * log(beta),
                       2 * tau1 + log(2) + log(beta),
                       2 * tau2 + log(beta)}));

  return {log_s_xxx, log_s_xxy, log_s_xyx, log_s_xyz};
}

class SpeciesTripletJC69NegLogLikelihood {
 public:
  SpeciesTripletJC69NegLogLikelihood(const std::array<double, 4>& site_patterns,
                               double theta,
                               double mu)
      : site_patterns_(site_patterns), theta_(theta), mu_(mu) {}

  double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
    const double tau1 = x[0];
    const double tau2 = x[0] + x[1];

    const auto site_pattern_weights =
        compute_jc69_site_pattern_log_probs(tau1, tau2, theta_, mu_);

    const double cxxx = site_patterns_[0];
    const double cxxy = site_patterns_[1];
    const double cxyx = site_patterns_[2];
    const double cxyz = site_patterns_[3];

    const double fx =
        -(site_pattern_weights[0] * cxxx + site_pattern_weights[1] * cxxy +
          site_pattern_weights[2] * cxyx + site_pattern_weights[3] * cxyz);

    grad.resize(2);

    const double x0_plus = std::min(0.001, (tau2 - tau1) * 0.1);
    const auto site_pattern_weights_tau1_plus =
        compute_jc69_site_pattern_log_probs(tau1 + x0_plus, tau2 + x0_plus, theta_, mu_);

    const double x0_minus = std::min(0.001, tau1 * 0.1);
    const auto site_pattern_weights_tau1_minus =
        compute_jc69_site_pattern_log_probs(tau1 - x0_minus, tau2 - x0_minus, theta_, mu_);
    grad[0] = -((site_pattern_weights_tau1_plus[0] -
                site_pattern_weights_tau1_minus[0]) *
                   cxxx +
               (site_pattern_weights_tau1_plus[1] -
                site_pattern_weights_tau1_minus[1]) *
                   cxxy +
               (site_pattern_weights_tau1_plus[2] -
                site_pattern_weights_tau1_minus[2]) *
                   cxyx +
               (site_pattern_weights_tau1_plus[3] -
                site_pattern_weights_tau1_minus[3]) *
                   cxyz) /
              (x0_plus + x0_minus);

    const double x1_plus = 0.001;
    const auto site_pattern_weights_tau2_plus =
        compute_jc69_site_pattern_log_probs(tau1, tau2 + x1_plus, theta_, mu_);

    const double x1_minus = std::min(0.001, (tau2 - tau1) * 0.1);
    const auto site_pattern_weights_tau2_minus =
        compute_jc69_site_pattern_log_probs(tau1, tau2 - x1_minus, theta_, mu_);
    grad[1] = -((site_pattern_weights_tau2_plus[0] -
                site_pattern_weights_tau2_minus[0]) *
                   cxxx +
               (site_pattern_weights_tau2_plus[1] -
                site_pattern_weights_tau2_minus[1]) *
                   cxxy +
               (site_pattern_weights_tau2_plus[2] -
                site_pattern_weights_tau2_minus[2]) *
                   cxyx +
               (site_pattern_weights_tau2_plus[3] -
                site_pattern_weights_tau2_minus[3]) *
                   cxyz) /
              (x1_plus + x1_minus);

    return fx;
  }

 private:
  const std::array<double, 4>& site_patterns_;
  const double theta_;
  const double mu_;
};

constexpr size_t site_pattern_index(size_t a, size_t b, size_t c) {
  return a * 16 + b * 4 + c;
}

/**
 * Compute JC69 triplet weights from pre-computed site patterns.
 *
 * The three weights correspond to topologies where a, b, c are the sorted
 * taxon names (matching site_patterns.taxa order):
 *   result[0] -- (a,b)|c
 *   result[1] -- (a,c)|b
 *   result[2] -- (b,c)|a
 */
inline std::array<double, 3> compute_jc69_triplet_weights(
    const TuplePatterns& site_patterns,
    const double theta,
    const double mu = 4.0 / 3.0) {
  const double n = site_patterns.n_counted;
  if (n == 0.0) {
    // TODO: warn when no sites are counted.
    return {0.0, 0.0, 0.0};
  }

  // cxxx, the count of sites where all three taxa have the same nucleotide
  double cxxx = 0.0;
  for (size_t i = 0; i < 4; ++i) {
    cxxx += site_patterns.site_pattern[site_pattern_index(i, i, i)];
  }

  // cxxy, the count of sites where the cherry taxa have same nucleotide and the
  // outgroup has a different nucleotide
  std::array<double, 3> cxxy = {0.0, 0.0, 0.0};
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      if (i == j)
        continue;
      cxxy[0] += site_patterns.site_pattern[site_pattern_index(i, i, j)];
      cxxy[1] += site_patterns.site_pattern[site_pattern_index(i, j, i)];
      cxxy[2] += site_patterns.site_pattern[site_pattern_index(j, i, i)];
    }
  }

  // cxyx, the count of sites where the cherry taxa different nucleotides and
  // the outgroup matches one of the cherry taxa
  std::array<double, 3> cxyx = {0.0, 0.0, 0.0};
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      if (i == j)
        continue;
      cxyx[0] += site_patterns.site_pattern[site_pattern_index(i, j, i)] +
                 site_patterns.site_pattern[site_pattern_index(j, i, i)];
      cxyx[1] += site_patterns.site_pattern[site_pattern_index(i, i, j)] +
                 site_patterns.site_pattern[site_pattern_index(j, i, i)];
      cxyx[2] += site_patterns.site_pattern[site_pattern_index(i, i, j)] +
                 site_patterns.site_pattern[site_pattern_index(i, j, i)];
    }
  }

  // cxyz, the count of sites where the cherry taxa and outgroup all have
  // different nucleotides
  double cxyz = 0.0;
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      for (size_t k = 0; k < 4; ++k) {
        if (i == j || i == k || j == k)
          continue;
        cxyz += site_patterns.site_pattern[site_pattern_index(i, j, k)];
      }
    }
  }

  const double pxxx = cxxx / n;
  const std::array<double, 3> pxxy = {cxxy[0] / n, cxxy[1] / n, cxxy[2] / n};
  const std::array<double, 3> pxyx = {cxyx[0] / n, cxyx[1] / n, cxyx[2] / n};
  const double pxyz = cxyz / n;

  const double log_3 = log(3.0);
  const double log_4 = log(4.0);
  const double log1p_43theta_term = log1p((4.0 / 3.0) * theta);

  // exp((-8/3)*(\tau2))
  //   = (4/3)*(pxxx - pxxy + 2*(pxyx - pxyz)) (3 + 4*\theta)
  //   = 4*(pxxx - pxxy + 2*(pxyx - pxyz)) * (1 + 4/3*\theta)
  std::array<double, 3> tau2 = {
      std::max(0.0,
               (-3.0 / 8.0) * (log(3 * pxxx - pxxy[0] + pxyx[0] - pxyz) -
                               log(n) - log_3 + log1p_43theta_term)),
      std::max(0.0,
               (-3.0 / 8.0) * (log(3 * pxxx - pxxy[1] + pxyx[1] - pxyz) -
                               log(n) - log_3 + log1p_43theta_term)),
      std::max(0.0,
               (-3.0 / 8.0) * (log(3 * pxxx - pxxy[2] + pxyx[2] - pxyz) -
                               log(n) - log_3 + log1p_43theta_term))};
  // exp((-8/3)*(\tau1))
  //   = (4/3)*(pxxx + 3*pxxy - 2*(pxyx + pxyz))(3 + 4*\theta)
  //   = 4*(pxxx + 3*pxxy - 2*(pxyx + pxyz)) * (1 + 4/3*\theta)
  std::array<double, 3> tau1 = {
      std::min(tau2[0],
               std::max(0.0,
                        (-3.0 / 8.0) *
                            (log(3 * pxxx + 3 * pxxy[0] - pxyx[0] + pxyz) -
                             log(n) - log_3 + log1p_43theta_term))),
      std::min(tau2[1],
               std::max(0.0,
                        (-3.0 / 8.0) *
                            (log(3 * pxxx + 3 * pxxy[1] - pxyx[1] + pxyz) -
                             log(n) - log_3 + log1p_43theta_term))),
      std::min(tau2[2],
               std::max(0.0,
                        (-3.0 / 8.0) *
                            (log(3 * pxxx + 3 * pxxy[2] - pxyx[2] + pxyz) -
                             log(n) - log_3 + log1p_43theta_term)))};

  std::array<double, 3> weights = {0.0, 0.0, 0.0};
  for (size_t i = 0; i < 3; ++i) {

    const double best_tau1_init = std::max(0.001, tau1[i]);
    const double best_tau2_init = std::max(best_tau1_init+0.001, tau2[i]);

    const auto site_pattern_weights =
        compute_jc69_site_pattern_log_probs(best_tau1_init, best_tau2_init, theta, mu);

    const double init_neg_ll =
        -(site_pattern_weights[0] * cxxx + site_pattern_weights[1] * cxxy[i] +
        site_pattern_weights[2] * cxyx[i] + site_pattern_weights[3] * cxyz);

    auto objective = SpeciesTripletJC69NegLogLikelihood(std::array<double,4>{cxxx, cxxy[i], cxyx[i], cxyz}, theta, mu);

    LBFGSpp::LBFGSBParam<double> param;
    param.epsilon = 1e-6;
    param.max_iterations = 200;
    param.past = 0;

    LBFGSpp::LBFGSBSolver<double> solver(param);

    Eigen::VectorXd lb(2), ub(2);
    lb << 1e-6, 1e-6;
    ub << 50.0, 50.0;  // TODO: determine reasonable upper bounds

    Eigen::VectorXd x(2);
    x << best_tau1_init, best_tau2_init;

    std::cout << best_tau1_init << " " << best_tau2_init << std::endl;

    double fx = init_neg_ll;
    try {
      solver.minimize(objective, x, fx, lb, ub);

      tau1[i] = x[0];
      tau2[i] = x[0] + x[1];
    } catch (const std::exception&) {
      // If L-BFGS-B fails, keep the initial guess.
      if (fx < init_neg_ll) {
        PHYLOTRIPLET_WARN("L-BFGS-B optimization failed, but found a better solution. Using the L-BFGS-B solution.");
        tau1[i] = x[0];
        tau2[i] = x[0] + x[1];
      } else {
        PHYLOTRIPLET_WARN("L-BFGS-B optimization failed, using the initial guess.");
        tau1[i] = best_tau1_init;
        tau2[i] = best_tau2_init;
        fx = init_neg_ll;
      }
    }

    weights[i] = -fx;

  }

  std::cout << tau1[0] << " " << tau2[0] << " " << tau1[1] << " " << tau2[1] << " "
            << tau1[2] << " " << tau2[2] << std::endl;

  return weights;
}

}  // namespace phylotriplet