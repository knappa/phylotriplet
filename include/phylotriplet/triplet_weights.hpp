#pragma once
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#include "sequence.hpp"
#include "warn.hpp"

namespace phylotriplet {

inline Eigen::Matrix<double, 4, 4> symmetrized_gtr4_q_matrix(
    double pi_a, double pi_c, double pi_g, double pi_t, double s_ac,
    double s_ag, double s_at, double s_cg, double s_ct, double s_gt) {
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
    const std::array<double, 4>& pi, const std::array<double, 6>& s) {
  return symmetrized_gtr4_q_matrix(pi[0], pi[1], pi[2], pi[3], s[0], s[1], s[2],
                                   s[3], s[4], s[5]);
}

inline double gtr4_rate(double pi_a, double pi_c, double pi_g, double pi_t,
                        double s_ac, double s_ag, double s_at, double s_cg,
                        double s_ct, double s_gt) {
  return 2.0 * (pi_a * pi_c * s_ac + pi_a * pi_g * s_ag + pi_a * pi_t * s_at +
                pi_c * pi_g * s_cg + pi_c * pi_t * s_ct + pi_g * pi_t * s_gt);
}

inline double gtr4_rate(const std::array<double, 4>& pi,
                        const std::array<double, 6>& s) {
  return gtr4_rate(pi[0], pi[1], pi[2], pi[3], s[0], s[1], s[2], s[3], s[4],
                   s[5]);
}

inline std::pair<double, std::array<double, 6>> normalize_gtr_params(
    const std::array<double, 4>& pi, const std::array<double, 6>& s) {
  const double rate = gtr4_rate(pi, s);
  if (rate == 0.0) throw std::invalid_argument("Zero rate GTR parameters");

  return std::make_pair(
      rate, std::array<double, 6>{s[0] / rate, s[1] / rate, s[2] / rate,
                                  s[3] / rate, s[4] / rate, s[5] / rate});
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
  auto pi_sqrt_inv = Eigen::Vector4d(1.0 / sqrt(pis[0]), 1.0 / sqrt(pis[1]),
                                     1.0 / sqrt(pis[2]), 1.0 / sqrt(pis[3]))
                         .asDiagonal();

  auto left = pi_sqrt_inv * evecs;
  auto right = evecs.transpose() * pi_sqrt;

  return std::make_tuple(left, right, evals);
}

inline std::array<Eigen::Matrix4d, 4> compute_transition_matrices(
    Eigen::Matrix4d left, Eigen::Matrix4d right, Eigen::Vector4d evals,
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
    const std::array<double, 4>& pis, const std::array<double, 6>& ss,
    const std::array<double, 4>& branch_lengths) {
  Eigen::Matrix4d left, right;
  Eigen::Vector4d evals;
  std::tie(left, right, evals) = compute_decomposition(pis, ss);
  return compute_transition_matrices(left, right, evals, branch_lengths);
}

inline double compute_log_likelihood(const std::array<double, 4>& pis,
                                     Eigen::Matrix4d p_t1, Eigen::Matrix4d p_t2,
                                     Eigen::Matrix4d p_t3, Eigen::Matrix4d p_t4,
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
        if (site_pattern_count == 0.0) continue;

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

        std::sort(lse_summands.begin(), lse_summands.end(),
                  std::less<double>());

        double accum = 0.0;
        for (double s : lse_summands) {
          accum = s + log1p(exp(accum - s));
        }
        log_likelihood += site_pattern_count * accum;
      }
    }
  }

  return log_likelihood;
}

inline double compute_log_likelihood(
    const std::array<double, 4>& pis, Eigen::Matrix4d left,
    Eigen::Matrix4d right, Eigen::Vector4d evals,
    const std::array<double, 4>& branch_lengths,
    const TuplePatterns& site_patterns, const std::array<int, 3>& taxon_order) {
  const auto mats =
      compute_transition_matrices(left, right, evals, branch_lengths);
  return compute_log_likelihood(pis, mats[0], mats[1], mats[2], mats[3],
                                site_patterns, taxon_order);
}

inline double compute_log_likelihood(
    const std::array<double, 4>& pis, const std::array<double, 6>& ss,
    const std::array<double, 4>& branch_lengths,
    const TuplePatterns& site_patterns, const std::array<int, 3>& taxon_order) {
  const auto mats = compute_transition_matrices(pis, ss, branch_lengths);
  return compute_log_likelihood(pis, mats[0], mats[1], mats[2], mats[3],
                                site_patterns, taxon_order);
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

double log_diff_exp(double log_a, double log_b) {
  if (log_a < log_b) {
    throw std::invalid_argument(
        "log_diff_exp: log_a must be >= log_b, got log_a=" +
        std::to_string(log_a) + " log_b=" + std::to_string(log_b));
  }
  return log_a + log1p(-exp(log_b - log_a));
}

inline std::array<double, 4> compute_jc69_site_pattern_log_probs(double tau1,
                                                                 double tau2,
                                                                 double theta) {
  if (tau1 < 0.0 || tau2 < 0.0 || theta < 0.0) {
    throw std::invalid_argument(
        "compute_jc69_site_pattern_weights: tau1, tau2, and theta must be "
        "non-negative");
  }
  if (tau1 > tau2) {
    throw std::invalid_argument(
        "compute_jc69_site_pattern_weights: must have tau1 <= tau2 but " +
        std::to_string(tau1) + " > " + std::to_string(tau2));
  }
  if (tau1 == tau2) {
    // Special case: star-tree topology. All patterns are equally likely, so
    // weights are just proportional to the number of patterns in each class.
    return {log(4.0 / 64.0), log(12.0 / 64.0), log(24.0 / 64.0),
            log(24.0 / 64.0)};
  }

  bool warn = false;

  /* log becomes constexpr in C++26? change when that becomes common */
  const double log_2 = log(2.0);
  const double log_3 = log(3.0);
  const double log_6 = log(6.0);
  const double log_8 = log(8.0);
  const double log_9 = log(9.0);
  const double log_16 = log(16.0);
  const double log_18 = log(18.0);
  const double log_54 = log(54.0);
  const double log_64 = log(64.0);

  const double log1p_2theta_term = log1p(2.0 * theta / 3.0);
  const double log1p_4theta_term = log1p(4.0 * theta / 3.0);

  const double log_s_xxx =
      (-log_16 - log1p_4theta_term +
       log_sum_exp(std::array<double, 4>{
           log1p_4theta_term, -8.0 * tau1 / 3.0 + log_3,
           -8.0 * tau2 / 3.0 + log_6,
           -4.0 * (tau1 + 2 * tau2) / 3.0 + log_6 - log1p_2theta_term}));

  // log_s_xxy
  const double positive_log_sum_exp_term_s_xxy =
      log_sum_exp(std::array<double, 2>{
          log1p_4theta_term + log_3,
          -8.0 * tau1 / 3.0 + log_9,
      });
  const double negative_log_sum_exp_term_s_xxy =
      log_sum_exp(std::array<double, 2>{
          -8.0 * tau2 / 3.0 + log_6,
          -4.0 * (tau1 + 2 * tau2) / 3.0 + log_6 - log1p_2theta_term});

  // guard against numerical errors
  if (negative_log_sum_exp_term_s_xxy >= positive_log_sum_exp_term_s_xxy) {
    warn = true;
  }

  const double log_s_xxy =
      positive_log_sum_exp_term_s_xxy > negative_log_sum_exp_term_s_xxy
          ? (-log_16 - log1p_4theta_term +
             log_diff_exp(positive_log_sum_exp_term_s_xxy,
                          negative_log_sum_exp_term_s_xxy))
          : -std::numeric_limits<double>::max();

  // log_s_xyx

  const double positive_log_sum_exp_term_s_xyx =
      log_sum_exp(std::array<double, 2>{
          log1p_4theta_term,
          -8.0 * tau2 / 3.0 + log_2,
      });
  const double negative_log_sum_exp_term_s_xyx = log_sum_exp(
      std::array<double, 2>{-8.0 * tau1 / 3.0, -4.0 * (tau1 + 2 * tau2) / 3.0 +
                                                   log_2 - log1p_2theta_term});

  // guard against numerical errors
  if (negative_log_sum_exp_term_s_xyx >= positive_log_sum_exp_term_s_xyx) {
    warn = true;
  }

  const double log_s_xyx =
      positive_log_sum_exp_term_s_xyx > negative_log_sum_exp_term_s_xyx
          ? (log_3 - log_8 - log1p_4theta_term +
             log_diff_exp(positive_log_sum_exp_term_s_xyx,
                          negative_log_sum_exp_term_s_xyx))
          : -std::numeric_limits<double>::max();

  // log_s_xyz

  const double positive_log_sum_exp_term_s_xyz =
      log_sum_exp(std::array<double, 2>{
          log1p_4theta_term,
          -4.0 * (tau1 + 2 * tau2) / 3.0 + log_2 - log1p_2theta_term,
      });
  const double negative_log_sum_exp_term_s_xyz =
      log_sum_exp(std::array<double, 2>{
          -8.0 * tau1 / 3.0,
          -8.0 * tau2 / 3.0 + log_2,
      });

  // guard against numerical errors
  if (negative_log_sum_exp_term_s_xyz >= positive_log_sum_exp_term_s_xyz) {
    warn = true;
  }

  const double log_s_xyz =
      positive_log_sum_exp_term_s_xyz > negative_log_sum_exp_term_s_xyz
          ? (log_3 - log_8 - log1p_4theta_term +
             log_diff_exp(positive_log_sum_exp_term_s_xyz,
                          negative_log_sum_exp_term_s_xyz))
          : -std::numeric_limits<double>::max();

  if (warn) {
    PHYLOTRIPLET_WARN(
        "compute_jc69_site_pattern_weights: numerical instability; some "
        "weights may be inaccurate. Examine tau1, tau2, theta values. ");
  }

  return {log_s_xxx, log_s_xxy, log_s_xyx, log_s_xyz};
}

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
inline std::array<double, 3> compute_jc69_triplet_log_weights(
    const TuplePatterns& site_patterns, double theta) {
  const double n = site_patterns.n_counted;
  if (n == 0.0) {
    // TODO: warn when no sites are counted.
    return {0.0, 0.0, 0.0};
  }

  double pxxx = 0.0;
  for (size_t i = 0; i < 4; ++i) {
    pxxx += site_patterns.site_pattern[site_pattern_index(i, i, i)] / n;
  }

  std::array<double, 3> pxxy = {0.0, 0.0, 0.0};
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      if (i == j) continue;
      pxxy[0] += site_patterns.site_pattern[site_pattern_index(i, i, j)] / n;
      pxxy[1] += site_patterns.site_pattern[site_pattern_index(i, j, i)] / n;
      pxxy[2] += site_patterns.site_pattern[site_pattern_index(j, i, i)] / n;
    }
  }

  std::array<double, 3> pxyx = {0.0, 0.0, 0.0};
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      if (i == j) continue;
      pxyx[0] += site_patterns.site_pattern[site_pattern_index(i, j, i)] / n +
                 site_patterns.site_pattern[site_pattern_index(j, i, i)] / n;
      pxyx[1] += site_patterns.site_pattern[site_pattern_index(i, i, j)] / n +
                 site_patterns.site_pattern[site_pattern_index(j, i, i)] / n;
      pxyx[2] += site_patterns.site_pattern[site_pattern_index(i, i, j)] / n +
                 site_patterns.site_pattern[site_pattern_index(i, j, i)] / n;
    }
  }

  double pxyz = 0.0;
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      for (size_t k = 0; k < 4; ++k) {
        if (i == j || i == k || j == k) continue;
        pxyz += site_patterns.site_pattern[site_pattern_index(i, j, k)] / n;
      }
    }
  }

  // exp((-8/3)*(\tau1))
  //   = (4/3)*(pxxx + 3*pxxy - 2*(pxyx + pxyz))(3 + 4*\theta)
  //   = 4*(pxxx + 3*pxxy - 2*(pxyx + pxyz)) * (1 + 4/3*\theta)
  const double log_3 = log(3.0);
  const double log_4 = log(4.0);
  const double log1p_43theta_term = log1p((4.0 / 3.0) * theta);
  std::array<double, 3> tau1 = {
      (-3.0 / 8.0) * (log(3 * pxxx + 3 * pxxy[0] - pxyx[0] + pxyz) - log_3 +
                      log1p_43theta_term),
      (-3.0 / 8.0) * (log(3 * pxxx + 3 * pxxy[1] - pxyx[1] + pxyz) - log_3 +
                      log1p_43theta_term),
      (-3.0 / 8.0) * (log(3 * pxxx + 3 * pxxy[2] - pxyx[2] + pxyz) - log_3 +
                      log1p_43theta_term)};
  // exp((-8/3)*(\tau2))
  //   = (4/3)*(pxxx - pxxy + 2*(pxyx - pxyz)) (3 + 4*\theta)
  //   = 4*(pxxx - pxxy + 2*(pxyx - pxyz)) * (1 + 4/3*\theta)
  std::array<double, 3> tau2 = {
      (-3.0 / 8.0) * (log(3 * pxxx - pxxy[0] + pxyx[0] - pxyz) - log_3 +
                      log1p_43theta_term),
      (-3.0 / 8.0) * (log(3 * pxxx - pxxy[1] + pxyx[1] - pxyz) - log_3 +
                      log1p_43theta_term),
      (-3.0 / 8.0) * (log(3 * pxxx - pxxy[2] + pxyx[2] - pxyz) - log_3 +
                      log1p_43theta_term)};

  std::array<double, 3> weights = {0.0, 0.0, 0.0};
  for (size_t i = 0; i < 3; ++i) {
    if (tau1[i] < 0 || tau2[i] < 0 || tau1[i] > tau2[i]) {
      // zero weight for this topology
      weights[i] = -std::numeric_limits<double>::infinity();  
    } else {
      auto site_pattern_weights =
          compute_jc69_site_pattern_log_probs(tau1[i], tau2[i], theta);
      weights[i] =
          site_pattern_weights[0] * pxxx + site_pattern_weights[1] * pxxy[i] +
          site_pattern_weights[2] * pxyx[i] + site_pattern_weights[3] * pxyz;
    }
  }

  return weights;
}

}  // namespace phylotriplet