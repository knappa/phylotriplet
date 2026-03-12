#pragma once
#include <Eigen/Core>
#include <array>
#include <stdexcept>

#include "sequence.hpp"

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
  auto [left, right, evals] = compute_decomposition(pis, ss);

  return compute_transition_matrices(left, right, evals, branch_lengths);
}

inline double compute_log_likelihood(
    const std::array<double, 4>& pis, Eigen::Matrix4d left,
    Eigen::Matrix4d right, Eigen::Vector4d evals,
    const std::array<double, 4>& branch_lengths,
    const TuplePatterns& site_patterns, const std::array<int, 3>& taxon_order) {
  // branch_lengths are (a,i), (b,i), (c,r), (i,r)
  // taxon_order identifies which taxon is a,b,c in site_patterns (which is
  // sorted by name)

  auto [p_t1, p_t2, p_t3, p_t4] =
      compute_transition_matrices(left, right, evals, branch_lengths);

  return compute_log_likelihood(pis, p_t1, p_t2, p_t3, p_t4, site_patterns,
                                taxon_order);
}

inline double compute_log_likelihood(
    const std::array<double, 4>& pis, const std::array<double, 6>& ss,
    const std::array<double, 4>& branch_lengths,
    const TuplePatterns& site_patterns, const std::array<int, 3>& taxon_order) {
  // branch_lengths are (a,i), (b,i), (c,r), (i,r)
  // taxon_order identifies which taxon is a,b,c in site_patterns (which is
  // sorted by name)
  auto [p_t1, p_t2, p_t3, p_t4] =
      compute_transition_matrices(pis, ss, branch_lengths);

  return compute_log_likelihood(pis, p_t1, p_t2, p_t3, p_t4, site_patterns,
                                taxon_order);
}

inline double compute_log_likelihood(const std::array<double, 4>& pis,
                                     Eigen::Matrix4d p_t1, Eigen::Matrix4d p_t2,
                                     Eigen::Matrix4d p_t3, Eigen::Matrix4d p_t4,
                                     const TuplePatterns& site_patterns,
                                     const std::array<int, 3>& taxon_order) {
  // branch_lengths are (a,i), (b,i), (c,r), (i,r)
  // taxon_order identifies which taxon is a,b,c in site_patterns (which is
  // sorted by name)
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

/**
 * STUB: Compute weights for the three possible rooted triplet topologies for
 * a specific taxon triplet, based on site-pattern counts across all genes.
 *
 * The three weights correspond to topologies where a, b, c are the sorted
 * taxon names:
 *   result[0] -- (a,b)|c
 *   result[1] -- (a,c)|b
 *   result[2] -- (b,c)|a
 *
 * Not yet implemented -- always throws std::runtime_error.
 */
inline std::array<double, 3> compute_triplet_weights(
    const MultiGeneAlignment& alignment, const std::string& taxon_a,
    const std::string& taxon_b, const std::string& taxon_c) {
  auto pis = alignment.nucleotide_frequencies();
  auto ss = std::array<double, 6>{1.0, 1.0, 1.0, 1.0, 1.0, 1.0};  // JC69 model

  return {0.0, 0.0, 0.0};
}

/**
 * STUB: Compute weights for the three possible rooted triplet topologies for
 * a specific taxon triplet, based on site-pattern counts in a single gene.
 *
 * The three weights correspond to topologies where a, b, c are the sorted
 * taxon names:
 *   result[0] -- (a,b)|c
 *   result[1] -- (a,c)|b
 *   result[2] -- (b,c)|a
 *
 * Not yet implemented -- always throws std::runtime_error.
 */
inline std::array<double, 3> compute_triplet_weights(
    const GeneAlignment& gene, const std::string& taxon_a,
    const std::string& taxon_b, const std::string& taxon_c) {
  throw std::runtime_error("compute_triplet_weights: not yet implemented");
}

}  // namespace phylotriplet