#pragma once

#include <LBFGSB.h>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <unsupported/Eigen/MatrixFunctions>

#include "sequence.hpp"
#include "util.hpp"
#include "warn.hpp"

namespace phylotriplet {

/* ********** ********** ********** ********** ********** ********** ********** */
/* GTR rate matrix computations */
/* ********** ********** ********** ********** ********** ********** ********** */

inline Eigen::Matrix4d symmetrized_gtr4_q_matrix(const double pi_a,
                                                 const double pi_c,
                                                 const double pi_g,
                                                 const double pi_t,
                                                 const double s_ac,
                                                 const double s_ag,
                                                 const double s_at,
                                                 const double s_cg,
                                                 const double s_ct,
                                                 const double s_gt) {
  Eigen::Matrix4d q;
  q << -(pi_c * s_ac + pi_g * s_ag + pi_t * s_at), sqrt(pi_a * pi_c) * s_ac,
      sqrt(pi_a * pi_g) * s_ag, sqrt(pi_a * pi_t) * s_at, sqrt(pi_a * pi_c) * s_ac,
      -(pi_a * s_ac + pi_g * s_cg + pi_t * s_ct), sqrt(pi_c * pi_g) * s_cg,
      sqrt(pi_c * pi_t) * s_ct, sqrt(pi_a * pi_g) * s_ag, sqrt(pi_c * pi_g) * s_cg,
      -(pi_a * s_ag + pi_c * s_cg + pi_t * s_gt), sqrt(pi_g * pi_t) * s_gt,
      sqrt(pi_a * pi_t) * s_at, sqrt(pi_c * pi_t) * s_ct, sqrt(pi_g * pi_t) * s_gt,
      -(pi_a * s_at + pi_c * s_ct + pi_g * s_gt);

  return q;
}

inline Eigen::Matrix4d symmetrized_gtr4_q_matrix(const std::array<double, 4>& pi,
                                                 const std::array<double, 6>& s) {
  return symmetrized_gtr4_q_matrix(pi[0], pi[1], pi[2], pi[3], s[0], s[1], s[2], s[3], s[4], s[5]);
}

inline Eigen::Matrix4d gtr4_q_matrix(const double pi_a,
                                     const double pi_c,
                                     const double pi_g,
                                     const double pi_t,
                                     const double s_ac,
                                     const double s_ag,
                                     const double s_at,
                                     const double s_cg,
                                     const double s_ct,
                                     const double s_gt) {
  Eigen::Matrix4d q;
  q << -(s_ac * pi_c + s_ag * pi_g + s_at * pi_t), s_ac * pi_c, s_ag * pi_g, s_at * pi_t,
      pi_a * s_ac, -(pi_a * s_ac + s_cg * pi_g + s_ct * pi_t), s_cg * pi_g, s_ct * pi_t,
      pi_a * s_ag, pi_c * s_cg, -(pi_a * s_ag + pi_c * s_cg + s_gt * pi_t), s_gt * pi_t,
      pi_a * s_at, pi_c * s_ct, pi_g * s_gt, -(pi_a * s_at + pi_c * s_ct + pi_g * s_gt);
  return q;
}

inline Eigen::Matrix4d gtr4_q_matrix(const std::array<double, 4>& pi,
                                     const std::array<double, 6>& s) {
  return gtr4_q_matrix(pi[0], pi[1], pi[2], pi[3], s[0], s[1], s[2], s[3], s[4], s[5]);
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
  return 2.0 * (pi_a * pi_c * s_ac + pi_a * pi_g * s_ag + pi_a * pi_t * s_at + pi_c * pi_g * s_cg +
                pi_c * pi_t * s_ct + pi_g * pi_t * s_gt);
}

inline double gtr4_rate(const std::array<double, 4>& pi, const std::array<double, 6>& s) {
  return gtr4_rate(pi[0], pi[1], pi[2], pi[3], s[0], s[1], s[2], s[3], s[4], s[5]);
}

inline std::pair<double, std::array<double, 6>> normalize_gtr_params(
    const std::array<double, 4>& pi,
    const std::array<double, 6>& s) {
  const double rate = gtr4_rate(pi, s);
  if (rate == 0.0)
    throw std::invalid_argument("Zero rate GTR parameters");

  return std::make_pair(
      rate,
      std::array<double, 6>{
          s[0] / rate, s[1] / rate, s[2] / rate, s[3] / rate, s[4] / rate, s[5] / rate});
}

inline std::tuple<Eigen::Matrix4d, Eigen::Matrix4d, Eigen::Vector4d> compute_decomposition(
    const std::array<double, 4>& pis,
    const std::array<double, 6>& ss) {
  Eigen::Matrix4d sym_q = symmetrized_gtr4_q_matrix(pis, ss);

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> eigen_solver(sym_q);

  auto evals = eigen_solver.eigenvalues();
  auto evecs = eigen_solver.eigenvectors();

  Eigen::Matrix4d pi_sqrt =
      Eigen::Vector4d(sqrt(pis[0]), sqrt(pis[1]), sqrt(pis[2]), sqrt(pis[3])).asDiagonal();
  Eigen::Matrix4d pi_sqrt_inv =
      Eigen::Vector4d(
          1.0 / sqrt(pis[0]), 1.0 / sqrt(pis[1]), 1.0 / sqrt(pis[2]), 1.0 / sqrt(pis[3]))
          .asDiagonal();

  Eigen::Matrix4d left = pi_sqrt_inv * evecs;
  Eigen::Matrix4d right = evecs.transpose() * pi_sqrt;

  return std::make_tuple(left, right, evals);
}

template <size_t N>
inline std::array<Eigen::Matrix4d, N> compute_transition_matrices(
    Eigen::Matrix4d left,
    Eigen::Matrix4d right,
    Eigen::Vector4d evals,
    const std::array<double, N>& branch_lengths) {
  std::array<Eigen::Matrix4d, N> transition_matrices;
  for (size_t i = 0; i < N; ++i) {
    transition_matrices[i] =
        left * (evals.array() * branch_lengths[i]).exp().matrix().asDiagonal() * right;
  }
  return transition_matrices;
}

template <size_t N>
inline std::array<Eigen::Matrix4d, N> compute_transition_matrices(
    const std::array<double, 4>& pis,
    const std::array<double, 6>& ss,
    const std::array<double, N>& branch_lengths) {
  std::array<Eigen::Matrix4d, N> transition_matrices;
  for (size_t i = 0; i < N; ++i) {
    transition_matrices[i] = (gtr4_q_matrix(pis, ss) * branch_lengths[i]).exp();
  }
  return transition_matrices;
}

/* ********** ********** ********** ********** ********** ********** ********** */
/* GTR gene log likelihood computations */
/* ********** ********** ********** ********** ********** ********** ********** */

inline double compute_gene_gtr_log_likelihood(const std::array<double, 4>& pis,
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

  auto log_pis = Eigen::Vector4d(pis[0], pis[1], pis[2], pis[3]).array().max(1e-10).log().matrix();

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
            lse_summands[i++] =
                log_p_t1(x, a) + log_p_t2(x, b) + log_p_t3(y, c) + log_p_t4(y, x) + log_pis(y);
          }
        }

        log_likelihood += site_pattern_count * log_sum_exp(lse_summands);
      }
    }
  }

  return log_likelihood;
}

inline double compute_gene_gtr_log_likelihood(const std::array<double, 4>& pis,
                                              Eigen::Matrix4d left,
                                              Eigen::Matrix4d right,
                                              Eigen::Vector4d evals,
                                              const std::array<double, 4>& branch_lengths,
                                              const TuplePatterns& site_patterns,
                                              const std::array<int, 3>& taxon_order) {
  const auto transition_mats = compute_transition_matrices(left, right, evals, branch_lengths);
  return compute_gene_gtr_log_likelihood(pis,
                                         transition_mats[0],
                                         transition_mats[1],
                                         transition_mats[2],
                                         transition_mats[3],
                                         site_patterns,
                                         taxon_order);
}

inline double compute_gene_gtr_log_likelihood(const std::array<double, 4>& pis,
                                              const std::array<double, 6>& ss,
                                              const std::array<double, 4>& branch_lengths,
                                              const TuplePatterns& site_patterns,
                                              const std::array<int, 3>& taxon_order) {
  const auto transition_mats = compute_transition_matrices(pis, ss, branch_lengths);
  return compute_gene_gtr_log_likelihood(pis,
                                         transition_mats[0],
                                         transition_mats[1],
                                         transition_mats[2],
                                         transition_mats[3],
                                         site_patterns,
                                         taxon_order);
}

inline double gtr_neg_log_likelihood(const TuplePatterns& site_patterns,
                                     const std::array<double, 4>& pis,
                                     const std::array<double, 6>& ss,
                                     const std::array<double, 4>& branch_lengths,
                                     const std::array<int, 3>& taxon_order) {
  return -compute_gene_gtr_log_likelihood(pis, ss, branch_lengths, site_patterns, taxon_order);
}

namespace detail {

inline std::array<double, 64> site_pattern_probs_from_transition_mats(
    const std::array<Eigen::Matrix4d, 2>& transition_mats) {
  std::array<double, 64> site_pattern_probs;
  site_pattern_probs.fill(0.0);
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      for (size_t k = 0; k < 4; ++k) {
        const size_t idx = site_pattern_index(i, j, k);
        for (size_t n = 0; n < 4; ++n) {
          site_pattern_probs[idx] +=
              transition_mats[0](n, i) * transition_mats[0](n, j) * transition_mats[1](n, k) / 4.0;
        }
      }
    }
  }
  return site_pattern_probs;
}

}  // namespace detail

inline std::array<double, 64> compute_gtr_triplet_gene_density(
    const std::array<double, 4>& pis,
    const std::array<double, 6>& ss,
    const std::array<double, 2>& t_branch_lengths) {
  // ultrametric tree. first branch length is from the cherries to their
  // bifurcation, the second is from the outgroup to the root.
  return detail::site_pattern_probs_from_transition_mats(
      compute_transition_matrices(pis, ss, t_branch_lengths));
}

inline std::array<double, 64> compute_gtr_triplet_gene_density(
    const Eigen::Matrix4d left,
    const Eigen::Matrix4d right,
    const Eigen::Vector4d evals,
    const std::array<double, 2>& t_branch_lengths) {
  // ultrametric tree. first branch length is from the cherries to their
  // bifurcation, the second is from the outgroup to the root.
  return detail::site_pattern_probs_from_transition_mats(
      compute_transition_matrices(left, right, evals, t_branch_lengths));
}

/* ********** ********** ********** ********** ********** ********** ********** */
/* GTR species log likelihood computations */
/* ********** ********** ********** ********** ********** ********** ********** */

// inline double f1(const double theta, double t1, double t2) {
//   return (4 / (theta * theta)) * exp((-2 / theta) * (t1 + t2));
// }

inline double f1(const double theta, double t1, double t2, double tau1, double tau2) {
  return (4 / (theta * theta)) * exp((-2 / theta) * (t1 - tau1 + t2 - tau2));
}

// inline double f2(const double theta, double t1, double t2, double tau1, double tau2) {
//   return (4 / (theta * theta)) * exp((-2 / theta) * (tau2 - tau1 + t2 + 2 * t1));
// }

inline double f2(const double theta, double t1, double t2, double tau1, double tau2) {
  return (4 / (theta * theta)) * exp((-2 / theta) * (2 * t1 + t2 - 2 * tau2 - tau1));
}

namespace detail {

template <typename DensityFn>
inline std::array<double, 64> first_density_integral_impl(
    DensityFn density_fn,
    const double theta,
    const std::array<double, 2>& tau_branch_lengths,
    const double max_branch_length) {
  std::array<double, 64> gene_site_pattern_probs;
  std::array<double, 64> compensation;  // for Kahan summation
  gene_site_pattern_probs.fill(0.0);
  compensation.fill(0.0);

  const double tau1 = tau_branch_lengths[0];
  const double tau2 = tau_branch_lengths[1];

  const double h1 = std::max(0.01, (tau2 - tau1) * 0.01);
  const double h2 = 0.01;

  for (double t2 = tau2; t2 < max_branch_length; t2 += h2) {
    for (double t1 = tau1; t1 < tau2; t1 += h1) {
      const auto gene_probs = density_fn(t1, t2);
      const double weight =
          f1(theta, t1, t2, tau1, tau2) * ((t1 + h1) > tau2 ? (tau2 - t1) : h1) * h2;
      for (size_t i = 0; i < 64; ++i) {
        const double term = gene_probs[i] * weight - compensation[i];
        const double new_sum = gene_site_pattern_probs[i] + term;
        compensation[i] = (new_sum - gene_site_pattern_probs[i]) - term;
        gene_site_pattern_probs[i] = new_sum;
      }
    }
  }
  return gene_site_pattern_probs;
}

}  // namespace detail

inline std::array<double, 64> first_density_integral(
    const std::array<double, 4>& pis,
    const std::array<double, 6>& ss,
    const double theta,
    const std::array<double, 2>& tau_branch_lengths,
    const double max_branch_length = 20.0) {
  return detail::first_density_integral_impl(
      [&](double t1, double t2) { return compute_gtr_triplet_gene_density(pis, ss, {t1, t2}); },
      theta,
      tau_branch_lengths,
      max_branch_length);
}

inline std::array<double, 64> first_density_integral(
    const Eigen::Matrix4d left,
    const Eigen::Matrix4d right,
    const Eigen::Vector4d evals,
    const double theta,
    const std::array<double, 2>& tau_branch_lengths,
    const double max_branch_length = 20.0) {
  return detail::first_density_integral_impl(
      [&](double t1, double t2) {
        return compute_gtr_triplet_gene_density(left, right, evals, {t1, t2});
      },
      theta,
      tau_branch_lengths,
      max_branch_length);
}

namespace detail {

template <typename DensityFn>
inline std::array<double, 64> first_density_integral_derivative_t1_fixed_cpt_impl(
    DensityFn density_fn,
    const double theta,
    const std::array<double, 2>& tau_branch_lengths,
    const double max_branch_length) {
  std::array<double, 64> partial_species_site_pattern_prob_deriv;
  std::array<double, 64> compensation;  // for Kahan summation
  partial_species_site_pattern_prob_deriv.fill(0.0);
  compensation.fill(0.0);

  const double tau1 = tau_branch_lengths[0];
  const double tau2 = tau_branch_lengths[1];

  const double h2 = 0.01;
  
  for (double t2 = tau2; t2 < max_branch_length; t2 += h2) {
    const auto gene_probs = density_fn(tau2, t2);
    const double weight = f1(theta, tau2, t2, tau1, tau2) * h2;
    for (size_t i = 0; i < 64; ++i) {
      const double term = gene_probs[i] * weight - compensation[i];
      const double new_sum = partial_species_site_pattern_prob_deriv[i] + term;
      compensation[i] = (new_sum - partial_species_site_pattern_prob_deriv[i]) - term;
      partial_species_site_pattern_prob_deriv[i] = new_sum;
    }
  }
  return partial_species_site_pattern_prob_deriv;
}

}  // namespace detail

inline std::array<double, 64> first_density_integral_derivative_t1_fixed_cpt(
    const std::array<double, 4>& pis,
    const std::array<double, 6>& ss,
    const double theta,
    const std::array<double, 2>& tau_branch_lengths,
    const double max_branch_length = 20.0) {
  return detail::first_density_integral_derivative_t1_fixed_cpt_impl(
      [&](double t1, double t2) { return compute_gtr_triplet_gene_density(pis, ss, {t1, t2}); },
      theta,
      tau_branch_lengths,
      max_branch_length);
}

inline std::array<double, 64> first_density_integral_derivative_t1_fixed_cpt(
    const Eigen::Matrix4d& left,
    const Eigen::Matrix4d& right,
    const Eigen::Vector4d& evals,
    const double theta,
    const std::array<double, 2>& tau_branch_lengths,
    const double max_branch_length = 20.0) {
  return detail::first_density_integral_derivative_t1_fixed_cpt_impl(
      [&](double t1, double t2) {
        return compute_gtr_triplet_gene_density(left, right, evals, {t1, t2});
      },
      theta,
      tau_branch_lengths,
      max_branch_length);
}


namespace detail {

template <typename DensityFn>
inline std::array<double, 64> first_density_integral_derivative_t2_fixed_cpt_impl(
    DensityFn density_fn,
    const double theta,
    const std::array<double, 2>& tau_branch_lengths,
    const double max_branch_length) {
  std::array<double, 64> partial_species_site_pattern_prob_deriv;
  std::array<double, 64> compensation;  // for Kahan summation
  partial_species_site_pattern_prob_deriv.fill(0.0);
  compensation.fill(0.0);

  const double tau1 = tau_branch_lengths[0];
  const double tau2 = tau_branch_lengths[1];

  const double h1 = 0.01;
  
  for (double t1 = tau1; t1 < tau2; t1 += h1) {
    const auto gene_probs = density_fn(t1, tau2);
    const double weight = f1(theta, t1, tau2, tau1, tau2) * std::min(h1, tau2 - t1);
    for (size_t i = 0; i < 64; ++i) {
      const double term = gene_probs[i] * weight - compensation[i];
      const double new_sum = partial_species_site_pattern_prob_deriv[i] + term;
      compensation[i] = (new_sum - partial_species_site_pattern_prob_deriv[i]) - term;
      partial_species_site_pattern_prob_deriv[i] = new_sum;
    }
  }
  return partial_species_site_pattern_prob_deriv;
}

}  // namespace detail

inline std::array<double, 64> first_density_integral_derivative_t2_fixed_cpt(
    const std::array<double, 4>& pis,
    const std::array<double, 6>& ss,
    const double theta,
    const std::array<double, 2>& tau_branch_lengths,
    const double max_branch_length = 20.0) {
  return detail::first_density_integral_derivative_t2_fixed_cpt_impl(
      [&](double t1, double t2) { return compute_gtr_triplet_gene_density(pis, ss, {t1, t2}); },
      theta,
      tau_branch_lengths,
      max_branch_length);
}

inline std::array<double, 64> first_density_integral_derivative_t2_fixed_cpt(
    const Eigen::Matrix4d& left,
    const Eigen::Matrix4d& right,
    const Eigen::Vector4d& evals,
    const double theta,
    const std::array<double, 2>& tau_branch_lengths,
    const double max_branch_length = 20.0) {
  return detail::first_density_integral_derivative_t2_fixed_cpt_impl(
      [&](double t1, double t2) {
        return compute_gtr_triplet_gene_density(left, right, evals, {t1, t2});
      },
      theta,
      tau_branch_lengths,
      max_branch_length);
}



namespace detail {

template <typename DensityFn>
inline std::array<double, 64> second_density_integral_impl(
    DensityFn density_fn,
    const double theta,
    const std::array<double, 2>& tau_branch_lengths,
    const double max_branch_length) {
  std::array<double, 64> partial_species_site_pattern_prob;
  std::array<double, 64> compensation;  // for Kahan summation
  partial_species_site_pattern_prob.fill(0.0);
  compensation.fill(0.0);

  const double h1 = 0.01;
  const double h2 = 0.01;
  const double tau1 = tau_branch_lengths[0];
  const double tau2 = tau_branch_lengths[1];

  for (double t2 = 0.0; t2 < max_branch_length; t2 += h2) {
    for (double t1 = 0.0; t1 < t2; t1 += h1) {
      const auto gene_probs = density_fn(t1, t2);
      const double weight = f2(theta, t1, t2, tau1, tau2) * std::min(h1, t2 - t1) * h2;
      for (size_t a = 0; a < 4; ++a) {
        for (size_t b = 0; b < 4; ++b) {
          for (size_t c = 0; c < 4; ++c) {
            const size_t i = site_pattern_index(a, b, c);
            const size_t j = site_pattern_index(a, c, b);
            const size_t k = site_pattern_index(b, c, a);

            const double term = (gene_probs[i] + gene_probs[j] + gene_probs[k]) * weight - compensation[i];
            const double new_sum = partial_species_site_pattern_prob[i] + term;
            compensation[i] = (new_sum - partial_species_site_pattern_prob[i]) - term;
            partial_species_site_pattern_prob[i] = new_sum;
          }
        } 
      }
    }
  }
  return partial_species_site_pattern_prob;
}

}  // namespace detail

inline std::array<double, 64> second_density_integral(
    const std::array<double, 4>& pis,
    const std::array<double, 6>& ss,
    const double theta,
    const std::array<double, 2>& tau_branch_lengths,
    const double max_branch_length = 20.0) {
  return detail::second_density_integral_impl(
      [&](double t1, double t2) { return compute_gtr_triplet_gene_density(pis, ss, {t1, t2}); },
      theta,
      tau_branch_lengths,
      max_branch_length);
}

inline std::array<double, 64> second_density_integral(
    const Eigen::Matrix4d& left,
    const Eigen::Matrix4d& right,
    const Eigen::Vector4d& evals,
    const double theta,
    const std::array<double, 2>& tau_branch_lengths,
    const double max_branch_length = 20.0) {
  return detail::second_density_integral_impl(
      [&](double t1, double t2) {
        return compute_gtr_triplet_gene_density(left, right, evals, {t1, t2});
      },
      theta,
      tau_branch_lengths,
      max_branch_length);
}



namespace detail {

template <typename DensityFn>
inline std::array<double, 64> second_density_integral_derivative_t1_fixed_cpt_impl(
    DensityFn density_fn,
    const double theta,
    const std::array<double, 2>& tau_branch_lengths,
    const double max_branch_length) {
  std::array<double, 64> partial_species_site_pattern_prob_deriv;
  std::array<double, 64> compensation;  // for Kahan summation
  partial_species_site_pattern_prob_deriv.fill(0.0);
  compensation.fill(0.0);

  const double h2 = 0.01;
  const double tau1 = tau_branch_lengths[0];
  const double tau2 = tau_branch_lengths[1];

  for (double t2 = 0.0; t2 < max_branch_length; t2 += h2) {
    const auto gene_probs = density_fn(tau2, t2);
    const double weight = f2(theta, tau2, t2, tau1, tau2) * h2;
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = 0; b < 4; ++b) {
        for (size_t c = 0; c < 4; ++c) {
          const size_t i = site_pattern_index(a, b, c);
          const size_t j = site_pattern_index(a, c, b);
          const size_t k = site_pattern_index(b, c, a);

          const double term = (gene_probs[i] + gene_probs[j] + gene_probs[k]) * weight - compensation[i];
          const double new_sum = partial_species_site_pattern_prob_deriv[i] + term;
          compensation[i] = (new_sum - partial_species_site_pattern_prob_deriv[i]) - term;
          partial_species_site_pattern_prob_deriv[i] = new_sum;
        }
      } 
    }
  }
  return partial_species_site_pattern_prob_deriv;
}

}  // namespace detail

inline std::array<double, 64> second_density_integral_derivative_t1_fixed_cpt(
    const std::array<double, 4>& pis,
    const std::array<double, 6>& ss,
    const double theta,
    const std::array<double, 2>& tau_branch_lengths,
    const double max_branch_length = 20.0) {
  return detail::second_density_integral_derivative_t1_fixed_cpt_impl(
      [&](double t1, double t2) { return compute_gtr_triplet_gene_density(pis, ss, {t1, t2}); },
      theta,
      tau_branch_lengths,
      max_branch_length);
}

inline std::array<double, 64> second_density_integral_derivative_t1_fixed_cpt(
    const Eigen::Matrix4d& left,
    const Eigen::Matrix4d& right,
    const Eigen::Vector4d& evals,
    const double theta,
    const std::array<double, 2>& tau_branch_lengths,
    const double max_branch_length = 20.0) {
  return detail::second_density_integral_derivative_t1_fixed_cpt_impl(
      [&](double t1, double t2) {
        return compute_gtr_triplet_gene_density(left, right, evals, {t1, t2});
      },
      theta,
      tau_branch_lengths,
      max_branch_length);
}



inline std::array<double, 64> compute_gtr_triplet_species_site_pattern_probs(
    const std::array<double, 4>& pis,
    const std::array<double, 6>& ss,
    const double theta,
    const std::array<double, 2>& tau_branch_lengths) {
  auto first_integral = first_density_integral(pis, ss, theta, tau_branch_lengths);
  auto second_integral = second_density_integral(pis, ss, theta, tau_branch_lengths);

  std::array<double, 64> species_site_pattern_probs;
  for (size_t i = 0; i < 64; ++i) {
    species_site_pattern_probs[i] = first_integral[i] + second_integral[i];
  }
  return species_site_pattern_probs;
}

namespace detail {

inline std::tuple<std::array<double, 64>, std::array<double, 64>, std::array<double, 64>>
combine_tau_deriv_integrals(const std::array<double, 64>& first_density_integral,
                            const std::array<double, 64>& first_di_deriv_t1_fixed_cpt,
                            const std::array<double, 64>& first_di_deriv_t2_fixed_cpt,
                            const std::array<double, 64>& second_density_integral,
                            const std::array<double, 64>& second_di_deriv_t1_fixed_cpt,
                            const double theta) {
  std::array<double, 64> probs;
  std::array<double, 64> tau1_deriv;
  std::array<double, 64> tau2_deriv;
  const double inv_theta_2 = 2.0 / theta;
  for (size_t i = 0; i < 64; ++i) {
    probs[i] = first_density_integral[i] + second_density_integral[i];
    tau1_deriv[i] = inv_theta_2 * first_density_integral[i] -
                    2 * inv_theta_2 * second_density_integral[i] - first_di_deriv_t1_fixed_cpt[i];
    tau2_deriv[i] = inv_theta_2 * first_density_integral[i] +
                    2 * inv_theta_2 * second_density_integral[i] + first_di_deriv_t1_fixed_cpt[i] -
                    first_di_deriv_t2_fixed_cpt[i] - second_di_deriv_t1_fixed_cpt[i];
  }
  return std::make_tuple(probs, tau1_deriv, tau2_deriv);
}

}  // namespace detail

inline std::tuple<std::array<double, 64>, std::array<double, 64>, std::array<double, 64>>
compute_gtr_triplet_species_site_pattern_probs_with_tau_derivs(
    const std::array<double, 4>& pis,
    const std::array<double, 6>& ss,
    const double theta,
    const std::array<double, 2>& tau_branch_lengths) {
  return detail::combine_tau_deriv_integrals(
      first_density_integral(pis, ss, theta, tau_branch_lengths),
      first_density_integral_derivative_t1_fixed_cpt(pis, ss, theta, tau_branch_lengths),
      first_density_integral_derivative_t2_fixed_cpt(pis, ss, theta, tau_branch_lengths),
      second_density_integral(pis, ss, theta, tau_branch_lengths),
      second_density_integral_derivative_t1_fixed_cpt(pis, ss, theta, tau_branch_lengths),
      theta);
}

inline std::tuple<std::array<double, 64>, std::array<double, 64>, std::array<double, 64>>
compute_gtr_triplet_species_site_pattern_probs_with_tau_derivs(
    const Eigen::Matrix4d& left,
    const Eigen::Matrix4d& right,
    const Eigen::Vector4d& evals,
    const double theta,
    const std::array<double, 2>& tau_branch_lengths) {
  return detail::combine_tau_deriv_integrals(
      first_density_integral(left, right, evals, theta, tau_branch_lengths),
      first_density_integral_derivative_t1_fixed_cpt(left, right, evals, theta, tau_branch_lengths),
      first_density_integral_derivative_t2_fixed_cpt(left, right, evals, theta, tau_branch_lengths),
      second_density_integral(left, right, evals, theta, tau_branch_lengths),
      second_density_integral_derivative_t1_fixed_cpt(left, right, evals, theta, tau_branch_lengths),
      theta);
}

inline double compute_gtr_triplet_log_likelihood(const std::array<double, 64>& site_patterns,
                                                 const std::array<double, 4>& pis,
                                                 const std::array<double, 6>& ss,
                                                 const double theta,
                                                 const std::array<double, 2>& tau_branch_lengths) {
  auto species_site_pattern_probs =
      compute_gtr_triplet_species_site_pattern_probs(pis, ss, theta, tau_branch_lengths);

  double log_likelihood = 0.0;
  for (size_t i = 0; i < 64; ++i) {
    if (site_patterns[i] > 0.0) {
      log_likelihood += site_patterns[i] * log(species_site_pattern_probs[i]);
    }
  }
  return log_likelihood;
}

namespace detail {

inline std::array<double, 3> accumulate_ll_with_tau_derivs(
    const std::array<double, 64>& site_patterns,
    const std::array<double, 64>& probs,
    const std::array<double, 64>& tau1_deriv_probs,
    const std::array<double, 64>& tau2_deriv_probs) {
  double log_likelihood = 0.0;
  double tau1_deriv = 0.0;
  double tau2_deriv = 0.0;
  for (size_t i = 0; i < 64; ++i) {
    if (site_patterns[i] > 0.0) {
      log_likelihood += site_patterns[i] * log(probs[i]);
      if (probs[i] > 0.0) {
        // logarithmic derivatives
        tau1_deriv += site_patterns[i] * tau1_deriv_probs[i] / probs[i];
        tau2_deriv += site_patterns[i] * tau2_deriv_probs[i] / probs[i];
      }
    }
  }
  return {log_likelihood, tau1_deriv, tau2_deriv};
}

}  // namespace detail

inline std::array<double, 3> compute_gtr_triplet_log_likelihood_with_tau_derivs(
    const std::array<double, 64>& site_patterns,
    const std::array<double, 4>& pis,
    const std::array<double, 6>& ss,
    const double theta,
    const std::array<double, 2>& tau_branch_lengths) {
  std::array<double, 64> probs, tau1_deriv_probs, tau2_deriv_probs;
  std::tie(probs, tau1_deriv_probs, tau2_deriv_probs) =
      compute_gtr_triplet_species_site_pattern_probs_with_tau_derivs(
          pis, ss, theta, tau_branch_lengths);
  return detail::accumulate_ll_with_tau_derivs(
      site_patterns, probs, tau1_deriv_probs, tau2_deriv_probs);
}

inline std::array<double, 3> compute_gtr_triplet_log_likelihood_with_tau_derivs(
    const std::array<double, 64>& site_patterns,
    const Eigen::Matrix4d& left,
    const Eigen::Matrix4d& right,
    const Eigen::Vector4d& evals,
    const double theta,
    const std::array<double, 2>& tau_branch_lengths) {
  std::array<double, 64> probs, tau1_deriv_probs, tau2_deriv_probs;
  std::tie(probs, tau1_deriv_probs, tau2_deriv_probs) =
      compute_gtr_triplet_species_site_pattern_probs_with_tau_derivs(
          left, right, evals, theta, tau_branch_lengths);
  return detail::accumulate_ll_with_tau_derivs(
      site_patterns, probs, tau1_deriv_probs, tau2_deriv_probs);
}

class SpeciesTripletGTRNegLogLikelihoodJointObjective {
 private:
  const std::array<double, 64>& site_patterns_;
  const std::array<double, 4>& pis_;
  const double theta_;
  const double mu_;

 public:
  SpeciesTripletGTRNegLogLikelihoodJointObjective(const std::array<double, 64>& site_patterns,
                                                  const std::array<double, 4>& pis,
                                                  const double theta,
                                                  const double mu)
      : site_patterns_(site_patterns), pis_(pis), theta_(theta), mu_(mu) {}

  double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
    const double tau1 = x[0];
    const double tau2 = x[0] + x[1];
    const std::array<double, 6> ss = {x[2], x[3], x[4], x[5], x[6], x[7]};

    double rate;
    std::array<double, 6> normalized_ss;
    std::tie(rate, normalized_ss) = normalize_gtr_params(pis_, ss);

    // const double fx = -compute_gtr_triplet_log_likelihood(
    //                       site_patterns_, pis_, normalized_ss, theta_, {tau1, tau2}) +
    //                   pow(rate - mu_, 2);

    auto fx_and_derivs = compute_gtr_triplet_log_likelihood_with_tau_derivs(
        site_patterns_, pis_, normalized_ss, theta_, {tau1, tau2});
    const double fx = -fx_and_derivs[0] + pow(rate - mu_, 2);
    const double fx_tau1_deriv = -fx_and_derivs[1];
    const double fx_tau2_deriv = -fx_and_derivs[2];

    std::cout << "Current parameters: tau1=" << tau1 << ", tau2=" << tau2
              << ", ss=" << normalized_ss[0] << "," << normalized_ss[1] << "," << normalized_ss[2]
              << "," << normalized_ss[3] << "," << normalized_ss[4] << "," << normalized_ss[5]
              << " with neg log likelihood ";
    std::cout.setf(std::ios::fixed);
    std::cout.precision(10);
    std::cout << (-fx) << std::endl;

    grad.resize(8);

    // as tau1 = x[0] and tau2 = x[0] + x[1]
    grad[0] = fx_tau1_deriv + fx_tau2_deriv;
    grad[1] = fx_tau2_deriv;

    std::cout << "Gradient w.r.t. tau1: " << fx_tau1_deriv << std::endl;
    std::cout << "Gradient w.r.t. x[0] (tau1&tau2): " << grad[0] << std::endl;
    std::cout << "Gradient w.r.t. x[1] (tau2): " << grad[1] << std::endl;

    for (size_t i = 0; i < 6; ++i) {
      std::array<double, 6> ss_plus = ss;
      ss_plus[i] += 0.001;

      double ss_plus_rate;
      std::array<double, 6> normalized_ss_plus;
      std::tie(ss_plus_rate, normalized_ss_plus) = normalize_gtr_params(pis_, ss_plus);

      const auto fx_ss_plus = -compute_gtr_triplet_log_likelihood(
                                  site_patterns_, pis_, normalized_ss_plus, theta_, {tau1, tau2}) +
                              pow(ss_plus_rate - mu_, 2);

      std::array<double, 6> ss_minus = ss;
      ss_minus[i] = std::max(0.0, ss_minus[i] - 0.001);

      double ss_minus_rate;
      std::array<double, 6> normalized_ss_minus;
      std::tie(ss_minus_rate, normalized_ss_minus) = normalize_gtr_params(pis_, ss_minus);

      const auto fx_ss_minus =
          -compute_gtr_triplet_log_likelihood(
              site_patterns_, pis_, normalized_ss_minus, theta_, {tau1, tau2}) +
          pow(ss_minus_rate - mu_, 2);

      grad[2 + i] = (fx_ss_plus - fx_ss_minus) / (ss_plus[i] - ss_minus[i]);
    }

    return fx;
  }
};

class SpeciesTripletGTRNegLogLikelihoodBranchLengthObjective {
 private:
  const std::array<double, 64>& site_patterns_;
  const std::array<double, 4>& pis_;
  const std::array<double, 6>& ss_;
  std::array<double, 6> normalized_ss_;
  Eigen::Matrix4d left_;
  Eigen::Matrix4d right_;
  Eigen::Vector4d evals_;
  const double theta_;
  const double mu_;

 public:
  SpeciesTripletGTRNegLogLikelihoodBranchLengthObjective(
      const std::array<double, 64>& site_patterns,
      const std::array<double, 4>& pis,
      const std::array<double, 6>& ss,
      const double theta,
      const double mu)
      : site_patterns_(site_patterns), pis_(pis), ss_(ss), theta_(theta), mu_(mu) {
    double rate;
    std::array<double, 6> normalized_ss;
    std::tie(rate, normalized_ss) = normalize_gtr_params(pis_, ss);
    normalized_ss_ = normalized_ss;

    std::tie(left_, right_, evals_) = compute_decomposition(pis_, normalized_ss_);
  }

  double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
    const double tau1 = x[0];
    const double tau2 = x[0] + x[1];

    auto fx_and_derivs = compute_gtr_triplet_log_likelihood_with_tau_derivs(
        site_patterns_, left_, right_, evals_, theta_, {tau1, tau2});
    const double fx = -fx_and_derivs[0];
    const double fx_tau1_deriv = -fx_and_derivs[1];
    const double fx_tau2_deriv = -fx_and_derivs[2];

    std::cout.setf(std::ios::fixed);
    std::cout.precision(10);
    std::cout << "Current parameters: tau1=" << tau1 << ", tau2=" << tau2
              << " with neg log likelihood ";
    std::cout << (-fx) << std::endl;

    grad.resize(2);

    // as tau1 = x[0] and tau2 = x[0] + x[1] (and it's a _negative_ log likelihood)
    grad[0] = -fx_tau1_deriv - fx_tau2_deriv;
    grad[1] = -fx_tau2_deriv;

    std::cout << "Gradient w.r.t. x[0] (tau1&tau2): " << grad[0] << std::endl;
    std::cout << "Gradient w.r.t. x[1] (tau2): " << grad[1] << std::endl;

    return fx;
  }
};

inline std::array<double, 64> permute_triplet_site_patterns(const TuplePatterns& site_patterns,
                                                            size_t i) {
  // i gives:
  // 0 -> (a,b)|c,
  // 1 -> (a,c)|b,
  // 2 -> (b,c)|a
  std::array<double, 64> permuted_site_patterns;
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = 0; b < 4; ++b) {
      for (size_t c = 0; c < 4; ++c) {
        size_t idx;
        if (i == 0) {
          idx = site_pattern_index(a, b, c);
        } else if (i == 1) {
          idx = site_pattern_index(a, c, b);
        } else if (i == 2) {
          idx = site_pattern_index(b, c, a);
        } else {
          throw std::invalid_argument("Invalid topology index");
        }
        permuted_site_patterns[site_pattern_index(a, b, c)] = site_patterns.site_pattern[idx];
      }
    }
  }
  return permuted_site_patterns;
}

/**
 * Compute GTR triplet weights from pre-computed site patterns.
 *
 * The three weights correspond to topologies where a, b, c are the sorted
 * taxon names (matching site_patterns.taxa order):
 *   result[0] -- (a,b)|c
 *   result[1] -- (a,c)|b
 *   result[2] -- (b,c)|a
 */
inline std::array<double, 3> compute_gtr_triplet_weights(const TuplePatterns& site_patterns,
                                                         const double theta,
                                                         const double mu = 1.0) {
  const double n = site_patterns.n_counted;
  if (n == 0.0) {
    // TODO: warn when no sites are counted.
    return {0.0, 0.0, 0.0};
  }

  auto calc_pi = [&site_patterns, n](size_t nuc) {
    double pi = 0.0;
    for (size_t i = 0; i < 4; ++i) {
      for (size_t j = 0; j < 4; ++j) {
        pi += site_patterns.site_pattern[site_pattern_index(nuc, i, j)];
        pi += site_patterns.site_pattern[site_pattern_index(i, nuc, j)];
        pi += site_patterns.site_pattern[site_pattern_index(i, j, nuc)];
      }
    }
    return pi / (3 * n);
  };
  const double pi_a = calc_pi(0);
  const double pi_c = calc_pi(1);
  const double pi_g = calc_pi(2);
  const double pi_t = calc_pi(3);
  const std::array<double, 4> pis = {pi_a, pi_c, pi_g, pi_t};

  std::array<double, 3> weights = {0.0, 0.0, 0.0};
  for (size_t i = 0; i < 3; ++i) {
    const std::array<double, 64> tree_permuted_site_patterns =
        permute_triplet_site_patterns(site_patterns, i);

    double tau1_init = 1.0;
    double tau2_init = 2.0;
    const std::array<double, 6> ss_init = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    double rate;
    std::array<double, 6> normalized_ss_init;
    std::tie(rate, normalized_ss_init) = normalize_gtr_params(pis, ss_init);

    const double init_neg_ll =
        -compute_gtr_triplet_log_likelihood(
            tree_permuted_site_patterns, pis, normalized_ss_init, theta, {tau1_init, tau2_init}) +
        pow(rate - mu, 2);
    double fx = init_neg_ll;

    /* ********** ********** ********** ********** */
    /* Initial Branch length optimization */
    /* ********** ********** ********** ********** */

    auto branch_length_objective = SpeciesTripletGTRNegLogLikelihoodBranchLengthObjective(
        tree_permuted_site_patterns, pis, normalized_ss_init, theta, mu);

    LBFGSpp::LBFGSBParam<double> param_branch_length;
    param_branch_length.epsilon = 1e-6;
    param_branch_length.max_iterations = 200;
    param_branch_length.past = 0;

    LBFGSpp::LBFGSBSolver<double> solver_branch_length(param_branch_length);

    Eigen::VectorXd lb_branch_length(2), ub_branch_length(2);
    lb_branch_length << 1e-6, 1e-6;
    ub_branch_length << 50.0, 50.0;  // TODO: determine reasonable upper bounds

    Eigen::VectorXd x_branch_length(2);
    x_branch_length << tau1_init, tau2_init - tau1_init;

    try {
      solver_branch_length.minimize(
          branch_length_objective, x_branch_length, fx, lb_branch_length, ub_branch_length);

      tau1_init = x_branch_length[0];
      tau2_init = x_branch_length[0] + x_branch_length[1];
    } catch (const std::exception&) {
      // If L-BFGS-B fails, keep the initial guess.
      if (fx < init_neg_ll) {
        PHYLOTRIPLET_WARN(
            "L-BFGS-B branch length optimization failed, but found a better solution. Using "
            "the L-BFGS-B solution.");
      } else {
        PHYLOTRIPLET_WARN("L-BFGS-B branch length optimization failed, using the initial guess.");
        fx = init_neg_ll;
      }
    }

    /* ********** ********** ********** ********** */
    /* Joint branch length and parameter optimization */
    /* ********** ********** ********** ********** */

    auto joint_objective = SpeciesTripletGTRNegLogLikelihoodJointObjective(
        tree_permuted_site_patterns, pis, theta, mu);

    LBFGSpp::LBFGSBParam<double> param_joint;
    param_joint.epsilon = 1e-6;
    param_joint.max_iterations = 200;
    param_joint.past = 0;

    LBFGSpp::LBFGSBSolver<double> solver_joint(param_joint);

    Eigen::VectorXd lb_joint(8),
        ub_joint(8);  // first two are for tau1 and tau2, the rest are for the ss parameters
    lb_joint << 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6;
    ub_joint << 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
        50.0;  // TODO: determine reasonable upper bounds

    Eigen::VectorXd x(8);
    x << tau1_init, tau2_init, ss_init[0], ss_init[1], ss_init[2], ss_init[3], ss_init[4],
        ss_init[5];

    try {
      solver_joint.minimize(joint_objective, x, fx, lb_joint, ub_joint);
    } catch (const std::exception&) {
      // If L-BFGS-B fails, keep the initial guess.
      if (fx < init_neg_ll) {
        PHYLOTRIPLET_WARN(
            "L-BFGS-B optimization failed, but found a better solution. Using "
            "the L-BFGS-B solution.");
      } else {
        PHYLOTRIPLET_WARN("L-BFGS-B optimization failed, using the initial guess.");
        fx = init_neg_ll;
      }
    }

    weights[i] = -fx;
  }

  return weights;
}

}  // namespace phylotriplet
