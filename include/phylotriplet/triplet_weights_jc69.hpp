#pragma once

#include <LBFGSB.h>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#include "sequence.hpp"
#include "util.hpp"
#include "warn.hpp"

namespace phylotriplet {

inline std::array<double, 4> compute_jc69_site_pattern_log_probs(const double tau1,
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

  const double log_s_xxx =
      -2 * tau1 - 2 * tau2 - 6 * log(2) - log(2 * beta + 1) - log(beta + 1) +
      log_sum_exp(std::array<double, 8>{2 * tau1 + 2 * tau2 + log(2) + 2 * log(beta),
                                        2 * tau1 + 2 * tau2 + log(3) + log(beta),
                                        2 * tau1 + log(12) + 2 * log(beta),
                                        2 * tau2 + log(6) + 2 * log(beta),
                                        2 * tau1 + 2 * tau2,
                                        tau1 + log(12) + 2 * log(beta),
                                        2 * tau1 + log(6) + log(beta),
                                        2 * tau2 + log(3) + log(beta)});

  const double log_s_xxy =
      -2 * tau1 - 2 * tau2 - 6 * log(2) - log(2 * beta + 1) - log(beta + 1) +
      log_diff_exp(log_sum_exp(std::array<double, 5>{2 * tau1 + 2 * tau2 + log(2) + 2 * log(beta),
                                                     2 * tau1 + 2 * tau2 + log(3) + log(beta),
                                                     2 * tau2 + log(6) + 2 * log(beta),
                                                     2 * tau1 + 2 * tau2,
                                                     2 * tau2 + log(3) + log(beta)}),
                   log_sum_exp(std::array<double, 3>{2 * tau1 + 2 * log(2) + 2 * log(beta),
                                                     tau1 + 2 * log(2) + 2 * log(beta),
                                                     2 * tau1 + log(2) + log(beta)}));

  const double log_s_xyx =
      -2 * tau1 - 2 * tau2 - 6 * log(2) - log(2 * beta + 1) - log(beta + 1) +
      log_diff_exp(log_sum_exp(std::array<double, 5>{2 * tau1 + 2 * tau2 + log(2) + 2 * log(beta),
                                                     2 * tau1 + 2 * tau2 + log(3) + log(beta),
                                                     2 * tau1 + 2 * log(2) + 2 * log(beta),
                                                     2 * tau1 + 2 * tau2,
                                                     2 * tau1 + log(2) + log(beta)}),
                   log_sum_exp(std::array<double, 3>{2 * tau2 + log(2) + 2 * log(beta),
                                                     tau1 + 2 * log(2) + 2 * log(beta),
                                                     2 * tau2 + log(beta)}));

  const double log_s_xyz =
      -2 * tau1 - 2 * tau2 - 6 * log(2) - log(2 * beta + 1) - log(beta + 1) +
      log_diff_exp(log_sum_exp(std::array<double, 4>{2 * tau1 + 2 * tau2 + log(2) + 2 * log(beta),
                                                     2 * tau1 + 2 * tau2 + log(3) + log(beta),
                                                     2 * tau1 + 2 * tau2,
                                                     tau1 + 2 * log(2) + 2 * log(beta)}),
                   log_sum_exp(std::array<double, 4>{2 * tau1 + 2 * log(2) + 2 * log(beta),
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

    const auto site_pattern_weights = compute_jc69_site_pattern_log_probs(tau1, tau2, theta_, mu_);

    const double cxxx = site_patterns_[0];
    const double cxxy = site_patterns_[1];
    const double cxyx = site_patterns_[2];
    const double cxyz = site_patterns_[3];

    const double fx = -(site_pattern_weights[0] * cxxx + site_pattern_weights[1] * cxxy +
                        site_pattern_weights[2] * cxyx + site_pattern_weights[3] * cxyz);

    grad.resize(2);

    const double x0_plus = std::min(0.001, (tau2 - tau1) * 0.1);
    const auto site_pattern_weights_tau1_plus =
        compute_jc69_site_pattern_log_probs(tau1 + x0_plus, tau2 + x0_plus, theta_, mu_);

    const double x0_minus = std::min(0.001, tau1 * 0.1);
    const auto site_pattern_weights_tau1_minus =
        compute_jc69_site_pattern_log_probs(tau1 - x0_minus, tau2 - x0_minus, theta_, mu_);
    grad[0] = -((site_pattern_weights_tau1_plus[0] - site_pattern_weights_tau1_minus[0]) * cxxx +
                (site_pattern_weights_tau1_plus[1] - site_pattern_weights_tau1_minus[1]) * cxxy +
                (site_pattern_weights_tau1_plus[2] - site_pattern_weights_tau1_minus[2]) * cxyx +
                (site_pattern_weights_tau1_plus[3] - site_pattern_weights_tau1_minus[3]) * cxyz) /
              (x0_plus + x0_minus);

    const double x1_plus = 0.001;
    const auto site_pattern_weights_tau2_plus =
        compute_jc69_site_pattern_log_probs(tau1, tau2 + x1_plus, theta_, mu_);

    const double x1_minus = std::min(0.001, (tau2 - tau1) * 0.1);
    const auto site_pattern_weights_tau2_minus =
        compute_jc69_site_pattern_log_probs(tau1, tau2 - x1_minus, theta_, mu_);
    grad[1] = -((site_pattern_weights_tau2_plus[0] - site_pattern_weights_tau2_minus[0]) * cxxx +
                (site_pattern_weights_tau2_plus[1] - site_pattern_weights_tau2_minus[1]) * cxxy +
                (site_pattern_weights_tau2_plus[2] - site_pattern_weights_tau2_minus[2]) * cxyx +
                (site_pattern_weights_tau2_plus[3] - site_pattern_weights_tau2_minus[3]) * cxyz) /
              (x1_plus + x1_minus);

    return fx;
  }

 private:
  const std::array<double, 4>& site_patterns_;
  const double theta_;
  const double mu_;
};

/**
 * Compute JC69 triplet weights from pre-computed site patterns.
 *
 * The three weights correspond to topologies where a, b, c are the sorted
 * taxon names (matching site_patterns.taxa order):
 *   result[0] -- (a,b)|c
 *   result[1] -- (a,c)|b
 *   result[2] -- (b,c)|a
 */
inline std::array<double, 3> compute_jc69_triplet_weights(const TuplePatterns& site_patterns,
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
  std::array<double, 3> tau2 = {std::max(0.0,
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
                        (-3.0 / 8.0) * (log(3 * pxxx + 3 * pxxy[0] - pxyx[0] + pxyz) - log(n) -
                                        log_3 + log1p_43theta_term))),
      std::min(tau2[1],
               std::max(0.0,
                        (-3.0 / 8.0) * (log(3 * pxxx + 3 * pxxy[1] - pxyx[1] + pxyz) - log(n) -
                                        log_3 + log1p_43theta_term))),
      std::min(tau2[2],
               std::max(0.0,
                        (-3.0 / 8.0) * (log(3 * pxxx + 3 * pxxy[2] - pxyx[2] + pxyz) - log(n) -
                                        log_3 + log1p_43theta_term)))};

  std::array<double, 3> weights = {0.0, 0.0, 0.0};
  for (size_t i = 0; i < 3; ++i) {
    const double best_tau1_init = std::max(0.001, tau1[i]);
    const double best_tau2_init = std::max(best_tau1_init + 0.001, tau2[i]);

    const auto init_site_pattern_weights =
        compute_jc69_site_pattern_log_probs(best_tau1_init, best_tau2_init, theta, mu);

    const double init_neg_ll =
        -(init_site_pattern_weights[0] * cxxx + init_site_pattern_weights[1] * cxxy[i] +
          init_site_pattern_weights[2] * cxyx[i] + init_site_pattern_weights[3] * cxyz);

    auto objective = SpeciesTripletJC69NegLogLikelihood(
        std::array<double, 4>{cxxx, cxxy[i], cxyx[i], cxyz}, theta, mu);

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
        PHYLOTRIPLET_WARN(
            "L-BFGS-B optimization failed, but found a better solution. Using "
            "the L-BFGS-B solution.");
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

  std::cout << tau1[0] << " " << tau2[0] << " " << tau1[1] << " " << tau2[1] << " " << tau1[2]
            << " " << tau2[2] << std::endl;

  return weights;
}

}  // namespace phylotriplet