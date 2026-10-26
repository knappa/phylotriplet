#pragma once

#include <LBFGSB.h>

#include <Eigen/Core>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "sequence.hpp"
#include "util.hpp"
#include "warn.hpp"

namespace phylotriplet {

// ---------------------------------------------------------------------------
// Site pattern extraction
// ---------------------------------------------------------------------------

// Five JC69-equivalent site pattern classes for a 3-taxon alignment.
// Taxa are ordered as in TuplePatterns::taxa (lexicographic).
// Topology encoding: 0 = (taxa[0],taxa[1])|taxa[2]
//                    1 = (taxa[0],taxa[2])|taxa[1]
//                    2 = (taxa[1],taxa[2])|taxa[0]
struct GeneSitePatterns {
  double n_xxx;  // all three taxa identical
  double n_xxy;  // taxa[0]=taxa[1] != taxa[2]
  double n_xyx;  // taxa[0]=taxa[2] != taxa[1]
  double n_yxx;  // taxa[1]=taxa[2] != taxa[0]
  double n_xyz;  // all three distinct
  double n_total;
};

inline GeneSitePatterns extract_jc69_gene_site_patterns(const TuplePatterns& tp) {
  if (tp.order() != 3) {
    throw std::invalid_argument("extract_jc69_gene_site_patterns: expected 3-taxon TuplePatterns");
  }
  GeneSitePatterns g{};
  g.n_total = tp.n_counted;
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = 0; b < 4; ++b) {
      for (size_t c = 0; c < 4; ++c) {
        const double count = tp.site_pattern[site_pattern_index(a, b, c)];
        if (count == 0.0)
          continue;
        if (a == b && b == c) {
          g.n_xxx += count;
        } else if (a == b && b != c) {
          g.n_xxy += count;
        } else if (a == c && a != b) {
          g.n_xyx += count;
        } else if (b == c && a != b) {
          g.n_yxx += count;
        } else {
          g.n_xyz += count;
        }
      }
    }
  }
  return g;
}

// ---------------------------------------------------------------------------
// JC69 unrooted 3-leaf likelihood
// ---------------------------------------------------------------------------

// JC69 same/diff probabilities for branch length e.
// s(e) = (1 + 3*exp(-4e/3)) / 4
// d(e) = (1 - exp(-4e/3)) / 4
// s + 3*d = 1
inline void jc69_sd(const double e, double& s, double& d) {
  const double exp_term = std::exp(-4.0 * e / 3.0);
  s = (1.0 + 3.0 * exp_term) / 4.0;
  d = (1.0 - exp_term) / 4.0;
}

// Raw (un-normalized by 1/4) site pattern probabilities for the unrooted
// 3-leaf JC69 star tree with edge lengths (e_A, e_B, e_C).
// The 1/4 stationarity factor is a constant in log-space and is omitted.
struct RawPatternProbs {
  double p_xxx, p_xxy, p_xyx, p_yxx, p_xyz;
};

inline RawPatternProbs jc69_unrooted_probs(const double e_A, const double e_B, const double e_C) {
  double sA, dA, sB, dB, sC, dC;
  jc69_sd(e_A, sA, dA);
  jc69_sd(e_B, sB, dB);
  jc69_sd(e_C, sC, dC);

  RawPatternProbs p;
  p.p_xxx = sA * sB * sC + 3.0 * dA * dB * dC;
  p.p_xxy = 3.0 * sA * sB * dC + 3.0 * dA * dB * sC + 6.0 * dA * dB * dC;
  p.p_xyx = 3.0 * sA * dB * sC + 3.0 * dA * sB * dC + 6.0 * dA * dB * dC;
  p.p_yxx = 3.0 * dA * sB * sC + 3.0 * sA * dB * dC + 6.0 * dA * dB * dC;
  p.p_xyz = 6.0 * sA * dB * dC + 6.0 * dA * sB * dC + 6.0 * dA * dB * sC + 6.0 * dA * dB * dC;
  return p;
}

inline double jc69_unrooted_neg_loglik(const double e_A,
                                       const double e_B,
                                       const double e_C,
                                       const GeneSitePatterns& g) {
  const auto p = jc69_unrooted_probs(e_A, e_B, e_C);
  const double eps = std::numeric_limits<double>::min();
  return -(g.n_xxx * std::log(p.p_xxx + eps) + g.n_xxy * std::log(p.p_xxy + eps) +
           g.n_xyx * std::log(p.p_xyx + eps) + g.n_yxx * std::log(p.p_yxx + eps) +
           g.n_xyz * std::log(p.p_xyz + eps));
}

inline double logp_xxx(const double mu, const double theta) {
  const double x = mu * theta;
  const double log_x = std::log(mu) + std::log(theta);
  return -5 * log(2) +
         log_sum_exp(std::array<double, 8>{7 * log_x + log(54),
                                           6 * log_x + log(414),
                                           6 * log_x + log(1407),
                                           4 * log_x + log(2927),
                                           3 * log_x + log(4804),
                                           2 * log_x + log(5796),
                                           log_x + log(3872),
                                           log(1024)}) -
         log1p(3 * x + 3) - 2 * log1p(3 * x + 1) - log1p(x + 1) - 3 * log1p(x);
}

inline double logp_xxy(const double mu, const double theta) {
  const double x = mu * theta;
  const double log_x = std::log(mu) + std::log(theta);
  return std::log(3.0 / 32.0) +
         log_sum_exp(std::array<double, 6>{5 * log_x + std::log(54),
                                           4 * log_x + std::log(360),
                                           3 * log_x + std::log(1047),
                                           2 * log_x + std::log(1664),
                                           log_x + std::log(1412),
                                           std::log(480)}) -
         std::log(3 * x + 4) - 2 * std::log(3 * x + 2) - std::log(x + 2) - 2 * std::log(x + 1) +
         std::log(x);
}

inline double logp_xyx(const double mu, const double theta) {
  const double x = mu * theta;
  const double log_x = std::log(mu) + std::log(theta);
  // TODO(reinspect): the seven terms below span degrees 6..0 in log_x and the
  // array size matches that count, so the sum has no zero-initialized slot
  // silently adding exp(0) = 1.  Check the derivation to confirm that seven is
  // the intended number of terms and none is missing.
  return std::log(3.0 / 32.0) +
         log_sum_exp(std::array<double, 7>{6 * log_x + std::log(54),
                                           5 * log_x + std::log(414),
                                           4 * log_x + std::log(1299),
                                           3 * log_x + std::log(2207),
                                           2 * log_x + std::log(2260),
                                           log_x + std::log(1348),
                                           std::log(352)}) -
         std::log1p(3 * x + 3) - 2 * std::log1p(3 * x + 1) - std::log1p(x + 1) - 3 * std::log1p(x) +
         std::log(x);
}

inline double logp_yxx(const double mu, const double theta) {
  return logp_xyx(mu, theta);
}

inline double logp_xyz(const double mu, const double theta) {
  const double x = mu * theta;
  const double log_x = std::log(mu) + std::log(theta);
  return std::log(3.0 / 16.0) +
         log_sum_exp(std::array<double, 6>{5 * log_x + std::log(54),
                                           4 * log_x + std::log(414),
                                           3 * log_x + std::log(1299),
                                           2 * log_x + std::log(2099),
                                           log_x + std::log(1684),
                                           std::log(516)}) -
         std::log1p(3 * x + 3) - 2 * std::log1p(3 * x + 1) - std::log1p(x + 1) - 3 * std::log1p(x) +
         2 * std::log1p(x);
}

inline std::array<double, 3> compute_jc69_species_triplet_weights(
    const TuplePatterns& site_patterns,
    const double theta,
    const double mu = 4.0 / 3.0) {
  auto sps = extract_jc69_gene_site_patterns(site_patterns);

  const double lp_xxx = logp_xxx(mu, theta);
  const double lp_xxy = logp_xxy(mu, theta);
  const double lp_xyx = logp_xyx(mu, theta);
  const double lp_yxx = logp_yxx(mu, theta);
  const double lp_xyz = logp_xyz(mu, theta);

  // clang-format off
  // species tree (a,b)|c
  const double weight_0 = lp_xxx * sps.n_xxx 
                        + lp_xxy * sps.n_xxy 
						+ lp_xyx * sps.n_xyx 
						+ lp_yxx * sps.n_yxx 
						+ lp_xyz * sps.n_xyz;

  // species tree (a,c)|b
  const double weight_1 = lp_xxx * sps.n_xxx 
  					    + lp_xxy * sps.n_xyx
						+ lp_xyx * sps.n_xxy 
						+ lp_yxx * sps.n_yxx
						+ lp_xyz * sps.n_xyz;

  // species tree (b,c)|a
  const double weight_2 = lp_xxx * sps.n_xxx
                        + lp_xxy * sps.n_yxx
						+ lp_xyx * sps.n_xyx 
						+ lp_yxx * sps.n_xxy 
						+ lp_xyz * sps.n_xyz;
  // clang-format on

  const double norm = log_sum_exp(std::array<double, 3>{weight_0, weight_1, weight_2});
  return {weight_0 - norm, weight_1 - norm, weight_2 - norm};
}

}  // namespace phylotriplet
