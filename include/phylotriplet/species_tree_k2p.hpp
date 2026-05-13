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
// Calculation for -log(P(site pattern)) under K2P + MSC
// depends on parameters, but not on the tree topology as this is integrated out
// ---------------------------------------------------------------------------
inline std::array<double, 64> k2p_pattern_probs(const double alpha,
                                                const double beta,
                                                const double theta) {
  // TODO: bound log probabilities away from -inf

  const double log_alpha = std::log(alpha);
  const double log_beta = std::log(beta);
  const double log_theta = std::log(theta);

  const double log_p000 =
      -6 * log(2) - 2 * log1p(2 * alpha * theta + 4 * beta * theta) -
      3 * log1p(2 * alpha * theta + 2 * beta * theta) - log1p(alpha * theta + 2 * beta * theta) -
      log1p(alpha * theta + beta * theta) - 3 * log1p(4 * beta * theta) - log1p(2 * beta * theta) +
      log_sum_exp(
          std::array<double, 68>{log(4096) + 7 * log_alpha + 4 * log_beta + 11 * log_theta,
                                 log(40960) + 6 * log_alpha + 5 * log_beta + 11 * log_theta,
                                 log(172032) + 5 * log_alpha + 6 * log_beta + 11 * log_theta,
                                 log(393216) + 4 * log_alpha + 7 * log_beta + 11 * log_theta,
                                 log(528384) + 3 * log_alpha + 8 * log_beta + 11 * log_theta,
                                 log(417792) + 2 * log_alpha + 9 * log_beta + 11 * log_theta,
                                 log(180224) + log_alpha + 10 * log_beta + 11 * log_theta,
                                 log(32768) + 11 * log_beta + 11 * log_theta,
                                 log(5120) + 7 * log_alpha + 3 * log_beta + 10 * log_theta,
                                 log(69632) + 6 * log_alpha + 4 * log_beta + 10 * log_theta,
                                 log(372736) + 5 * log_alpha + 5 * log_beta + 10 * log_theta,
                                 log(1042432) + 4 * log_alpha + 6 * log_beta + 10 * log_theta,
                                 log(1666048) + 3 * log_alpha + 7 * log_beta + 10 * log_theta,
                                 log(1533952) + 2 * log_alpha + 8 * log_beta + 10 * log_theta,
                                 log(757760) + log_alpha + 9 * log_beta + 10 * log_theta,
                                 log(155648) + 10 * log_beta + 10 * log_theta,
                                 log(2432) + 7 * log_alpha + 2 * log_beta + 9 * log_theta,
                                 log(47360) + 6 * log_alpha + 3 * log_beta + 9 * log_theta,
                                 log(335104) + 5 * log_alpha + 4 * log_beta + 9 * log_theta,
                                 log(1178112) + 4 * log_alpha + 5 * log_beta + 9 * log_theta,
                                 log(2287488) + 3 * log_alpha + 6 * log_beta + 9 * log_theta,
                                 log(2495744) + 2 * log_alpha + 7 * log_beta + 9 * log_theta,
                                 log(1433088) + log_alpha + 8 * log_beta + 9 * log_theta,
                                 log(336896) + 9 * log_beta + 9 * log_theta,
                                 log(608) + 7 * log_alpha + log_beta + 8 * log_theta,
                                 log(17024) + 6 * log_alpha + 2 * log_beta + 8 * log_theta,
                                 log(163968) + 5 * log_alpha + 3 * log_beta + 8 * log_theta,
                                 log(745920) + 4 * log_alpha + 4 * log_beta + 8 * log_theta,
                                 log(1805472) + 3 * log_alpha + 5 * log_beta + 8 * log_theta,
                                 log(2389696) + 2 * log_alpha + 6 * log_beta + 8 * log_theta,
                                 log(1630848) + log_alpha + 7 * log_beta + 8 * log_theta,
                                 log(448256) + 8 * log_beta + 8 * log_theta,
                                 log(128) + 7 * log_alpha + 7 * log_theta,
                                 log(4016) + 6 * log_alpha + log_beta + 7 * log_theta,
                                 log(50096) + 5 * log_alpha + 2 * log_beta + 7 * log_theta,
                                 log(297136) + 4 * log_alpha + 3 * log_beta + 7 * log_theta,
                                 log(918768) + 3 * log_alpha + 4 * log_beta + 7 * log_theta,
                                 log(1520224) + 2 * log_alpha + 5 * log_beta + 7 * log_theta,
                                 log(1276128) + log_alpha + 6 * log_beta + 7 * log_theta,
                                 log(427328) + 7 * log_beta + 7 * log_theta,
                                 log(576) + 6 * log_alpha + 6 * log_theta,
                                 log(10256) + 5 * log_alpha + log_beta + 6 * log_theta,
                                 log(79632) + 4 * log_alpha + 2 * log_beta + 6 * log_theta,
                                 log(320504) + 3 * log_alpha + 3 * log_beta + 6 * log_theta,
                                 log(684224) + 2 * log_alpha + 4 * log_beta + 6 * log_theta,
                                 log(732792) + log_alpha + 5 * log_beta + 6 * log_theta,
                                 log(310768) + 6 * log_beta + 6 * log_theta,
                                 log(1104) + 5 * log_alpha + 5 * log_theta,
                                 log(13944) + 4 * log_alpha + log_beta + 5 * log_theta,
                                 log(77224) + 3 * log_alpha + 2 * log_beta + 5 * log_theta,
                                 log(220844) + 2 * log_alpha + 3 * log_beta + 5 * log_theta,
                                 log(311372) + log_alpha + 4 * log_beta + 5 * log_theta,
                                 log(170640) + 5 * log_beta + 5 * log_theta,
                                 log(1192) + 4 * log_alpha + 4 * log_theta,
                                 log(11726) + 3 * log_alpha + log_beta + 4 * log_theta,
                                 log(48844) + 2 * log_alpha + 2 * log_beta + 4 * log_theta,
                                 log(94402) + log_alpha + 3 * log_beta + 4 * log_theta,
                                 log(67896) + 4 * log_beta + 4 * log_theta,
                                 log(828) + 3 * log_alpha + 3 * log_theta,
                                 log(6547) + 2 * log_alpha + log_beta + 3 * log_theta,
                                 log(19027) + log_alpha + 2 * log_beta + 3 * log_theta,
                                 log(18672) + 3 * log_beta + 3 * log_theta,
                                 log(394) + 2 * log_alpha + 2 * log_theta,
                                 log(2256) + log_alpha + log_beta + 2 * log_theta,
                                 log(3339) + 2 * log_beta + 2 * log_theta,
                                 log(118) + log_alpha + log_theta,
                                 log(348) + log_beta + log_theta,
                                 log(16)});

  const double log_p001 = -log(192) - 2 * log1p(2 * alpha * theta + 4 * beta * theta) -
                          3 * log1p(2 * alpha * theta + 2 * beta * theta) -
                          log1p(alpha * theta + 2 * beta * theta) -
                          log1p(alpha * theta + beta * theta) - 3 * log1p(4 * beta * theta) -
                          log1p(2 * beta * theta) + log_beta + log_theta +
                          log_sum_exp(std::array<double, 60>{
                              log(12288) + 7 * log_alpha + 3 * log_beta + 10 * log_theta,
                              log(122880) + 6 * log_alpha + 4 * log_beta + 10 * log_theta,
                              log(516096) + 5 * log_alpha + 5 * log_beta + 10 * log_theta,
                              log(1179648) + 4 * log_alpha + 6 * log_beta + 10 * log_theta,
                              log(1585152) + 3 * log_alpha + 7 * log_beta + 10 * log_theta,
                              log(1253376) + 2 * log_alpha + 8 * log_beta + 10 * log_theta,
                              log(540672) + log_alpha + 9 * log_beta + 10 * log_theta,
                              log(98304) + 10 * log_beta + 10 * log_theta,
                              log(15360) + 7 * log_alpha + 2 * log_beta + 9 * log_theta,
                              log(208896) + 6 * log_alpha + 3 * log_beta + 9 * log_theta,
                              log(1118208) + 5 * log_alpha + 4 * log_beta + 9 * log_theta,
                              log(3127296) + 4 * log_alpha + 5 * log_beta + 9 * log_theta,
                              log(4998144) + 3 * log_alpha + 6 * log_beta + 9 * log_theta,
                              log(4601856) + 2 * log_alpha + 7 * log_beta + 9 * log_theta,
                              log(2273280) + log_alpha + 8 * log_beta + 9 * log_theta,
                              log(466944) + 9 * log_beta + 9 * log_theta,
                              log(6784) + 7 * log_alpha + log_beta + 8 * log_theta,
                              log(136960) + 6 * log_alpha + 2 * log_beta + 8 * log_theta,
                              log(981760) + 5 * log_alpha + 3 * log_beta + 8 * log_theta,
                              log(3468800) + 4 * log_alpha + 4 * log_beta + 8 * log_theta,
                              log(6745216) + 3 * log_alpha + 5 * log_beta + 8 * log_theta,
                              log(7357184) + 2 * log_alpha + 6 * log_beta + 8 * log_theta,
                              log(4219392) + log_alpha + 7 * log_beta + 8 * log_theta,
                              log(990208) + 8 * log_beta + 8 * log_theta,
                              log(1184) + 7 * log_alpha + 7 * log_theta,
                              log(42368) + 6 * log_alpha + log_beta + 7 * log_theta,
                              log(442752) + 5 * log_alpha + 2 * log_beta + 7 * log_theta,
                              log(2077760) + 4 * log_alpha + 3 * log_beta + 7 * log_theta,
                              log(5083744) + 3 * log_alpha + 4 * log_beta + 7 * log_theta,
                              log(6734656) + 2 * log_alpha + 5 * log_beta + 7 * log_theta,
                              log(4574592) + log_alpha + 6 * log_beta + 7 * log_theta,
                              log(1247488) + 7 * log_beta + 7 * log_theta,
                              log(5328) + 6 * log_alpha + 6 * log_theta,
                              log(103760) + 5 * log_alpha + log_beta + 6 * log_theta,
                              log(713680) + 4 * log_alpha + 2 * log_beta + 6 * log_theta,
                              log(2322064) + 3 * log_alpha + 3 * log_beta + 6 * log_theta,
                              log(3881632) + 2 * log_alpha + 4 * log_beta + 6 * log_theta,
                              log(3222048) + log_alpha + 5 * log_beta + 6 * log_theta,
                              log(1051840) + 6 * log_beta + 6 * log_theta,
                              log(10128) + 5 * log_alpha + 5 * log_theta,
                              log(133616) + 4 * log_alpha + log_beta + 5 * log_theta,
                              log(643880) + 3 * log_alpha + 2 * log_beta + 5 * log_theta,
                              log(1444736) + 2 * log_alpha + 3 * log_beta + 5 * log_theta,
                              log(1533896) + log_alpha + 4 * log_beta + 5 * log_theta,
                              log(624528) + 5 * log_beta + 5 * log_theta,
                              log(10616) + 4 * log_alpha + 4 * log_theta,
                              log(100152) + 3 * log_alpha + log_beta + 4 * log_theta,
                              log(338164) + 2 * log_alpha + 2 * log_beta + 4 * log_theta,
                              log(491540) + log_alpha + 3 * log_beta + 4 * log_theta,
                              log(262512) + 4 * log_beta + 4 * log_theta,
                              log(6706) + 3 * log_alpha + 3 * log_theta,
                              log(45252) + 2 * log_alpha + log_beta + 3 * log_theta,
                              log(101350) + log_alpha + 2 * log_beta + 3 * log_theta,
                              log(75960) + 3 * log_beta + 3 * log_theta,
                              log(2633) + 2 * log_alpha + 2 * log_theta,
                              log(12085) + log_alpha + log_beta + 2 * log_theta,
                              log(14264) + 2 * log_beta + 2 * log_theta,
                              log(630) + log_alpha + log_theta,
                              log(1553) + log_beta + log_theta,
                              log(74)});

  const double log_p002 = -log(192) - 2 * log1p(2 * alpha * theta + 4 * beta * theta) -
                          3 * log1p(2 * alpha * theta + 2 * beta * theta) -
                          log1p(alpha * theta + 2 * beta * theta) -
                          log1p(alpha * theta + beta * theta) - 3 * log1p(4 * beta * theta) -
                          log1p(2 * beta * theta) + log_theta +
                          log_sum_exp(std::array<double, 66>{
                              log(12288) + 7 * log_alpha + 4 * log_beta + 10 * log_theta,
                              log(122880) + 6 * log_alpha + 5 * log_beta + 10 * log_theta,
                              log(516096) + 5 * log_alpha + 6 * log_beta + 10 * log_theta,
                              log(1179648) + 4 * log_alpha + 7 * log_beta + 10 * log_theta,
                              log(1585152) + 3 * log_alpha + 8 * log_beta + 10 * log_theta,
                              log(1253376) + 2 * log_alpha + 9 * log_beta + 10 * log_theta,
                              log(540672) + log_alpha + 10 * log_beta + 10 * log_theta,
                              log(98304) + 11 * log_beta + 10 * log_theta,
                              log(15360) + 7 * log_alpha + 3 * log_beta + 9 * log_theta,
                              log(208896) + 6 * log_alpha + 4 * log_beta + 9 * log_theta,
                              log(1118208) + 5 * log_alpha + 5 * log_beta + 9 * log_theta,
                              log(3127296) + 4 * log_alpha + 6 * log_beta + 9 * log_theta,
                              log(4998144) + 3 * log_alpha + 7 * log_beta + 9 * log_theta,
                              log(4601856) + 2 * log_alpha + 8 * log_beta + 9 * log_theta,
                              log(2273280) + log_alpha + 9 * log_beta + 9 * log_theta,
                              log(466944) + 10 * log_beta + 9 * log_theta,
                              log(7296) + 7 * log_alpha + 2 * log_beta + 8 * log_theta,
                              log(142080) + 6 * log_alpha + 3 * log_beta + 8 * log_theta,
                              log(1001216) + 5 * log_alpha + 4 * log_beta + 8 * log_theta,
                              log(3501568) + 4 * log_alpha + 5 * log_beta + 8 * log_theta,
                              log(6760064) + 3 * log_alpha + 6 * log_beta + 8 * log_theta,
                              log(7331584) + 2 * log_alpha + 7 * log_beta + 8 * log_theta,
                              log(4184576) + log_alpha + 8 * log_beta + 8 * log_theta,
                              log(977920) + 9 * log_beta + 8 * log_theta,
                              log(1824) + 7 * log_alpha + log_beta + 7 * log_theta,
                              log(51072) + 6 * log_alpha + 2 * log_beta + 7 * log_theta,
                              log(486784) + 5 * log_alpha + 3 * log_beta + 7 * log_theta,
                              log(2178368) + 4 * log_alpha + 4 * log_beta + 7 * log_theta,
                              log(5167584) + 3 * log_alpha + 5 * log_beta + 7 * log_theta,
                              log(6683712) + 2 * log_alpha + 6 * log_beta + 7 * log_theta,
                              log(4446080) + log_alpha + 7 * log_beta + 7 * log_theta,
                              log(1189120) + 8 * log_beta + 7 * log_theta,
                              log(384) + 7 * log_alpha + 6 * log_theta,
                              log(12048) + 6 * log_alpha + log_beta + 6 * log_theta,
                              log(147728) + 5 * log_alpha + 2 * log_beta + 6 * log_theta,
                              log(846864) + 4 * log_alpha + 3 * log_beta + 6 * log_theta,
                              log(2491856) + 3 * log_alpha + 4 * log_beta + 6 * log_theta,
                              log(3873056) + 2 * log_alpha + 5 * log_beta + 6 * log_theta,
                              log(3010976) + log_alpha + 6 * log_beta + 6 * log_theta,
                              log(917440) + 7 * log_beta + 6 * log_theta,
                              log(1728) + 6 * log_alpha + 5 * log_theta,
                              log(30128) + 5 * log_alpha + log_beta + 5 * log_theta,
                              log(221488) + 4 * log_alpha + 2 * log_beta + 5 * log_theta,
                              log(808360) + 3 * log_alpha + 3 * log_beta + 5 * log_theta,
                              log(1507904) + 2 * log_alpha + 4 * log_beta + 5 * log_theta,
                              log(1359400) + log_alpha + 5 * log_beta + 5 * log_theta,
                              log(461776) + 6 * log_beta + 5 * log_theta,
                              log(3248) + 5 * log_alpha + 4 * log_theta,
                              log(38312) + 4 * log_alpha + log_beta + 4 * log_theta,
                              log(181848) + 3 * log_alpha + 2 * log_beta + 4 * log_theta,
                              log(408932) + 2 * log_alpha + 3 * log_beta + 4 * log_theta,
                              log(419812) + log_alpha + 4 * log_beta + 4 * log_theta,
                              log(150832) + 5 * log_beta + 4 * log_theta,
                              log(3288) + 4 * log_alpha + 3 * log_theta,
                              log(26666) + 3 * log_alpha + log_beta + 3 * log_theta,
                              log(77908) + 2 * log_alpha + 2 * log_beta + 3 * log_theta,
                              log(90646) + log_alpha + 3 * log_beta + 3 * log_theta,
                              log(30760) + 4 * log_beta + 3 * log_theta,
                              log(1892) + 3 * log_alpha + 2 * log_theta,
                              log(9673) + 2 * log_alpha + log_beta + 2 * log_theta,
                              log(13865) + log_alpha + 2 * log_beta + 2 * log_theta,
                              log(3552) + 3 * log_beta + 2 * log_theta,
                              log(582) + 2 * log_alpha + log_theta,
                              log(1424) + log_alpha + log_beta + log_theta,
                              log(177) + 2 * log_beta + log_theta,
                              log(74) + log_alpha});

  const double log_p012 = -log(192) - 2 * log1p(2 * alpha * theta + 4 * beta * theta) -
                          3 * log1p(2 * alpha * theta + 2 * beta * theta) -
                          log1p(alpha * theta + 2 * beta * theta) -
                          log1p(alpha * theta + beta * theta) - 3 * log1p(4 * beta * theta) -
                          log1p(2 * beta * theta) + log_beta + 2 * log_theta +
                          log_sum_exp(std::array<double, 59>{
                              log(12288) + 7 * log_alpha + 3 * log_beta + 9 * log_theta,
                              log(122880) + 6 * log_alpha + 4 * log_beta + 9 * log_theta,
                              log(516096) + 5 * log_alpha + 5 * log_beta + 9 * log_theta,
                              log(1179648) + 4 * log_alpha + 6 * log_beta + 9 * log_theta,
                              log(1585152) + 3 * log_alpha + 7 * log_beta + 9 * log_theta,
                              log(1253376) + 2 * log_alpha + 8 * log_beta + 9 * log_theta,
                              log(540672) + log_alpha + 9 * log_beta + 9 * log_theta,
                              log(98304) + 10 * log_beta + 9 * log_theta,
                              log(15360) + 7 * log_alpha + 2 * log_beta + 8 * log_theta,
                              log(208896) + 6 * log_alpha + 3 * log_beta + 8 * log_theta,
                              log(1118208) + 5 * log_alpha + 4 * log_beta + 8 * log_theta,
                              log(3127296) + 4 * log_alpha + 5 * log_beta + 8 * log_theta,
                              log(4998144) + 3 * log_alpha + 6 * log_beta + 8 * log_theta,
                              log(4601856) + 2 * log_alpha + 7 * log_beta + 8 * log_theta,
                              log(2273280) + log_alpha + 8 * log_beta + 8 * log_theta,
                              log(466944) + 9 * log_beta + 8 * log_theta,
                              log(6784) + 7 * log_alpha + log_beta + 7 * log_theta,
                              log(136960) + 6 * log_alpha + 2 * log_beta + 7 * log_theta,
                              log(979712) + 5 * log_alpha + 3 * log_beta + 7 * log_theta,
                              log(3452416) + 4 * log_alpha + 4 * log_beta + 7 * log_theta,
                              log(6694016) + 3 * log_alpha + 5 * log_beta + 7 * log_theta,
                              log(7279360) + 2 * log_alpha + 6 * log_beta + 7 * log_theta,
                              log(4162048) + log_alpha + 7 * log_beta + 7 * log_theta,
                              log(973824) + 8 * log_beta + 7 * log_theta,
                              log(1184) + 7 * log_alpha + 6 * log_theta,
                              log(42368) + 6 * log_alpha + log_beta + 6 * log_theta,
                              log(440192) + 5 * log_alpha + 2 * log_beta + 6 * log_theta,
                              log(2048064) + 4 * log_alpha + 3 * log_beta + 6 * log_theta,
                              log(4959328) + 3 * log_alpha + 4 * log_beta + 6 * log_theta,
                              log(6491968) + 2 * log_alpha + 5 * log_beta + 6 * log_theta,
                              log(4351360) + log_alpha + 6 * log_beta + 6 * log_theta,
                              log(1169664) + 7 * log_beta + 6 * log_theta,
                              log(5328) + 6 * log_alpha + 5 * log_theta,
                              log(102736) + 5 * log_alpha + log_beta + 5 * log_theta,
                              log(694480) + 4 * log_alpha + 2 * log_beta + 5 * log_theta,
                              log(2207760) + 3 * log_alpha + 3 * log_beta + 5 * log_theta,
                              log(3577760) + 2 * log_alpha + 4 * log_beta + 5 * log_theta,
                              log(2850464) + log_alpha + 5 * log_beta + 5 * log_theta,
                              log(881856) + 6 * log_beta + 5 * log_theta,
                              log(10000) + 5 * log_alpha + 4 * log_theta,
                              log(128368) + 4 * log_alpha + log_beta + 4 * log_theta,
                              log(593928) + 3 * log_alpha + 2 * log_beta + 4 * log_theta,
                              log(1252352) + 2 * log_alpha + 3 * log_beta + 4 * log_theta,
                              log(1210472) + log_alpha + 4 * log_beta + 4 * log_theta,
                              log(428688) + 5 * log_beta + 4 * log_theta,
                              log(10104) + 4 * log_alpha + 3 * log_theta,
                              log(89736) + 3 * log_alpha + log_beta + 3 * log_theta,
                              log(273348) + 2 * log_alpha + 2 * log_beta + 3 * log_theta,
                              log(333780) + log_alpha + 3 * log_beta + 3 * log_theta,
                              log(133680) + 4 * log_beta + 3 * log_theta,
                              log(5874) + 3 * log_alpha + 2 * log_theta,
                              log(34188) + 2 * log_alpha + log_beta + 2 * log_theta,
                              log(57966) + log_alpha + 2 * log_beta + 2 * log_theta,
                              log(25752) + 3 * log_beta + 2 * log_theta,
                              log(1881) + 2 * log_alpha + log_theta,
                              log(5805) + log_alpha + log_beta + log_theta,
                              log(2784) + 2 * log_beta + log_theta,
                              log(258) + log_alpha,
                              log(129) + log_beta});

  const std::array<double, 64> p = {log_p000,
                                    log_p001,
                                    log_p002,
                                    log_p001,
                                    log_p001,
                                    log_p001,
                                    log_p012,
                                    log_p012,
                                    log_p002,
                                    log_p012,
                                    log_p002,
                                    log_p012,
                                    log_p001,
                                    log_p012,
                                    log_p012,
                                    log_p001,
                                    log_p001,
                                    log_p001,
                                    log_p012,
                                    log_p012,
                                    log_p001,
                                    log_p000,
                                    log_p001,
                                    log_p002,
                                    log_p012,
                                    log_p001,
                                    log_p001,
                                    log_p012,
                                    log_p012,
                                    log_p002,
                                    log_p012,
                                    log_p002,
                                    log_p002,
                                    log_p012,
                                    log_p002,
                                    log_p012,
                                    log_p012,
                                    log_p001,
                                    log_p001,
                                    log_p012,
                                    log_p002,
                                    log_p001,
                                    log_p000,
                                    log_p001,
                                    log_p012,
                                    log_p012,
                                    log_p001,
                                    log_p001,
                                    log_p001,
                                    log_p012,
                                    log_p012,
                                    log_p001,
                                    log_p012,
                                    log_p002,
                                    log_p012,
                                    log_p002,
                                    log_p012,
                                    log_p012,
                                    log_p001,
                                    log_p001,
                                    log_p001,
                                    log_p002,
                                    log_p001,
                                    log_p000};

  return p;
}

inline double k2p_sequence_nll(const double alpha, const double beta, const double theta, const TuplePatterns& tp) {
  if (tp.order() != 3) {
    throw std::invalid_argument("k2p_neg_loglikelihood: expected 3-taxon TuplePatterns");
  }
  const auto log_ps = k2p_pattern_probs(alpha, beta, theta);
  
  double neg_loglikelihood = 0.0;
  for (size_t i = 0; i < 64; ++i) {
    neg_loglikelihood -= tp.site_pattern[i] * log_ps[i];
  }
  return neg_loglikelihood;

}


struct K2PMLEResult {
	double alpha_ml;
	double beta_ml;
	double neg_logL;
};

// Objective for L-BFGS-B: single free parameter alpha, with beta = (mu - alpha)/2.
class K2PObjective {
public:
	K2PObjective(const TuplePatterns& tp, double theta, double mu)
	    : tp_(tp), theta_(theta), mu_(mu) {}

	double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
		const double alpha = x[0];
		const double beta  = (mu_ - alpha) / 2.0;
		const double fx    = k2p_sequence_nll(alpha, beta, theta_, tp_);

		// central differences, step sizes chosen to keep both params positive
		const double h_plus  = std::min(1e-4, (mu_ - alpha) * 0.1);
		const double h_minus = std::min(1e-4, alpha * 0.1);
		const double fx_p = k2p_sequence_nll(alpha + h_plus,  (mu_ - alpha - h_plus)  / 2.0, theta_, tp_);
		const double fx_m = k2p_sequence_nll(alpha - h_minus, (mu_ - alpha + h_minus) / 2.0, theta_, tp_);

		grad.resize(1);
		grad[0] = (fx_p - fx_m) / (h_plus + h_minus);
		return fx;
	}

private:
	const TuplePatterns& tp_;
	const double theta_;
	const double mu_;
};

// Find the MLE of alpha and beta subject to alpha + 2*beta = mu.
// theta and mu are held fixed; alpha is the single free parameter.
// Initial value: alpha = beta = mu/3 (JC69-equivalent point).
inline K2PMLEResult k2p_mle(const TuplePatterns& tp,
                             const double theta,
                             const double mu) {
	if (tp.order() != 3)
		throw std::invalid_argument("k2p_mle: expected 3-taxon TuplePatterns");

	const double eps = 1e-6;
	const double alpha_init = std::max(eps, std::min(mu - 2.0 * eps, mu / 3.0));

	K2PObjective objective(tp, theta, mu);

	LBFGSpp::LBFGSBParam<double> param;
	param.epsilon        = 1e-6;
	param.max_iterations = 200;
	param.past           = 0;

	LBFGSpp::LBFGSBSolver<double> solver(param);

	Eigen::VectorXd lb(1), ub(1);
	lb << eps;
	ub << mu - 2.0 * eps;

	Eigen::VectorXd x(1);
	x << alpha_init;

	double fx = k2p_sequence_nll(alpha_init, (mu - alpha_init) / 2.0, theta, tp);
	const double init_fx = fx;

	try {
		solver.minimize(objective, x, fx, lb, ub);
	} catch (const std::exception&) {
		if (fx < init_fx) {
			PHYLOTRIPLET_WARN("k2p_mle: L-BFGS-B failed, but found a better solution; using it.");
		} else {
			PHYLOTRIPLET_WARN("k2p_mle: L-BFGS-B failed, using initial guess.");
			x[0] = alpha_init;
			fx   = init_fx;
		}
	}

	K2PMLEResult result;
	result.alpha_ml = x[0];
	result.beta_ml  = (mu - x[0]) / 2.0;
	result.neg_logL = fx;
	return result;
}




// ---------------------------------------------------------------------------
// Calculation for -log(P(species tree))) under K2P + MSC
// depends on parameters, but not on the branch lengths as this is integrated out
// ---------------------------------------------------------------------------
inline std::array<double, 64> k2p_species_tree_probs(const double alpha,
                                                const double beta,
                                                const double theta) {
  // TODO: bound log probabilities away from -inf

  const double log_alpha = std::log(alpha);
  const double log_beta = std::log(beta);
  const double log_theta = std::log(theta);


const double log_p000 =
  -6*log(2)
  - 2*log1p(2*alpha*theta + 4*beta*theta)
  - 3*log1p(2*alpha*theta + 2*beta*theta)
  - log1p(alpha*theta + 2*beta*theta)
  - log1p(alpha*theta + beta*theta)
  - 3*log1p(4*beta*theta)
  - log1p(2*beta*theta)
  + log_sum_exp(std::array<double,68>{
      12*log(2) + 7*log_alpha + 4*log_beta + 11*log_theta,
      log(40960) + 6*log_alpha + 5*log_beta + 11*log_theta,
      log(172032) + 5*log_alpha + 6*log_beta + 11*log_theta,
      log(393216) + 4*log_alpha + 7*log_beta + 11*log_theta,
      log(528384) + 3*log_alpha + 8*log_beta + 11*log_theta,
      log(417792) + 2*log_alpha + 9*log_beta + 11*log_theta,
      log(180224) + log_alpha + 10*log_beta + 11*log_theta,
      15*log(2) + 11*log_beta + 11*log_theta,
      log(5120) + 7*log_alpha + 3*log_beta + 10*log_theta,
      log(69632) + 6*log_alpha + 4*log_beta + 10*log_theta,
      log(372736) + 5*log_alpha + 5*log_beta + 10*log_theta,
      log(1042432) + 4*log_alpha + 6*log_beta + 10*log_theta,
      log(1666048) + 3*log_alpha + 7*log_beta + 10*log_theta,
      log(1533952) + 2*log_alpha + 8*log_beta + 10*log_theta,
      log(757760) + log_alpha + 9*log_beta + 10*log_theta,
      log(155648) + 10*log_beta + 10*log_theta,
      log(2432) + 7*log_alpha + 2*log_beta + 9*log_theta,
      log(47360) + 6*log_alpha + 3*log_beta + 9*log_theta,
      log(335104) + 5*log_alpha + 4*log_beta + 9*log_theta,
      log(1178112) + 4*log_alpha + 5*log_beta + 9*log_theta,
      log(2287488) + 3*log_alpha + 6*log_beta + 9*log_theta,
      log(2495744) + 2*log_alpha + 7*log_beta + 9*log_theta,
      log(1433088) + log_alpha + 8*log_beta + 9*log_theta,
      log(336896) + 9*log_beta + 9*log_theta,
      log(608) + 7*log_alpha + log_beta + 8*log_theta,
      log(17024) + 6*log_alpha + 2*log_beta + 8*log_theta,
      log(163968) + 5*log_alpha + 3*log_beta + 8*log_theta,
      log(745920) + 4*log_alpha + 4*log_beta + 8*log_theta,
      log(1805472) + 3*log_alpha + 5*log_beta + 8*log_theta,
      log(2389696) + 2*log_alpha + 6*log_beta + 8*log_theta,
      log(1630848) + log_alpha + 7*log_beta + 8*log_theta,
      log(448256) + 8*log_beta + 8*log_theta,
      7*log(2) + 7*log_alpha + 7*log_theta,
      log(4016) + 6*log_alpha + log_beta + 7*log_theta,
      log(50096) + 5*log_alpha + 2*log_beta + 7*log_theta,
      log(297136) + 4*log_alpha + 3*log_beta + 7*log_theta,
      log(918768) + 3*log_alpha + 4*log_beta + 7*log_theta,
      log(1520224) + 2*log_alpha + 5*log_beta + 7*log_theta,
      log(1276128) + log_alpha + 6*log_beta + 7*log_theta,
      log(427328) + 7*log_beta + 7*log_theta,
      log(576) + 6*log_alpha + 6*log_theta,
      log(10256) + 5*log_alpha + log_beta + 6*log_theta,
      log(79632) + 4*log_alpha + 2*log_beta + 6*log_theta,
      log(320504) + 3*log_alpha + 3*log_beta + 6*log_theta,
      log(684224) + 2*log_alpha + 4*log_beta + 6*log_theta,
      log(732792) + log_alpha + 5*log_beta + 6*log_theta,
      log(310768) + 6*log_beta + 6*log_theta,
      log(1104) + 5*log_alpha + 5*log_theta,
      log(13944) + 4*log_alpha + log_beta + 5*log_theta,
      log(77224) + 3*log_alpha + 2*log_beta + 5*log_theta,
      log(220844) + 2*log_alpha + 3*log_beta + 5*log_theta,
      log(311372) + log_alpha + 4*log_beta + 5*log_theta,
      log(170640) + 5*log_beta + 5*log_theta,
      log(1192) + 4*log_alpha + 4*log_theta,
      log(11726) + 3*log_alpha + log_beta + 4*log_theta,
      log(48844) + 2*log_alpha + 2*log_beta + 4*log_theta,
      log(94402) + log_alpha + 3*log_beta + 4*log_theta,
      log(67896) + 4*log_beta + 4*log_theta,
      log(828) + 3*log_alpha + 3*log_theta,
      log(6547) + 2*log_alpha + log_beta + 3*log_theta,
      log(19027) + log_alpha + 2*log_beta + 3*log_theta,
      log(18672) + 3*log_beta + 3*log_theta,
      log(394) + 2*log_alpha + 2*log_theta,
      log(2256) + log_alpha + log_beta + 2*log_theta,
      log(3339) + 2*log_beta + 2*log_theta,
      log(118) + log_alpha + log_theta,
      log(348) + log_beta + log_theta,
      4*log(2)});


const double log_p001 =
  -6*log(2)
  - 2*log1p(2*alpha*theta + 4*beta*theta)
  - 3*log1p(2*alpha*theta + 2*beta*theta)
  - log1p(alpha*theta + 2*beta*theta)
  - log1p(alpha*theta + beta*theta)
  - 3*log1p(4*beta*theta)
  - log1p(2*beta*theta)
  + log_beta
  + log_theta
  + log_sum_exp(std::array<double,60>{
      12*log(2) + 7*log_alpha + 3*log_beta + 10*log_theta,
      log(40960) + 6*log_alpha + 4*log_beta + 10*log_theta,
      log(172032) + 5*log_alpha + 5*log_beta + 10*log_theta,
      log(393216) + 4*log_alpha + 6*log_beta + 10*log_theta,
      log(528384) + 3*log_alpha + 7*log_beta + 10*log_theta,
      log(417792) + 2*log_alpha + 8*log_beta + 10*log_theta,
      log(180224) + log_alpha + 9*log_beta + 10*log_theta,
      15*log(2) + 10*log_beta + 10*log_theta,
      log(5120) + 7*log_alpha + 2*log_beta + 9*log_theta,
      log(69632) + 6*log_alpha + 3*log_beta + 9*log_theta,
      log(372736) + 5*log_alpha + 4*log_beta + 9*log_theta,
      log(1042432) + 4*log_alpha + 5*log_beta + 9*log_theta,
      log(1666048) + 3*log_alpha + 6*log_beta + 9*log_theta,
      log(1533952) + 2*log_alpha + 7*log_beta + 9*log_theta,
      log(757760) + log_alpha + 8*log_beta + 9*log_theta,
      log(155648) + 9*log_beta + 9*log_theta,
      log(2432) + 7*log_alpha + log_beta + 8*log_theta,
      log(47360) + 6*log_alpha + 2*log_beta + 8*log_theta,
      log(335104) + 5*log_alpha + 3*log_beta + 8*log_theta,
      log(1178112) + 4*log_alpha + 4*log_beta + 8*log_theta,
      log(2287488) + 3*log_alpha + 5*log_beta + 8*log_theta,
      log(2495744) + 2*log_alpha + 6*log_beta + 8*log_theta,
      log(1433088) + log_alpha + 7*log_beta + 8*log_theta,
      log(336896) + 8*log_beta + 8*log_theta,
      log(480) + 7*log_alpha + 7*log_theta,
      log(15744) + 6*log_alpha + log_beta + 7*log_theta,
      log(158592) + 5*log_alpha + 2*log_beta + 7*log_theta,
      log(732608) + 4*log_alpha + 3*log_beta + 7*log_theta,
      log(1781792) + 3*log_alpha + 4*log_beta + 7*log_theta,
      log(2358208) + 2*log_alpha + 5*log_beta + 7*log_theta,
      log(1604736) + log_alpha + 6*log_beta + 7*log_theta,
      log(439040) + 7*log_beta + 7*log_theta,
      log(2160) + 6*log_alpha + 6*log_theta,
      log(39664) + 5*log_alpha + log_beta + 6*log_theta,
      log(265072) + 4*log_alpha + 2*log_beta + 6*log_theta,
      log(851376) + 3*log_alpha + 3*log_beta + 6*log_theta,
      log(1416672) + 2*log_alpha + 4*log_beta + 6*log_theta,
      log(1175392) + log_alpha + 5*log_beta + 6*log_theta,
      log(384064) + 6*log_beta + 6*log_theta,
      log(4144) + 5*log_alpha + 5*log_theta,
      log(52624) + 4*log_alpha + log_beta + 5*log_theta,
      log(248600) + 3*log_alpha + 2*log_beta + 5*log_theta,
      log(552128) + 2*log_alpha + 3*log_beta + 5*log_theta,
      log(582328) + log_alpha + 4*log_beta + 5*log_theta,
      log(235312) + 5*log_beta + 5*log_theta,
      log(4424) + 4*log_alpha + 4*log_theta,
      log(40744) + 3*log_alpha + log_beta + 4*log_theta,
      log(135340) + 2*log_alpha + 2*log_beta + 4*log_theta,
      log(193596) + log_alpha + 3*log_beta + 4*log_theta,
      log(101216) + 4*log_beta + 4*log_theta,
      log(2862) + 3*log_alpha + 3*log_theta,
      log(18908) + 2*log_alpha + log_beta + 3*log_theta,
      log(41226) + log_alpha + 2*log_beta + 3*log_theta,
      log(29768) + 3*log_beta + 3*log_theta,
      log(1143) + 2*log_alpha + 2*log_theta,
      log(5055) + log_alpha + log_beta + 2*log_theta,
      log(5660) + 2*log_beta + 2*log_theta,
      log(270) + log_alpha + log_theta,
      log(623) + log_beta + log_theta,
      log(30)});


const double log_p002 =
  -6*log(2)
  - 2*log1p(2*alpha*theta + 4*beta*theta)
  - 3*log1p(2*alpha*theta + 2*beta*theta)
  - log1p(alpha*theta + 2*beta*theta)
  - log1p(alpha*theta + beta*theta)
  - 3*log1p(4*beta*theta)
  - log1p(2*beta*theta)
  + log_theta
  + log_sum_exp(std::array<double,66>{
      12*log(2) + 7*log_alpha + 4*log_beta + 10*log_theta,
      log(40960) + 6*log_alpha + 5*log_beta + 10*log_theta,
      log(172032) + 5*log_alpha + 6*log_beta + 10*log_theta,
      log(393216) + 4*log_alpha + 7*log_beta + 10*log_theta,
      log(528384) + 3*log_alpha + 8*log_beta + 10*log_theta,
      log(417792) + 2*log_alpha + 9*log_beta + 10*log_theta,
      log(180224) + log_alpha + 10*log_beta + 10*log_theta,
      15*log(2) + 11*log_beta + 10*log_theta,
      log(5120) + 7*log_alpha + 3*log_beta + 9*log_theta,
      log(69632) + 6*log_alpha + 4*log_beta + 9*log_theta,
      log(372736) + 5*log_alpha + 5*log_beta + 9*log_theta,
      log(1042432) + 4*log_alpha + 6*log_beta + 9*log_theta,
      log(1666048) + 3*log_alpha + 7*log_beta + 9*log_theta,
      log(1533952) + 2*log_alpha + 8*log_beta + 9*log_theta,
      log(757760) + log_alpha + 9*log_beta + 9*log_theta,
      log(155648) + 10*log_beta + 9*log_theta,
      log(2432) + 7*log_alpha + 2*log_beta + 8*log_theta,
      log(47360) + 6*log_alpha + 3*log_beta + 8*log_theta,
      log(335104) + 5*log_alpha + 4*log_beta + 8*log_theta,
      log(1178112) + 4*log_alpha + 5*log_beta + 8*log_theta,
      log(2287488) + 3*log_alpha + 6*log_beta + 8*log_theta,
      log(2495744) + 2*log_alpha + 7*log_beta + 8*log_theta,
      log(1433088) + log_alpha + 8*log_beta + 8*log_theta,
      log(336896) + 9*log_beta + 8*log_theta,
      log(608) + 7*log_alpha + log_beta + 7*log_theta,
      log(17024) + 6*log_alpha + 2*log_beta + 7*log_theta,
      log(163968) + 5*log_alpha + 3*log_beta + 7*log_theta,
      log(743872) + 4*log_alpha + 4*log_beta + 7*log_theta,
      log(1791136) + 3*log_alpha + 5*log_beta + 7*log_theta,
      log(2352832) + 2*log_alpha + 6*log_beta + 7*log_theta,
      log(1589888) + log_alpha + 7*log_beta + 7*log_theta,
      log(431872) + 8*log_beta + 7*log_theta,
      7*log(2) + 7*log_alpha + 6*log_theta,
      log(4016) + 6*log_alpha + log_beta + 6*log_theta,
      log(50096) + 5*log_alpha + 2*log_beta + 6*log_theta,
      log(294064) + 4*log_alpha + 3*log_beta + 6*log_theta,
      log(887536) + 3*log_alpha + 4*log_beta + 6*log_theta,
      log(1415776) + 2*log_alpha + 5*log_beta + 6*log_theta,
      log(1129696) + log_alpha + 6*log_beta + 6*log_theta,
      log(353088) + 7*log_beta + 6*log_theta,
      log(576) + 6*log_alpha + 5*log_theta,
      log(10256) + 5*log_alpha + log_beta + 5*log_theta,
      log(77968) + 4*log_alpha + 2*log_beta + 5*log_theta,
      log(295032) + 3*log_alpha + 3*log_beta + 5*log_theta,
      log(570560) + 2*log_alpha + 4*log_beta + 5*log_theta,
      log(533112) + log_alpha + 5*log_beta + 5*log_theta,
      log(187632) + 6*log_beta + 5*log_theta,
      log(1104) + 5*log_alpha + 4*log_theta,
      log(13560) + 4*log_alpha + log_beta + 4*log_theta,
      log(67400) + 3*log_alpha + 2*log_beta + 4*log_theta,
      log(158636) + 2*log_alpha + 3*log_beta + 4*log_theta,
      log(170380) + log_alpha + 4*log_beta + 4*log_theta,
      log(64240) + 5*log_beta + 4*log_theta,
      log(1160) + 4*log_alpha + 3*log_theta,
      log(9918) + 3*log_alpha + log_beta + 3*log_theta,
      log(30556) + 2*log_alpha + 2*log_beta + 3*log_theta,
      log(37506) + log_alpha + 3*log_beta + 3*log_theta,
      log(13624) + 4*log_beta + 3*log_theta,
      log(700) + 3*log_alpha + 2*log_theta,
      log(3787) + 2*log_alpha + log_beta + 2*log_theta,
      log(5747) + log_alpha + 2*log_beta + 2*log_theta,
      log(1624) + 3*log_beta + 2*log_theta,
      log(226) + 2*log_alpha + log_theta,
      log(584) + log_alpha + log_beta + log_theta,
      log(83) + 2*log_beta + log_theta,
      log(30) + log_alpha});


const double log_p010 =
  -6*log(2)
  - 2*log1p(2*alpha*theta + 4*beta*theta)
  - 3*log1p(2*alpha*theta + 2*beta*theta)
  - log1p(alpha*theta + 2*beta*theta)
  - log1p(alpha*theta + beta*theta)
  - 3*log1p(4*beta*theta)
  - log1p(2*beta*theta)
  + log_beta
  + log_theta
  + log_sum_exp(std::array<double,60>{
      12*log(2) + 7*log_alpha + 3*log_beta + 10*log_theta,
      log(40960) + 6*log_alpha + 4*log_beta + 10*log_theta,
      log(172032) + 5*log_alpha + 5*log_beta + 10*log_theta,
      log(393216) + 4*log_alpha + 6*log_beta + 10*log_theta,
      log(528384) + 3*log_alpha + 7*log_beta + 10*log_theta,
      log(417792) + 2*log_alpha + 8*log_beta + 10*log_theta,
      log(180224) + log_alpha + 9*log_beta + 10*log_theta,
      15*log(2) + 10*log_beta + 10*log_theta,
      log(5120) + 7*log_alpha + 2*log_beta + 9*log_theta,
      log(69632) + 6*log_alpha + 3*log_beta + 9*log_theta,
      log(372736) + 5*log_alpha + 4*log_beta + 9*log_theta,
      log(1042432) + 4*log_alpha + 5*log_beta + 9*log_theta,
      log(1666048) + 3*log_alpha + 6*log_beta + 9*log_theta,
      log(1533952) + 2*log_alpha + 7*log_beta + 9*log_theta,
      log(757760) + log_alpha + 8*log_beta + 9*log_theta,
      log(155648) + 9*log_beta + 9*log_theta,
      log(2176) + 7*log_alpha + log_beta + 8*log_theta,
      log(44800) + 6*log_alpha + 2*log_beta + 8*log_theta,
      log(323328) + 5*log_alpha + 3*log_beta + 8*log_theta,
      log(1145344) + 4*log_alpha + 4*log_beta + 8*log_theta,
      log(2228864) + 3*log_alpha + 5*log_beta + 8*log_theta,
      log(2430720) + 2*log_alpha + 6*log_beta + 8*log_theta,
      log(1393152) + log_alpha + 7*log_beta + 8*log_theta,
      log(326656) + 8*log_beta + 8*log_theta,
      log(352) + 7*log_alpha + 7*log_theta,
      log(13312) + 6*log_alpha + log_beta + 7*log_theta,
      log(142080) + 5*log_alpha + 2*log_beta + 7*log_theta,
      log(672576) + 4*log_alpha + 3*log_beta + 7*log_theta,
      log(1650976) + 3*log_alpha + 4*log_beta + 7*log_theta,
      log(2188224) + 2*log_alpha + 5*log_beta + 7*log_theta,
      log(1484928) + log_alpha + 6*log_beta + 7*log_theta,
      log(404224) + 7*log_beta + 7*log_theta,
      log(1584) + 6*log_alpha + 6*log_theta,
      log(32048) + 5*log_alpha + log_beta + 6*log_theta,
      log(224304) + 4*log_alpha + 2*log_beta + 6*log_theta,
      log(735344) + 3*log_alpha + 3*log_beta + 6*log_theta,
      log(1232480) + 2*log_alpha + 4*log_beta + 6*log_theta,
      log(1023328) + log_alpha + 5*log_beta + 6*log_theta,
      log(333888) + 6*log_beta + 6*log_theta,
      log(2992) + 5*log_alpha + 5*log_theta,
      log(40496) + 4*log_alpha + log_beta + 5*log_theta,
      log(197640) + 3*log_alpha + 2*log_beta + 5*log_theta,
      log(446304) + 2*log_alpha + 3*log_beta + 5*log_theta,
      log(475784) + log_alpha + 4*log_beta + 5*log_theta,
      log(194608) + 5*log_beta + 5*log_theta,
      log(3096) + 4*log_alpha + 4*log_theta,
      log(29704) + 3*log_alpha + log_beta + 4*log_theta,
      log(101412) + 2*log_alpha + 2*log_beta + 4*log_theta,
      log(148972) + log_alpha + 3*log_beta + 4*log_theta,
      log(80648) + 4*log_beta + 4*log_theta,
      log(1922) + 3*log_alpha + 3*log_theta,
      log(13172) + 2*log_alpha + log_beta + 3*log_theta,
      log(30062) + log_alpha + 2*log_beta + 3*log_theta,
      log(23096) + 3*log_beta + 3*log_theta,
      log(745) + 2*log_alpha + 2*log_theta,
      log(3515) + log_alpha + log_beta + 2*log_theta,
      log(4302) + 2*log_beta + 2*log_theta,
      log(180) + log_alpha + log_theta,
      log(465) + log_beta + log_theta,
      log(22)});


const double log_p012 =
  -6*log(2)
  - 2*log1p(2*alpha*theta + 4*beta*theta)
  - 3*log1p(2*alpha*theta + 2*beta*theta)
  - log1p(alpha*theta + 2*beta*theta)
  - log1p(alpha*theta + beta*theta)
  - 3*log1p(4*beta*theta)
  - log1p(2*beta*theta)
  + log_beta
  + 2*log_theta
  + log_sum_exp(std::array<double,59>{
      12*log(2) + 7*log_alpha + 3*log_beta + 9*log_theta,
      log(40960) + 6*log_alpha + 4*log_beta + 9*log_theta,
      log(172032) + 5*log_alpha + 5*log_beta + 9*log_theta,
      log(393216) + 4*log_alpha + 6*log_beta + 9*log_theta,
      log(528384) + 3*log_alpha + 7*log_beta + 9*log_theta,
      log(417792) + 2*log_alpha + 8*log_beta + 9*log_theta,
      log(180224) + log_alpha + 9*log_beta + 9*log_theta,
      15*log(2) + 10*log_beta + 9*log_theta,
      log(5120) + 7*log_alpha + 2*log_beta + 8*log_theta,
      log(69632) + 6*log_alpha + 3*log_beta + 8*log_theta,
      log(372736) + 5*log_alpha + 4*log_beta + 8*log_theta,
      log(1042432) + 4*log_alpha + 5*log_beta + 8*log_theta,
      log(1666048) + 3*log_alpha + 6*log_beta + 8*log_theta,
      log(1533952) + 2*log_alpha + 7*log_beta + 8*log_theta,
      log(757760) + log_alpha + 8*log_beta + 8*log_theta,
      log(155648) + 9*log_beta + 8*log_theta,
      log(2176) + 7*log_alpha + log_beta + 7*log_theta,
      log(44800) + 6*log_alpha + 2*log_beta + 7*log_theta,
      log(323328) + 5*log_alpha + 3*log_beta + 7*log_theta,
      log(1145344) + 4*log_alpha + 4*log_beta + 7*log_theta,
      log(2228864) + 3*log_alpha + 5*log_beta + 7*log_theta,
      log(2430720) + 2*log_alpha + 6*log_beta + 7*log_theta,
      log(1393152) + log_alpha + 7*log_beta + 7*log_theta,
      log(326656) + 8*log_beta + 7*log_theta,
      log(352) + 7*log_alpha + 6*log_theta,
      log(13312) + 6*log_alpha + log_beta + 6*log_theta,
      log(142080) + 5*log_alpha + 2*log_beta + 6*log_theta,
      log(671552) + 4*log_alpha + 3*log_beta + 6*log_theta,
      log(1643808) + 3*log_alpha + 4*log_beta + 6*log_theta,
      log(2169792) + 2*log_alpha + 5*log_beta + 6*log_theta,
      log(1464448) + log_alpha + 6*log_beta + 6*log_theta,
      log(396032) + 7*log_beta + 6*log_theta,
      log(1584) + 6*log_alpha + 5*log_theta,
      log(32048) + 5*log_alpha + log_beta + 5*log_theta,
      log(223280) + 4*log_alpha + 2*log_beta + 5*log_theta,
      log(724848) + 3*log_alpha + 3*log_beta + 5*log_theta,
      log(1193568) + 2*log_alpha + 4*log_beta + 5*log_theta,
      log(963424) + log_alpha + 5*log_beta + 5*log_theta,
      log(301376) + 6*log_beta + 5*log_theta,
      log(2992) + 5*log_alpha + 4*log_theta,
      log(40176) + 4*log_alpha + log_beta + 4*log_theta,
      log(192072) + 3*log_alpha + 2*log_beta + 4*log_theta,
      log(415328) + 2*log_alpha + 3*log_beta + 4*log_theta,
      log(409736) + log_alpha + 4*log_beta + 4*log_theta,
      log(147632) + 5*log_beta + 4*log_theta,
      log(3064) + 4*log_alpha + 3*log_theta,
      log(28440) + 3*log_alpha + log_beta + 3*log_theta,
      log(89700) + 2*log_alpha + 2*log_beta + 3*log_theta,
      log(112716) + log_alpha + 3*log_beta + 3*log_theta,
      log(46296) + 4*log_beta + 3*log_theta,
      log(1818) + 3*log_alpha + 2*log_theta,
      log(11052) + 2*log_alpha + log_beta + 2*log_theta,
      log(19438) + log_alpha + 2*log_beta + 2*log_theta,
      log(8952) + 3*log_beta + 2*log_theta,
      log(597) + 2*log_alpha + log_theta,
      log(1923) + log_alpha + log_beta + log_theta,
      log(970) + 2*log_beta + log_theta,
      log(84) + log_alpha,
      log(45) + log_beta});


const double log_p020 =
  -6*log(2)
  - 2*log1p(2*alpha*theta + 4*beta*theta)
  - 3*log1p(2*alpha*theta + 2*beta*theta)
  - log1p(alpha*theta + 2*beta*theta)
  - log1p(alpha*theta + beta*theta)
  - 3*log1p(4*beta*theta)
  - log1p(2*beta*theta)
  + log_theta
  + log_sum_exp(std::array<double,66>{
      12*log(2) + 7*log_alpha + 4*log_beta + 10*log_theta,
      log(40960) + 6*log_alpha + 5*log_beta + 10*log_theta,
      log(172032) + 5*log_alpha + 6*log_beta + 10*log_theta,
      log(393216) + 4*log_alpha + 7*log_beta + 10*log_theta,
      log(528384) + 3*log_alpha + 8*log_beta + 10*log_theta,
      log(417792) + 2*log_alpha + 9*log_beta + 10*log_theta,
      log(180224) + log_alpha + 10*log_beta + 10*log_theta,
      15*log(2) + 11*log_beta + 10*log_theta,
      log(5120) + 7*log_alpha + 3*log_beta + 9*log_theta,
      log(69632) + 6*log_alpha + 4*log_beta + 9*log_theta,
      log(372736) + 5*log_alpha + 5*log_beta + 9*log_theta,
      log(1042432) + 4*log_alpha + 6*log_beta + 9*log_theta,
      log(1666048) + 3*log_alpha + 7*log_beta + 9*log_theta,
      log(1533952) + 2*log_alpha + 8*log_beta + 9*log_theta,
      log(757760) + log_alpha + 9*log_beta + 9*log_theta,
      log(155648) + 10*log_beta + 9*log_theta,
      log(2432) + 7*log_alpha + 2*log_beta + 8*log_theta,
      log(47360) + 6*log_alpha + 3*log_beta + 8*log_theta,
      log(333056) + 5*log_alpha + 4*log_beta + 8*log_theta,
      log(1161728) + 4*log_alpha + 5*log_beta + 8*log_theta,
      log(2236288) + 3*log_alpha + 6*log_beta + 8*log_theta,
      log(2417920) + 2*log_alpha + 7*log_beta + 8*log_theta,
      log(1375744) + log_alpha + 8*log_beta + 8*log_theta,
      log(320512) + 9*log_beta + 8*log_theta,
      log(608) + 7*log_alpha + log_beta + 7*log_theta,
      log(17024) + 6*log_alpha + 2*log_beta + 7*log_theta,
      log(161408) + 5*log_alpha + 3*log_beta + 7*log_theta,
      log(717248) + 4*log_alpha + 4*log_beta + 7*log_theta,
      log(1688224) + 3*log_alpha + 5*log_beta + 7*log_theta,
      log(2165440) + 2*log_alpha + 6*log_beta + 7*log_theta,
      log(1428096) + log_alpha + 7*log_beta + 7*log_theta,
      log(378624) + 8*log_beta + 7*log_theta,
      7*log(2) + 7*log_alpha + 6*log_theta,
      log(4016) + 6*log_alpha + log_beta + 6*log_theta,
      log(48816) + 5*log_alpha + 2*log_beta + 6*log_theta,
      log(276400) + 4*log_alpha + 3*log_beta + 6*log_theta,
      log(802160) + 3*log_alpha + 4*log_beta + 6*log_theta,
      log(1228640) + 2*log_alpha + 5*log_beta + 6*log_theta,
      log(940640) + log_alpha + 6*log_beta + 6*log_theta,
      log(282176) + 7*log_beta + 6*log_theta,
      log(576) + 6*log_alpha + 5*log_theta,
      log(9936) + 5*log_alpha + log_beta + 5*log_theta,
      log(71760) + 4*log_alpha + 2*log_beta + 5*log_theta,
      log(256664) + 3*log_alpha + 3*log_beta + 5*log_theta,
      log(468672) + 2*log_alpha + 4*log_beta + 5*log_theta,
      log(413144) + log_alpha + 5*log_beta + 5*log_theta,
      log(137072) + 6*log_beta + 5*log_theta,
      log(1072) + 5*log_alpha + 4*log_theta,
      log(12376) + 4*log_alpha + log_beta + 4*log_theta,
      log(57224) + 3*log_alpha + 2*log_beta + 4*log_theta,
      log(125148) + 2*log_alpha + 3*log_beta + 4*log_theta,
      log(124716) + log_alpha + 4*log_beta + 4*log_theta,
      log(43296) + 5*log_beta + 4*log_theta,
      log(1064) + 4*log_alpha + 3*log_theta,
      log(8374) + 3*log_alpha + log_beta + 3*log_theta,
      log(23676) + 2*log_alpha + 2*log_beta + 3*log_theta,
      log(26570) + log_alpha + 3*log_beta + 3*log_theta,
      log(8568) + 4*log_beta + 3*log_theta,
      log(596) + 3*log_alpha + 2*log_theta,
      log(2943) + 2*log_alpha + log_beta + 2*log_theta,
      log(4059) + log_alpha + 2*log_beta + 2*log_theta,
      log(964) + 3*log_beta + 2*log_theta,
      log(178) + 2*log_alpha + log_theta,
      log(420) + log_alpha + log_beta + log_theta,
      log(47) + 2*log_beta + log_theta,
      log(22) + log_alpha});

const double log_p021 =
  -6*log(2)
  - 2*log(2*alpha*theta + 4*beta*theta + 1)
  - 3*log(2*alpha*theta + 2*beta*theta + 1)
  - log(alpha*theta + 2*beta*theta + 1)
  - log(alpha*theta + beta*theta + 1)
  - 3*log(4*beta*theta + 1)
  - log(2*beta*theta + 1)
  + log_beta
  + 2*log_theta
  + log_sum_exp(std::array<double,59>{
      12*log(2) + 7*log_alpha + 3*log_beta + 9*log_theta,
      log(40960) + 6*log_alpha + 4*log_beta + 9*log_theta,
      log(172032) + 5*log_alpha + 5*log_beta + 9*log_theta,
      log(393216) + 4*log_alpha + 6*log_beta + 9*log_theta,
      log(528384) + 3*log_alpha + 7*log_beta + 9*log_theta,
      log(417792) + 2*log_alpha + 8*log_beta + 9*log_theta,
      log(180224) + log_alpha + 9*log_beta + 9*log_theta,
      15*log(2) + 10*log_beta + 9*log_theta,
      log(5120) + 7*log_alpha + 2*log_beta + 8*log_theta,
      log(69632) + 6*log_alpha + 3*log_beta + 8*log_theta,
      log(372736) + 5*log_alpha + 4*log_beta + 8*log_theta,
      log(1042432) + 4*log_alpha + 5*log_beta + 8*log_theta,
      log(1666048) + 3*log_alpha + 6*log_beta + 8*log_theta,
      log(1533952) + 2*log_alpha + 7*log_beta + 8*log_theta,
      log(757760) + log_alpha + 8*log_beta + 8*log_theta,
      log(155648) + 9*log_beta + 8*log_theta,
      log(2432) + 7*log_alpha + log_beta + 7*log_theta,
      log(47360) + 6*log_alpha + 2*log_beta + 7*log_theta,
      log(333056) + 5*log_alpha + 3*log_beta + 7*log_theta,
      log(1161728) + 4*log_alpha + 4*log_beta + 7*log_theta,
      log(2236288) + 3*log_alpha + 5*log_beta + 7*log_theta,
      log(2417920) + 2*log_alpha + 6*log_beta + 7*log_theta,
      log(1375744) + log_alpha + 7*log_beta + 7*log_theta,
      log(320512) + 8*log_beta + 7*log_theta,
      log(480) + 7*log_alpha + 6*log_theta,
      log(15744) + 6*log_alpha + log_beta + 6*log_theta,
      log(156032) + 5*log_alpha + 2*log_beta + 6*log_theta,
      log(704960) + 4*log_alpha + 3*log_beta + 6*log_theta,
      log(1671712) + 3*log_alpha + 4*log_beta + 6*log_theta,
      log(2152384) + 2*log_alpha + 5*log_beta + 6*log_theta,
      log(1422464) + log_alpha + 6*log_beta + 6*log_theta,
      log(377600) + 7*log_beta + 6*log_theta,
      log(2160) + 6*log_alpha + 5*log_theta,
      log(38640) + 5*log_alpha + log_beta + 5*log_theta,
      log(247920) + 4*log_alpha + 2*log_beta + 5*log_theta,
      log(758064) + 3*log_alpha + 3*log_beta + 5*log_theta,
      log(1190624) + 2*log_alpha + 4*log_beta + 5*log_theta,
      log(923616) + log_alpha + 5*log_beta + 5*log_theta,
      log(279104) + 6*log_beta + 5*log_theta,
      log(4016) + 5*log_alpha + 4*log_theta,
      log(48016) + 4*log_alpha + log_beta + 4*log_theta,
      log(209784) + 3*log_alpha + 2*log_beta + 4*log_theta,
      log(421696) + 2*log_alpha + 3*log_beta + 4*log_theta,
      log(391000) + log_alpha + 4*log_beta + 4*log_theta,
      log(133424) + 5*log_beta + 4*log_theta,
      log(3976) + 4*log_alpha + 3*log_theta,
      log(32856) + 3*log_alpha + log_beta + 3*log_theta,
      log(93948) + 2*log_alpha + 2*log_beta + 3*log_theta,
      log(108348) + log_alpha + 3*log_beta + 3*log_theta,
      log(41088) + 4*log_beta + 3*log_theta,
      log(2238) + 3*log_alpha + 2*log_theta,
      log(12084) + 2*log_alpha + log_beta + 2*log_theta,
      log(19090) + log_alpha + 2*log_beta + 2*log_theta,
      log(7848) + 3*log_beta + 2*log_theta,
      log(687) + 2*log_alpha + log_theta,
      log(1959) + log_alpha + log_beta + log_theta,
      log(844) + 2*log_beta + log_theta,
      log(90) + log_alpha,
      log(39) + log_beta});

  const std::array<double, 64> p = {log_p000, log_p001, log_p002, log_p001, log_p010, log_p010, log_p012, log_p012, log_p020, log_p021, log_p020, log_p021, log_p010, log_p012, log_p012, log_p010, log_p010, log_p010, log_p012, log_p012, log_p001, log_p000, log_p001, log_p002, log_p012, log_p010, log_p010, log_p012, log_p021, log_p020, log_p021, log_p020, log_p020, log_p021, log_p020, log_p021, log_p012, log_p010, log_p010, log_p012, log_p002, log_p001, log_p000, log_p001, log_p012, log_p012, log_p010, log_p010, log_p010, log_p012, log_p012, log_p010, log_p021, log_p020, log_p021, log_p020, log_p012, log_p012, log_p010, log_p010, log_p001, log_p002, log_p001, log_p000};

  return p;
}


inline double k2p_species_tree_nll(const double alpha, const double beta, const double theta, const TuplePatterns& tp) {
  if (tp.order() != 3) {
    throw std::invalid_argument("k2p_neg_loglikelihood: expected 3-taxon TuplePatterns");
  }
  const auto log_ps = k2p_species_tree_probs(alpha, beta, theta);
  
  double neg_loglikelihood = 0.0;
  for (size_t i = 0; i < 64; ++i) {
    neg_loglikelihood -= tp.site_pattern[i] * log_ps[i];
  }
  return neg_loglikelihood;

}




// Permute a 3-taxon TuplePatterns so that new position i draws from old position pos_i.
// Used to place the cherry taxa in positions 0,1 for each of the 3 topologies.
inline TuplePatterns permute_taxa(const TuplePatterns& tp,
                                    const size_t pos0,
                                    const size_t pos1,
                                    const size_t pos2) {
	TuplePatterns result(std::vector<std::string>{tp.taxa[pos0], tp.taxa[pos1], tp.taxa[pos2]});
	result.n_total   = tp.n_total;
	result.n_counted = tp.n_counted;
	size_t bases[3];
	for (size_t a = 0; a < 4; ++a)
		for (size_t b = 0; b < 4; ++b)
			for (size_t c = 0; c < 4; ++c) {
				bases[pos0] = a; bases[pos1] = b; bases[pos2] = c;
				result.site_pattern[site_pattern_index(a, b, c)] =
				    tp.site_pattern[site_pattern_index(bases[0], bases[1], bases[2])];
			}
	return result;
}

// Estimate alpha and beta via k2p_mle, then evaluate k2p_species_tree_nll for
// each of the 3 topologies.  Returns normalized log-weights:
//   result[0] -- log P((taxa[0],taxa[1])|taxa[2] | data)
//   result[1] -- log P((taxa[0],taxa[2])|taxa[1] | data)
//   result[2] -- log P((taxa[1],taxa[2])|taxa[0] | data)
inline std::array<double, 3> compute_k2p_species_triplet_weights(
    const TuplePatterns& tp,
    const double theta,
    const double mu) {
	const auto mle   = k2p_mle(tp, theta, mu);
	const double alpha = mle.alpha_ml;
	const double beta  = mle.beta_ml;

  std::cout << "tp: " << tp.taxa[0] << ", " << tp.taxa[1] << ", " << tp.taxa[2] << std::endl;
  
	// topology 0: (taxa[0],taxa[1])|taxa[2] — cherry already in positions 0,1
	const double nll_0 = k2p_species_tree_nll(alpha, beta, theta, tp);

	// topology 1: (taxa[0],taxa[2])|taxa[1] — swap positions 1 and 2
	const double nll_1 = k2p_species_tree_nll(alpha, beta, theta, permute_taxa(tp, 0, 2, 1));

	// topology 2: (taxa[1],taxa[2])|taxa[0] — rotate: old 1->0, 2->1, 0->2
	const double nll_2 = k2p_species_tree_nll(alpha, beta, theta, permute_taxa(tp, 1, 2, 0));

	const std::array<double, 3> log_weights = {-nll_0, -nll_1, -nll_2};
	const double norm = log_sum_exp(log_weights);
	return {log_weights[0] - norm, log_weights[1] - norm, log_weights[2] - norm};
}

}  // namespace phylotriplet
