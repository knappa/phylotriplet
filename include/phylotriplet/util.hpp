#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

template <typename T, size_t N>
constexpr inline T log_sum_exp(std::array<T, N> values) {
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

constexpr inline double log_diff_exp(double log_a, double log_b) {
  if (log_a < log_b) {
    throw std::invalid_argument("log_diff_exp: log_a must be >= log_b, got log_a=" +
                                std::to_string(log_a) + " log_b=" + std::to_string(log_b));
  }
  return log_a + log1p(-exp(log_b - log_a));
}

constexpr size_t site_pattern_index(size_t a, size_t b, size_t c) {
  // a,b,c are ints in [0,3]
  return (a << 4) + (b << 2) + c;
}
