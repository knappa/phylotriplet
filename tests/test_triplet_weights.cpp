#define CATCH_CONFIG_MAIN
#include <Eigen/Eigenvalues>
#include <cmath>
#include <limits>

#include "catch.hpp"
#include "phylotriplet/triplet_weights.hpp"

using namespace phylotriplet;

// --- log_sum_exp -------------------------------------------------------------

TEST_CASE("log_sum_exp: all -infinity returns -infinity") {
  const double neginf = -std::numeric_limits<double>::infinity();
  std::array<double, 3> vals = {neginf, neginf, neginf};
  double result = log_sum_exp(vals);
  REQUIRE(result == neginf);
}

TEST_CASE("log_sum_exp: single -infinity among finite values") {
  const double neginf = -std::numeric_limits<double>::infinity();
  // log_sum_exp({-inf, 0, 0}) = log(0 + 1 + 1) = log(2)
  std::array<double, 3> vals = {neginf, 0.0, 0.0};
  double result = log_sum_exp(vals);
  REQUIRE(result == Approx(std::log(2.0)));
}

TEST_CASE("log_sum_exp: multiple -infinity among finite values") {
  const double neginf = -std::numeric_limits<double>::infinity();
  // log_sum_exp({-inf, -inf, 1.0}) = 1.0
  std::array<double, 3> vals = {neginf, neginf, 1.0};
  double result = log_sum_exp(vals);
  REQUIRE(result == Approx(1.0));
}

TEST_CASE("log_sum_exp: equal finite values") {
  // log_sum_exp({a, a, a}) = a + log(3)
  std::array<double, 3> vals = {2.0, 2.0, 2.0};
  double result = log_sum_exp(vals);
  REQUIRE(result == Approx(2.0 + std::log(3.0)));
}

TEST_CASE("log_sum_exp: known values") {
  // log(exp(1) + exp(2)) = log(e + e^2) = 2 + log(1 + exp(-1))
  std::array<double, 2> vals = {1.0, 2.0};
  double result = log_sum_exp(vals);
  REQUIRE(result == Approx(2.0 + std::log1p(std::exp(-1.0))));
}

TEST_CASE("log_sum_exp: size 1 returns the value") {
  std::array<double, 1> vals = {42.0};
  double result = log_sum_exp(vals);
  REQUIRE(result == Approx(42.0));
}

TEST_CASE("log_sum_exp: size 1 -infinity returns -infinity") {
  const double neginf = -std::numeric_limits<double>::infinity();
  std::array<double, 1> vals = {neginf};
  double result = log_sum_exp(vals);
  REQUIRE(result == neginf);
}

TEST_CASE("log_sum_exp: large magnitude values do not overflow") {
  // log(exp(1000) + exp(1001)) = 1001 + log(1 + exp(-1))
  std::array<double, 2> vals = {1000.0, 1001.0};
  double result = log_sum_exp(vals);
  REQUIRE(result == Approx(1001.0 + std::log1p(std::exp(-1.0))));
}

TEST_CASE("log_sum_exp: very negative values") {
  // log(exp(-1000) + exp(-999)) = -999 + log(1 + exp(-1))
  std::array<double, 2> vals = {-1000.0, -999.0};
  double result = log_sum_exp(vals);
  REQUIRE(result == Approx(-999.0 + std::log1p(std::exp(-1.0))));
}
