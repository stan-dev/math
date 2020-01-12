#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

template <typename T_N, typename T_n>
void test_binom_coefficient(const T_N& N, const T_n& n) {
  using stan::math::binomial_coefficient_log;
  EXPECT_FLOAT_EQ(lgamma(N + 1) - lgamma(n + 1) - lgamma(N - n + 1),
                  binomial_coefficient_log(N, n));
}

TEST(MathFunctions, binomial_coefficient_log) {
  using stan::math::binomial_coefficient_log;
  EXPECT_FLOAT_EQ(1.0, exp(binomial_coefficient_log(2.0, 2.0)));
  EXPECT_FLOAT_EQ(2.0, exp(binomial_coefficient_log(2.0, 1.0)));
  EXPECT_FLOAT_EQ(3.0, exp(binomial_coefficient_log(3.0, 1.0)));
  EXPECT_NEAR(3.0, exp(binomial_coefficient_log(3.0, 2.0)), 0.0001);

  EXPECT_FLOAT_EQ(29979.16, binomial_coefficient_log(100000, 91116));

  for (int n = 0; n < 1010; ++n) {
    test_binom_coefficient(1010, n);
    test_binom_coefficient(1010.0, n);
    test_binom_coefficient(1010, static_cast<double>(n));
    test_binom_coefficient(1010.0, static_cast<double>(n));
  }

  test_binom_coefficient(1e9, 1e5);
  test_binom_coefficient(1e50, 1e45);
  test_binom_coefficient(1e20, 1e15);
}

TEST(MathFunctions, binomial_coefficient_log_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::binomial_coefficient_log(2.0, nan)));
  EXPECT_TRUE(std::isnan(stan::math::binomial_coefficient_log(nan, 2.0)));
  EXPECT_TRUE(std::isnan(stan::math::binomial_coefficient_log(nan, nan)));
}
