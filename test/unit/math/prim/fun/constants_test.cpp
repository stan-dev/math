#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathsConstants, sqrt2) {
  EXPECT_FLOAT_EQ(std::sqrt(2.0), stan::math::SQRT_TWO);
}
TEST(MathsConstants, log2) {
  EXPECT_FLOAT_EQ(std::log(2.0), stan::math::LOG_TWO);
}
TEST(MathsConstants, log10) {
  EXPECT_FLOAT_EQ(std::log(10.0), stan::math::LOG_TEN);
}

TEST(MathsConstants, infty) {
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), stan::math::INFTY);
}
TEST(MathsConstants, neg_infty) {
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(),
                  stan::math::NEGATIVE_INFTY);
}
TEST(MathsConstants, nan) { EXPECT_TRUE(std::isnan(stan::math::NOT_A_NUMBER)); }

TEST(MathsConstants, epsilon) { EXPECT_TRUE(stan::math::EPSILON > 0.0); }
TEST(MathsConstants, pi_fun) {
  EXPECT_FLOAT_EQ(4.0 * std::atan(1.0), stan::math::pi());
}
TEST(MathsConstants, e_fun) { EXPECT_FLOAT_EQ(std::exp(1.0), stan::math::e()); }
TEST(MathsConstants, sqrt2_fun) {
  EXPECT_FLOAT_EQ(std::sqrt(2.0), stan::math::sqrt2());
}
TEST(MathsConstants, log10_fun) {
  EXPECT_FLOAT_EQ(std::log(10.0), stan::math::log10());
}

TEST(MathsConstants, infty_fun) {
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(),
                  stan::math::positive_infinity());
}
TEST(MathsConstants, neg_infty_fun) {
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(),
                  stan::math::negative_infinity());
}
TEST(MathsConstants, not_a_number_fun) {
  EXPECT_TRUE(std::isnan(stan::math::not_a_number()));
}
TEST(MathsConstants, machine_precision_fun) {
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::epsilon(),
                  stan::math::machine_precision());
}
TEST(MathsConstants, two_over_sqrt_pi) {
  EXPECT_FLOAT_EQ(1.128379167, stan::math::TWO_OVER_SQRT_PI);
}
TEST(MathsConstants, inv_sqrt_two_pi) {
  EXPECT_FLOAT_EQ(.39894228040, stan::math::INV_SQRT_TWO_PI);
}
TEST(MathsConstants, egamma) {
  EXPECT_FLOAT_EQ(-stan::math::digamma(1), stan::math::egamma());
}

TEST(MathsConstants, derivedConstantsCorrectlyRounded) {
  // Correctly rounded doubles of log(2 pi), log(2 pi) / 2 and sqrt(2)/sqrt(pi),
  // computed with mpmath at 60 decimal digits.
  EXPECT_EQ(1.8378770664093456, stan::math::LOG_TWO_PI);
  EXPECT_EQ(0.9189385332046728, stan::math::HALF_LOG_TWO_PI);
  EXPECT_EQ(0.7978845608028654, stan::math::SQRT_TWO_OVER_SQRT_PI);
  // Negation and doubling are exact, so the three log constants agree bitwise.
  EXPECT_EQ(-stan::math::HALF_LOG_TWO_PI, stan::math::NEG_LOG_SQRT_TWO_PI);
  EXPECT_EQ(2 * stan::math::HALF_LOG_TWO_PI, stan::math::LOG_TWO_PI);
}
