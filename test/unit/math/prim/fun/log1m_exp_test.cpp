#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, log1m_exp) {
  using stan::math::log1m_exp;

  // exp(10000.0) overflows
  EXPECT_FLOAT_EQ(0, log1m_exp(-1e10));
  EXPECT_FLOAT_EQ(0, log1m_exp(-1000));
  EXPECT_FLOAT_EQ(-3.720076e-44, log1m_exp(-100));
  EXPECT_FLOAT_EQ(-4.540096e-05, log1m_exp(-10));
  EXPECT_FLOAT_EQ(-0.4586751, log1m_exp(-1));
  EXPECT_FLOAT_EQ(-2.352168, log1m_exp(-0.1));
  EXPECT_FLOAT_EQ(-11.51293, log1m_exp(-1e-5));
  EXPECT_FLOAT_EQ(-23.02585, log1m_exp(-1e-10));
  EXPECT_FLOAT_EQ(-46.0517, log1m_exp(-1e-20));
  EXPECT_FLOAT_EQ(-92.1034, log1m_exp(-1e-40));
  EXPECT_NO_THROW(log1m_exp(1));
  EXPECT_FLOAT_EQ(log1m_exp(0), -stan::math::INFTY);
  EXPECT_TRUE(std::isnan(log1m_exp(0.001)));
}

TEST(MathFunctions, log1m_exp_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::log1m_exp(nan)));
}
