#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, log1p_exp) {
  using stan::math::log1p_exp;

  // exp(10000.0) overflows
  EXPECT_FLOAT_EQ(10000.0, log1p_exp(10000.0));
  EXPECT_FLOAT_EQ(0.0, log1p_exp(-10000.0));
}

TEST(MathFunctions, log1p_exp_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::log1p_exp(nan)));
}

TEST(MathFunctions, log1p_exp_works_with_other_functions) {
  Eigen::VectorXd a(5);
  a << 1.1, 1.2, 1.3, 1.4, 1.5;
  Eigen::RowVectorXd b(5);
  b << 1.1, 1.2, 1.3, 1.4, 1.5;
  stan::math::multiply(a, stan::math::log1p_exp(b));
}
