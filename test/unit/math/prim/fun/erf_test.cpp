#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>

TEST(MathFunctions, erf) {
  using stan::math::erf;
  EXPECT_FLOAT_EQ(-0.3286267594591274, erf(-0.3));
  EXPECT_FLOAT_EQ(0, erf(0));
  EXPECT_FLOAT_EQ(0.9999939742388482, erf(3.2));
}

TEST(MathFunctions, erfOverflow) {
  EXPECT_FLOAT_EQ(-1, erf(-100));
  EXPECT_FLOAT_EQ(1, erf(100));
}

TEST(MathFunctions, erfNan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(std::isnan(stan::math::erf(nan)));
}

TEST(MathFunctions, erf_works_with_other_functions) {
  Eigen::VectorXd a(5);
  a << 1.1, 1.2, 1.3, 1.4, 1.5;
  Eigen::RowVectorXd b(5);
  b << 1.1, 1.2, 1.3, 1.4, 1.5;
  stan::math::multiply(a, stan::math::erf(b));
}
