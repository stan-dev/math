#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathsSpecialFunctions, inv) {
  double y = 2.0;
  EXPECT_FLOAT_EQ(1 / y, stan::math::inv(y));

  y = 0.0;
  EXPECT_FLOAT_EQ(stan::math::positive_infinity(), stan::math::inv(y));

  y = -32.7;
  EXPECT_FLOAT_EQ(1 / y, stan::math::inv(y));
}

TEST(MathFunctions, inv_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::inv(nan)));
}

TEST(MathFunctions, inv_works_with_other_functions) {
  Eigen::VectorXd a(5);
  a << 1.1, 1.2, 1.3, 1.4, 1.5;
  Eigen::RowVectorXd b(5);
  b << 1.1, 1.2, 1.3, 1.4, 1.5;
  stan::math::multiply(a, stan::math::inv(b));
}
