#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, inv_sqrt) {
  double y = 4.0;
  EXPECT_FLOAT_EQ(1 / 2.0, stan::math::inv_sqrt(y));

  y = 25.0;
  EXPECT_FLOAT_EQ(1 / 5.0, stan::math::inv_sqrt(y));

  y = 0.0;
  EXPECT_FLOAT_EQ(stan::math::positive_infinity(), stan::math::inv_sqrt(y));

  y = -50.0;
  std::isnan(stan::math::inv_sqrt(y));
}

TEST(MathFunctions, inv_sqrt_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::inv_sqrt(nan)));
}

TEST(MathFunctions, inv_sqrt_works_with_other_functions) {
  Eigen::VectorXd a(5);
  a << 1.1, 1.2, 1.3, 1.4, 1.5;
  Eigen::RowVectorXd b(5);
  b << 1.1, 1.2, 1.3, 1.4, 1.5;
  stan::math::multiply(a, stan::math::inv_sqrt(b));
}
