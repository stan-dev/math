#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, square) {
  double y = 2.0;
  EXPECT_FLOAT_EQ(y * y, stan::math::square(y));

  y = 0.0;
  EXPECT_FLOAT_EQ(y * y, stan::math::square(y));

  y = -32.7;
  EXPECT_FLOAT_EQ(y * y, stan::math::square(y));
}

TEST(MathFunctions, square_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::square(nan)));
}

TEST(MathFunctions, square_works_with_other_functions) {
  Eigen::VectorXd a(5);
  Eigen::RowVectorXd b(5);
  stan::math::multiply(a, stan::math::square(b));
}
