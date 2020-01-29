#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, round) {
  using stan::math::round;
  EXPECT_FLOAT_EQ(-27, round(-27.3239));
  EXPECT_FLOAT_EQ(-1, round(-0.5));
  EXPECT_FLOAT_EQ(0, round(0));
  EXPECT_FLOAT_EQ(0, round(0.0));
  EXPECT_FLOAT_EQ(1, round(0.5));
  EXPECT_FLOAT_EQ(27, round(27.3239));
}

TEST(MathFunctions, roundNaN) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(std::isnan(stan::math::round(nan)));
}

TEST(MathFunctions, round_works_with_other_functions) {
  Eigen::VectorXd a(5);
  a << 1.1, 1.2, 1.3, 1.4, 1.5;
  Eigen::RowVectorXd b(5);
  b << 1.1, 1.2, 1.3, 1.4, 1.5;
  stan::math::multiply(a, stan::math::round(b));
}
