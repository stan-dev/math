#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, trunc) {
  using stan::math::trunc;
  EXPECT_FLOAT_EQ(-27, trunc(-27.8239));
  EXPECT_FLOAT_EQ(-1, trunc(-1.5));
  EXPECT_FLOAT_EQ(0, trunc(0));
  EXPECT_FLOAT_EQ(0, trunc(0.0));
  EXPECT_FLOAT_EQ(0, trunc(0.5));
  EXPECT_FLOAT_EQ(27, trunc(27.3239));
}

TEST(MathFunctions, truncNaN) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(std::isnan(stan::math::trunc(nan)));
}

TEST(MathFunctions, trunc_works_with_other_functions) {
  Eigen::VectorXd a(5);
  a << 1.1, 1.2, 1.3, 1.4, 1.5;
  Eigen::RowVectorXd b(5);
  b << 1.1, 1.2, 1.3, 1.4, 1.5;
  stan::math::multiply(a, stan::math::trunc(b));
}
