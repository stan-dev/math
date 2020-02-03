#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>

TEST(MathFunctions, cbrt) {
  using stan::math::cbrt;
  EXPECT_FLOAT_EQ(-2.0, cbrt(-8.0));
  EXPECT_FLOAT_EQ(-1.392476650083834, cbrt(-2.7));
  EXPECT_FLOAT_EQ(0, cbrt(0));
  EXPECT_FLOAT_EQ(2.0, cbrt(8.0));
}

TEST(MathFunctions, cbrt_inf_return) {
  EXPECT_EQ(-std::numeric_limits<double>::infinity(),
            stan::math::cbrt(-std::numeric_limits<double>::infinity()));
  EXPECT_EQ(std::numeric_limits<double>::infinity(),
            stan::math::cbrt(std::numeric_limits<double>::infinity()));
}

TEST(MathFunctions, cbrt_nan) {
  using stan::math::cbrt;
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(std::isnan(stan::math::cbrt(nan)));
}

TEST(MathFunctions, cbrt_works_with_other_functions) {
  Eigen::VectorXd a(5);
  a << 1.1, 1.2, 1.3, 1.4, 1.5;
  Eigen::RowVectorXd b(5);
  b << 1.1, 1.2, 1.3, 1.4, 1.5;
  stan::math::multiply(a, stan::math::cbrt(b));
}
