#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, log2_fun) {
  using stan::math::log2;

  EXPECT_FLOAT_EQ(std::log(2.0), log2());
}

TEST(MathFunctions, log2) {
  using stan::math::log2;

  EXPECT_FLOAT_EQ(1.0, log2(2.0));
  EXPECT_FLOAT_EQ(2.0, log2(4.0));
  EXPECT_FLOAT_EQ(3.0, log2(8.0));
}

TEST(MathFunctions, log2_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::log2(nan)));
}

TEST(MathFunctions, log2_works_with_other_functions) {
  Eigen::VectorXd a(5);
  a << 1.1, 1.2, 1.3, 1.4, 1.5;
  Eigen::RowVectorXd b(5);
  b << 1.1, 1.2, 1.3, 1.4, 1.5;
  stan::math::multiply(a, stan::math::log2(b));
}
