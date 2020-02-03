#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, inv_cloglog) {
  EXPECT_EQ(1 - std::exp(-std::exp(3.7)), stan::math::inv_cloglog(3.7));
  EXPECT_EQ(1 - std::exp(-std::exp(0.0)), stan::math::inv_cloglog(0.0));
  EXPECT_EQ(1 - std::exp(-std::exp(-2.93)), stan::math::inv_cloglog(-2.93));
}

TEST(MathFunctions, inv_cloglog_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::inv_cloglog(nan)));
}

TEST(MathFunctions, inv_cloglog_works_with_other_functions) {
  Eigen::VectorXd a(5);
  a << 1.1, 1.2, 1.3, 1.4, 1.5;
  Eigen::RowVectorXd b(5);
  b << 1.1, 1.2, 1.3, 1.4, 1.5;
  stan::math::multiply(a, stan::math::inv_cloglog(b));
}
